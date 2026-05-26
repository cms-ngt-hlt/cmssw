// Analyzer for validation histograms for tau objects at HLT/RECO
// E. Vernazza Apr. 10, 2026

#include "Validation/RecoTau/interface/TauValidationRECO.h"
#include "FWCore/Framework/interface/MakerMacros.h"

using namespace edm;
using namespace reco;
using namespace std;

std::string TauValidationRECO::convertId(double cut) {
  if (cut == 0.0) return "0p0";
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(1) << cut;
  std::string result = oss.str();
  for (char& c : result) {
    if (c == '.') c = 'p';
  }
  return result;
}

bool TauValidationRECO::passIdCut(const reco::PFTauRef& tauRef, const std::vector<double> &cutIDs,
  const reco::TauDiscriminatorContainer* recoTauIDjet,
  const reco::TauDiscriminatorContainer* recoTauIDe,
  const reco::TauDiscriminatorContainer* recoTauIDmu) {

  if (cutIDs[0] == 0.0 && cutIDs[1] == 0.0 && cutIDs[2] == 0.0) return true;

  if (cutIDs[0] != 0.0) { // vs Jet ID cut
    if (!recoTauIDjet) return false;
    const auto &jetDisc = (*recoTauIDjet)[tauRef];
    if (jetDisc.rawValues.empty() || jetDisc.rawValues[0] <= cutIDs[0]) return false;
  }
  if (cutIDs[1] != 0.0) { // vs E ID cut
    if (!recoTauIDe) return false;
    const auto &eDisc = (*recoTauIDe)[tauRef];
    if (eDisc.rawValues.empty() || eDisc.rawValues[0] <= cutIDs[1]) return false;
  }
  if (cutIDs[2] != 0.0) { // vs Mu ID cut
    if (!recoTauIDmu) return false;
    const auto &muDisc = (*recoTauIDmu)[tauRef];
    if (muDisc.rawValues.empty() || muDisc.rawValues[0] <= cutIDs[2]) return false;
  }

  return true;
}

TauValidationRECO::TauValidationRECO(const edm::ParameterSet& iConfig) {
  genTauToken_ = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genTauCollection"));
  recoTauToken_ = consumes<reco::PFTauCollection>(iConfig.getParameter<edm::InputTag>("recoTauCollection"));
  recoTauIDjetToken_ = consumes<reco::TauDiscriminatorContainer>(iConfig.getParameter<edm::InputTag>("recoTauIDCollectionVsJet"));
  recoTauIDeToken_ = consumes<reco::TauDiscriminatorContainer>(iConfig.getParameter<edm::InputTag>("recoTauIDCollectionVsE"));
  recoTauIDmuToken_ = consumes<reco::TauDiscriminatorContainer>(iConfig.getParameter<edm::InputTag>("recoTauIDCollectionVsMu"));
  cutIDs = iConfig.getParameter<std::vector<double>>("cutIDs"); // expects 3 entries: vsJet, vsE, vsMu
  matchingDeltaR = iConfig.getParameter<double>("minDeltaR");
  outFolder = iConfig.getParameter<std::string>("outFolder");
  isHLT = iConfig.getUntrackedParameter<bool>("isHLT");

  if (cutIDs.size() < 3) cutIDs.resize(3, 0.0);  // fills missing entries with 0.0

}

void TauValidationRECO::bookHistograms(DQMStore::IBooker& ibooker, edm::Run const& iRun, edm::EventSetup const&) {

  // ---------------------------- Book Summary Histograms -------------------------------
  applyIdCuts = !cutIDs.empty() && !(cutIDs[0] == 0.0 && cutIDs[1] == 0.0 && cutIDs[2] == 0.0);
  if (applyIdCuts) {
    outFolder += "_IdCuts";
    outFolder += cutIDs[0] != 0.0 ? "_jet" + convertId(cutIDs[0]) : "";
    outFolder += cutIDs[1] != 0.0 ? "_e" + convertId(cutIDs[1]) : "";
    outFolder += cutIDs[2] != 0.0 ? "_mu" + convertId(cutIDs[2]) : "";
  }
  ibooker.setCurrentFolder(outFolder);
  
  // Book 1D histograms for gen and reco tau kinematics
  for (auto& hVar : histoVars) {
    auto [nBins, hMin, hMax] = hVar.second;
    h_recoTau_[hVar.first] = ibooker.book1D("recoTau_" + hVar.first, ";#tau^{reco};" + hVar.first, nBins, hMin, hMax);
    h_recoTauMatched_[hVar.first] = ibooker.book1D("recoTauMatched_" + hVar.first, ";#tau^{reco} (Matched);" + hVar.first, nBins, hMin, hMax);
    h_recoTauMultiMatched_[hVar.first] = ibooker.book1D("recoTauMultiMatched_" + hVar.first, ";#tau^{reco} (Multi-Matched);" + hVar.first, nBins, hMin, hMax);
    h_genTau_[hVar.first] = ibooker.book1D("genTau_" + hVar.first, "#tau^{gen};" + hVar.first, nBins, hMin, hMax);
    h_genTauMatched_[hVar.first] = ibooker.book1D("genTauMatched_" + hVar.first, ";#tau^{gen} (Matched);" + hVar.first, nBins, hMin, hMax);
    h_genTauMultiMatched_[hVar.first] = ibooker.book1D("genTauMultiMatched_" + hVar.first, ";#tau^{gen} (Multi-Matched);" + hVar.first, nBins, hMin, hMax);
    h2d_responsePt_[hVar.first] = ibooker.book2D("responsePt_" + hVar.first, ";#tau Pt Response;" + hVar.first, nBins, hMin, hMax, 50, 0., 2.);
    h2d_responseMass_[hVar.first] = ibooker.book2D("responseMass_" + hVar.first, ";#tau Mass Response;" + hVar.first, nBins, hMin, hMax, 50, 0., 2.);
  }

  // Book 1D histograms for reco tau ID discriminators
  for (auto& hVar : histoVarsReco) {
    auto [nBins, hMin, hMax] = hVar.second;
    h_recoTau_[hVar.first] = ibooker.book1D("recoTau_" + hVar.first, ";#tau^{reco};" + hVar.first, nBins, hMin, hMax);
    h_recoTauMatched_[hVar.first] = ibooker.book1D("recoTauMatched_" + hVar.first, ";#tau^{reco} (Matched);" + hVar.first, nBins, hMin, hMax);
    h_recoTauMultiMatched_[hVar.first] = ibooker.book1D("recoTauMultiMatched_" + hVar.first, ";#tau^{reco} (Multi-Matched);" + hVar.first, nBins, hMin, hMax);
  }

  // Book 2D histograms for gen and reco tau kinematics
  for (auto& h2dVar : histoVars2D) {
    auto [nBinsX, hMinX, hMaxX, nBinsY, hMinY, hMaxY] = h2dVar.second;
    auto x_title = h2dVar.first.substr(0, h2dVar.first.find("_"));
    auto y_title = h2dVar.first.substr(h2dVar.first.find("_") + 1);
    h2d_recoTau_[h2dVar.first] = ibooker.book2D("recoTau_" + h2dVar.first, ";#tau^{reco}" + x_title + ";" + y_title, nBinsX, hMinX, hMaxX, nBinsY, hMinY, hMaxY);
    h2d_recoTauMatched_[h2dVar.first] = ibooker.book2D("recoTauMatched_" + h2dVar.first, ";#tau^{reco} (Matched);" + x_title + ";" + y_title, nBinsX, hMinX, hMaxX, nBinsY, hMinY, hMaxY);
    h2d_recoTauMultiMatched_[h2dVar.first] = ibooker.book2D("recoTauMultiMatched_" + h2dVar.first, ";#tau^{reco} (Multi-Matched);" + x_title + ";" + y_title, nBinsX, hMinX, hMaxX, nBinsY, hMinY, hMaxY);
    h2d_genTau_[h2dVar.first] = ibooker.book2D("genTau_" + h2dVar.first, ";#tau^{gen}" + x_title + ";" + y_title, nBinsX, hMinX, hMaxX, nBinsY, hMinY, hMaxY);
    h2d_genTauMatched_[h2dVar.first] = ibooker.book2D("genTauMatched_" + h2dVar.first, ";#tau^{gen} (Matched);" + x_title + ";" + y_title, nBinsX, hMinX, hMaxX, nBinsY, hMinY, hMaxY);
    h2d_genTauMultiMatched_[h2dVar.first] = ibooker.book2D("genTauMultiMatched_" + h2dVar.first, ";#tau^{gen} (Multi-Matched);" + x_title + ";" + y_title, nBinsX, hMinX, hMaxX, nBinsY, hMinY, hMaxY);
  }

  // Book 2D histograms for reco tau ID discriminators vs kinematics
  for (auto& h2dVar : histoVarsReco2D) {
    auto [nBinsX, hMinX, hMaxX, nBinsY, hMinY, hMaxY] = h2dVar.second;
    auto x_title = h2dVar.first.substr(0, h2dVar.first.find("_"));
    auto y_title = h2dVar.first.substr(h2dVar.first.find("_") + 1);
    h2d_recoTau_[h2dVar.first] = ibooker.book2D("recoTau_" + h2dVar.first, ";#tau^{reco}" + x_title + ";" + y_title, nBinsX, hMinX, hMaxX, nBinsY, hMinY, hMaxY);
    h2d_recoTauMatched_[h2dVar.first] = ibooker.book2D("recoTauMatched_" + h2dVar.first, ";#tau^{reco} (Matched);" + x_title + ";" + y_title, nBinsX, hMinX, hMaxX, nBinsY, hMinY, hMaxY);
    h2d_recoTauMultiMatched_[h2dVar.first] = ibooker.book2D("recoTauMultiMatched_" + h2dVar.first, ";#tau^{reco} (Multi-Matched);" + x_title + ";" + y_title, nBinsX, hMinX, hMaxX, nBinsY, hMinY, hMaxY);
  }

}

//------------------------------------------------------------------------------
// ~TauValidationRECO
//------------------------------------------------------------------------------
TauValidationRECO::~TauValidationRECO() {}

//------------------------------------------------------------------------------
// analyze
//------------------------------------------------------------------------------
void TauValidationRECO::analyze(const edm::Event& mEvent, const edm::EventSetup& mSetup) {

  edm::Handle<reco::PFTauCollection> recoTaus;
  mEvent.getByToken(recoTauToken_, recoTaus);
  if (!recoTaus.isValid()) {
    edm::LogPrint("TauValidationRECO") << " Reco Tau collection not found while running TauValidationRECO.cc ";
    return;
  }
  // std::cout << "Number of reco taus: " << recoTaus->size() << std::endl; // [DEBUG]

  edm::Handle<reco::GenJetCollection> genTaus;
  mEvent.getByToken(genTauToken_, genTaus);
  if (!genTaus.isValid()) {
    edm::LogPrint("TauValidationRECO") << " Gen Tau collection not found while running TauValidationRECO.cc ";
    return;
  }
  // std::cout << "Number of gen taus: " << genTaus->size() << std::endl; // [DEBUG]

  bool useIdvsJet = true;
  edm::Handle<reco::TauDiscriminatorContainer> recoTauIDjet;
  mEvent.getByToken(recoTauIDjetToken_, recoTauIDjet);
  if (!recoTauIDjet.isValid()) {
    edm::LogPrint("TauValidationRECO") << "Reco Tau Identifier vs Jet collection not found while running TauValidationRECO.cc ";
    useIdvsJet = false;
  }

  bool useIdvsE = true;
  edm::Handle<reco::TauDiscriminatorContainer> recoTauIDe;
  mEvent.getByToken(recoTauIDeToken_, recoTauIDe);
  if (!recoTauIDe.isValid()) {
    edm::LogPrint("TauValidationRECO") << "Reco Tau Identifier vs Electrons collection not found while running TauValidationRECO.cc ";
    useIdvsE = false;
  }

  bool useIdvsMu = true;
  edm::Handle<reco::TauDiscriminatorContainer> recoTauIDmu;
  mEvent.getByToken(recoTauIDmuToken_, recoTauIDmu);
  if (!recoTauIDmu.isValid()) {
    edm::LogPrint("TauValidationRECO") << "Reco Tau Identifier vs Muons collection not found while running TauValidationRECO.cc ";
    useIdvsMu = false;
  }

  // if (useIdvsJet && useIdvsE && useIdvsMu) { // [DEBUG]
  //   for (size_t i = 0; i < recoTaus->size(); ++i) {
  //     reco::PFTauRef tauRef(recoTaus, i);
  //     const reco::SingleTauDiscriminatorContainer& discContainerJet = (*recoTauIDjet)[tauRef];
  //     const reco::SingleTauDiscriminatorContainer& discContainerE = (*recoTauIDe)[tauRef];
  //     const reco::SingleTauDiscriminatorContainer& discContainerMu = (*recoTauIDmu)[tauRef];
  //     std::cout << "Tau " << i << " ID vs Jet = " << discContainerJet.rawValues[0] <<
  //                   ", ID vs E = " << discContainerE.rawValues[0] <<
  //                   ", ID vs Mu = " << discContainerMu.rawValues[0] << std::endl;
  //     // for (size_t j = 0; j < discContainerJet.workingPoints.size(); ++j) {
  //     //   std::cout << "Tau " << i << " WP[" << j << "] = " << discContainerJet.workingPoints[j] << std::endl;
  //     // }
  //   }
  // }

  // Loop for efficiency 
  for (uint itau = 0; itau < genTaus->size(); ++itau) {
    
    h_genTau_["pt"]->Fill(genTaus->at(itau).pt());
    h_genTau_["eta"]->Fill(genTaus->at(itau).eta());
    h_genTau_["phi"]->Fill(genTaus->at(itau).phi());
    h_genTau_["mass"]->Fill(genTaus->at(itau).mass());
    h2d_genTau_["pt_eta"]->Fill(genTaus->at(itau).pt(), genTaus->at(itau).eta());
    h2d_genTau_["pt_phi"]->Fill(genTaus->at(itau).pt(), genTaus->at(itau).phi());
    h2d_genTau_["pt_mass"]->Fill(genTaus->at(itau).pt(), genTaus->at(itau).mass());
    h2d_genTau_["mass_eta"]->Fill(genTaus->at(itau).mass(), genTaus->at(itau).eta());
    h2d_genTau_["mass_phi"]->Fill(genTaus->at(itau).mass(), genTaus->at(itau).phi());

    // Count how many reco taus are matched to the gen tau
    int nRecoMatchedToOneGen = 0;
    float bestDeltaR = 999.;
    float ResponsePt_bestDeltaR = 0.;
    float ResponseMass_bestDeltaR = 0.;
    for (uint jtau = 0; jtau < recoTaus->size(); ++jtau) {
      reco::PFTauRef tauRef(recoTaus, jtau);
      if (applyIdCuts && !passIdCut(tauRef, cutIDs, recoTauIDjet.product(), recoTauIDe.product(), recoTauIDmu.product())) {
        continue; // skip if tau doesn't pass ID cuts
      }
      float deltaRValue = deltaR(genTaus->at(itau), recoTaus->at(jtau));
      if (deltaRValue < matchingDeltaR) {
        nRecoMatchedToOneGen++;
        if (deltaRValue < bestDeltaR) {
          bestDeltaR = deltaRValue;
          ResponsePt_bestDeltaR = recoTaus->at(jtau).pt() / genTaus->at(itau).pt();
          ResponseMass_bestDeltaR = recoTaus->at(jtau).mass() / genTaus->at(itau).mass();
        }
      }
    }

    // Fill histograms for gen taus matched to at least one reco tau
    if (nRecoMatchedToOneGen > 0) {
      // Fill gen tau histograms for matched taus
      h_genTauMatched_["pt"]->Fill(genTaus->at(itau).pt());
      h_genTauMatched_["eta"]->Fill(genTaus->at(itau).eta());
      h_genTauMatched_["phi"]->Fill(genTaus->at(itau).phi());
      h_genTauMatched_["mass"]->Fill(genTaus->at(itau).mass());
      h2d_genTauMatched_["pt_eta"]->Fill(genTaus->at(itau).pt(), genTaus->at(itau).eta());
      h2d_genTauMatched_["pt_phi"]->Fill(genTaus->at(itau).pt(), genTaus->at(itau).phi());
      h2d_genTauMatched_["pt_mass"]->Fill(genTaus->at(itau).pt(), genTaus->at(itau).mass());
      h2d_genTauMatched_["mass_eta"]->Fill(genTaus->at(itau).mass(), genTaus->at(itau).eta());
      h2d_genTauMatched_["mass_phi"]->Fill(genTaus->at(itau).mass(), genTaus->at(itau).phi());
      // Fill response histograms for matched taus
      h2d_responsePt_["pt"]->Fill(genTaus->at(itau).pt(), ResponsePt_bestDeltaR);
      h2d_responsePt_["eta"]->Fill(genTaus->at(itau).eta(), ResponsePt_bestDeltaR);
      h2d_responsePt_["phi"]->Fill(genTaus->at(itau).phi(), ResponsePt_bestDeltaR);
      h2d_responsePt_["mass"]->Fill(genTaus->at(itau).mass(), ResponsePt_bestDeltaR);
      h2d_responseMass_["pt"]->Fill(genTaus->at(itau).pt(), ResponseMass_bestDeltaR);
      h2d_responseMass_["eta"]->Fill(genTaus->at(itau).eta(), ResponseMass_bestDeltaR);
      h2d_responseMass_["phi"]->Fill(genTaus->at(itau).phi(), ResponseMass_bestDeltaR);
      h2d_responseMass_["mass"]->Fill(genTaus->at(itau).mass(), ResponseMass_bestDeltaR);

      if (nRecoMatchedToOneGen > 1) {
        // Fill gen tau histograms for multi-matched taus
        h_genTauMultiMatched_["pt"]->Fill(genTaus->at(itau).pt());
        h_genTauMultiMatched_["eta"]->Fill(genTaus->at(itau).eta());
        h_genTauMultiMatched_["phi"]->Fill(genTaus->at(itau).phi());
        h_genTauMultiMatched_["mass"]->Fill(genTaus->at(itau).mass());
        h2d_genTauMultiMatched_["pt_eta"]->Fill(genTaus->at(itau).pt(), genTaus->at(itau).eta());
        h2d_genTauMultiMatched_["pt_phi"]->Fill(genTaus->at(itau).pt(), genTaus->at(itau).phi());
        h2d_genTauMultiMatched_["pt_mass"]->Fill(genTaus->at(itau).pt(), genTaus->at(itau).mass());
        h2d_genTauMultiMatched_["mass_eta"]->Fill(genTaus->at(itau).mass(), genTaus->at(itau).eta());
        h2d_genTauMultiMatched_["mass_phi"]->Fill(genTaus->at(itau).mass(), genTaus->at(itau).phi());
      }
    }
  }

  // Define tau ID structres for cleaner code in the loop below (vs Jet, vs E, vs Mu)
  const std::array<std::tuple<const reco::TauDiscriminatorContainer*, bool, const char*>, 3> idInfo = {{
    {recoTauIDjet.isValid() ? recoTauIDjet.product() : nullptr, useIdvsJet, "Jet"},
    {recoTauIDe.isValid() ? recoTauIDe.product() : nullptr, useIdvsE, "E"},
    {recoTauIDmu.isValid() ? recoTauIDmu.product() : nullptr, useIdvsMu, "Mu"}
  }};

  // Lambda function to get ID value or return -1 if ID collection is not valid
  auto idValue = [&](const reco::TauDiscriminatorContainer* container, bool valid, const reco::PFTauRef& tauRef) {
    return valid ? (*container)[tauRef].rawValues[0] : -1.;
  };

  // Loop for fake rate 
  for (uint itau = 0; itau < recoTaus->size(); ++itau) {
    reco::PFTauRef tauRef(recoTaus, itau);

    if (applyIdCuts && !passIdCut(tauRef, cutIDs, recoTauIDjet.product(), recoTauIDe.product(), recoTauIDmu.product())) {
      continue; // skip if tau doesn't pass ID cuts
    }

    h_recoTau_["pt"]->Fill(recoTaus->at(itau).pt());
    h_recoTau_["eta"]->Fill(recoTaus->at(itau).eta());
    h_recoTau_["phi"]->Fill(recoTaus->at(itau).phi());
    h_recoTau_["mass"]->Fill(recoTaus->at(itau).mass());
    h2d_recoTau_["pt_eta"]->Fill(recoTaus->at(itau).pt(), recoTaus->at(itau).eta());
    h2d_recoTau_["pt_phi"]->Fill(recoTaus->at(itau).pt(), recoTaus->at(itau).phi());
    h2d_recoTau_["pt_mass"]->Fill(recoTaus->at(itau).pt(), recoTaus->at(itau).mass());
    h2d_recoTau_["mass_eta"]->Fill(recoTaus->at(itau).mass(), recoTaus->at(itau).eta());
    h2d_recoTau_["mass_phi"]->Fill(recoTaus->at(itau).mass(), recoTaus->at(itau).phi());

    for (auto const& [container, valid, label] : idInfo) {
      const double idRawValue = idValue(container, valid, tauRef);
      const std::string idName = std::string("id") + label;
      h_recoTau_[idName]->Fill(idRawValue);
      h2d_recoTau_[idName + "_pt"]->Fill(idRawValue, recoTaus->at(itau).pt());
      h2d_recoTau_[idName + "_eta"]->Fill(idRawValue, recoTaus->at(itau).eta());
      h2d_recoTau_[idName + "_phi"]->Fill(idRawValue, recoTaus->at(itau).phi());
      h2d_recoTau_[idName + "_mass"]->Fill(idRawValue, recoTaus->at(itau).mass());
    }

    // Count how many gen taus are matched to the reco tau
    int nGenMatchedToOneReco = 0;
    for (uint jtau = 0; jtau < genTaus->size(); ++jtau) {
      if (deltaR(genTaus->at(jtau), recoTaus->at(itau)) < matchingDeltaR) {
        nGenMatchedToOneReco++;
      }
    }

    // Fill histograms for reco taus matched to at least one gen tau
    if (nGenMatchedToOneReco > 0) {
      // Fill reco tau histograms for matched taus
      h_recoTauMatched_["pt"]->Fill(recoTaus->at(itau).pt());
      h_recoTauMatched_["eta"]->Fill(recoTaus->at(itau).eta());
      h_recoTauMatched_["phi"]->Fill(recoTaus->at(itau).phi());
      h_recoTauMatched_["mass"]->Fill(recoTaus->at(itau).mass());
      h2d_recoTauMatched_["pt_eta"]->Fill(recoTaus->at(itau).pt(), recoTaus->at(itau).eta());
      h2d_recoTauMatched_["pt_phi"]->Fill(recoTaus->at(itau).pt(), recoTaus->at(itau).phi());
      h2d_recoTauMatched_["pt_mass"]->Fill(recoTaus->at(itau).pt(), recoTaus->at(itau).mass());
      h2d_recoTauMatched_["mass_eta"]->Fill(recoTaus->at(itau).mass(), recoTaus->at(itau).eta());
      h2d_recoTauMatched_["mass_phi"]->Fill(recoTaus->at(itau).mass(), recoTaus->at(itau).phi());

      for (auto const& [container, valid, label] : idInfo) {
        const double idRawValue = idValue(container, valid, tauRef);
        const std::string idName = std::string("id") + label;
        h_recoTauMatched_[idName]->Fill(idRawValue);
        h2d_recoTauMatched_[idName + "_pt"]->Fill(idRawValue, recoTaus->at(itau).pt());
        h2d_recoTauMatched_[idName + "_eta"]->Fill(idRawValue, recoTaus->at(itau).eta());
        h2d_recoTauMatched_[idName + "_phi"]->Fill(idRawValue, recoTaus->at(itau).phi());
        h2d_recoTauMatched_[idName + "_mass"]->Fill(idRawValue, recoTaus->at(itau).mass());
      }

      if (nGenMatchedToOneReco > 1) {
        // Fill reco tau histograms for multi-matched taus
        h_recoTauMultiMatched_["pt"]->Fill(recoTaus->at(itau).pt());
        h_recoTauMultiMatched_["eta"]->Fill(recoTaus->at(itau).eta());
        h_recoTauMultiMatched_["phi"]->Fill(recoTaus->at(itau).phi());
        h_recoTauMultiMatched_["mass"]->Fill(recoTaus->at(itau).mass());
        h2d_recoTauMultiMatched_["pt_eta"]->Fill(recoTaus->at(itau).pt(), recoTaus->at(itau).eta());
        h2d_recoTauMultiMatched_["pt_phi"]->Fill(recoTaus->at(itau).pt(), recoTaus->at(itau).phi());
        h2d_recoTauMultiMatched_["pt_mass"]->Fill(recoTaus->at(itau).pt(), recoTaus->at(itau).mass());
        h2d_recoTauMultiMatched_["mass_eta"]->Fill(recoTaus->at(itau).mass(), recoTaus->at(itau).eta());
        h2d_recoTauMultiMatched_["mass_phi"]->Fill(recoTaus->at(itau).mass(), recoTaus->at(itau).phi());

        for (auto const& [container, valid, label] : idInfo) {
          const double idRawValue = idValue(container, valid, tauRef);
          const std::string idName = std::string("id") + label;
          h_recoTauMultiMatched_[idName]->Fill(idRawValue);
          h2d_recoTauMultiMatched_[idName + "_pt"]->Fill(idRawValue, recoTaus->at(itau).pt());
          h2d_recoTauMultiMatched_[idName + "_eta"]->Fill(idRawValue, recoTaus->at(itau).eta());
          h2d_recoTauMultiMatched_[idName + "_phi"]->Fill(idRawValue, recoTaus->at(itau).phi());
          h2d_recoTauMultiMatched_[idName + "_mass"]->Fill(idRawValue, recoTaus->at(itau).mass());
        }
        }
    }
  }
}

//------------------------------------------------------------------------------
// fill description
//------------------------------------------------------------------------------
void TauValidationRECO::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  // Default tau validation HLT
  desc.add<edm::InputTag>("genTauCollection", edm::InputTag("tauGenJets"));
  desc.add<edm::InputTag>("recoTauCollection", edm::InputTag("hltHpsPFTauProducer"));
  desc.add<edm::InputTag>("recoTauIDCollectionVsJet", edm::InputTag("hltHpsPFTauDeepTauProducer:VSjet"));
  desc.add<edm::InputTag>("recoTauIDCollectionVsE", edm::InputTag("hltHpsPFTauDeepTauProducer:VSe"));
  desc.add<edm::InputTag>("recoTauIDCollectionVsMu", edm::InputTag("hltHpsPFTauDeepTauProducer:VSmu"));
  desc.add<std::vector<double>>("cutIDs", std::vector<double>{0.0, 0.0, 0.0});
  desc.add<double>("minDeltaR", 0.3);
  desc.add<std::string>("outFolder", "HLT/Tau/TauValidation");
  desc.addUntracked<bool>("isHLT", true);
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(TauValidationRECO);