// Analyzer for validation histograms for tau objects at HLT/RECO
// E. Vernazza Apr. 10, 2026

#include "Validation/RecoTau/interface/TauValidationRECO.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include <cstddef>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <array>

using namespace edm;
using namespace reco;
using namespace std;

std::string TauValidationRECO::convertId(double cut) {
  if (cut == 0.0) return "0p0";
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(2) << cut;
  std::string result = oss.str();
  for (char& c : result) {
    if (c == '.') c = 'p';
  }
  return result;
}

bool TauValidationRECO::passIdCut(const reco::PFTauRef& tauRef,
                                  const std::vector<const reco::TauDiscriminatorContainer*>& validRecoTauIDs,
                                  const std::vector<double>& validCutIDs_raw, const std::vector<int>& validCutIDs_wp,
                                  bool use_raw, bool use_wp) {

  if (use_raw) {
    for (size_t i = 0; i < validRecoTauIDs.size(); ++i) {
      const auto& disc = (*validRecoTauIDs[i])[tauRef];
      if (validCutIDs_raw[i] != 0.0) {
        if (disc.rawValues.empty() || disc.rawValues[0] <= validCutIDs_raw[i]) {
          return false;
        }
      }
    }
    return true; // Passes all raw cuts
  }

  if (use_wp) {
    for (size_t i = 0; i < validRecoTauIDs.size(); ++i) {
      const auto& disc = (*validRecoTauIDs[i])[tauRef];
      if (validCutIDs_wp[i] >= 0) {
        if (disc.workingPoints.size() <= static_cast<size_t>(validCutIDs_wp[i]) || !disc.workingPoints[validCutIDs_wp[i]]) {
          return false;
        }
      }
    }
    return true; // Passes all WP cuts
  }

  return true; // If no cuts specified, pass by default

}

TauValidationRECO::TauValidationRECO(const edm::ParameterSet& iConfig) {
  genTauToken_ = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genTauCollection"));
  recoTauToken_ = consumes<reco::PFTauCollection>(iConfig.getParameter<edm::InputTag>("recoTauCollection"));
  matchingDeltaR = iConfig.getParameter<double>("minDeltaR");
  outFolder = iConfig.getParameter<std::string>("outFolder");
  isHLT = iConfig.getUntrackedParameter<bool>("isHLT");

  std::vector<edm::InputTag> idTags = iConfig.getParameter<std::vector<edm::InputTag>>("recoTauIDCollections");
  for (const auto& tag : idTags) {
    // Extract label from the instance part of the InputTag (e.g., "VSjet" from "module:VSjet")
    recoTauIDTokens_.push_back(consumes<reco::TauDiscriminatorContainer>(tag));
    recoTauIDLabels_.push_back(tag.instance());
  }

  cutIDs_wp = iConfig.getParameter<std::vector<int>>("cutIDs_wp");
  cutIDs_raw = iConfig.getParameter<std::vector<double>>("cutIDs_raw");

  use_wp = std::any_of(cutIDs_wp.begin(), cutIDs_wp.end(), [](int x){ return x >= 0; });
  use_raw = std::any_of(cutIDs_raw.begin(), cutIDs_raw.end(), [](double x){ return x != 0.0; });

  if (use_wp && use_raw) {
      throw cms::Exception("Configuration") << "Specify either cutIDs_wp OR cutIDs_raw, not both";
  }

  if (!recoTauIDLabels_.empty()) {
    if ((use_wp) && (cutIDs_wp.size() != recoTauIDLabels_.size())) {
      cutIDs_wp.resize(recoTauIDLabels_.size(), -1);
      edm::LogPrint("TauValidationRECO") << "Warning: cutIDs_wp size (" << cutIDs_wp.size() 
        << ") adjusted to match idLabels size (" << recoTauIDLabels_.size() << ")";
    }
    if ((use_raw) && (cutIDs_raw.size() != recoTauIDLabels_.size())) {
      cutIDs_raw.resize(recoTauIDLabels_.size(), 0.0);
      edm::LogPrint("TauValidationRECO") << "Warning: cutIDs_raw size (" << cutIDs_raw.size() 
        << ") adjusted to match idLabels size (" << recoTauIDLabels_.size() << ")";
    }
  }

}

void TauValidationRECO::bookHistograms(DQMStore::IBooker& ibooker, edm::Run const& iRun, edm::EventSetup const&) {

  // ---------------------------- Book Summary Histograms -------------------------------

  if (use_wp) {
    outFolder += "/CutWP";
    for (size_t i = 0; i < recoTauIDLabels_.size(); ++i) {
      if (cutIDs_wp[i] >= 0) {
        outFolder += "_" + recoTauIDLabels_[i] + std::to_string(cutIDs_wp[i]);
      }
    }
  }
  if (use_raw) {
    outFolder += "/CutID";
    for (size_t i = 0; i < recoTauIDLabels_.size(); ++i) {
      if (cutIDs_raw[i] > 0) {
        outFolder += "_" + recoTauIDLabels_[i] + convertId(cutIDs_raw[i]);
      }
    }
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

    // Book 2D histograms for reco tau ID discriminators vs kinematics (dynamic based on idLabels)
    for (const auto& label : recoTauIDLabels_) {
      std::string idName = "id" + label + "_" + hVar.first;
      h2d_recoTau_[idName] = ibooker.book2D("recoTau_" + idName, ";#tau^{reco}; ID" + label + ";" + hVar.first, 50, 0., 1., nBins, hMin, hMax);
      h2d_recoTauMatched_[idName] = ibooker.book2D("recoTauMatched_" + idName, ";#tau^{reco} (Matched); ID" + label + ";" + hVar.first, 50, 0., 1., nBins, hMin, hMax);
      h2d_recoTauMultiMatched_[idName] = ibooker.book2D("recoTauMultiMatched_" + idName, ";#tau^{reco} (Multi-Matched); ID" + label + ";" + hVar.first, 50, 0., 1., nBins, hMin, hMax);
    }
  }

  // Book 1D histograms for reco tau ID discriminators (dynamic based on idLabels)
  for (const auto& label : recoTauIDLabels_) {
    std::string idName = "id" + label;
    h_recoTau_[idName] = ibooker.book1D("recoTau_" + idName, ";#tau^{reco};" + idName, 50, 0., 1.);
    h_recoTauMatched_[idName] = ibooker.book1D("recoTauMatched_" + idName, ";#tau^{reco} (Matched);" + idName, 50, 0., 1.);
    h_recoTauMultiMatched_[idName] = ibooker.book1D("recoTauMultiMatched_" + idName, ";#tau^{reco} (Multi-Matched);" + idName, 50, 0., 1.);
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

  std::vector<const reco::TauDiscriminatorContainer*> validRecoTauIDs;
  std::vector<std::string> validRecoTauIDLabels;
  std::vector<double> validCutIDs_raw;
  std::vector<int> validCutIDs_wp;
  for (size_t i = 0; i < recoTauIDTokens_.size(); ++i) {
    edm::Handle<reco::TauDiscriminatorContainer> recoTauID;
    mEvent.getByToken(recoTauIDTokens_[i], recoTauID);
    if (!recoTauID.isValid()) {
      edm::LogPrint("TauValidationRECO") << "Reco Tau Identifier " << recoTauIDLabels_[i] << " collection not found while running TauValidationRECO.cc ";
      continue;
    }
    validRecoTauIDs.push_back(recoTauID.product());
    validRecoTauIDLabels.push_back(recoTauIDLabels_[i]);
    validCutIDs_raw.push_back(cutIDs_raw[i]);
    validCutIDs_wp.push_back(cutIDs_wp[i]);
  }

  bool plotId = validRecoTauIDs.size() > 0;
  bool applyIdCuts = plotId && (use_wp || use_raw);

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
      if (applyIdCuts && !passIdCut(tauRef, validRecoTauIDs, validCutIDs_raw, validCutIDs_wp, use_raw, use_wp)) {
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

  // Loop for fake rate 
  for (uint itau = 0; itau < recoTaus->size(); ++itau) {
    reco::PFTauRef tauRef(recoTaus, itau);

    if (applyIdCuts && !passIdCut(tauRef, validRecoTauIDs, validCutIDs_raw, validCutIDs_wp, use_raw, use_wp)) {
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

    if (plotId) {
      for (size_t i = 0; i < validRecoTauIDLabels.size(); ++i) {
        const auto& disc = (*validRecoTauIDs[i])[tauRef];
        const double idRawValue = disc.rawValues.empty() ? -1. : disc.rawValues[0];
        const std::string idName = "id" + validRecoTauIDLabels[i];
        h_recoTau_[idName]->Fill(idRawValue);
        h2d_recoTau_[idName + "_pt"]->Fill(idRawValue, recoTaus->at(itau).pt());
        h2d_recoTau_[idName + "_eta"]->Fill(idRawValue, recoTaus->at(itau).eta());
        h2d_recoTau_[idName + "_phi"]->Fill(idRawValue, recoTaus->at(itau).phi());
        h2d_recoTau_[idName + "_mass"]->Fill(idRawValue, recoTaus->at(itau).mass());
      }
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

      if (plotId) {
        for (size_t i = 0; i < validRecoTauIDLabels.size(); ++i) {
          const auto& disc = (*validRecoTauIDs[i])[tauRef];
          const double idRawValue = disc.rawValues.empty() ? -1. : disc.rawValues[0];
          const std::string idName = "id" + validRecoTauIDLabels[i];
          h_recoTauMatched_[idName]->Fill(idRawValue);
          h2d_recoTauMatched_[idName + "_pt"]->Fill(idRawValue, recoTaus->at(itau).pt());
          h2d_recoTauMatched_[idName + "_eta"]->Fill(idRawValue, recoTaus->at(itau).eta());
          h2d_recoTauMatched_[idName + "_phi"]->Fill(idRawValue, recoTaus->at(itau).phi());
          h2d_recoTauMatched_[idName + "_mass"]->Fill(idRawValue, recoTaus->at(itau).mass());
        }
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

        if (plotId) {
          for (size_t i = 0; i < validRecoTauIDLabels.size(); ++i) {
            const auto& disc = (*validRecoTauIDs[i])[tauRef];
            const double idRawValue = disc.rawValues.empty() ? -1. : disc.rawValues[0];
            const std::string idName = "id" + validRecoTauIDLabels[i];
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
}

//------------------------------------------------------------------------------
// fill description
//------------------------------------------------------------------------------
void TauValidationRECO::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  // Default tau validation HLT
  desc.add<edm::InputTag>("genTauCollection", edm::InputTag("tauGenJets"));
  desc.add<edm::InputTag>("recoTauCollection", edm::InputTag("hltHpsPFTauProducer"));
  desc.add<std::vector<edm::InputTag>>("recoTauIDCollections", 
    std::vector<edm::InputTag>{
      edm::InputTag("hltHpsPFTauDeepTauProducer:VSjet"),
      edm::InputTag("hltHpsPFTauDeepTauProducer:VSe"),
      edm::InputTag("hltHpsPFTauDeepTauProducer:VSmu")
    });
  desc.add<std::vector<double>>("cutIDs_raw", std::vector<double>{0.0, 0.0, 0.0});
  desc.add<std::vector<int>>("cutIDs_wp", std::vector<int>{-1, -1, -1});
  desc.add<double>("minDeltaR", 0.3);
  desc.add<std::string>("outFolder", "HLT/Tau/TauValidation");
  desc.addUntracked<bool>("isHLT", true);
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(TauValidationRECO);