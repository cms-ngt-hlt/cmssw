// Analyzer for validation histograms for tau objects at HLT/RECO
// E. Vernazza Apr. 10, 2026

#include "Validation/RecoTau/interface/TauValidation.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include <vector>
#include <format>

using namespace edm;
using namespace reco;
using namespace std;

TauValidation::TauValidation(const edm::ParameterSet& iConfig) {
  recoTauToken_ = consumes<reco::PFTauCollection>(iConfig.getParameter<edm::InputTag>("recoTauCollection"));
  genTauToken_ = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genTauCollection"));
  decayModes = iConfig.getParameter<std::vector<std::string>>("decayModes");
  matchingDeltaR = iConfig.getParameter<double>("minDeltaR");
  outFolder_ = iConfig.getParameter<std::string>("outFolder");
  isHLT = iConfig.getUntrackedParameter<bool>("isHLT");
}

bool TauValidation::isSelectedDecayMode(const reco::GenJet& genTau, const std::vector<std::string>& decayModes) const {
  std::string tauGenJetDecayMode = JetMCTagUtils::genTauDecayMode(genTau);
  for (std::vector<std::string>::const_iterator i_dm = decayModes.begin(); i_dm != decayModes.end(); ++i_dm) {
    if (tauGenJetDecayMode == (*i_dm)) {
      // edm::LogPrint("TauValidation") << "Selected gen tau decay mode: " << tauGenJetDecayMode;
      return true;
    }
  }
  return false;
}


void TauValidation::bookHistograms(DQMStore::IBooker& ibooker, edm::Run const& iRun, edm::EventSetup const&) {

  // ---------------------------- Book Summary Histograms -------------------------------
  ibooker.setCurrentFolder(outFolder_);
  
  for (auto& hVar : histoVars) {
    auto [nBins, hMin, hMax] = hVar.second;
    h_recoTau_[hVar.first] = ibooker.book1D("recoTau_" + hVar.first, ";#tau^{reco};" + hVar.first, nBins, hMin, hMax);
    h_recoTauMatched_[hVar.first] = ibooker.book1D("recoTauMatched_" + hVar.first, ";#tau^{reco} (Matched);" + hVar.first, nBins, hMin, hMax);
    h_recoTauMultiMatched_[hVar.first] = ibooker.book1D("recoTauMultiMatched_" + hVar.first, ";#tau^{reco} (Multi-Matched);" + hVar.first, nBins, hMin, hMax);
    h_genTau_[hVar.first] = ibooker.book1D("genTau_" + hVar.first, "#tau^{gen};" + hVar.first, nBins, hMin, hMax);
    h_genTauMatched_[hVar.first] = ibooker.book1D("genTauMatched_" + hVar.first, ";#tau^{gen} (Matched);" + hVar.first, nBins, hMin, hMax);
    h_genTauMultiMatched_[hVar.first] = ibooker.book1D("genTauMultiMatched_" + hVar.first, ";#tau^{gen} (Multi-Matched);" + hVar.first, nBins, hMin, hMax);
  }

}

//------------------------------------------------------------------------------
// ~TauValidation
//------------------------------------------------------------------------------
TauValidation::~TauValidation() {}

//------------------------------------------------------------------------------
// analyze
//------------------------------------------------------------------------------
void TauValidation::analyze(const edm::Event& mEvent, const edm::EventSetup& mSetup) {

  edm::Handle<reco::PFTauCollection> recoTaus;
  mEvent.getByToken(recoTauToken_, recoTaus);
  if (!recoTaus.isValid()) {
    edm::LogPrint("TauValidation") << " Reco Tau collection not found while running TauValidation.cc ";
    return;
  }
  // std::cout << "Number of reco taus: " << recoTaus->size() << std::endl; // [DEBUG]

  edm::Handle<reco::GenJetCollection> genTaus;
  mEvent.getByToken(genTauToken_, genTaus);
  if (!genTaus.isValid()) {
    edm::LogPrint("TauValidation") << " Gen Tau collection not found while running TauValidation.cc ";
    return;
  }
  // std::cout << "Number of gen taus: " << genTaus->size() << std::endl; // [DEBUG]

  // Loop for efficiency 
  for (uint itau = 0; itau < genTaus->size(); ++itau) {

    if (!isSelectedDecayMode(genTaus->at(itau), decayModes)){
      continue;
    }
    
    h_genTau_["pt"]->Fill(genTaus->at(itau).pt());
    h_genTau_["eta"]->Fill(genTaus->at(itau).eta());
    h_genTau_["phi"]->Fill(genTaus->at(itau).phi());
    h_genTau_["mass"]->Fill(genTaus->at(itau).mass());
    
    // Count how many reco taus are matched to the gen tau
    int nRecoMatchedToOneGen = 0;
    for (uint jtau = 0; jtau < recoTaus->size(); ++jtau) {
      if (deltaR(genTaus->at(itau), recoTaus->at(jtau)) < matchingDeltaR) {
        nRecoMatchedToOneGen++;
      }
    }

    // Fill histograms for gen taus matched to at least one reco tau
    if (nRecoMatchedToOneGen > 0) {
      h_genTauMatched_["pt"]->Fill(genTaus->at(itau).pt());
      h_genTauMatched_["eta"]->Fill(genTaus->at(itau).eta());
      h_genTauMatched_["phi"]->Fill(genTaus->at(itau).phi());
      h_genTauMatched_["mass"]->Fill(genTaus->at(itau).mass());
      if (nRecoMatchedToOneGen > 1) {
        h_genTauMultiMatched_["pt"]->Fill(genTaus->at(itau).pt());
        h_genTauMultiMatched_["eta"]->Fill(genTaus->at(itau).eta());
        h_genTauMultiMatched_["phi"]->Fill(genTaus->at(itau).phi());
        h_genTauMultiMatched_["mass"]->Fill(genTaus->at(itau).mass());
      }
    }
  }

  // Loop for fake rate
  for (uint itau = 0; itau < recoTaus->size(); ++itau) {
    h_recoTau_["pt"]->Fill(recoTaus->at(itau).pt());
    h_recoTau_["eta"]->Fill(recoTaus->at(itau).eta());
    h_recoTau_["phi"]->Fill(recoTaus->at(itau).phi());
    h_recoTau_["mass"]->Fill(recoTaus->at(itau).mass());

    // Count how many gen taus are matched to the reco tau
    int nGenMatchedToOneReco = 0;
    for (uint jtau = 0; jtau < genTaus->size(); ++jtau) {
      if (deltaR(genTaus->at(jtau), recoTaus->at(itau)) < matchingDeltaR) {
        nGenMatchedToOneReco++;
      }

    // Fill histograms for reco taus matched to at least one gen tau
      if (nGenMatchedToOneReco > 0) {
        h_recoTauMatched_["pt"]->Fill(recoTaus->at(itau).pt());
        h_recoTauMatched_["eta"]->Fill(recoTaus->at(itau).eta());
        h_recoTauMatched_["phi"]->Fill(recoTaus->at(itau).phi());
        h_recoTauMatched_["mass"]->Fill(recoTaus->at(itau).mass());
        if (nGenMatchedToOneReco > 1) {
          h_recoTauMultiMatched_["pt"]->Fill(recoTaus->at(itau).pt());
          h_recoTauMultiMatched_["eta"]->Fill(recoTaus->at(itau).eta());
          h_recoTauMultiMatched_["phi"]->Fill(recoTaus->at(itau).phi());
          h_recoTauMultiMatched_["mass"]->Fill(recoTaus->at(itau).mass());
        }
      }
    }
  }

}

//------------------------------------------------------------------------------
// fill description
//------------------------------------------------------------------------------
void TauValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  // Default tau validation HLT
  desc.add<edm::InputTag>("recoTauCollection", edm::InputTag("hltHpsPFTauProducer"));
  desc.add<edm::InputTag>("genTauCollection", edm::InputTag("tauGenJets"));
  desc.add<std::vector<std::string>>("decayModes", {"oneProng0Pi0", "oneProng1Pi0", "oneProng2Pi0", "oneProngOther", "threeProng0Pi0", "threeProng1Pi0", "threeProngOther", "rare"});
  desc.add<double>("minDeltaR", 0.3);
  desc.add<std::string>("outFolder", "HLT/Tau/TauValidation");
  desc.addUntracked<bool>("isHLT", true);
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(TauValidation);