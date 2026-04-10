// Analyzer for validation histograms for tau objects at HLT/RECO
// E. Vernazza Apr. 10, 2026

#include "Validation/RecoTau/interface/TauValidationRECO.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include <vector>
#include <format>

using namespace edm;
using namespace reco;
using namespace std;

TauValidationRECO::TauValidationRECO(const edm::ParameterSet& iConfig) {
  recoTauToken_ = consumes<reco::PFTauCollection>(iConfig.getParameter<edm::InputTag>("recoTauCollection"));
  genTauToken_ = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genTauCollection"));
  isHLT = iConfig.getUntrackedParameter<bool>("isHLT");
}

void TauValidationRECO::bookHistograms(DQMStore::IBooker& ibooker, edm::Run const& iRun, edm::EventSetup const&) {

  // ---------------------------- Book Summary Histograms -------------------------------
  if (isHLT) {
    ibooker.setCurrentFolder("HLT/Tau/TauValidation");
  } else {
    ibooker.setCurrentFolder("Tau/TauValidation");
  }
  
  for (auto& hVar : histoVars) {
    auto [nBins, hMin, hMax] = hVar.second;
    h_recoTau_[hVar.first] = ibooker.book1D("recoTau_" + hVar.first, "#tau^{reco};" + hVar.first, nBins, hMin, hMax);
    h_genTau_[hVar.first] = ibooker.book1D("genTau_" + hVar.first, "#tau^{gen};" + hVar.first, nBins, hMin, hMax);
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
  std::cout << "Number of reco taus: " << recoTaus->size() << std::endl; // [DEBUG]

  edm::Handle<reco::GenJetCollection> genTaus;
  mEvent.getByToken(genTauToken_, genTaus);
  if (!genTaus.isValid()) {
    edm::LogPrint("TauValidationRECO") << " Gen Tau collection not found while running TauValidationRECO.cc ";
    return;
  }
  std::cout << "Number of gen taus: " << genTaus->size() << std::endl; // [DEBUG]

  // Loop for efficiency 
  for (uint itau = 0; itau < genTaus->size(); ++itau) {
    h_genTau_["pt"]->Fill(genTaus->at(itau).pt());
    h_genTau_["eta"]->Fill(genTaus->at(itau).eta());
    h_genTau_["phi"]->Fill(genTaus->at(itau).phi());
    h_genTau_["mass"]->Fill(genTaus->at(itau).mass());
  }

  // Loop for fake rate
  for (uint itau = 0; itau < recoTaus->size(); ++itau) {
    h_recoTau_["pt"]->Fill(recoTaus->at(itau).pt());
    h_recoTau_["eta"]->Fill(recoTaus->at(itau).eta());
    h_recoTau_["phi"]->Fill(recoTaus->at(itau).phi());
    h_recoTau_["mass"]->Fill(recoTaus->at(itau).mass());
  }

}

//------------------------------------------------------------------------------
// fill description
//------------------------------------------------------------------------------
void TauValidationRECO::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  // Default tau validation HLT
  desc.add<edm::InputTag>("recoTauCollection", edm::InputTag("hltHpsPFTauProducer"));
  desc.add<edm::InputTag>("genTauCollection", edm::InputTag("tauGenJets"));
  desc.addUntracked<bool>("isHLT", true);
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(TauValidationRECO);