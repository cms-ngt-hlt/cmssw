#ifndef TauValidationMiniAOD_h
#define TauValidationMiniAOD_h

// -*- C++ -*-
//
// Package:    TauValidationMiniAOD
// Class:      TauValidationMiniAOD
//
/* *\class TauValidationMiniAOD TauValidationMiniAOD.cc

 Description: EDAnalyzer to validate tau collection in miniAOD
 Implementation:

*/
// Original Author: Aniello Spiezia On August 13, 2019
// Updated April, 2020 by Ece Asilar and Gage DeZoort
// Updated July, 2023 by Gourab Saha

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/Math/interface/deltaR.h"

// Include DQM core
#include <DQMServices/Core/interface/DQMStore.h>
#include <DQMServices/Core/interface/MonitorElement.h>
#include <DQMServices/Core/interface/DQMEDAnalyzer.h>

struct histoInfo {
  int nbins;
  double min;
  double max;
  histoInfo(int n, double m, double M) {
    nbins = n;
    min = m;
    max = M;
  }
  histoInfo(const edm::ParameterSet &config) {
    nbins = config.getParameter<int>("nbins");
    min = config.getParameter<double>("min");
    max = config.getParameter<double>("max");
  }
};

// class declaration
class TauValidationMiniAOD : public DQMEDAnalyzer {
public:
  explicit TauValidationMiniAOD(const edm::ParameterSet &);
  ~TauValidationMiniAOD() override;

  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  void analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) override;
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2#Decay_Mode_Reconstruction
  int findDecayMode(int Nc, int Np) { return 5 * (Nc - 1) + Np; };

private:
  edm::EDGetTokenT<std::vector<pat::Tau>> patTauToken_;
  edm::EDGetTokenT<reco::PFTauCollection> pfTauToken_;
  edm::EDGetTokenT<edm::View<reco::Candidate>> genRefToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> pvToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>> prunedGenToken_;

  const std::unordered_map<std::string, std::tuple<unsigned, float, float>> histoVars = {
    {"pt", std::make_tuple(200, 0., 1000.)},
    {"eta", std::make_tuple(60, -4.0, 4.0)},
    {"phi", std::make_tuple(50, -3.5, 3.5)},
    {"mass", std::make_tuple(200, 0, 10.)},
    {"pu", std::make_tuple(100, 0., 100.)},
  };

  const std::vector<std::pair<std::vector<int>, std::string>> dm_list = {
    {{0}, "dm0"},
    {{1, 2}, "dm1p2"},
    {{5}, "dm5"},
    {{6}, "dm6"},
    {{10}, "dm10"},
    {{11}, "dm11"}
  };

  using UMap = std::unordered_map<std::string, MonitorElement*>;
  UMap h_tausMatchedToRef_;
  UMap h_tausMatchedToRef_TightvsJet;
  UMap h_tausMatchedToRef_TightvsEle;
  UMap h_tausMatchedToRef_TightvsMuo;
  UMap h_tausMatchedToRef_MediumvsJet;
  UMap h_tausMatchedToRef_MediumvsEle;
  UMap h_tausMatchedToRef_MediumvsMuo;
  UMap h_tausMatchedToRef_LoosevsJet;
  UMap h_tausMatchedToRef_LoosevsEle;
  UMap h_tausMatchedToRef_LoosevsMuo;

  MonitorElement* DeepTau2018v2p5VSele;
  MonitorElement* DeepTau2018v2p5VSjet;
  MonitorElement* DeepTau2018v2p5VSmuo;
  MonitorElement* decayModeFinding;
  MonitorElement* decayModeTauReco;
  MonitorElement* decayModeTauGen;
  MonitorElement* dmMigration;
  MonitorElement* nTau_vs_dm;

  std::vector<MonitorElement*> h_pTOverProng_dm_;
  std::vector<MonitorElement*> h_TauMass_dm_;

  std::map<std::string, MonitorElement *> summaryMap;

  edm::ParameterSet histoSettings_;
  std::string extensionName_;
  std::vector<edm::ParameterSet> discriminators_;

  // Parameters
  bool isMini;
  bool isReco;
  bool isHLT;
  std::string TauType;
};

#endif
