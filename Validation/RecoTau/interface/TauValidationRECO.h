#ifndef ValidationRecoTau_TauValidationRECO_h
#define ValidationRecoTau_TauValidationRECO_h

// Analyzer for validation histograms for tau objects at HLT/RECO
// E. Vernazza Apr. 10, 2026

#include <cmath>
#include <string>
#include <string_view>

#include "DQMServices/Core/interface/DQMStore.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include <DQMServices/Core/interface/DQMEDAnalyzer.h>
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
#include "DataFormats/Math/interface/deltaR.h"

class TauValidationRECO : public DQMEDAnalyzer {

public:
  TauValidationRECO(const edm::ParameterSet &);
  ~TauValidationRECO() override;

  void analyze(const edm::Event &, const edm::EventSetup &) override;
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

  bool isSelectedDecayMode(const reco::GenJet &genTau, const std::vector<std::string> &decayModes) const;

private:

  edm::EDGetTokenT<reco::PFTauCollection> recoTauToken_;
  edm::EDGetTokenT<reco::GenJetCollection> genTauToken_;

  const std::unordered_map<std::string, std::tuple<unsigned, float, float>> histoVars = {
    {"pt", std::make_tuple(200, 0., 1000.)},
    {"eta", std::make_tuple(60, -4.0, 4.0)},
    {"phi", std::make_tuple(50, -3.5, 3.5)},
    {"mass", std::make_tuple(200, 0, 10.)},
  };

  using UMap = std::unordered_map<std::string, MonitorElement*>;
  UMap h_recoTau_;
  UMap h_recoTauMatched_;
  UMap h_recoTauMultiMatched_;
  UMap h_genTau_;
  UMap h_genTauMatched_;
  UMap h_genTauMultiMatched_;

  bool isHLT;
  float matchingDeltaR;
  std::vector<std::string> decayModes;

};

#endif