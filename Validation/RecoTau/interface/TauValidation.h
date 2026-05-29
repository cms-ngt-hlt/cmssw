#ifndef ValidationRecoTau_TauValidation_h
#define ValidationRecoTau_TauValidation_h

// Analyzer for validation histograms for tau objects at HLT/RECO
// E. Vernazza Apr. 10, 2026

#include <string>
// #include <string_view>
#include <vector>
#include <tuple>
#include <unordered_map>
// #include <array>
// #include <sstream>
// #include <iomanip>

#include "DQMServices/Core/interface/DQMStore.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include <DQMServices/Core/interface/DQMEDAnalyzer.h>
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
#include "DataFormats/TauReco/interface/TauDiscriminatorContainer.h"
#include "DataFormats/Math/interface/deltaR.h"

class TauValidation : public DQMEDAnalyzer {

public:
  TauValidation(const edm::ParameterSet &);
  ~TauValidation() override;

  void analyze(const edm::Event &, const edm::EventSetup &) override;
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  std::string convertId(double cut);
  bool passIdCut(const std::vector<double> idValuesForTau, const std::vector<std::vector<bool>> wpValuesForTau,
                const std::vector<double>& validCutIDs_raw, const std::vector<int>& validCutIDs_wp,
                bool use_raw, bool use_wp);

private:

  edm::EDGetTokenT<reco::GenJetCollection> genTauToken_;
  edm::EDGetTokenT<reco::PFTauCollection> recoTauToken_;
  edm::EDGetTokenT<pat::TauCollection> patTauToken_;
  std::vector<edm::EDGetTokenT<reco::TauDiscriminatorContainer>> recoTauIDTokens_;
  std::vector<std::string> recoTauIDLabels_;
  edm::InputTag recoTauCollection;

  const std::unordered_map<std::string, std::tuple<unsigned, float, float>> histoVars = {
    {"pt", std::make_tuple(200, 0., 1000.)},
    {"eta", std::make_tuple(60, -4.0, 4.0)},
    {"phi", std::make_tuple(50, -3.5, 3.5)},
    {"mass", std::make_tuple(200, 0, 10.)},
  };

  const std::unordered_map<std::string, std::tuple<unsigned, float, float, unsigned, float, float>> histoVars2D = {
    {"pt_eta", std::make_tuple(200, 0., 1000., 60, -4.0, 4.0)},
    {"pt_phi", std::make_tuple(200, 0., 1000., 50, -3.5, 3.5)},
    {"pt_mass", std::make_tuple(200, 0., 1000., 200, 0., 10.)},
    {"mass_eta", std::make_tuple(200, 0., 10., 60, -4.0, 4.0)},
    {"mass_phi", std::make_tuple(200, 0., 10., 50, -3.5, 3.5)},
  };


  using UMap = std::unordered_map<std::string, MonitorElement*>;
  UMap h_recoTau_;
  UMap h_recoTauMatched_;
  UMap h_recoTauMultiMatched_;
  UMap h_genTau_;
  UMap h_genTauMatched_;
  UMap h_genTauMultiMatched_;
  UMap h2d_recoTau_;
  UMap h2d_recoTauMatched_;
  UMap h2d_recoTauMultiMatched_;
  UMap h2d_genTau_;
  UMap h2d_genTauMatched_;
  UMap h2d_genTauMultiMatched_;
  UMap h2d_responsePt_;
  UMap h2d_responseMass_;

  std::vector<int> cutIDs_wp;  // Working-point indices (WP mode)
  bool use_wp;
  std::vector<double> cutIDs_raw;    // Raw discriminator value cuts (raw mode)
  bool use_raw;

  bool isPatTaus;
  float matchingDeltaR;
  std::string outFolder;

};

#endif