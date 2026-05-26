#ifndef ValidationRecoTau_TauValidationRECO_h
#define ValidationRecoTau_TauValidationRECO_h

// Analyzer for validation histograms for tau objects at HLT/RECO
// E. Vernazza Apr. 10, 2026

#include <cmath>
#include <string>
#include <string_view>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <array>
#include <sstream>
#include <iomanip>

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
#include "DataFormats/TauReco/interface/TauDiscriminatorContainer.h"
#include "DataFormats/Math/interface/deltaR.h"

class TauValidationRECO : public DQMEDAnalyzer {

public:
  TauValidationRECO(const edm::ParameterSet &);
  ~TauValidationRECO() override;

  void analyze(const edm::Event &, const edm::EventSetup &) override;
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  std::string convertId(double cut);
  bool passIdCut(const reco::PFTauRef& tauRef, const std::vector<double> &cutIDs,
                 const reco::TauDiscriminatorContainer* recoTauIDjet,
                 const reco::TauDiscriminatorContainer* recoTauIDe,
                 const reco::TauDiscriminatorContainer* recoTauIDmu);

private:

  edm::EDGetTokenT<reco::GenJetCollection> genTauToken_;
  edm::EDGetTokenT<reco::PFTauCollection> recoTauToken_;
  edm::EDGetTokenT<reco::TauDiscriminatorContainer> recoTauIDjetToken_;
  edm::EDGetTokenT<reco::TauDiscriminatorContainer> recoTauIDeToken_;
  edm::EDGetTokenT<reco::TauDiscriminatorContainer> recoTauIDmuToken_;

  const std::unordered_map<std::string, std::tuple<unsigned, float, float>> histoVars = {
    {"pt", std::make_tuple(200, 0., 1000.)},
    {"eta", std::make_tuple(60, -4.0, 4.0)},
    {"phi", std::make_tuple(50, -3.5, 3.5)},
    {"mass", std::make_tuple(200, 0, 10.)},
  };

  const std::unordered_map<std::string, std::tuple<unsigned, float, float>> histoVarsReco = {
    {"idJet", std::make_tuple(50, 0., 1.)},
    {"idE", std::make_tuple(50, 0., 1.)},
    {"idMu", std::make_tuple(50, 0., 1.)},
  };

  const std::unordered_map<std::string, std::tuple<unsigned, float, float, unsigned, float, float>> histoVars2D = {
    {"pt_eta", std::make_tuple(200, 0., 1000., 60, -4.0, 4.0)},
    {"pt_phi", std::make_tuple(200, 0., 1000., 50, -3.5, 3.5)},
    {"pt_mass", std::make_tuple(200, 0., 1000., 200, 0., 10.)},
    {"mass_eta", std::make_tuple(200, 0., 10., 60, -4.0, 4.0)},
    {"mass_phi", std::make_tuple(200, 0., 10., 50, -3.5, 3.5)},
  };

  const std::unordered_map<std::string, std::tuple<unsigned, float, float, unsigned, float, float>> histoVarsReco2D = {
    {"idJet_pt", std::make_tuple(50, 0., 1., 200, 0., 1000.)},
    {"idE_pt", std::make_tuple(50, 0., 1., 200, 0., 1000.)},
    {"idMu_pt", std::make_tuple(50, 0., 1., 200, 0., 1000.)},
    {"idJet_mass", std::make_tuple(50, 0., 1., 200, 0., 10.)},
    {"idE_mass", std::make_tuple(50, 0., 1., 200, 0., 10.)},
    {"idMu_mass", std::make_tuple(50, 0., 1., 200, 0., 10.)},
    {"idJet_eta", std::make_tuple(50, 0., 1., 60, -4.0, 4.0)},
    {"idE_eta", std::make_tuple(50, 0., 1., 60, -4.0, 4.0)},
    {"idMu_eta", std::make_tuple(50, 0., 1., 60, -4.0, 4.0)},
    {"idJet_phi", std::make_tuple(50, 0., 1., 50, -3.5, 3.5)},
    {"idE_phi", std::make_tuple(50, 0., 1., 50, -3.5, 3.5)},
    {"idMu_phi", std::make_tuple(50, 0., 1., 50, -3.5, 3.5)},
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

  std::vector<double> cutIDs;
  bool applyIdCuts;
  bool isHLT;
  float matchingDeltaR;
  std::string outFolder;

};

#endif