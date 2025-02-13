// -*- C++ -*-
//
// Package:    Validation/TrackingMCTruth
// Class:      SimDoubletsAnalyzer
//
/**\class SimDoubletsAnalyzer SimDoubletsAnalyzer.cc Validation/TrackingMCTruth/plugins/SimDoubletsAnalyzer.cc

 Description: DQM analyzer for true RecHit doublets (SimDoublets) of the inner tracker in Phase-2 HLT

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Luca Ferragina, Elena Vernazza, Jan Schulz
//         Created:  Thu, 16 Jan 2025 13:46:21 GMT
//
//

#include <string>
#include <map>
#include <vector>

// user include files
#include "DataFormats/Histograms/interface/MonitorElementCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/approx_atan2.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/SimDoublets.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "Geometry/CommonTopologies/interface/SimplePixelTopology.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DQMServices/Core/interface/DQMEDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"

#define CLUSTER_DEBUG

// -------------------------------------------------------------------------------------------------------------
// class declaration
// -------------------------------------------------------------------------------------------------------------

class SimDoubletsAnalyzer : public DQMEDAnalyzer {
public:
  explicit SimDoubletsAnalyzer(const edm::ParameterSet&);
  ~SimDoubletsAnalyzer() override;

  void dqmBeginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;

  void analyze(const edm::Event&, const edm::EventSetup&) override;

  // ------------ member data ------------

  const TrackerGeometry* trackerGeometry_ = nullptr;
  const TrackerTopology* trackerTopology_ = nullptr;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geometry_getToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> topology_getToken_;
  const edm::EDGetTokenT<SimDoubletsCollection> simDoublets_getToken_;
  const edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster>> siPixelClusters_getToken_;

  // cutting parameters
  std::vector<double> cellMinz_;
  std::vector<double> cellMaxz_;
  std::vector<int> cellPhiCuts_;
  std::vector<double> cellMaxr_;
  int cellMinYSizeB1_;
  int cellMinYSizeB2_;
  int cellMaxDYSize12_;
  int cellMaxDYSize_;
  int cellMaxDYPred_;
  double cellZ0Cut_;
  double cellPtCut_;

  std::string folder_;  // main folder in the DQM file
  int eventCount_ = 0;  // event counter

  // monitor elements (histograms) to be filled
  MonitorElement* h_layerPairs_;
  MonitorElement* h_numSkippedLayers_;
  MonitorElement* h_numSimDoubletsPerTrackingParticle_;
  MonitorElement* h_numLayersPerTrackingParticle_;
  MonitorElement* h_numMissedLayersPerTrackingParticle_;
  MonitorElement* h_numTPVsPt_;
  MonitorElement* h_pass_numTPVsPt_;
  MonitorElement* h_numTPVsEta_;
  MonitorElement* h_pass_numTPVsEta_;
  MonitorElement* h_numVsPt_;
  MonitorElement* h_pass_numVsPt_;
  MonitorElement* h_numVsEta_;
  MonitorElement* h_pass_numVsEta_;
  MonitorElement* h_z0_;
  MonitorElement* h_curvatureR_;
  MonitorElement* h_pTFromR_;
  MonitorElement* h_YsizeB1_;
  MonitorElement* h_YsizeB2_;
  MonitorElement* h_DYsize12_;
  MonitorElement* h_DYsize_;
  MonitorElement* h_DYPred_;
  MonitorElement* h_pass_layerPairs_;
  MonitorElement* h_pass_z0_;
  MonitorElement* h_pass_pTFromR_;
  MonitorElement* h_pass_YsizeB1_;
  MonitorElement* h_pass_YsizeB2_;
  MonitorElement* h_pass_DYsize12_;
  MonitorElement* h_pass_DYsize_;
  MonitorElement* h_pass_DYPred_;
  std::vector<MonitorElement*> hVector_dr_;
  std::vector<MonitorElement*> hVector_dphi_;
  std::vector<MonitorElement*> hVector_idphi_;
  std::vector<MonitorElement*> hVector_innerZ_;
  std::vector<MonitorElement*> hVector_Ysize_;
  std::vector<MonitorElement*> hVector_DYsize_;
  std::vector<MonitorElement*> hVector_DYPred_;
  std::vector<MonitorElement*> hVector_pass_dr_;
  std::vector<MonitorElement*> hVector_pass_idphi_;
  std::vector<MonitorElement*> hVector_pass_innerZ_;
};

namespace simdoublets {
  std::pair<std::string, std::string> getInnerOuterLayerNames(int const layerPairId) {
    // make a string from the Id (int)
    std::string index = std::to_string(layerPairId);
    // determine inner and outer layer name
    std::string innerLayerName;
    std::string outerLayerName;
    if (index.size() < 3) {
      innerLayerName = "0";
      outerLayerName = index;
    } else if (index.size() == 3) {
      innerLayerName = index.substr(0, 1);
      outerLayerName = index.substr(1, 3);
    } else {
      innerLayerName = index.substr(0, 2);
      outerLayerName = index.substr(2, 4);
    }
    if (outerLayerName[0] == '0') {
      outerLayerName = outerLayerName.substr(1, 2);
    }

    return {innerLayerName, outerLayerName};
  }

  // make bins logarithmic
  void BinLogX(TH1* h) {
    TAxis* axis = h->GetXaxis();
    int bins = axis->GetNbins();

    float from = axis->GetXmin();
    float to = axis->GetXmax();
    float width = (to - from) / bins;
    std::vector<float> new_bins(bins + 1, 0);

    for (int i = 0; i <= bins; i++) {
      new_bins[i] = TMath::Power(10, from + i * width);
    }
    axis->Set(bins, new_bins.data());
  }

  // function to produce histogram with log scale on x (taken from MultiTrackValidator)
  template <typename... Args>
  dqm::reco::MonitorElement* make1DLogX(dqm::reco::DQMStore::IBooker& ibook, Args&&... args) {
    auto h = std::make_unique<TH1F>(std::forward<Args>(args)...);
    BinLogX(h.get());
    const auto& name = h->GetName();
    return ibook.book1D(name, h.release());
  }

  // function that checks if two vector share a common element
  template <typename T>
  bool haveCommonElement(std::vector<T> const& v1, std::vector<T> const& v2) {
    return std::find_first_of(v1.begin(), v1.end(), v2.begin(), v2.end()) != v1.end();
  }

  std::vector<uint8_t> getBarrelMissingLayerIds(const std::vector<uint8_t>& layerIds) {
    std::vector<uint8_t> missingLayerIds;
    for (uint8_t i = 0; i < 4; ++i) {
      if (std::find(layerIds.begin(), layerIds.end(), i) == layerIds.end()) {
        missingLayerIds.push_back(i);
      }
    }
    return missingLayerIds;
  }

  std::vector<DetId> getRecHitsDetIds(const SiPixelRecHitRefVector& recHits) {
    std::vector<DetId> detIds;
    for (const auto& recHit : recHits) {
      detIds.push_back(recHit->geographicalId());
    }
    return detIds;
  }

}  // namespace simdoublets

//
// static data member definitions
//
// map that takes the layerPairId as defined in the SimDoublets
// and gives the position of the histogram in the histogram vector
// NOTE: It is absolutely necessary that the map is sorted here,
// otherwise the histograms will not be labeled corresponding to the correct layer pair but are mixed up
static const std::map<int, int> layerPairId2Index{
    {1, 0},     {4, 1},     {16, 2},    {102, 3},   {104, 4},   {116, 5},   {203, 6},   {204, 7},   {216, 8},
    {405, 9},   {506, 10},  {607, 11},  {708, 12},  {809, 13},  {910, 14},  {1011, 15}, {1617, 16}, {1718, 17},
    {1819, 18}, {1920, 19}, {2021, 20}, {2122, 21}, {2223, 22}, {2, 23},    {5, 24},    {17, 25},   {6, 26},
    {18, 27},   {103, 28},  {105, 29},  {117, 30},  {106, 31},  {118, 32},  {1112, 33}, {1213, 34}, {1314, 35},
    {1415, 36}, {2324, 37}, {2425, 38}, {2526, 39}, {2627, 40}, {406, 41},  {507, 42},  {608, 43},  {709, 44},
    {810, 45},  {911, 46},  {1012, 47}, {1618, 48}, {1719, 49}, {1820, 50}, {1921, 51}, {2022, 52}, {2123, 53},
    {2224, 54}, {1315, 55}, {217, 56},  {205, 57},  {2325, 58}, {1113, 59}, {2426, 60}, {1214, 61}, {2527, 62}};

static const size_t numLayerPairs = layerPairId2Index.size();

// -------------------------------
// constructors and destructor
// -------------------------------
SimDoubletsAnalyzer::SimDoubletsAnalyzer(const edm::ParameterSet& iConfig)
    : geometry_getToken_(esConsumes<TrackerGeometry, TrackerDigiGeometryRecord, edm::Transition::BeginRun>()),
      topology_getToken_(esConsumes<TrackerTopology, TrackerTopologyRcd, edm::Transition::BeginRun>()),
      simDoublets_getToken_(consumes(iConfig.getParameter<edm::InputTag>("simDoubletsSrc"))),
      siPixelClusters_getToken_(consumes(iConfig.getParameter<edm::InputTag>("siPixelClustersSrc"))),
      cellMinz_(iConfig.getParameter<std::vector<double>>("cellMinz")),
      cellMaxz_(iConfig.getParameter<std::vector<double>>("cellMaxz")),
      cellPhiCuts_(iConfig.getParameter<std::vector<int>>("cellPhiCuts")),
      cellMaxr_(iConfig.getParameter<std::vector<double>>("cellMaxr")),
      cellMinYSizeB1_(iConfig.getParameter<int>("cellMinYSizeB1")),
      cellMinYSizeB2_(iConfig.getParameter<int>("cellMinYSizeB2")),
      cellMaxDYSize12_(iConfig.getParameter<int>("cellMaxDYSize12")),
      cellMaxDYSize_(iConfig.getParameter<int>("cellMaxDYSize")),
      cellMaxDYPred_(iConfig.getParameter<int>("cellMaxDYPred")),
      cellZ0Cut_(iConfig.getParameter<double>("cellZ0Cut")),
      cellPtCut_(iConfig.getParameter<double>("cellPtCut")),
      folder_(iConfig.getParameter<std::string>("folder")) {
  hVector_dr_.resize(numLayerPairs);
  hVector_dphi_.resize(numLayerPairs);
  hVector_idphi_.resize(numLayerPairs);
  hVector_innerZ_.resize(numLayerPairs);
  hVector_Ysize_.resize(numLayerPairs);
  hVector_DYsize_.resize(numLayerPairs);
  hVector_DYPred_.resize(numLayerPairs);
  hVector_pass_dr_.resize(numLayerPairs);
  hVector_pass_idphi_.resize(numLayerPairs);
  hVector_pass_innerZ_.resize(numLayerPairs);
}

SimDoubletsAnalyzer::~SimDoubletsAnalyzer() {}

// -----------------------
// member functions
// -----------------------

void SimDoubletsAnalyzer::dqmBeginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
  trackerGeometry_ = &iSetup.getData(geometry_getToken_);
  trackerTopology_ = &iSetup.getData(topology_getToken_);
}

void SimDoubletsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  eventCount_++;

  // get simDoublets
  SimDoubletsCollection const& simDoubletsCollection = iEvent.get(simDoublets_getToken_);
  edmNew::DetSetVector<SiPixelCluster> const& siPixelClusters = iEvent.get(siPixelClusters_getToken_);

  // create vectors for inner and outer RecHits of SimDoublets passing all cuts
  std::vector<SiPixelRecHitRef> innerRecHitsPassing;
  std::vector<SiPixelRecHitRef> outerRecHitsPassing;

  // loop over SimDoublets (= loop over TrackingParticles)
  for (auto const& simDoublets : simDoubletsCollection) {
    // get true pT of the TrackingParticle
    auto true_pT = simDoublets.trackingParticle()->pt();
    auto true_eta = simDoublets.trackingParticle()->eta();

    // create the true RecHit doublets of the TrackingParticle
    auto doublets = simDoublets.getSimDoublets(trackerGeometry_);

    int numSimDoublets = doublets.size();
    float weight = 1. / float(numSimDoublets);

    // fill histograms for number of SimDoublets
    h_numSimDoubletsPerTrackingParticle_->Fill(numSimDoublets);
    h_numLayersPerTrackingParticle_->Fill(simDoublets.numLayers());

    // get all detIds of the RecHits of the SimDoublets
    auto recHitsDetIds = simdoublets::getRecHitsDetIds(simDoublets.recHits());
    // get the missing layers for the TrackingParticle
    auto missingLayers = simdoublets::getBarrelMissingLayerIds(simDoublets.layerIds());

#ifdef CLUSTER_DEBUG
    // Find clusters that have no corresponding RecHits
    bool printRecHits = false;
    edmNew::DetSetVector<SiPixelCluster>::const_iterator it;
    for (it = siPixelClusters.begin(); it != siPixelClusters.end(); ++it) {
      auto clusterDetId = it->detId();
      const PixelGeomDetUnit* geomDetUnit =
          dynamic_cast<const PixelGeomDetUnit*>(trackerGeometry_->idToDetUnit(clusterDetId));
      const unsigned int layer = trackerTopology_->pxbLayer(geomDetUnit->geographicalId());
      const unsigned int ladder = trackerTopology_->pxbLadder(geomDetUnit->geographicalId());
      const unsigned int module = trackerTopology_->pxbModule(geomDetUnit->geographicalId());
      for (const SiPixelCluster& cluster : *it) {
        bool missingLayer = std::find(missingLayers.begin(), missingLayers.end(), layer) != missingLayers.end();
        bool missingDetId = std::find(recHitsDetIds.begin(), recHitsDetIds.end(), clusterDetId) == recHitsDetIds.end();
        if (missingLayer && missingDetId) {
          printRecHits = true;
          LocalPoint clusterLocalPos =
              geomDetUnit->specificTopology().localPosition(MeasurementPoint(cluster.x(), cluster.y()));
          GlobalPoint clusterGlobalPos = geomDetUnit->surface().toGlobal(clusterLocalPos);
          std::cout << "Found cluster in (layer, ladder, module, detId): (" << layer << ", " << ladder << ", " << module
                    << ", " << clusterDetId << ") where no TP recHits were produced!" << std::endl;
          std::cout << "Cluster size: " << cluster.size() << ", local position:" << clusterLocalPos
                    << ", global position:" << clusterGlobalPos << std::endl;
          h_numMissedLayersPerTrackingParticle_->Fill(layer);
        }
      }
    }
    // print RecHits for TPs with at least a cluster without a recHit
    if (printRecHits) {
      std::cout << "RecHits for this TrackingParticle:" << std::endl;
      int i = 1;
      for (const auto& recHit : simDoublets.recHits()) {
        auto id = recHit->geographicalId();
        const unsigned int layer = trackerTopology_->pxbLayer(id);
        const unsigned int ladder = trackerTopology_->pxbLadder(id);
        const unsigned int module = trackerTopology_->pxbModule(id);
        std::cout << "recHit " << i << " in (layer, ladder, module, detId): (" << layer << ", " << ladder << ", "
                  << module << ", " << id << ")" << std::endl;
        std::cout << "recHit local position:" << recHit->localPosition()
                  << ", recHit global position:" << recHit->globalPosition() << std::endl;
        std::cout << "---------------------------------------------------------" << std::endl;
        ++i;
      }
    }
#endif

    // fill histograms for number of TrackingParticles
    h_numTPVsPt_->Fill(true_pT);
    h_numTPVsEta_->Fill(true_eta);

    // clear passing inner and outer RecHits
    innerRecHitsPassing.clear();
    outerRecHitsPassing.clear();

    // loop over those doublets
    for (auto const& doublet : doublets) {
      // RecHit properties
      auto inner_r = doublet.innerGlobalPos().perp();
      auto inner_z = doublet.innerGlobalPos().z();
      auto inner_phi = doublet.innerGlobalPos().barePhi();  // returns float, whereas .phi() returns phi object
      auto inner_iphi = unsafe_atan2s<7>(doublet.innerGlobalPos().y(), doublet.innerGlobalPos().x());
      auto outer_r = doublet.outerGlobalPos().perp();
      auto outer_z = doublet.outerGlobalPos().z();
      auto outer_phi = doublet.outerGlobalPos().barePhi();
      auto outer_iphi = unsafe_atan2s<7>(doublet.outerGlobalPos().y(), doublet.outerGlobalPos().x());

      auto dz = outer_z - inner_z;
      auto dr = outer_r - inner_r;
      auto dphi = reco::deltaPhi(inner_phi, outer_phi);
      auto idphi = std::min(std::abs(int16_t(outer_iphi - inner_iphi)), std::abs(int16_t(inner_iphi - outer_iphi)));

      // ----------------------------------------------------------
      // layer pair independent plots (main folder)
      // ----------------------------------------------------------

      // outer layer vs inner layer of SimDoublets
      h_layerPairs_->Fill(doublet.innerLayerId(), doublet.outerLayerId());

      // number of skipped layers by SimDoublets
      h_numSkippedLayers_->Fill(doublet.numSkippedLayers());

      // longitudinal impact parameter with respect to the beamspot
      double z0 = std::abs(inner_r * outer_z - inner_z * outer_r) / dr;
      h_z0_->Fill(z0);

      // radius of the circle defined by the two RecHits and the beamspot
      auto curvature = 1.f / 2.f * std::sqrt((dr / dphi) * (dr / dphi) + (inner_r * outer_r));
      h_curvatureR_->Fill(curvature);

      // pT that this curvature radius corresponds to
      auto pT = curvature / 87.78f;
      h_pTFromR_->Fill(pT);

      // ----------------------------------------------------------
      // layer pair dependent plots (sub-folders for layer pairs)
      // ----------------------------------------------------------

      // first, get layer pair Id and exclude layer pairs that are not considered
      int layerPairId = doublet.layerPairId();
      if (layerPairId2Index.find(layerPairId) == layerPairId2Index.end()) {
        continue;
      }

      // get the position of the layer pair in the vectors of histograms
      int layerPairIdIndex = layerPairId2Index.at(layerPairId);

      // dr = (outer_r - inner_r) histogram
      hVector_dr_[layerPairIdIndex]->Fill(dr);

      // dphi histogram
      hVector_dphi_[layerPairIdIndex]->Fill(dphi);
      hVector_idphi_[layerPairIdIndex]->Fill(idphi);

      // z of the inner RecHit histogram
      hVector_innerZ_[layerPairIdIndex]->Fill(inner_z);

      // ----------------------------------------------------------
      // cluster size plots (main + sub-folders for layer pairs)
      // ----------------------------------------------------------

      // cluster size in local y histogram
      auto innerClusterSizeY = doublet.innerRecHit()->cluster()->sizeY();
      hVector_Ysize_[layerPairIdIndex]->Fill(innerClusterSizeY);

      // create bool that indicates if the doublet gets cut
      bool doubletGetsCut = false;
      // create bools that trace if doublet is subject to any clsuter size cut
      bool subjectToYsizeB1 = false;
      bool subjectToYsizeB2 = false;
      bool subjectToDYsize = false;
      bool subjectToDYsize12 = false;
      bool subjectToDYPred = false;

      // apply all cuts that do not depend on the cluster size
      // z window cut
      if (inner_z < cellMinz_[layerPairIdIndex] || inner_z > cellMaxz_[layerPairIdIndex]) {
        doubletGetsCut = true;
      }
      // z0cutoff
      if (dr > cellMaxr_[layerPairIdIndex] || dr < 0 || z0 > cellZ0Cut_) {
        doubletGetsCut = true;
      }
      // ptcut
      if (pT < cellPtCut_) {
        doubletGetsCut = true;
      }
      // iphicut
      if (idphi > cellPhiCuts_[layerPairIdIndex]) {
        doubletGetsCut = true;
      }

      // determine the moduleId
      DetId detIdObject(doublet.innerRecHit()->geographicalId());
      const GeomDetUnit* geomDetUnit = trackerGeometry_->idToDetUnit(detIdObject);
      const uint32_t moduleId = geomDetUnit->index();

      // define bools needed to decide on cutting parameters
      const bool innerInB1 = (doublet.innerLayerId() == 0);
      const bool innerInB2 = (doublet.innerLayerId() == 1);
      const bool isOuterLadder = (0 == (moduleId / 8) % 2);  // check if this even makes sense in Phase-2
      const bool innerInBarrel = (doublet.innerLayerId() < 4);
      const bool outerInBarrel = (doublet.outerLayerId() < 4);

      // histograms for clusterCut
      // cluster size in local y
      if (!outerInBarrel) {
        if (innerInB1 && isOuterLadder) {
          subjectToYsizeB1 = true;
          h_YsizeB1_->Fill(innerClusterSizeY);
          // apply the cut
          if (innerClusterSizeY < cellMinYSizeB1_) {
            doubletGetsCut = true;
          }
        }
        if (innerInB2) {
          subjectToYsizeB2 = true;
          h_YsizeB2_->Fill(innerClusterSizeY);
          // apply the cut
          if (innerClusterSizeY < cellMinYSizeB2_) {
            doubletGetsCut = true;
          }
        }
      }

      // histograms for zSizeCut
      int DYsize{0}, DYPred{0};
      if (innerInBarrel) {
        if (outerInBarrel) {  // onlyBarrel
          DYsize = std::abs(innerClusterSizeY - doublet.outerRecHit()->cluster()->sizeY());
          if (innerInB1 && isOuterLadder) {
            subjectToDYsize12 = true;
            hVector_DYsize_[layerPairIdIndex]->Fill(DYsize);
            h_DYsize12_->Fill(DYsize);
            // apply the cut
            if (DYsize > cellMaxDYSize12_) {
              doubletGetsCut = true;
            }
          } else if (!innerInB1) {
            subjectToDYsize = true;
            hVector_DYsize_[layerPairIdIndex]->Fill(DYsize);
            h_DYsize_->Fill(DYsize);
            // apply the cut
            if (DYsize > cellMaxDYSize_) {
              doubletGetsCut = true;
            }
          }
        } else {  // not onlyBarrel
          subjectToDYPred = true;
          DYPred = std::abs(innerClusterSizeY - int(std::abs(dz / dr) * pixelTopology::Phase2::dzdrFact + 0.5f));
          hVector_DYPred_[layerPairIdIndex]->Fill(DYPred);
          h_DYPred_->Fill(DYPred);
          // apply the cut
          if (DYPred > cellMaxDYPred_) {
            doubletGetsCut = true;
          }
        }
      }

      // fill the number histograms
      // histogram of all valid doublets
      h_numVsPt_->Fill(true_pT);
      h_numVsEta_->Fill(true_eta);
      // fill histogram of doublets that pass all cuts
      if (!doubletGetsCut) {
        h_pass_layerPairs_->Fill(doublet.innerLayerId(), doublet.outerLayerId());
        h_pass_numVsPt_->Fill(true_pT, weight);
        h_pass_numVsEta_->Fill(true_eta, weight);

        // also put the inner/outer RecHit in the respective vector
        innerRecHitsPassing.push_back(doublet.innerRecHit());
        outerRecHitsPassing.push_back(doublet.outerRecHit());

        // fill pass_ histograms
        h_pass_z0_->Fill(z0);
        h_pass_pTFromR_->Fill(pT);
        hVector_pass_dr_[layerPairIdIndex]->Fill(dr);
        hVector_pass_idphi_[layerPairIdIndex]->Fill(idphi);
        hVector_pass_innerZ_[layerPairIdIndex]->Fill(inner_z);
        if (subjectToDYPred) {
          h_pass_DYPred_->Fill(DYPred);
        }
        if (subjectToDYsize) {
          h_pass_DYsize_->Fill(DYsize);
        }
        if (subjectToDYsize12) {
          h_pass_DYsize12_->Fill(DYsize);
        }
        if (subjectToYsizeB1) {
          h_pass_YsizeB1_->Fill(innerClusterSizeY);
        }
        if (subjectToYsizeB2) {
          h_pass_YsizeB2_->Fill(innerClusterSizeY);
        }
      }
    }  // end loop over those doublets

    // Now check if the TrackingParticle is reconstructable by at least two conencted SimDoublets surviving the cuts
    if (simdoublets::haveCommonElement<SiPixelRecHitRef>(innerRecHitsPassing, outerRecHitsPassing)) {
      h_pass_numTPVsPt_->Fill(true_pT);
      h_pass_numTPVsEta_->Fill(true_eta);
    }
  }  // end loop over SimDoublets (= loop over TrackingParticles)
}

void SimDoubletsAnalyzer::bookHistograms(DQMStore::IBooker& ibook, edm::Run const& run, edm::EventSetup const& iSetup) {
  // set some common parameters
  int pTNBins = 50;
  double pTmin = log10(0.01);
  double pTmax = log10(1000);
  int etaNBins = 80;
  double etamin = -4.;
  double etamax = 4.;

  // ----------------------------------------------------------
  // booking general histograms (general folder)
  // ----------------------------------------------------------

  ibook.setCurrentFolder(folder_ + "/general");

  // overview histograms
  h_layerPairs_ = ibook.book2D(
      "layerPairs", "Layer pairs in SimDoublets; Inner layer ID; Outer layer ID", 28, -0.5, 27.5, 28, -0.5, 27.5);
  h_pass_layerPairs_ = ibook.book2D(
      "pass_layerPairs", "Layer pairs in SimDoublets passing all cuts; Inner layer ID; Outer layer ID", 28, -0.5, 27.5, 28, -0.5, 27.5);
  h_numSkippedLayers_ = ibook.book1D(
      "numSkippedLayers", "Number of skipped layers; Number of skipped layers; Number of SimDoublets", 16, -1.5, 14.5);
  h_numSimDoubletsPerTrackingParticle_ =
      ibook.book1D("numSimDoubletsPerTrackingParticle",
                   "Number of SimDoublets per Tracking Particle; Number of SimDoublets; Number of Tracking Particles",
                   31,
                   -0.5,
                   30.5);
  h_numLayersPerTrackingParticle_ =
      ibook.book1D("numLayersPerTrackingParticle",
                   "Number of layers hit by Tracking Particle; Number of layers; Number of Tracking Particles",
                   29,
                   -0.5,
                   28.5);
  h_numMissedLayersPerTrackingParticle_ =
      ibook.book1D("numMissedLayersPerTrackingParticle",
                   "Number of missed layers by Tracking Particle; Layer ID; Number of "
                   "Tracking Particles",
                   5,
                   0,
                   5);
  h_numTPVsPt_ = simdoublets::make1DLogX(
      ibook,
      "numTPVsPt",
      "Total number of TrackingParticles; True transverse momentum p_{T} [GeV]; Total number of TrackingParticles",
      pTNBins,
      pTmin,
      pTmax);
  h_pass_numTPVsPt_ = simdoublets::make1DLogX(ibook,
                                             "pass_numTPVsPt",
                                             "Reconstructable TrackingParticles (two or more connected SimDoublets "
                                             "pass cuts); True transverse momentum p_{T} [GeV]; "
                                             "Number of reconstructable TrackingParticles",
                                             pTNBins,
                                             pTmin,
                                             pTmax);
  h_numTPVsEta_ =
      ibook.book1D("numTPVsEta",
                   "Total number of TrackingParticles; True pseudorapidity #eta; Total number of TrackingParticles",
                   etaNBins,
                   etamin,
                   etamax);
  h_pass_numTPVsEta_ = ibook.book1D("pass_numTPVsEta",
                                   "Reconstructable TrackingParticles (two or more connected SimDoublets "
                                   "pass cuts); True pseudorapidity #eta; Number of reconstructable TrackingParticles",
                                   etaNBins,
                                   etamin,
                                   etamax);
  h_numVsPt_ = simdoublets::make1DLogX(
      ibook,
      "numVsPt",
      "Total number of SimDoublets; True transverse momentum p_{T} [GeV]; Total number of SimDoublets",
      pTNBins,
      pTmin,
      pTmax);
  h_pass_numVsPt_ =
      simdoublets::make1DLogX(ibook,
                              "pass_numVsPt",
                              "Weighted number of passing SimDoublets; True transverse momentum p_{T} [GeV]; "
                              "Number of SimDoublets passing all cuts",
                              pTNBins,
                              pTmin,
                              pTmax);
  h_numVsEta_ = ibook.book1D("numVsEta",
                                "Total number of SimDoublets; True pseudorapidity #eta; Total number of SimDoublets",
                                etaNBins,
                                etamin,
                                etamax);
  h_pass_numVsEta_ =
      ibook.book1D("pass_numVsEta",
                   "Weighted number of SimDoublets; True pseudorapidity #eta; Number of SimDoublets passing all cuts",
                   etaNBins,
                   etamin,
                   etamax);

  // ----------------------------------------------------------
  // booking layer pair independent cut histograms (global folder)
  // ----------------------------------------------------------

  ibook.setCurrentFolder(folder_ + "/cutParameters/global");

  // histogram for z0cutoff  (z0Cut)
  h_z0_ = ibook.book1D("z0", "z_{0}; Longitudinal impact parameter z_{0} [cm]; Number of SimDoublets", 51, -1, 50);
  h_pass_z0_ = ibook.book1D(
      "pass_z0",
      "z_{0} of SimDoublets passing all cuts; Longitudinal impact parameter z_{0} [cm]; Number of SimDoublets",
      51,
      -1,
      50);

  // histograms for ptcut  (ptCut)
  h_curvatureR_ = ibook.book1D(
      "curvatureR", "Curvature from SimDoublet+beamspot; Curvature radius [cm] ; Number of SimDoublets", 100, 0, 1000);
  h_pTFromR_ = simdoublets::make1DLogX(
      ibook,
      "pTFromR",
      "Transverse momentum from curvature; Transverse momentum p_{T} [GeV]; Number of SimDoublets",
      pTNBins,
      pTmin,
      pTmax);
  h_pass_pTFromR_ = simdoublets::make1DLogX(ibook,
                                            "pass_pTFromR",
                                            "Transverse momentum from curvature of SimDoublets passing all cuts; "
                                            "Transverse momentum p_{T} [GeV]; Number of SimDoublets",
                                            pTNBins,
                                            pTmin,
                                            pTmax);

  // histograms for clusterCut  (minYsizeB1 and minYsizeB2)
  h_YsizeB1_ = ibook.book1D(
      "YsizeB1",
      "Cluster size along z (inner from B1); Size along z of inner cluster [num of pixels]; Number of SimDoublets",
      51,
      -1,
      50);
  h_YsizeB2_ = ibook.book1D(
      "YsizeB2",
      "Cluster size along z (inner not from B1); Size along z of inner cluster [num of pixels]; Number of SimDoublets",
      51,
      -1,
      50);
  h_pass_YsizeB1_ = ibook.book1D("pass_YsizeB1",
                                 "Cluster size along z of SimDoublets passing all cuts (inner from B1); Size along z "
                                 "of inner cluster [num of pixels]; Number of SimDoublets",
                                 51,
                                 -1,
                                 50);
  h_pass_YsizeB2_ = ibook.book1D("pass_YsizeB2",
                                 "Cluster size along z of SimDoublets passing all cuts (inner not from B1); Size along "
                                 "z of inner cluster [num of pixels]; Number of SimDoublets",
                                 51,
                                 -1,
                                 50);

  // histograms for zSizeCut  (maxDYsize12, maxDYsize and maxDYPred)
  h_DYsize12_ =
      ibook.book1D("DYsize12",
                   "Difference in cluster size along z (inner from B1); Absolute difference in cluster size along z of "
                   "the two RecHits [num of pixels]; Number of SimDoublets",
                   31,
                   -1,
                   30);
  h_DYsize_ = ibook.book1D("DYsize",
                           "Difference in cluster size along z; Absolute difference in cluster size along z of the two "
                           "RecHits [num of pixels]; Number of SimDoublets",
                           31,
                           -1,
                           30);
  h_DYPred_ = ibook.book1D("DYPred",
                           "Difference between actual and predicted cluster size along z of inner cluster; Absolute "
                           "difference [num of pixels]; Number of SimDoublets",
                           201,
                           -1,
                           200);
  h_pass_DYsize12_ = ibook.book1D("pass_DYsize12",
                                  "Difference in cluster size along z of SimDoublets passing all cuts (inner from B1); "
                                  "Absolute difference in cluster size along z of "
                                  "the two RecHits [num of pixels]; Number of SimDoublets",
                                  31,
                                  -1,
                                  30);
  h_pass_DYsize_ =
      ibook.book1D("pass_DYsize",
                   "Difference in cluster size along z of SimDoublets passing all cuts; Absolute difference in "
                   "cluster size along z of the two RecHits [num of pixels]; Number of SimDoublets",
                   31,
                   -1,
                   30);
  h_pass_DYPred_ =
      ibook.book1D("pass_DYPred",
                   "Difference between actual and predicted cluster size along z of inner cluster of SimDoublets "
                   "passing all cuts; Absolute difference [num of pixels]; Number of SimDoublets",
                   201,
                   -1,
                   200);

  // -----------------------------------------------------------------------
  // booking layer pair dependent histograms (sub-folders for layer pairs)
  // -----------------------------------------------------------------------

  // loop through valid layer pairs and add for each one booked hist per vector
  for (auto id = layerPairId2Index.begin(); id != layerPairId2Index.end(); ++id) {
    // get the position of the layer pair in the histogram vectors
    int layerPairIdIndex = id->second;

    // get layer names from the layer pair Id
    auto layerNames = simdoublets::getInnerOuterLayerNames(id->first);
    std::string innerLayerName = layerNames.first;
    std::string outerLayerName = layerNames.second;

    // name the sub-folder for the layer pair "lp_${innerLayerId}_${outerLayerId}"
    std::string subFolderName = "/cutParameters/lp_" + innerLayerName + "_" + outerLayerName;

    // layer mentioning in histogram titles
    std::string layerTitle = "(layers (" + innerLayerName + "," + outerLayerName + "))";

    // set folder to the sub-folder for the layer pair
    ibook.setCurrentFolder(folder_ + subFolderName);

    // histogram for z0cutoff  (maxr)
    hVector_dr_.at(layerPairIdIndex) = ibook.book1D(
        "dr",
        "dr of RecHit pair " + layerTitle + "; dr between outer and inner RecHit [cm]; Number of SimDoublets",
        31,
        -1,
        30);
    hVector_pass_dr_.at(layerPairIdIndex) = ibook.book1D(
        "pass_dr",
        "dr of RecHit pair " + layerTitle +
            " for SimDoublets passing all cuts; dr between outer and inner RecHit [cm]; Number of SimDoublets",
        31,
        -1,
        30);

    // histograms for iphicut  (phiCuts)
    hVector_dphi_.at(layerPairIdIndex) = ibook.book1D(
        "dphi",
        "dphi of RecHit pair " + layerTitle + "; d#phi between outer and inner RecHit [rad]; Number of SimDoublets",
        50,
        -M_PI,
        M_PI);
    hVector_idphi_.at(layerPairIdIndex) =
        ibook.book1D("idphi",
                     "idphi of RecHit pair " + layerTitle +
                         "; Absolute int d#phi between outer and inner RecHit; Number of SimDoublets",
                     50,
                     0,
                     1000);
    hVector_pass_idphi_.at(layerPairIdIndex) = ibook.book1D("pass_idphi",
                                                            "idphi of RecHit pair " + layerTitle +
                                                                " for SimDoublets passing all cuts; Absolute int d#phi "
                                                                "between outer and inner RecHit; Number of SimDoublets",
                                                            50,
                                                            0,
                                                            1000);

    // histogram for z window  (minz and maxz)
    hVector_innerZ_.at(layerPairIdIndex) =
        ibook.book1D("innerZ",
                     "z of the inner RecHit " + layerTitle + "; z of inner RecHit [cm]; Number of SimDoublets",
                     100,
                     -300,
                     300);
    hVector_pass_innerZ_.at(layerPairIdIndex) =
        ibook.book1D("pass_innerZ",
                     "z of the inner RecHit " + layerTitle +
                         " for SimDoublets passing all cuts; z of inner RecHit [cm]; Number of SimDoublets",
                     100,
                     -300,
                     300);

    // histograms for cluster size and size differences
    hVector_DYsize_.at(layerPairIdIndex) =
        ibook.book1D("DYsize",
                     "Difference in cluster size along z between outer and inner RecHit " + layerTitle +
                         "; Absolute difference in cluster size along z of the two "
                         "RecHits [num of pixels]; Number of SimDoublets",
                     51,
                     -1,
                     50);
    hVector_DYPred_.at(layerPairIdIndex) =
        ibook.book1D("DYPred",
                     "Difference between actual and predicted cluster size along z of inner cluster " + layerTitle +
                         "; Absolute difference [num of pixels]; Number of SimDoublets",
                     51,
                     -1,
                     50);
    hVector_Ysize_.at(layerPairIdIndex) = ibook.book1D(
        "Ysize",
        "Cluster size along z " + layerTitle + "; Size along z of inner cluster [num of pixels]; Number of SimDoublets",
        51,
        -1,
        50);
  }
}

void SimDoubletsAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("folder", "Tracking/TrackingMCTruth/SimDoublets");
  desc.add<edm::InputTag>("simDoubletsSrc", edm::InputTag("simDoubletsProducer"));
  desc.add<edm::InputTag>("siPixelClustersSrc", edm::InputTag("hltSiPixelClusters"));

  // cutting parameters
  desc.add<std::vector<double>>("cellMinz", std::vector<double>(55, -20.))->setComment("Minimum z for each layer pair");
  desc.add<std::vector<double>>("cellMaxz", std::vector<double>(55, 20.))->setComment("Maximum z for each layer pair");
  desc.add<std::vector<int>>("cellPhiCuts", std::vector<int>(55, 20))->setComment("Cuts in phi for cells");
  desc.add<std::vector<double>>("cellMaxr", std::vector<double>(55, 20.))->setComment("Cut for dr of cells");
  desc.add<int>("cellMinYSizeB1", 25)->setComment("Minimum cluster size for B1");
  desc.add<int>("cellMinYSizeB2", 15)->setComment("Minimum cluster size for B2");
  desc.add<int>("cellMaxDYSize12", 12)->setComment("Maximum cluster size difference for B1/B2");
  desc.add<int>("cellMaxDYSize", 10)->setComment("Maximum cluster size difference");
  desc.add<int>("cellMaxDYPred", 20)->setComment("Maximum cluster size difference prediction");
  desc.add<double>("cellZ0Cut", 7.5)->setComment("Maximum longitudinal impact parameter");
  desc.add<double>("cellPtCut", 0.85)->setComment("Minimum tranverse momentum");

  descriptions.addWithDefaultLabel(desc);
}

// define this as a plug-in
DEFINE_FWK_MODULE(SimDoubletsAnalyzer);
