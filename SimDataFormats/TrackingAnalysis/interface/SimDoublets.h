#ifndef SimDataFormats_TrackingAnalysis_SimDoublets_h
#define SimDataFormats_TrackingAnalysis_SimDoublets_h

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

// FIXME better move those definitions to the appropriate package for SiPixelRecHits
// persistent reference to a SiPixelRecHit in a SiPixelRecHitCollection
typedef edm::Ref<SiPixelRecHitCollection, SiPixelRecHit> SiPixelRecHitRef;
// persistent vector of references to SiPixelRecHits in a SiPixelRecHitCollection
typedef edm::RefVector<SiPixelRecHitCollection, SiPixelRecHit> SiPixelRecHitRefVector;

/** @brief Semi-Monte Carlo truth information used for pixel-tracking opimization.
 *
 * SimDoublets holds references to all pixel RecHits of a MC TrackingParticle.
 * Ones those RecHits are sorted according to their global position by the method
 * sortRecHits(), you can create the true doublets of RecHits that the TrackingParticle
 * left in the detector. These SimDoublets::Doublet objects can be used to optimize the 
 * doublet creation in the reconstruction.
 *
 * The Doublets are generated as the RecHit pairs between two consecutively hit layers.
 * I.e., if a TrackingParticle produces
 *  - 1 hit (A) in 1st layer
 *  - 2 hits (B, C) in 3rd layer
 *  - 1 hit (D) in 4th layer
 * then, the true Doublets are:
 *  (A-B), (A-C), (B-D) and (C-D).
 * So, neither does it matter that the 2nd layer got "skipped" as there are no hits,
 * nor is the Doublet of (A-D) formed since there is a layer with hits in between.
 * Doublets are not created between hits within the same layer.
 *
 * @author Jan Schulz (jan.gerrit.schulz@cern.ch)
 * @date January 2025
 */
class SimDoublets {
public:
  /**
     * Sub-class for true doublets of RecHits
     *  - first hit = inner RecHit
     *  - second hit = outer RecHit
     */
  class Doublet {
  public:
    enum layer { inner = 0, outer = 1 };  // not used anywhere yet...

    // default constructor
    Doublet() {}

    // constructor with explicit setting of useClusterLocalPosition_
    Doublet(SimDoublets const&, size_t const, size_t const, const TrackerGeometry*, const TrackerTopology*, bool);
    // constrcutor with automatic setting of useClusterLocalPosition_
    // checks if RecHit localPosition is meaningful: if yes, prefer RecHit localPosition over cluster
    Doublet(SimDoublets const&, size_t const, size_t const, const TrackerGeometry*, const TrackerTopology*);

    // method to access the layers pair
    std::pair<uint8_t, uint8_t> layerIds() const { return layerIds_; }

    // method to access the RecHit pair
    std::pair<SiPixelRecHitRef, SiPixelRecHitRef> recHits() const { return recHitRefs_; }

    // method to access the number of skipped layers
    int8_t numSkippedLayers() const { return numSkippedLayers_; }

    // method to access the layer pair ID
    int16_t layerPairId() const { return layerPairId_; }

    // method to access the inner layerId
    uint8_t innerLayerId() const { return layerIds_.first; }

    // method to access the outer layerId
    uint8_t outerLayerId() const { return layerIds_.second; }

    // method to access the inner RecHit
    SiPixelRecHitRef innerRecHit() const { return recHitRefs_.first; }

    // method to access the outer RecHit
    SiPixelRecHitRef outerRecHit() const { return recHitRefs_.second; }

    // method to access the global position of the inner RecHit
    GlobalPoint innerGlobalPos() const;

    // method to access the global position of the outer RecHit
    GlobalPoint outerGlobalPos() const;

  private:
    const TrackerGeometry* trackerGeometry_ = nullptr;  // pointer to the tracker geometry
    bool useClusterLocalPosition_;  // bool that decides whether to use the local position of the RecHit or the cluster
    TrackingParticleRef trackingParticleRef_;                   // reference to the TrackingParticle
    std::pair<SiPixelRecHitRef, SiPixelRecHitRef> recHitRefs_;  // reference pair to RecHits of the Doublet
    std::pair<uint8_t, uint8_t> layerIds_;                      // pair of layer IDs corresponding to the RecHits
    int8_t numSkippedLayers_;                                   // number of layers skipped by the Doublet
    int16_t layerPairId_;            // ID of the layer pair as defined in the reconstruction for the doublets
    GlobalVector beamSpotPosition_;  // global position of the beam spot (needed to correct the global RecHit position)
  };

  // default contructor
  SimDoublets() {}

  // constructor
  SimDoublets(TrackingParticleRef const trackingParticleRef, reco::BeamSpot const& beamSpot)
      : trackingParticleRef_(trackingParticleRef), beamSpotPosition_(beamSpot.x0(), beamSpot.y0(), beamSpot.z0()) {}

  // method to add a RecHitRef with its layer
  void addRecHit(SiPixelRecHitRef const recHitRef, uint8_t const layerId) {
    recHitsAreSorted_ = false;  // set sorted to false again
    // check if the layerId is not present in the layerIdVector yet
    if (std::find(layerIdVector_.begin(), layerIdVector_.end(), layerId) == layerIdVector_.end()) {
      // if it does not exist, increment number of layers
      numLayers_++;
    }
    // add recHit and layerId to the vectors
    recHitRefVector_.push_back(recHitRef);
    layerIdVector_.push_back(layerId);
  }

  // method to access the reference to the TrackingParticle
  TrackingParticleRef trackingParticle() const { return trackingParticleRef_; }

  // method to access the reference vector to the RecHits
  SiPixelRecHitRefVector recHits() const { return recHitRefVector_; }

  // method to access the RecHit at index i
  SiPixelRecHitRef recHits(size_t i) const { return recHitRefVector_[i]; }

  // method to access the layer id vector
  std::vector<uint8_t> layerIds() const { return layerIdVector_; }

  // method to access the layer id at index i
  uint8_t layerIds(size_t i) const { return layerIdVector_[i]; }

  // method to access the beam spot position
  GlobalVector beamSpotPosition() const { return beamSpotPosition_; }

  // method to get the number of layers
  int numLayers() const { return numLayers_; }

  // method to get number of RecHits in the SimDoublets
  int numRecHits() const { return layerIdVector_.size(); }

  // method to sort the RecHits according to the position
  void sortRecHits(const TrackerGeometry*);

  // method to produce the SimDoublets from the RecHits
  std::vector<Doublet> getSimDoublets(const TrackerGeometry* trackerGeometry = nullptr, const TrackerTopology* trackerTopology = nullptr) const;

private:
  TrackingParticleRef trackingParticleRef_;  // reference to the TrackingParticle
  SiPixelRecHitRefVector recHitRefVector_;   // reference vector to RecHits associated to the TP (sorted afer building)
  std::vector<uint8_t> layerIdVector_;       // vector of layer IDs corresponding to the RecHits
  GlobalVector beamSpotPosition_;  // global position of the beam spot (needed to correct the global RecHit position)
  bool recHitsAreSorted_{false};   // true if RecHits were sorted
  int numLayers_{0};               // number of layers hit by the TrackingParticle
};

// collection of SimDoublets
typedef std::vector<SimDoublets> SimDoubletsCollection;

#endif
