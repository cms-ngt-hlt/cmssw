import FWCore.ParameterSet.Config as cms

from RecoTracker.FinalTrackSelectors.TrackOnlineDnnSelector import TrackOnlineDnnSelector as _TrackOnlineDnnSelector
TrackOnlineDnnSelector = _TrackOnlineDnnSelector(
    maxRecHits = 16,
    probThreshold = 0.095,
    tracksSrc = "hltPhase2PixelTracksCAExtension",
    beamSpot = "hltOnlineBeamSpot"
)
