import FWCore.ParameterSet.Config as cms

hltTrackOnlineDnnSelector = cms.EDProducer("TrackOnlineDnnSelector",
    beamSpot = cms.InputTag("hltOnlineBeamSpot"),
    inputNames = cms.vstring('hits', 'tracks'),
    maxRecHits = cms.uint32(16),
    onnxModelPath = cms.FileInPath('RecoTracker/FinalTrackSelectors/data/TrackOnlineDnnSelector/track_classifier.onnx'),
    output_en = cms.vstring('probabilities'),
    probThreshold = cms.double(0.095),
    tracksSrc = cms.InputTag("hltPhase2PixelTracksCAExtension"),
)
