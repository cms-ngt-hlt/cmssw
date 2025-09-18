import FWCore.ParameterSet.Config as cms

from RecoVertex.PixelVertexFinding.PixelVertexProducerAlpakaPhase2_alpaka import PixelVertexProducerAlpakaPhase2_alpaka as _PixelVertexProducerAlpakaPhase2_alpaka

hltPhase2PixelVerticesSoA = _PixelVertexProducerAlpakaPhase2_alpaka(
    pixelTrackSrc = "hltPhase2PixelTracksSoA",
    #eps = 0.07,
    #errmax = 0.01,
    maxVertices = 512,
    useDBSCAN = cms.bool(False),
    useDensity = cms.bool(True),
    useIterative = cms.bool(False),
)
