import FWCore.ParameterSet.Config as cms

hltSiPixelRecHits = cms.EDProducer('SiPixelRecHitFromSoAAlpaka',
    pixelRecHitSrc = cms.InputTag('hltPhase2SiPixelRecHitsSoA'),
    src = cms.InputTag('hltSiPixelClusters'),
)

_hltSiPixelRecHits = cms.EDProducer("SiPixelRecHitConverter",
    CPE = cms.string('PixelCPEGeneric'),
    src = cms.InputTag("hltSiPixelClusters")
)

from Configuration.ProcessModifiers.phase2LegacyPixelTracks_cff import phase2LegacyPixelTracks
phase2LegacyPixelTracks.toReplaceWith(hltSiPixelRecHits, _hltSiPixelRecHits)
