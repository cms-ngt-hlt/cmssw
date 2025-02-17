"""
This script runs the SimDoubletsProducer and SimDoubletsAnalyzer.
It is just meant for testing and development.

!!! NOTE !!!
You have to change at least the input file for now.
"""

import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process("SIMDOUBLETS",Phase2C17I13M9)

# maximum number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    fileNames = cms.untracked.vstring('file:../1_performance_comparison_legacy_alpaka/legacy/Phase2_L1P2GT_HLT_debug.root'),
    inputCommands = cms.untracked.vstring(
        'keep *',
        'drop *_hltSiPixelRecHits_*_HLTX'   # we will reproduce them to have their local position available
    ),
    secondaryFileNames = cms.untracked.vstring()
)


### conditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

### standard includes
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

### load hltTPClusterProducer
process.load("Validation.RecoTrack.associators_cff")
### load hltSiPixelRecHit Producer
process.load("HLTrigger.Configuration.HLT_75e33.modules.hltSiPixelRecHits_cfi")
### load the new EDProducer "SimDoubletsProducer"
process.load("SimTracker.TrackerHitAssociation.simDoubletsProducerPhase2_cfi")
### load the new DQM EDAnalyzer "SimDoubletsAnalyzer"
process.load("Validation.TrackingMCTruth.simDoubletsAnalyzerPhase2_cfi")

####  set up the paths
process.simDoubletProduction = cms.Path(
    process.hltTPClusterProducer *
    process.hltSiPixelRecHits *
    process.simDoubletsProducerPhase2 *
    process.simDoubletsAnalyzerPhase2
)

# Output definition
process.MyTestoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_*_*_SIMDOUBLETS'  # just keep the newly produced branches
    ),
    fileName = cms.untracked.string('file:simDoublets_TEST.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('TEST')
    )
)


process.endjob_step = cms.EndPath(process.endOfProcess)
process.MyTestoutput_step = cms.EndPath(process.MyTestoutput)


process.schedule = cms.Schedule(
      process.simDoubletProduction,process.endjob_step,process.MyTestoutput_step
)

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(1),
    wantSummary = cms.untracked.bool(True)
)
