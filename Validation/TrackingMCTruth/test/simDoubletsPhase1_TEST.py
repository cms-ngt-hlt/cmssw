"""
This script runs the SimDoubletsProducer and SimDoubletsAnalyzer.
It is just meant for testing and development.

!!! NOTE !!!
You have to change at least the input file for now.
"""

import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_2025_cff import Run3_2025

process = cms.Process("SIMDOUBLETS",Run3_2025)

# maximum number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    fileNames = cms.untracked.vstring('file:step2.root'),
    inputCommands = cms.untracked.vstring(
        'keep *',
        'drop *_hltSiPixelRecHits_*_HLTX'   # we will reproduce them to have their local position available
    ),
    secondaryFileNames = cms.untracked.vstring()
)


### conditions
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '142X_mcRun3_2025_realistic_v4', '')

### standard includes
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
# process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

### load hltTPClusterProducer
process.load("Validation.RecoTrack.associators_cff")
### load the new EDProducer "SimDoubletsProducerPhase1"
process.load("SimTracker.TrackerHitAssociation.simDoubletsProducerPhase1_cfi")
### load the new DQM EDAnalyzer "SimDoubletsAnalyzerPhase1"
process.load("Validation.TrackingMCTruth.simDoubletsAnalyzerPhase1_cfi")

####  set up the paths
process.simDoubletProduction = cms.Path(
    process.HLTDoLocalPixelSequence *
    process.hltTPClusterProducer *
    process.simDoubletsProducerPhase1 *
    process.simDoubletsAnalyzerPhase1
)

# Output definition
process.MyTestoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_simDoubletsProducerPhase1_*_SIMDOUBLETS'  # just keep the newly produced branches
    ),
    fileName = cms.untracked.string('file:simDoublets_TEST.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('TEST')
    )
)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:simDoublets_TEST_DQMIO.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)


process.endjob_step = cms.EndPath(process.endOfProcess)
process.MyTestoutput_step = cms.EndPath(process.MyTestoutput)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)


process.schedule = cms.Schedule(
      process.simDoubletProduction,process.endjob_step,process.MyTestoutput_step, process.DQMoutput_step
)

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(1),
    wantSummary = cms.untracked.bool(True)
)