import FWCore.ParameterSet.Config as cms

hltPhase2PixelTracksSoA = cms.EDProducer('CAHitNtupletAlpakaPhase2@alpaka',
    pixelRecHitSrc = cms.InputTag('hltPhase2SiPixelRecHitsSoA'),
    CPE = cms.string('PixelCPEFastParamsPhase2'),
    ptmin = cms.double(0.9),
    CAThetaCutBarrel = cms.double(0.002),
    CAThetaCutForward = cms.double(0.003),
    hardCurvCut = cms.double(0.0328407225),
    dcaCutInnerTriplet = cms.double(0.15),
    dcaCutOuterTriplet = cms.double(0.25),
    earlyFishbone = cms.bool(True),
    lateFishbone = cms.bool(False),
    fillStatistics = cms.bool(False),
    minHitsPerNtuplet = cms.uint32(4),
    cellMinz = cms.vdouble(
        -16, 4, -22, -17, 6,-22, -18, 11, -22, 23, 30, 39, 50,
        65, 82, 109, -28, -35, -44, -55, -70, -87, -113, -16, 
        7, -22, 11,-22, -17, 9,-22, 13, -22, 137, 173, 199, 229,
        -142, -177, -203, -233, 23, 30, 39, 50, 65, 82, 109, -28,
        -35, -44, -55, -70, -87, -113
    ),
    cellMaxz = cms.vdouble(
        17, 22, -4, 17, 22, -6, 18, 22, -11, 28, 35, 44, 55, 70,
        87, 113, -23, -30, -39, -50, -65, -82, -109, 17, 22, -7,
        22, -10, 17, 22, -9, 22, -13, 142, 177, 203, 233, -137,
        -173, -199, -229, 28, 35, 44, 55, 70, 87, 113, -23, -30,
        -39, -50, -65, -82, -109
    ),
    cellMaxr = cms.vdouble(
        5.0, 5.0, 5.0, 7.0, 8.0, 8.0,  7.0, 7.0, 7.0, 6.0, 6.0, 6.0, 6.0, 5.0,
        6.0, 5.0, 6.0, 6.0, 6.0, 6.0,  5.0, 6.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
        5.0, 8.0, 8.0, 8.0, 8.0, 6.0,  5.0, 5.0, 5.0, 6.0, 5.0, 5.0, 5.0, 9.0,
        9.0, 9.0, 8.0, 8.0, 8.0, 11.0, 9.0, 9.0, 9.0, 8.0, 8.0, 8.0, 11.0
    ),
    cellMinYSizeB1 = cms.int32(25),
    cellMinYSizeB2 = cms.int32(15),
    cellZ0Cut = cms.double(7.5),
    phiCuts = cms.vint32(
        522, 522, 522, 626, 730, 730, 626, 730, 730, 522, 522,
        522, 522, 522, 522, 522, 522, 522, 522, 522, 522, 522,
        522, 522, 522, 522, 522, 522, 522, 730, 730, 730, 730,
        730, 730, 730, 730, 730, 730, 730, 730, 730, 730, 730,
        730, 730, 730, 522, 522, 522, 522, 522, 522, 522, 522
    ),
    cellMaxDYSize12 = cms.int32(12),
    cellMaxDYSize = cms.int32(10),
    cellMaxDYPred = cms.int32(20),
    cellPtCut = cms.double(0.85),
    maxNumberOfDoublets = cms.uint32(5*512*1024),
    minHitsForSharingCut = cms.uint32(10),
    fitNas4 = cms.bool(False),
    doClusterCut = cms.bool(True),
    doZ0Cut = cms.bool(True),
    doPtCut = cms.bool(True),
    useRiemannFit = cms.bool(False),
    doSharedHitCut = cms.bool(True),
    dupPassThrough = cms.bool(False),
    useSimpleTripletCleaner = cms.bool(True),
    idealConditions = cms.bool(False),
    includeJumpingForwardDoublets = cms.bool(True),
    trackQualityCuts = cms.PSet(
        maxChi2 = cms.double(5.0),
        minPt   = cms.double(0.9),
        maxTip  = cms.double(0.3),
        maxZip  = cms.double(12.),
    ),
    # autoselect the alpaka backend
    alpaka = cms.untracked.PSet(backend = cms.untracked.string(''))
)

_hltPhase2PixelTracksSoASingleIterPatatrack = hltPhase2PixelTracksSoA.clone( minHitsPerNtuplet = 3 )

from Configuration.ProcessModifiers.singleIterPatatrack_cff import singleIterPatatrack
singleIterPatatrack.toReplaceWith(hltPhase2PixelTracksSoA, _hltPhase2PixelTracksSoASingleIterPatatrack)
