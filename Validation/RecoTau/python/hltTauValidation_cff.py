import FWCore.ParameterSet.Config as cms

from Validation.RecoTau.TauValidationRECO import TauValidationRECO as _TauValidationRECO

hltTauValidation = _TauValidationRECO(
    recoTauCollection = "hltHpsPFTauProducer",
    genTauCollection = "tauGenJetsSelectorAllHadrons", # only GenTaus decaying hadronically
    minDeltaR = 0.3,
    outFolder = "HLT/Tau/TauValidation",
    isHLT = True
)

hltTauValidation_deltaR0p3 = hltTauValidation.clone(
    minDeltaR = 0.3,
    outFolder = "HLT/Tau/TauValidation_DeltaR/DeltaR0p3",
)

hltTauValidation_deltaR0p25 = hltTauValidation.clone(
    minDeltaR = 0.25,
    outFolder = "HLT/Tau/TauValidation_DeltaR/DeltaR0p25",
)

hltTauValidation_deltaR0p2 = hltTauValidation.clone(
    minDeltaR = 0.2,
    outFolder = "HLT/Tau/TauValidation_DeltaR/DeltaR0p2",
)

hltTauValidation_deltaR0p15 = hltTauValidation.clone(
    minDeltaR = 0.15,
    outFolder = "HLT/Tau/TauValidation_DeltaR/DeltaR0p15",
)

hltTauValidation_deltaR0p1 = hltTauValidation.clone(
    minDeltaR = 0.1,
    outFolder = "HLT/Tau/TauValidation_DeltaR/DeltaR0p1",
)

hltTauValidation_deltaR = cms.Sequence(
    hltTauValidation_deltaR0p3 +
    hltTauValidation_deltaR0p25 +
    hltTauValidation_deltaR0p2 +
    hltTauValidation_deltaR0p15 +
    hltTauValidation_deltaR0p1
)

hltTauValidationSequence = cms.Sequence(
    hltTauValidation +
    # only for testing DeltaR matching
    hltTauValidation_deltaR
)