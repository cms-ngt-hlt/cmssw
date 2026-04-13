import FWCore.ParameterSet.Config as cms

from Validation.RecoTau.TauValidationRECO import TauValidationRECO as _TauValidationRECO

hltTauValidation = _TauValidationRECO(
    recoTauCollection = "hltHpsPFTauProducer",
    genTauCollection = "tauGenJets",
    isHLT = True
)