import FWCore.ParameterSet.Config as cms

from Validation.RecoTau.TauValidation import TauValidation as _TauValidation

hltTauValidation = _TauValidation(
    recoTauCollection = "hltHpsPFTauProducer",
    genTauCollection = "tauGenJets",
    isHLT = True
)