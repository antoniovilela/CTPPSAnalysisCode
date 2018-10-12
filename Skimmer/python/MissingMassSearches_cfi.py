import FWCore.ParameterSet.Config as cms

missing_aod = cms.EDAnalyzer('MissingMassSearches',
    triggerResults = cms.InputTag('TriggerResults', '', 'HLT'),
)
