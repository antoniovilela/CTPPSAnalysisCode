import FWCore.ParameterSet.Config as cms

missing_mass = cms.EDAnalyzer('MissingMassSearches',
    triggerResults = cms.InputTag('TriggerResults', '', 'HLT'),
)
