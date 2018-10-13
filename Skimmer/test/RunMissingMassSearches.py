#flake8: noqa

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# Setting Input Parameters from Line Command
options = VarParsing ('analysis')
options.register('Lepton','Muon',VarParsing.multiplicity.singleton, VarParsing.varType.string,"Option to Run: Muon or Electron")
options.parseArguments()

process = cms.Process("MissMass")

if options.Lepton == "Muon":
  print("")
  print("#########")
  print("Data Muon")
  print("#########")
  print("")
  triggerlist = 'HLT_IsoMu27_v*', 'HLT_DoubleMu43NoFiltersNoVtx_v*', 'HLT_DoubleMu48NoFiltersNoVtx_v*'
  testfile = '/store/data/Run2017B/DoubleMuon/MINIAOD/31Mar2018-v1/90000/20A919B7-1737-E811-B9D6-1866DAEB296C.root'
elif options.Lepton == "Electron":
  print("")
  print("#############")
  print("Data Electron")
  print("#############")
  print("")
  triggerlist = 'HLT_DoubleEle33_CaloIdL_MW_v*','HLT_Ele27_WPTight_Gsf_v*'
  testfile = '/store/data/Run2017B/DoubleEG/MINIAOD/31Mar2018-v1/00000/00AC968A-8437-E811-A3AD-90B11CBCFF8F.root'
else:
  print("")
  print("##########################################")
  print("Please, try lepton=Muon or lepton=Electron")
  print("##########################################")
  print("")

#########################
#    General options    #
#########################

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options   = cms.untracked.PSet(
    #wantSummary = cms.untracked.bool(True),
    #SkipEvent = cms.untracked.vstring('ProductNotFound'),
    allowUnscheduled = cms.untracked.bool(True),
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#########################
#      Input files      #
#########################

process.source = cms.Source("PoolSource",
				fileNames = cms.untracked.vstring(testfile),
				#firstEvent = cms.untracked.uint32(0)
)

#########################
#        Triggers       #
#########################
process.load("CTPPSAnalysisCode.Skimmer.HLTFilter_cfi")
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilter.HLTPaths = (triggerlist)

#########################
#      Preskimming      #
#########################
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

#########################
#      Electron ID      #
#########################

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runVID=False, #saves CPU time by not needlessly re-running VID
                       era='2017-Nov17ReReco')  

#########################
#     Proton Filter     #
#########################

process.protontagFilter = cms.EDFilter("ProtonTagFilter",
					debugging = cms.bool(False),
					protonTag = cms.InputTag('ctppsLocalTrackLiteProducer'),
					pixelsppsTag = cms.InputTag('ctppsPixelLocalTracks'),
					timingppsTag = cms.InputTag('ctppsDiamondLocalTracks'),
					stripstotemTag = cms.InputTag('totemRPLocalTrackFitter')
					)

#########################
#       Analysis        #
#########################

process.load("CTPPSAnalysisCode.Skimmer.countsAnalyzer_cfi")
process.CounterBeforePPSTagging = process.countsAnalyzer.clone()
process.CounterAfterPPSTagging = process.countsAnalyzer.clone()

process.load("CTPPSAnalysisCode.Skimmer.MissingMassSearches_cfi")
process.missing_mass.mode = cms.string(options.Lepton) #Muons or Electrons
process.missing_mass.debugging = cms.bool(False)
process.missing_mass.includeMuons = cms.bool(True)
process.missing_mass.includeElectrons = cms.bool(True)
process.missing_mass.includeJets = cms.bool(True)
process.missing_mass.includePhotons = cms.bool(True)
process.missing_mass.includeMET = cms.bool(True)
process.missing_mass.includeVertices = cms.bool(True)
process.missing_mass.includeParticleFlow = cms.bool(True)
process.missing_mass.includeProtons = cms.bool(True)
process.missing_mass.triggersList = process.hltFilter.HLTPaths
process.missing_mass.JetAlgoA = cms.InputTag('slimmedJets')
process.missing_mass.JetAlgoB = cms.InputTag('slimmedJetsAK8')
process.missing_mass.electronTag = cms.InputTag("slimmedElectrons")
process.missing_mass.muonTag = cms.InputTag("slimmedMuons")
process.missing_mass.pfTag = cms.InputTag('particleFlow')
process.missing_mass.packedTag = cms.InputTag('packedPFCandidates')
process.missing_mass.photonTag = cms.InputTag('slimmedPhotons')
process.missing_mass.metTag = cms.InputTag('slimmedMETs')
process.missing_mass.vertexTag = cms.InputTag('offlineSlimmedPrimaryVertices')
process.missing_mass.protonTag = cms.InputTag('ctppsLocalTrackLiteProducer')
process.missing_mass.pixelsppsTag = cms.InputTag('ctppsPixelLocalTracks')
process.missing_mass.timingppsTag = cms.InputTag('ctppsDiamondLocalTracks')
process.missing_mass.stripstotemTag = cms.InputTag('totemRPLocalTrackFitter')

# E/gamma identification
process.missing_mass.eleIdLabels = cms.PSet(
    mediumLabel = cms.InputTag('mvaEleID-Spring16-GeneralPurpose-V1-wp90'),
    tightLabel = cms.InputTag('mvaEleID-Spring16-GeneralPurpose-V1-wp80'),
)
process.missing_mass.phoIdLabels = cms.PSet(
    mediumLabel = cms.InputTag('mvaPhoID-Spring16-nonTrig-V1-wp90'),
    tightLabel = cms.InputTag('mvaPhoID-Spring16-nonTrig-V1-wp80'),
)

# prepare the output file
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('output.root'),
    closeFileFast = cms.untracked.bool(True)
)

process.p = cms.Path(
    process.hltFilter*
    process.CounterBeforePPSTagging*
    process.protontagFilter*
    process.CounterAfterPPSTagging*
    process.egammaPostRecoSeq*
    process.missing_mass
)
