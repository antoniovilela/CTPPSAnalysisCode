#!/usr/bin/env python

##############
# Multi CRAB #
##############

import time 

if __name__ == '__main__':
  from CRABAPI.RawCommand import crabCommand

def submit(config):
  res = crabCommand('submit', config = config)

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

dataset = {
  'MuonB-MINIAOD' : '/DoubleMuon/Run2017B-31Mar2018-v1/MINIAOD',
  'MuonC-MINIAOD' : '/DoubleMuon/Run2017C-31Mar2018-v1/MINIAOD',
  'MuonD-MINIAOD' : '/DoubleMuon/Run2017D-31Mar2018-v1/MINIAOD',
  'MuonE-MINIAOD' : '/DoubleMuon/Run2017E-31Mar2018-v1/MINIAOD',
  'MuonF-MINIAOD' : '/DoubleMuon/Run2017F-31Mar2018-v1/MINIAOD',
  'ElectronB-MINIAOD' : '/DoubleEG/Run2017B-31Mar2018-v1/MINIAOD',
  'ElectronC-MINIAOD' : '/DoubleEG/Run2017C-31Mar2018-v1/MINIAOD',
  'ElectronD-MINIAOD' : '/DoubleEG/Run2017D-31Mar2018-v1/MINIAOD',
  'ElectronE-MINIAOD' : '/DoubleEG/Run2017E-31Mar2018-v1/MINIAOD',
  'ElectronF-MINIAOD' : '/DoubleEG/Run2017F-31Mar2018-v1/MINIAOD',
}

mode = "" 
filesPerJob = 1

config.General.transferLogs = True
config.General.transferOutputs = True
config.JobType.pluginName = 'Analysis'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.publication = False
config.Site.storageSite = 'T2_CH_CERN'
config.JobType.psetName = 'RunGammaGammaLeptonLepton_cfg.py'
#config.Data.ignoreLocality = True 
#config.JobType.disableAutomaticOutputCollection = False

def doSubmit(listOfSamples):
  for sample in listOfSamples:
    config.JobType.outputFiles = ['output.root']
    config.General.workArea = 'crab_'+ sample
    config.General.requestName = sample
    config.Data.inputDataset = dataset[sample]
    config.Data.unitsPerJob = filesPerJob
    config.Data.lumiMask = 'combined_RPIN_CMS_2017.json'
    config.Data.outputDatasetTag = sample
    config.Data.outLFNDirBase = '/store/group/phys_pps/MissingMassSearch/' + sample
    config.Site.storageSite = 'T2_CH_CERN'
    config.JobType.pyCfgParams = ["Lepton="+mode]
    submit(config)

# MuonB
mode = "Muon"
filesPerJob = 200
listOfSamples = ['MuonB-MINIAOD']
doSubmit(listOfSamples)

# MuonC
mode = "Muon"
filesPerJob = 200
listOfSamples = ['MuonC-MINIAOD']
doSubmit(listOfSamples)

# MuonD
mode = "Muon"
filesPerJob = 200
listOfSamples = ['MuonD-MINIAOD']
doSubmit(listOfSamples)

# MuonE
mode = "Muon"
filesPerJob = 200
listOfSamples = ['MuonE-MINIAOD']
doSubmit(listOfSamples)

# MuonF
mode = "Muon"
filesPerJob = 200
listOfSamples = ['MuonF-MINIAOD']
doSubmit(listOfSamples)

# ElectronB
mode = "Electron"
filesPerJob = 200
listOfSamples = ['ElectronB-MINIAOD']
doSubmit(listOfSamples)

# ElectronC
mode = "Electron"
filesPerJob = 200
listOfSamples = ['ElectronC-MINIAOD']
doSubmit(listOfSamples)

# ElectronD
mode = "Electron"
filesPerJob = 200
listOfSamples = ['ElectronD-MINIAOD']
doSubmit(listOfSamples)

# ElectronE
mode = "Electron"
filesPerJob = 200
listOfSamples = ['ElectronE-MINIAOD']
doSubmit(listOfSamples)

# ElectronF
mode = "Electron"
filesPerJob = 200
listOfSamples = ['ElectronF-MINIAOD']
doSubmit(listOfSamples)
