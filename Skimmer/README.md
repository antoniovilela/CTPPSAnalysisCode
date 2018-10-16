# Missing Mass Searches Skimmer
Skimmer code includes CTPPS pixels, strips and timing containers. Furthermore, it includes jets, leptons, particle flow collections used for Missing Mass Searches.
Instructions how to produce your ntuples:

## Running Locally

```sh
cmsrel CMSSW_9_4_9_cand2
cd CMSSW_9_4_9_cand2/src
cmsenv
cd CTPPSAnalysisCode/Skimmer/test
cmsRun RunMissingMassSearches.py Lepton=Muon (or Electron)
```

## Running on Grid

```sh
cmsrel CMSSW_9_4_9_cand2
cd CMSSW_9_4_9_cand2/src
cmsenv
cd CTPPSAnalysisCode/Skimmer/test
python GridSubmission.py 
```

This script submits crab3 jobs for all the samples (muons, electrons) for the eras B, C, D, E and F.
