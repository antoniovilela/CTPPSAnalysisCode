# Missing Mass Searches Skimmer
Skimmer code which includes CTPPS pixels, strips and timing containers. Furthermore, it includes jets, leptons, particle flow collections used for Missing Mass Searches.
Instructions how to install:

```sh
cmsrel CMSSW_9_4_9_cand2
cd CMSSW_9_4_9_cand2/src
cmsenv
git cms-init
git cms-merge-topic cms-egamma:EgammaPostRecoTools_940
git clone https://github.com/dfigueiredo/CTPPSAnalysisCode.git
cd CTPPSAnalysisCode/
git checkout MissingMassSearches_v01
cd CMSSW_9_4_9_cand2/src
scram b -j 8
```
