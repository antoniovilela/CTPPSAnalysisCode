# Missing Mass Searches Tools

Such set of software tools have been made to be used for Missing Search Analysis. Basically, it consists in two steps: skimmer and Analyzer. 
The skimmer is used to read the CMS MiniAOD samples and it produces files which will be needed for the Analyzer. Finally, the outputs of the Analyzer
are the final files used for Physics. Detailed instructions for each step can be found into the tool folder ("Skimmer" or "Analyzer").

## Instalation

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
