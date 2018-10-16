# Missing Mass Searches Analyzer
Tool which produces the final samples used for Missing Mass Search analysis.

### Compilation

```sh
cmsrel CMSSW_9_4_9_cand2
cd CMSSW_9_4_9_cand2/src
cmsenv
cd CTPPSAnalysisCode/Analyzer
make
```

### Running Locally

You can use test the following command:

```sh
./MissingMassNtupleAnalyzer --f testfile_era_b_mode_muons.root --era B --mode Muon --jobid 1
```

Where the file 'testfile_era_b_mode_muons.root' has been produced by the Skimmer.

### Running on Lxbatch

In order to produce your own sample for the final physics analysis, you can use a script which automatically produces lxbatch jobs per each crab output file (produced by the Skimmer).

```sh
python RunTTreeProductionLxbatch.py
```

The output files will be saved into 'OutputJobs' folder created by the script.
