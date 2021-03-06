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

The output files will be saved into 'Output' folder created by the script, for each sample.

#### Lxbatch Commands 

```sh
LSB_JOB_REPORT_MAIL=N bsub -q 8nm -J job_id_0 -e jobs_error -o jobs_output < job.sh (Submit)
bkill 0 -u $USER (Kill all jobs)
```

**More information here:** [https://twiki.cern.ch/twiki/bin/view/Main/BatchJobs](https://twiki.cern.ch/twiki/bin/view/Main/BatchJobs)

### Running on HTC Condor

In order to produce your own sample for the final physics analysis, you can use a script which automatically produces condor jobs per each crab output file (produced by the Skimmer). @CERN, HTC condor is very recommended. 

```sh
python RunTTreeProductionCondor.py
```

#### Condor Commands

```sh
condor_submit job_condor.sub (Submit with automatic output transfer)
watch condor_q (Check job)
condor_q (Check job)
condor_submit -spool job_condor.sub (Submit, but not automatic output transfer)
condor_transfer_data $LOGNAME -const 'JobStatus == 4' (Retrieve data when submitted with -spool option)
condor_rm -all (kill all jobs)
```

**More information here:** [https://twiki.cern.ch/twiki/bin/view/ABPComputing/LxbatchHTCondor](https://twiki.cern.ch/twiki/bin/view/ABPComputing/LxbatchHTCondor)

**Gui for Monitoring:** [https://monit-grafana.cern.ch/d/000000869/user-batch-jobs?orgId=5&refresh=5m&var-cluster=cernprod](https://monit-grafana.cern.ch/d/000000869/user-batch-jobs?orgId=5&refresh=5m&var-cluster=cernprod)

