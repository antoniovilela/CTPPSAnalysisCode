#!/usr/bin/env python
import os, re
import commands
import math, time
import sys

from os import listdir
from os.path import isfile, join

queue = "1nh" # give bsub queue -- 8nm (8 minutes), 1nh (1 hour), 8nh, 1nd (1day), 2nd, 1nw (1 week), 2nw 

# Creating Output Folder
path = os.getcwd()
folderout = str(path)+"/OutputJobs"

if not os.path.exists(folderout):
	os.makedirs(folderout)

# Configure here the folders where your jobs are located. Please, fill era (B,C,D,E and F) and mode (Muon or Electron) correctly!
# format ['Crab output folder', 'Era', 'Mode']
mytask=[['/eos/cms/store/group/phys_pps/MissingMassSearch/MuonB-MiniAOD-v3/DoubleMuon/MuonB-MiniAOD-v3/181015_165808/0000/','B','Muon'],
	['/eos/cms/store/group/phys_pps/MissingMassSearch/MuonC-MiniAOD-v3/DoubleMuon/MuonC-MiniAOD-v3/181015_165821/0000/','C','Muon'],
	['/eos/cms/store/group/phys_pps/MissingMassSearch/MuonD-MiniAOD-v3/DoubleMuon/MuonD-MiniAOD-v3/181015_165835/0000/','D','Muon'],
	['/eos/cms/store/group/phys_pps/MissingMassSearch/MuonE-MiniAOD-v3/DoubleMuon/MuonE-MiniAOD-v3/181015_165848/0000/','E','Muon'],
	['/eos/cms/store/group/phys_pps/MissingMassSearch/MuonF-MiniAOD-v3/DoubleMuon/MuonF-MiniAOD-v3/181015_165902/0000/','F','Muon'],
	['/eos/cms/store/group/phys_pps/MissingMassSearch/ElectronB-MiniAOD-v3/DoubleEG/ElectronB-MiniAOD-v3/181015_165916/0000/','B','Electron'],
	['/eos/cms/store/group/phys_pps/MissingMassSearch/ElectronC-MiniAOD-v3/DoubleEG/ElectronC-MiniAOD-v3/181015_165933/0000/','C','Electron'],
	['/eos/cms/store/group/phys_pps/MissingMassSearch/ElectronD-MiniAOD-v3/DoubleEG/ElectronD-MiniAOD-v3/181015_165949/0000/','D','Electron'],
	['/eos/cms/store/group/phys_pps/MissingMassSearch/ElectronE-MiniAOD-v3/DoubleEG/ElectronE-MiniAOD-v3/181015_170002/0000/','E','Electron'],
	['/eos/cms/store/group/phys_pps/MissingMassSearch/ElectronF-MiniAOD-v3/DoubleEG/ElectronF-MiniAOD-v3/181015_170020/0000/','F','Electron']]

# Creating list of commands for all files inside each directory defined under mytask.
# i.e: ./MissingMassNtupleAnalyzer --f /eos/cms/store/group/phys_pps/MissingMassSearch/MuonF-MiniAOD-v3/DoubleMuon/MuonF-MiniAOD-v3/181015_165902/0000/output_99.root --era F --mode Muon --jobid 893

i = 0
jobid = 0
command = []
while i < len(mytask):
	onlyfiles = [f for f in listdir(mytask[i][0]) if isfile(join(mytask[i][0], f))]
	j = 0
	while j < len(onlyfiles):
		command.append("./MissingMassNtupleAnalyzer --f " + mytask[i][0] + onlyfiles[j] + " --era " + mytask[i][1]+" --mode "+mytask[i][2]+" --jobid id"+str(jobid))
		jobid += 1
		j += 1
	i += 1

# Sending Jobs @ lxbatch! Party is just in the beginning!
print '\nSending Lxbatch jobs to produce NTuples for Missing Mass CTPPS Analysis\n\n'

path = os.getcwd()
print path

i = 0
while i < len(command):
        with open('job.sh', 'w') as fout:
                fout.write("#!/bin/sh\n")
                fout.write("echo\n")
                fout.write("echo 'START---------------'\n")
                fout.write("cd "+str(path)+"\n")
                fout.write(command[i]+"\n")
                fout.write("cp *.root "+ str(folderout) +"/.\n")
                fout.write("echo 'STOP---------------'\n")
                fout.write("echo\n")
        os.system("chmod 755 job.sh")
        os.system("bsub -q "+queue+" -J job_"+str(i)+" < job.sh")
        print "job nr " + str(i) + " submitted"
        i += 1

print
print "your jobs:"
os.system("bjobs")
print
print 'END'
print
