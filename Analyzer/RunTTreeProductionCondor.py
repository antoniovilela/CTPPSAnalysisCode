#!/usr/bin/env python
import os, re
import commands
import math, time
import sys

from os import listdir
from os.path import isfile, join

# Creating Output Folder
path = os.getcwd()

# Configure here the folders where your jobs are located. Please, fill era (B,C,D,E and F) and mode (Muon or Electron) correctly!
# format ['Crab output folder', 'Era', 'Mode']
mytask= [['/eos/cms/store/group/phys_pps/MissingMassSearch/MuonB-MiniAOD-v3/DoubleMuon/MuonB-MiniAOD-v3/181015_165808/0000/','B','Muon','OutputMuonB'],
	['/eos/cms/store/group/phys_pps/MissingMassSearch/MuonC-MiniAOD-v3/DoubleMuon/MuonC-MiniAOD-v3/181015_165821/0000/','C','Muon','OutputMuonC'],
	['/eos/cms/store/group/phys_pps/MissingMassSearch/MuonD-MiniAOD-v3/DoubleMuon/MuonD-MiniAOD-v3/181015_165835/0000/','D','Muon','OutputMuonD'],
	['/eos/cms/store/group/phys_pps/MissingMassSearch/MuonE-MiniAOD-v3/DoubleMuon/MuonE-MiniAOD-v3/181015_165848/0000/','E','Muon','OutputMuonE'],
	['/eos/cms/store/group/phys_pps/MissingMassSearch/MuonF-MiniAOD-v3/DoubleMuon/MuonF-MiniAOD-v3/181015_165902/0000/','F','Muon','OutputMuonF'],
	['/eos/cms/store/group/phys_pps/MissingMassSearch/ElectronB-MiniAOD-v3/DoubleEG/ElectronB-MiniAOD-v3/181015_165916/0000/','B','Electron','OutputElectronB'],
	['/eos/cms/store/group/phys_pps/MissingMassSearch/ElectronC-MiniAOD-v3/DoubleEG/ElectronC-MiniAOD-v3/181015_165933/0000/','C','Electron','OutputElectronC'],
	['/eos/cms/store/group/phys_pps/MissingMassSearch/ElectronD-MiniAOD-v3/DoubleEG/ElectronD-MiniAOD-v3/181015_165949/0000/','D','Electron','OutputElectronD'],
	['/eos/cms/store/group/phys_pps/MissingMassSearch/ElectronE-MiniAOD-v3/DoubleEG/ElectronE-MiniAOD-v3/181015_170002/0000/','E','Electron','OutputElectronE'],
	['/eos/cms/store/group/phys_pps/MissingMassSearch/ElectronF-MiniAOD-v3/DoubleEG/ElectronF-MiniAOD-v3/181015_170020/0000/','F','Electron','OutputElectronF']]

# Creating list of commands for all files inside each directory defined under mytask.
# i.e: ./MissingMassNtupleAnalyzer --f /eos/cms/store/group/phys_pps/MissingMassSearch/MuonF-MiniAOD-v3/DoubleMuon/MuonF-MiniAOD-v3/181015_165902/0000/output_99.root --era F --mode Muon --jobid 893
i = 0
jobid = 0
parameters = []
element = []
while i < len(mytask):
	onlyfiles = [f for f in listdir(mytask[i][0]) if isfile(join(mytask[i][0], f))]
	j = 0
        
	# Creating Output Folder
	folderout = str(path)+"/"+mytask[i][3]

	if not os.path.exists(folderout):
		os.makedirs(folderout)
		os.system("cp *.csv "+folderout+"/.")
	else:
		os.system("rm -rf "+folderout)

	while j < len(onlyfiles):
		# path + file name, era, mode, jobid, output folder
		element = [mytask[i][0]+onlyfiles[j], mytask[i][1], mytask[i][2], str(jobid), mytask[i][3]]
                parameters.append(element)
		jobid += 1
		j += 1
	i += 1

# Sending Jobs @ Condor! Party is just in the beginning!
print '\nSending Condor jobs to produce NTuples for Missing Mass CTPPS Analysis\n'

i = 0
while i < len(parameters):
        with open('job_condor_tmp.sub', 'w') as fout:
                fout.write("initialdir = "+parameters[i][4]+"\n")
                fout.write("executable = ./MissingMassNtupleAnalyzer\n")
                fout.write("arguments = --f "+parameters[i][0]+" --era "+parameters[i][1]+" --mode "+parameters[i][2]+" --jobid "+parameters[i][3]+"\n")
                fout.write("transfer_input_files = ../xangle_afterTS2_STABLEBEAMS_CLEANUP.csv, ../xangle_tillTS2_STABLEBEAMS_CLEANUP.csv\n")
                fout.write("output = missing_mass_"+parameters[i][3]+".out\n")
                fout.write("error = missing_mass_"+parameters[i][3]+".err\n")
                fout.write("log = missing_mass_"+parameters[i][3]+".log\n")
                fout.write("getenv = True\n")

		###########################
		# espresso     = 20 minutes
		# microcentury = 1 hour
		# longlunch    = 2 hours
		# workday      = 8 hours
		# tomorrow     = 1 day
		# testmatch    = 3 days
		# nextweek     = 1 week
		##########################

                fout.write("+JobFlavour = microcentury\n")
                fout.write("queue\n")
	os.system('condor_submit -spool job_condor_tmp.sub')
	os.system('rm job_condor_tmp.sub')
        i += 1

print 'END\n'
