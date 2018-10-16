#define MissingMassNtupleAnalyzer_cxx

// ROOT header
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

// Header file for the classes stored in the TTree if any.
#include <iostream>
#include <fstream>
#include <algorithm>

// User Header
#include "statusbar.h"
#include "MissingMassNtupleAnalyzer.h"

// Code Constants
#define MASS_B 4.2 // GeV
#define MASS_MU 0.1057 // GeV
#define MASS_E  0.000511 // GeV
#define MASS_P  0.938272029 // GeV
#define pi 3.14159265359
#define ECM 13000.0 // GeV
#define ns_to_s_ 1e-9;
#define m_to_cm_ 1e2;

void MissingMassNtupleAnalyzer::Loop(char * era, char * mode, bool preTS2, char * jobid)
{

  // debugging variables (printout values)
  bool debug = false;
  TString filenameout;

  // saving TTree
  TFile* fileout;
  TTree* tout;

  if(jobid && jobid!=NULL){
    filenameout = "Ntuple_data_mode_"+std::string(mode)+"_"+std::string(jobid)+".root";
  }else{
    filenameout = "Ntuple_data_mode_"+std::string(mode)+".root";
  }
  fileout = new TFile(filenameout, "RECREATE");

  std::cout << "\t = Options =" << std::endl;
  std::cout << "\t\t Mode: " << mode << std::endl;
  std::cout << "\t\t Era: " << era << std::endl;
  std::cout << "\t\t preTS2: " << preTS2 << std::endl;
  if(jobid!=0) std::cout << "\t\t jobid: " << jobid << std::endl;
  std::cout << "\t\t output: " << filenameout << "\n" << std::endl;

  fileout->cd();

  int run_;
  int event_;
  int lumiblock_;

  bool triggerSingle = false;
  bool triggerDouble = false;

  int nLeptons = 0;
  float leadingLeptonEnergy = 0.;
  float leadingLeptonPt = 0.;
  float leadingLeptonEta = 0.;
  float leadingLeptonPhi = 0.;
  int leadingLeptonCharge = 0.;
  float leadingLeptonVz = 0.;
  bool leadingLeptonLooseId = false;
  bool leadingLeptonMediumId = false;
  bool leadingLeptonTightId = false;
  bool leadingLeptonPfIsoMedium = false;
  bool leadingLeptonMiniIsoTight = false;
  bool leadingLeptonPfIsoVeryTight = false;

  float secondLeptonEnergy = 0.;
  float secondLeptonPt = 0.;
  float secondLeptonEta = 0.;
  float secondLeptonPhi = 0.;
  int secondLeptonCharge = 0.;
  float secondLeptonVz = 0.;
  bool secondLeptonLooseId = false;
  bool secondLeptonMediumId = false;
  bool secondLeptonTightId = false;
  bool secondLeptonPfIsoMedium = false;
  bool secondLeptonMiniIsoTight = false;
  bool secondLeptonPfIsoVeryTight = false;

  int nJets = 0;
  float leadingJetEnergy = 0.;
  float leadingJetPt = 0.;
  float leadingJetEta = 0.;
  float leadingJetPhi = 0.;
  float leadingJetVz = 0.;
  float leadingJetBtag = 0.;
  float secondJetEnergy = 0.;
  float secondJetPt = 0.;
  float secondJetEta = 0.;
  float secondJetPhi = 0.;
  float secondJetVz = 0.;
  float secondJetBtag = 0.;

  int dileptonCharge = 0.;
  float dileptonMass = 0.;
  float dileptonEta = 0.;
  float dileptonPhi = 0.;
  float dileptonPt = 0.;
  float dileptonRapidity = 0.;

  float dijetMass = 0.;
  float dijetEta = 0.;
  float dijetPhi = 0.;
  float dijetPt = 0.;
  float dijetRapidity = 0.;

  float missingMass = 0.;
  float missingEta = 0.;
  float missingPhi = 0.;
  float missingPt = 0.;
  float missingRapidity = 0.;

  float missEt_ = 0.;
  float missEt_phi_ = 0.;
  int nPF_ = 0;
  float sumEnergyPF_ = 0;
  float dphi_ = 0;

  tout = new TTree("Events", "Events");
  tout->Branch("run",&run_,"run/I");
  tout->Branch("event",&event_,"event/I");
  tout->Branch("lumiblock",&lumiblock,"lumiblock/I");
  tout->Branch("triggerSingle",&triggerSingle,"triggerSingle/B");
  tout->Branch("triggerDouble",&triggerDouble,"triggerDouble/B");
  tout->Branch("nLeptons",&nLeptons,"nLeptons/I");
  tout->Branch("leadingLeptonEnergy",&leadingLeptonEnergy,"leadingLeptonEnergy/F");
  tout->Branch("leadingLeptonPt",&leadingLeptonPt,"leadingLeptonPt/F");
  tout->Branch("leadingLeptonEta",&leadingLeptonEta,"leadingLeptonEta/F");
  tout->Branch("leadingLeptonPhi",&leadingLeptonPhi,"leadingLeptonPhi/F");
  tout->Branch("leadingLeptonCharge",&leadingLeptonCharge,"leadingLeptonCharge/I");
  tout->Branch("leadingLeptonVz",&leadingLeptonVz,"leadingLeptonVz/F");
  tout->Branch("leadingLeptonLooseId",&leadingLeptonLooseId,"leadingLeptonLooseId/B");
  tout->Branch("leadingLeptonMediumId",&leadingLeptonMediumId,"leadingLeptonMediumId/B");
  tout->Branch("leadingLeptonTightId",&leadingLeptonTightId,"leadingLeptonTightId/B");
  tout->Branch("leadingLeptonPfIsoMedium",&leadingLeptonPfIsoMedium,"leadingLeptonPfIsoMedium/B");
  tout->Branch("leadingLeptonMiniIsoTight",&leadingLeptonMiniIsoTight,"leadingLeptonMiniIsoTight/B");
  tout->Branch("leadingLeptonPfIsoVeryTight",&leadingLeptonPfIsoVeryTight,"leadingLeptonPfIsoVeryTight/B");
  tout->Branch("secondLeptonEnergy",&secondLeptonEnergy,"secondLeptonEnergy/F");
  tout->Branch("secondLeptonPt",&secondLeptonPt,"secondLeptonPt/F");
  tout->Branch("secondLeptonEta",&secondLeptonEta,"secondLeptonEta/F");
  tout->Branch("secondLeptonPhi",&secondLeptonPhi,"secondLeptonPhi/F");
  tout->Branch("secondLeptonCharge",&secondLeptonCharge,"secondLeptonCharge/I");
  tout->Branch("secondLeptonVz",&secondLeptonVz,"secondLeptonVz/F");
  tout->Branch("secondLeptonLooseId",&secondLeptonLooseId,"secondLeptonLooseId/B");
  tout->Branch("secondLeptonMediumId",&secondLeptonMediumId,"secondLeptonMediumId/B");
  tout->Branch("secondLeptonTightId",&secondLeptonTightId,"secondLeptonTightId/B");
  tout->Branch("secondLeptonPfIsoMedium",&secondLeptonPfIsoMedium,"secondLeptonPfIsoMedium/B");
  tout->Branch("secondLeptonMiniIsoTight",&secondLeptonMiniIsoTight,"secondLeptonMiniIsoTight/B");
  tout->Branch("secondLeptonPfIsoVeryTight",&secondLeptonPfIsoVeryTight,"secondLeptonPfIsoVeryTight/B");
  tout->Branch("nJets",&nJets,"nJets/I");
  tout->Branch("leadingJetEnergy",&leadingJetEnergy,"leadingJetEnergy/F");
  tout->Branch("leadingJetPt",&leadingJetPt,"leadingJetPt/F");
  tout->Branch("leadingJetEta",&leadingJetEta,"leadingJetEta/F");
  tout->Branch("leadingJetPhi",&leadingJetPhi,"leadingJetPhi/F");
  tout->Branch("leadingJetVz",&leadingJetVz,"leadingJetVz/F");
  tout->Branch("leadingJetBtag",&leadingJetBtag,"leadingJetBtag/F");
  tout->Branch("secondJetEnergy",&secondJetEnergy,"secondJetEnergy/F");
  tout->Branch("secondJetPt",&secondJetPt,"secondJetPt/F");
  tout->Branch("secondJetEta",&secondJetEta,"secondJetEta/F");
  tout->Branch("secondJetPhi",&secondJetPhi,"secondJetPhi/F");
  tout->Branch("secondJetVz",&secondJetVz,"secondJetVz/F");
  tout->Branch("secondJetBtag",&secondJetBtag,"secondJetBtag/F");
  tout->Branch("dileptonCharge",&dileptonCharge,"dileptonCharge/I");
  tout->Branch("dileptonMass",&dileptonMass,"dileptonMass/F");
  tout->Branch("dileptonEta",&dileptonEta,"dileptonEta/F");
  tout->Branch("dileptonPhi",&dileptonPhi,"dileptonPhi/F");
  tout->Branch("dileptonPt",&dileptonPt,"dileptonPt/F");
  tout->Branch("dileptonRapidity",&dileptonRapidity,"dileptonRapidity/F");
  tout->Branch("dijetMass",&dijetMass,"dijetMass/F");
  tout->Branch("dijetEta",&dijetEta,"dijetEta/F");
  tout->Branch("dijetPhi",&dijetPhi,"dijetPhi/F");
  tout->Branch("dijetPt",&dijetPt,"dijetPt/F");
  tout->Branch("dijetRapidity",&dijetRapidity,"dijetRapidity/F");
  tout->Branch("missingMass",&missingMass,"missingMass/F");
  tout->Branch("missingEta",&missingEta,"missingEta/F");
  tout->Branch("missingPhi",&missingPhi,"missingPhi/F");
  tout->Branch("missingPt",&missingPt,"missingPt/F");
  tout->Branch("missingRapidity",&missingRapidity,"missingRapidity/F");
  tout->Branch("missEt",&missEt_,"missEt/F");
  tout->Branch("missEt_phi",&missEt_phi_,"missEt_phi/F");
  tout->Branch("nPF",&nPF_,"nPF/I");
  tout->Branch("sumEnergyPF",&sumEnergyPF_,"sumEnergyPF/F");
  tout->Branch("dphi",&dphi_,"dphi/F");

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    run_ = 0;
    event_ = 0;
    lumiblock_ = 0;

    triggerSingle = false;
    triggerDouble = false;

    nLeptons = 0;
    leadingLeptonEnergy = 0.;
    leadingLeptonPt = 0.;
    leadingLeptonEta = 0.;
    leadingLeptonPhi = 0.;
    leadingLeptonCharge = 0.;
    leadingLeptonVz = 0.;
    leadingLeptonLooseId = false;
    leadingLeptonMediumId = false;
    leadingLeptonTightId = false;
    leadingLeptonPfIsoMedium = false;
    leadingLeptonMiniIsoTight = false;
    leadingLeptonPfIsoVeryTight = false;

    secondLeptonEnergy = 0.;
    secondLeptonPt = 0.;
    secondLeptonEta = 0.;
    secondLeptonPhi = 0.;
    secondLeptonCharge = 0.;
    secondLeptonVz = 0.;
    secondLeptonLooseId = false;
    secondLeptonMediumId = false;
    secondLeptonTightId = false;
    secondLeptonPfIsoMedium = false;
    secondLeptonMiniIsoTight = false;
    secondLeptonPfIsoVeryTight = false;

    nJets = 0;
    leadingJetEnergy = 0.;
    leadingJetPt = 0.;
    leadingJetEta = 0.;
    leadingJetPhi = 0.;
    leadingJetVz = 0.;
    leadingJetBtag = 0.;
    secondJetEnergy = 0.;
    secondJetPt = 0.;
    secondJetEta = 0.;
    secondJetPhi = 0.;
    secondJetVz = 0.;
    secondJetBtag = 0.;

    dileptonCharge = 0.;
    dileptonMass = 0.;
    dileptonEta = 0.;
    dileptonPhi = 0.;
    dileptonPt = 0.;
    dileptonRapidity = 0.;

    dijetMass = 0.;
    dijetEta = 0.;
    dijetPhi = 0.;
    dijetPt = 0.;
    dijetRapidity = 0.;

    missingMass = 0.;
    missingEta = 0.;
    missingPhi = 0.;
    missingPt = 0.;
    missingRapidity = 0.;

    missEt_ = 0.;
    missEt_phi_ = 0.;
    nPF_ = 0;
    sumEnergyPF_ = 0.;
    dphi_ = 0.;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;

    if(!debug) loadBar(jentry, nentries);

    nb = fChain->GetEntry(jentry);   nbytes += nb;

    int protonArm45 = 0;
    int protonArm56 = 0;

    float protonArm45_x = 0;
    float protonArm45_y = 0;

    float protonArm56_x = 0;
    float protonArm56_y = 0;

    for(std::vector<int>::size_type i = 0; i != protonsArm->size(); i++){
      if(protonsArm->at(i)==0&&(protonsStation->at(i)==0||protonsStation->at(i)==2)&&protonsRP->at(i)==3){
	protonArm45++;
	protonArm45_x = protonsX->at(i);
	protonArm45_y = protonsY->at(i);
      }
      if(protonsArm->at(i)==1&&(protonsStation->at(i)==0||protonsStation->at(i)==2)&&protonsRP->at(i)==3){
	protonArm56++;
	protonArm56_x = protonsX->at(i);
	protonArm56_y = protonsY->at(i);
      }
    }

    bool OneProtonPerArm = false;
    if(protonArm45==1&&protonArm56==1) OneProtonPerArm = true;

    if(OneProtonPerArm){

      TLorentzVector p1;
      TLorentzVector p2;

      TLorentzVector lepton1;
      TLorentzVector lepton2;
      TLorentzVector dilepton;

      TLorentzVector jet1;
      TLorentzVector jet2;
      TLorentzVector dijet;

      TLorentzVector X;

      int angle = -1;
      float dispersionL = -1;
      float dispersionR = -1;

      float Dispersion[2][200];
      Dispersion[0][120]=9.145*10;
      Dispersion[0][130]=8.591*10;
      Dispersion[0][140]=8.226*10;
      Dispersion[1][120]=7.291*10;
      Dispersion[1][130]=6.621*10;
      Dispersion[1][140]=6.191*10;

      float relcorr=0;
      preTS2 ? relcorr=(42.05):relcorr=(42.2);

      float xi1, xi2 = -1;
      float diffmass = -1;
      float proton1_pz, proton2_pz = -1;

      if(preTS2){
	angle=getXangle(run, lumiblock, "./inp/new/xangle_tillTS2_STABLEBEAMS_CLEANUP.csv");
      }else{
	angle=getXangle(run, lumiblock, "./inp/xangle_afterTS2_cleanup.csv");
      }

      if(angle !=120 && angle != 130 && angle != 140){
	dispersionL = Dispersion[0][140] - 4.6*((angle-140.)/10.) ;
	dispersionR = Dispersion[1][140] - 4.6*((angle-140.)/10.) ;
      }

      if(angle == 120 || angle == 130 || angle == 140){
	dispersionL=Dispersion[0][angle];
	dispersionR=Dispersion[1][angle];
      }

      // CMS and TOTEM have inverted z direction!
      // Left = 45 = Totem z(-) = CMS z(+).
      // Right = 56 = Totem z(+) = CMS z(-).
      xi1 = (protonArm45_x-relcorr)/dispersionL;
      proton1_pz = (1.-(xi1))*(ECM/2.);

      xi2 = (protonArm56_x-relcorr)/dispersionR;
      proton2_pz = -(1.-(xi2))*(ECM/2.);

      p1.SetPxPyPzE(0.,0., ECM*xi1/2., ECM*xi1/2.);
      p2.SetPxPyPzE(0.,0., -ECM*xi2/2., ECM*xi2/2.);
      diffmass=ECM*sqrt((  (xi1)*(xi2) ));

      // Not associated protons... not physical result!
      if(!(xi1>0. && xi2>0.)) continue;

      // HLT Muons: 'HLT_IsoMu27_v*', 'HLT_DoubleMu43NoFiltersNoVtx_v*', 'HLT_DoubleMu48NoFiltersNoVtx_v*' 
      // HLT Electrons: 'HLT_DoubleEle33_CaloIdL_MW_v*','HLT_Ele27_WPTight_Gsf_v*'
      // If -1: not present in the path, 0: not triggered and 1 triggered.
      if(strcmp(mode, "Muon")==0){
	if(trigger->at(0)==1) triggerSingle = true;
	if(trigger->at(1)==1 || trigger->at(2)==1) triggerDouble = true;
      }else if(strcmp(mode, "Electron")==0){
	if(trigger->at(0)==1) triggerDouble = true;
	if(trigger->at(1)==1) triggerSingle = true;
      }else{
	std::cout << "\nNo Mode option!\n" << std::endl;
	exit(EXIT_FAILURE);
      }

      if(!(leptons_pt->size()>1&&jetsak4_pt->size()>1)) continue;

      // Leptons and Jets are already sorted by pT 
      lepton1.SetPtEtaPhiE(leptons_pt->at(0), leptons_eta->at(0), leptons_phi->at(0), leptons_energy->at(0));
      lepton2.SetPtEtaPhiE(leptons_pt->at(1), leptons_eta->at(1), leptons_phi->at(1), leptons_energy->at(1));

      jet1.SetPtEtaPhiE(jetsak4_pt->at(0), jetsak4_eta->at(0), jetsak4_phi->at(0), jetsak4_energy->at(0));
      jet2.SetPtEtaPhiE(jetsak4_pt->at(1), jetsak4_eta->at(1), jetsak4_phi->at(1), jetsak4_energy->at(1));

      dilepton=lepton1+lepton2;
      dijet=jet1+jet2;

      // Invisible particle associated with dilepton
      X = p1 + p2 - dilepton;

      // By chance, some protons not associated with the event...
      if(!(X.M()>0)) continue;

      if(debug){
	std::cout << "\n\t == Event " << jentry << " ==" << std::endl;
	for(std::vector<int>::size_type i = 0; i != trigger->size(); i++){
	  std::cout << "\t\t --> Trigger (" << i << "): " << trigger->at(i) << std::endl;
	}
	std::cout << "\t\t --> # Jets: " << jetsak4_pt->size() << std::endl;
	std::cout << "\t\t --> # Fat Jets: " << jetsak8_phi->size() << std::endl;
	std::cout << "\t\t --> # Leptons: " << leptons_pt->size() << std::endl;
	std::cout << "\t\t --> # Protons (strips && pixels && timing): " << protonsArm->size() << std::endl;
	std::cout << "\t\t --> @Pixel SEC45: Proton1, pz [GeV]: " << proton1_pz << ", xi1: " << xi1 << std::endl;
	std::cout << "\t\t --> @Pixel SEC56: Proton2, pz [GeV]: " << proton2_pz << ", xi2: " << xi2 << std::endl;
	std::cout << "\t\t --> Dilepton Mass [GeV]: " << dilepton.M() << std::endl;
	std::cout << "\t\t --> Dijet Mass [GeV]: " << dijet.M() << std::endl;
	std::cout << "\t\t --> X Missing Mass [GeV]: " << X.M() << std::endl;
      }

      // Filling TTree
      run_ = run;
      event_ = event;
      lumiblock_ = lumiblock;

      nLeptons = leptons_pt->size();
      leadingLeptonEnergy = leptons_energy->at(0);
      leadingLeptonPt = leptons_pt->at(0);
      leadingLeptonEta = leptons_eta->at(0);
      leadingLeptonPhi = leptons_phi->at(0);
      leadingLeptonCharge = leptons_charge->at(0);
      leadingLeptonVz = leptons_vz->at(0);
      leadingLeptonLooseId = leptons_looseId->at(0);
      leadingLeptonMediumId = leptons_mediumId->at(0);
      leadingLeptonTightId = leptons_tightId->at(0);
      leadingLeptonPfIsoMedium = leptons_pfIsoMedium_->at(0);
      leadingLeptonMiniIsoTight = leptons_miniIsoTight_->at(0);
      leadingLeptonPfIsoVeryTight = leptons_pfIsoVeryTight_->at(0);

      secondLeptonEnergy = leptons_energy->at(1);
      secondLeptonPt = leptons_pt->at(1);
      secondLeptonEta = leptons_eta->at(1);
      secondLeptonPhi = leptons_phi->at(1);
      secondLeptonCharge = leptons_charge->at(1);
      secondLeptonVz = leptons_vz->at(1);
      secondLeptonLooseId = leptons_looseId->at(1);
      secondLeptonMediumId = leptons_mediumId->at(1);
      secondLeptonTightId = leptons_tightId->at(1);
      secondLeptonPfIsoMedium = leptons_pfIsoMedium_->at(1);
      secondLeptonMiniIsoTight = leptons_miniIsoTight_->at(1);
      secondLeptonPfIsoVeryTight = leptons_pfIsoVeryTight_->at(1);

      nJets = jetsak4_pt->size();
      leadingJetEnergy = jetsak4_energy->at(0);
      leadingJetPt = jetsak4_pt->at(0);
      leadingJetEta = jetsak4_eta->at(0);
      leadingJetPhi = jetsak4_phi->at(0);
      leadingJetVz = jetsak4_vz->at(0);
      leadingJetBtag = jetsak4_bdis->at(0);
      secondJetEnergy = jetsak4_energy->at(1);
      secondJetPt = jetsak4_pt->at(1);
      secondJetEta = jetsak4_eta->at(1);
      secondJetPhi = jetsak4_phi->at(1);
      secondJetVz = jetsak4_vz->at(1);
      secondJetBtag = jetsak4_bdis->at(1);

      dileptonCharge = leptons_charge->at(0) + leptons_charge->at(1);
      dileptonMass = dilepton.M();
      dileptonEta = dilepton.Eta();
      dileptonPhi = dilepton.Phi();
      dileptonPt = dilepton.Pt();
      dileptonRapidity = dilepton.Rapidity();
      dijetMass = dijet.M();
      dijetEta = dijet.Eta();
      dijetPhi = dijet.Phi();
      dijetPt = dijet.Pt();
      dijetRapidity = dijet.Rapidity();
      missingMass = X.M();
      missingEta = X.Eta();
      missingPhi = X.Phi();
      missingPt = X.Pt();
      missingRapidity = X.Rapidity();

      missEt_ = missEt;
      missEt_phi_ = missEt_phi;
      nPF_ = nPF;
      sumEnergyPF_ = sumEnergyPF;

      double dphi = fabs(dilepton.Phi()-dijet.Phi());
      dphi = (dphi<pi) ? dphi : 2.*pi-dphi;
      dphi_ = dphi;

      tout->Fill();

    } // proton loop

  } // end of the event loop

  std::cout << "\n\n<END> good luck!\n" << std::endl;
  tout->Write();
  fileout->Close();

}

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)
  {
    return *itr;
  }
  return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
  return std::find(begin, end, option) != end;
}

int getXangle(int run, int lumi, const char* filename)
{

  TString drun;
  TString dfill;
  TString dlumi;
  TString temp;
  int Xangle=-1;

  TString runs;runs.Form("%d", run);
  TString lumis;lumis.Form("%d", lumi);
  std::ifstream F;

  F.open((const char*)filename);
  int counter=0;
  if(F){
    while (!F.eof())
    {
      F>>drun;
      F>>dfill;
      F>>dlumi;
      F>>temp;
      F>>Xangle;
      if( runs == drun &&  lumis==dlumi )
      {
	break;
      }
    }
  }//endif
  else std::cout << "[!] gerXangle():Error reading from file" << std::endl;
  if(F.eof())
  {
    std::cout << "[getXangle() warning:] No Xangle data for this run found!" << std::endl;
    F.close();
    return -1;
  }
  else  {F.close();
    return Xangle;
  }
}

int main(int argc, char * argv[])
{

  std::cout << "\n===================" << std::endl;
  std::cout << "Missing Mass Search" << std::endl;
  std::cout << "===================" << std::endl;

  TTree* tree = NULL;
  if(cmdOptionExists(argv, argv+argc, "--h")||cmdOptionExists(argv, argv+argc, "--help"))
  {
    std::cout << "\n== Help ==\n" << std::endl;
    std::cout << "\t --f filename.root (input file)" << std::endl;
    std::cout << "\t --era B (B, C, D, E or F)" << std::endl;
    std::cout << "\t --mode Muon (or Electron)" << std::endl;
    std::cout << "\t --jobid job1 (tag to be added in the outputfile)\n" << std::endl;
    return 0;
  }

  char * filename = getCmdOption(argv, argv + argc, "--f");
  char * era = getCmdOption(argv, argv + argc, "--era");
  char * mode = getCmdOption(argv, argv + argc, "--mode");
  char * jobid = getCmdOption(argv, argv + argc, "--jobid");

  if (filename && era && mode)
  {
    TTree* tree = NULL;
    TFile filecheck(filename);
    if(filecheck.IsZombie()){
      std::cout << "\n\t ---> Corrupted file! Please try another file!\n" << std::endl;
      return 0;
    }

    if (!filecheck.GetDirectory("missing_mass")){
      std::cout << "\n---------------------------------------------------" << std::endl;
      std::cout << " There is no directory/path "<< std::endl;
      std::cout << " in the file." << std::endl;
      std::cout << "---------------------------------------------------\n" << std::endl;
      return 0;
    }

    std::cout << "Reading file " << filename << std::endl;

    TFile* file = TFile::Open(filename,"READ");
    tree = (TTree*) file->Get( "missing_mass/analyzer" );

    bool preTS2 = false;
    if(strcmp(era, "B")==0 || strcmp(era, "C")==0 || strcmp(era, "D")==0 || strcmp(era, "b")==0 || strcmp(era, "c")==0 || strcmp(era, "d")==0) preTS2 = true;
    else if(strcmp(era, "E")==0 || strcmp(era, "F")==0 || strcmp(era, "e")==0 || strcmp(era, "f")==0){
      preTS2 = false;
    }else{
      std::cout << "\nNo --era option! It should be b, c, d, e or f.\n" << std::endl;
      exit(EXIT_FAILURE);
    }

    // Accessing Missing Mass Object
    MissingMassNtupleAnalyzer m(tree); 
    m.Loop(era, mode, preTS2, jobid);
  }else{
    std::cout << "\n\t --> Please, insert --f filename.root and --era B (or C, D, E, F) --mode Muon (or Electron)\n" << std::endl;
  }

  return 0;
}
