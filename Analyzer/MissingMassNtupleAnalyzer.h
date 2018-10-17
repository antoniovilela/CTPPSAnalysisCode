//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Oct 14 14:42:49 2018 by ROOT version 6.10/09
// from TTree analyzer/analyzer
// found on file: /eos/cms/store/group/phys_pps/MissingMassSearch/MuonB-MiniAOD-v2/DoubleMuon/MuonB-MiniAOD-v2/181013_221704/0000/output_30.root
//////////////////////////////////////////////////////////

#ifndef MissingMassNtupleAnalyzer_h
#define MissingMassNtupleAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

int getXangle(int , int, const char*);

class MissingMassNtupleAnalyzer {
  public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Fixed size dimensions of array or collections stored in the TTree if any.
    static constexpr Int_t kMaxleptons_pfIsoMedium = 1;
    static constexpr Int_t kMaxleptons_miniIsoTight = 1;
    static constexpr Int_t kMaxleptons_pfIsoVeryTight = 1;
    static constexpr Int_t kMaxleptons_pfIso = 1;
    static constexpr Int_t kMaxleptons_tkIso = 1;

    // Declaration of leaf types
    Int_t           run;
    Long64_t        event;
    Int_t           lumiblock;
    std::vector<int>     *trigger;
    std::vector<float>   *vertex_x;
    std::vector<float>   *vertex_y;
    std::vector<float>   *vertex_z;
    std::vector<float>   *vertex_ntrack;
    std::vector<float>   *vertex_chi2;
    std::vector<float>   *vertex_ndof;
    std::vector<float>   *jetsak4_pt;
    std::vector<float>   *jetsak4_energy;
    std::vector<float>   *jetsak4_phi;
    std::vector<float>   *jetsak4_eta;
    std::vector<float>   *jetsak4_vz;
    std::vector<float>   *jetsak4_bdis;
    std::vector<bool>    *jetsak4_looseId;
    std::vector<bool>    *jetsak4_tightId;
    std::vector<bool>    *jetsak4_lepVeto;
    std::vector<float>   *jetsak8_pt;
    std::vector<float>   *jetsak8_energy;
    std::vector<float>   *jetsak8_phi;
    std::vector<float>   *jetsak8_eta;
    std::vector<float>   *jetsak8_vz;
    std::vector<float>   *jetsak8_bdis;
    std::vector<bool>    *jetsak8_looseId;
    std::vector<bool>    *jetsak8_tightId;
    std::vector<bool>    *jetsak8_lepVeto;
    std::vector<float>   *leptons_energy;
    std::vector<float>   *leptons_pt;
    std::vector<float>   *leptons_eta;
    std::vector<float>   *leptons_phi;
    std::vector<float>   *leptons_px;
    std::vector<float>   *leptons_py;
    std::vector<float>   *leptons_pz;
    std::vector<float>   *leptons_charge;
    std::vector<float>   *leptons_vx;
    std::vector<float>   *leptons_vy;
    std::vector<float>   *leptons_vz;
    std::vector<bool>    *leptons_looseId;
    std::vector<bool>    *leptons_mediumId;
    std::vector<bool>    *leptons_tightId;
    std::vector<bool>    *leptons_pfIsoMedium_;
    std::vector<bool>    *leptons_miniIsoTight_;
    std::vector<bool>    *leptons_pfIsoVeryTight_;
    std::vector<float>   *leptons_pfIso_;
    std::vector<float>   *leptons_tkIso_;
    Float_t         missEt;
    Float_t         missEt_phi;
    Int_t           nPF;
    Float_t         sumEnergyPF;
    std::vector<int>     *protonsArm;
    std::vector<int>     *protonsStation;
    std::vector<int>     *protonsRP;
    std::vector<float>   *protonsX;
    std::vector<float>   *protonsXUnc;
    std::vector<float>   *protonsY;
    std::vector<float>   *protonsYUnc;
    std::vector<float>   *protonsTime;
    std::vector<float>   *protonsTimeUnc;

    // List of branches
    TBranch        *b_run;   //!
    TBranch        *b_event;   //!
    TBranch        *b_lumiblock;   //!
    TBranch        *b_trigger;   //!
    TBranch        *b_vertex_x;   //!
    TBranch        *b_vertex_y;   //!
    TBranch        *b_vertex_z;   //!
    TBranch        *b_vertex_ntrack;   //!
    TBranch        *b_vertex_chi2;   //!
    TBranch        *b_vertex_ndof;   //!
    TBranch        *b_jetsak4_pt;   //!
    TBranch        *b_jetsak4_energy;   //!
    TBranch        *b_jetsak4_phi;   //!
    TBranch        *b_jetsak4_eta;   //!
    TBranch        *b_jetsak4_vz;   //!
    TBranch        *b_jetsak4_bdis;   //!
    TBranch        *b_jetsak4_looseId;   //!
    TBranch        *b_jetsak4_tightId;   //!
    TBranch        *b_jetsak4_lepVeto;   //!
    TBranch        *b_jetsak8_pt;   //!
    TBranch        *b_jetsak8_energy;   //!
    TBranch        *b_jetsak8_phi;   //!
    TBranch        *b_jetsak8_eta;   //!
    TBranch        *b_jetsak8_vz;   //!
    TBranch        *b_jetsak8_bdis;   //!
    TBranch        *b_jetsak8_looseId;   //!
    TBranch        *b_jetsak8_tightId;   //!
    TBranch        *b_jetsak8_lepVeto;   //!
    TBranch        *b_leptons_energy;   //!
    TBranch        *b_leptons_pt;   //!
    TBranch        *b_leptons_eta;   //!
    TBranch        *b_leptons_phi;   //!
    TBranch        *b_leptons_px;   //!
    TBranch        *b_leptons_py;   //!
    TBranch        *b_leptons_pz;   //!
    TBranch        *b_leptons_charge;   //!
    TBranch        *b_leptons_vx;   //!
    TBranch        *b_leptons_vy;   //!
    TBranch        *b_leptons_vz;   //!
    TBranch        *b_leptons_looseId;   //!
    TBranch        *b_leptons_mediumId;   //!
    TBranch        *b_leptons_tightId;   //!
    TBranch        *b_leptons_pfIsoMedium_;   //!
    TBranch        *b_leptons_miniIsoTight_;   //!
    TBranch        *b_leptons_pfIsoVeryTight_;   //!
    TBranch        *b_leptons_pfIso_;   //!
    TBranch        *b_leptons_tkIso_;   //!
    TBranch        *b_missEt;   //!
    TBranch        *b_missEt_phi;   //!
    TBranch        *b_nPF;   //!
    TBranch        *b_SumEnergyPF;   //!
    TBranch        *b_protonsArm;   //!
    TBranch        *b_protonsStation;   //!
    TBranch        *b_protonsRP;   //!
    TBranch        *b_protonsX;   //!
    TBranch        *b_protonsXUnc;   //!
    TBranch        *b_protonsY;   //!
    TBranch        *b_protonsYUnc;   //!
    TBranch        *b_protonsTime;   //!
    TBranch        *b_protonsTimeUnc;   //!

    MissingMassNtupleAnalyzer(TTree *tree=0);
    virtual ~MissingMassNtupleAnalyzer();
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     Loop(char*, char*, bool, char*, char*);
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);

    char* getCmdOption(char**, char**, const std::string&);
    bool cmdOptionExists(char**, char**, const std::string&);
};

#endif

#ifdef MissingMassNtupleAnalyzer_cxx
MissingMassNtupleAnalyzer::MissingMassNtupleAnalyzer(TTree *tree) : fChain(0) 
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/cms/store/group/phys_pps/MissingMassSearch/MuonB-MiniAOD-v2/DoubleMuon/MuonB-MiniAOD-v2/181013_221704/0000/output_30.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("/eos/cms/store/group/phys_pps/MissingMassSearch/MuonB-MiniAOD-v2/DoubleMuon/MuonB-MiniAOD-v2/181013_221704/0000/output_30.root");
    }
    TDirectory * dir = (TDirectory*)f->Get("/eos/cms/store/group/phys_pps/MissingMassSearch/MuonB-MiniAOD-v2/DoubleMuon/MuonB-MiniAOD-v2/181013_221704/0000/output_30.root:/missing_mass");
    dir->GetObject("analyzer",tree);

  }
  Init(tree);
}

MissingMassNtupleAnalyzer::~MissingMassNtupleAnalyzer()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t MissingMassNtupleAnalyzer::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t MissingMassNtupleAnalyzer::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void MissingMassNtupleAnalyzer::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
  trigger = 0;
  vertex_x = 0;
  vertex_y = 0;
  vertex_z = 0;
  vertex_ntrack = 0;
  vertex_chi2 = 0;
  vertex_ndof = 0;
  jetsak4_pt = 0;
  jetsak4_energy = 0;
  jetsak4_phi = 0;
  jetsak4_eta = 0;
  jetsak4_vz = 0;
  jetsak4_bdis = 0;
  jetsak4_looseId = 0;
  jetsak4_tightId = 0;
  jetsak4_lepVeto = 0;
  jetsak8_pt = 0;
  jetsak8_energy = 0;
  jetsak8_phi = 0;
  jetsak8_eta = 0;
  jetsak8_vz = 0;
  jetsak8_bdis = 0;
  jetsak8_looseId = 0;
  jetsak8_tightId = 0;
  jetsak8_lepVeto = 0;
  leptons_energy = 0;
  leptons_pt = 0;
  leptons_eta = 0;
  leptons_phi = 0;
  leptons_px = 0;
  leptons_py = 0;
  leptons_pz = 0;
  leptons_charge = 0;
  leptons_vx = 0;
  leptons_vy = 0;
  leptons_vz = 0;
  leptons_looseId = 0;
  leptons_mediumId = 0;
  leptons_tightId = 0;
  leptons_pfIsoMedium_ = 0;
  leptons_miniIsoTight_ = 0;
  leptons_pfIsoVeryTight_ = 0;
  leptons_pfIso_ = 0;
  leptons_tkIso_ = 0;
  protonsArm = 0;
  protonsStation = 0;
  protonsRP = 0;
  protonsX = 0;
  protonsXUnc = 0;
  protonsY = 0;
  protonsYUnc = 0;
  protonsTime = 0;
  protonsTimeUnc = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("run", &run, &b_run);
  fChain->SetBranchAddress("event", &event, &b_event);
  fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
  fChain->SetBranchAddress("trigger", &trigger, &b_trigger);
  fChain->SetBranchAddress("vertex_x", &vertex_x, &b_vertex_x);
  fChain->SetBranchAddress("vertex_y", &vertex_y, &b_vertex_y);
  fChain->SetBranchAddress("vertex_z", &vertex_z, &b_vertex_z);
  fChain->SetBranchAddress("vertex_ntrack", &vertex_ntrack, &b_vertex_ntrack);
  fChain->SetBranchAddress("vertex_chi2", &vertex_chi2, &b_vertex_chi2);
  fChain->SetBranchAddress("vertex_ndof", &vertex_ndof, &b_vertex_ndof);
  fChain->SetBranchAddress("jetsak4_pt", &jetsak4_pt, &b_jetsak4_pt);
  fChain->SetBranchAddress("jetsak4_energy", &jetsak4_energy, &b_jetsak4_energy);
  fChain->SetBranchAddress("jetsak4_phi", &jetsak4_phi, &b_jetsak4_phi);
  fChain->SetBranchAddress("jetsak4_eta", &jetsak4_eta, &b_jetsak4_eta);
  fChain->SetBranchAddress("jetsak4_vz", &jetsak4_vz, &b_jetsak4_vz);
  fChain->SetBranchAddress("jetsak4_bdis", &jetsak4_bdis, &b_jetsak4_bdis);
  fChain->SetBranchAddress("jetsak4_looseId", &jetsak4_looseId, &b_jetsak4_looseId);
  fChain->SetBranchAddress("jetsak4_tightId", &jetsak4_tightId, &b_jetsak4_tightId);
  fChain->SetBranchAddress("jetsak4_lepVeto", &jetsak4_lepVeto, &b_jetsak4_lepVeto);
  fChain->SetBranchAddress("jetsak8_pt", &jetsak8_pt, &b_jetsak8_pt);
  fChain->SetBranchAddress("jetsak8_energy", &jetsak8_energy, &b_jetsak8_energy);
  fChain->SetBranchAddress("jetsak8_phi", &jetsak8_phi, &b_jetsak8_phi);
  fChain->SetBranchAddress("jetsak8_eta", &jetsak8_eta, &b_jetsak8_eta);
  fChain->SetBranchAddress("jetsak8_vz", &jetsak8_vz, &b_jetsak8_vz);
  fChain->SetBranchAddress("jetsak8_bdis", &jetsak8_bdis, &b_jetsak8_bdis);
  fChain->SetBranchAddress("jetsak8_looseId", &jetsak8_looseId, &b_jetsak8_looseId);
  fChain->SetBranchAddress("jetsak8_tightId", &jetsak8_tightId, &b_jetsak8_tightId);
  fChain->SetBranchAddress("jetsak8_lepVeto", &jetsak8_lepVeto, &b_jetsak8_lepVeto);
  fChain->SetBranchAddress("leptons_energy", &leptons_energy, &b_leptons_energy);
  fChain->SetBranchAddress("leptons_pt", &leptons_pt, &b_leptons_pt);
  fChain->SetBranchAddress("leptons_eta", &leptons_eta, &b_leptons_eta);
  fChain->SetBranchAddress("leptons_phi", &leptons_phi, &b_leptons_phi);
  fChain->SetBranchAddress("leptons_px", &leptons_px, &b_leptons_px);
  fChain->SetBranchAddress("leptons_py", &leptons_py, &b_leptons_py);
  fChain->SetBranchAddress("leptons_pz", &leptons_pz, &b_leptons_pz);
  fChain->SetBranchAddress("leptons_charge", &leptons_charge, &b_leptons_charge);
  fChain->SetBranchAddress("leptons_vx", &leptons_vx, &b_leptons_vx);
  fChain->SetBranchAddress("leptons_vy", &leptons_vy, &b_leptons_vy);
  fChain->SetBranchAddress("leptons_vz", &leptons_vz, &b_leptons_vz);
  fChain->SetBranchAddress("leptons_looseId", &leptons_looseId, &b_leptons_looseId);
  fChain->SetBranchAddress("leptons_mediumId", &leptons_mediumId, &b_leptons_mediumId);
  fChain->SetBranchAddress("leptons_tightId", &leptons_tightId, &b_leptons_tightId);
  fChain->SetBranchAddress("leptons_pfIsoMedium_", &leptons_pfIsoMedium_, &b_leptons_pfIsoMedium_);
  fChain->SetBranchAddress("leptons_miniIsoTight_", &leptons_miniIsoTight_, &b_leptons_miniIsoTight_);
  fChain->SetBranchAddress("leptons_pfIsoVeryTight_", &leptons_pfIsoVeryTight_, &b_leptons_pfIsoVeryTight_);
  fChain->SetBranchAddress("leptons_pfIso_", &leptons_pfIso_, &b_leptons_pfIso_);
  fChain->SetBranchAddress("leptons_tkIso_", &leptons_tkIso_, &b_leptons_tkIso_);
  fChain->SetBranchAddress("missEt", &missEt, &b_missEt);
  fChain->SetBranchAddress("missEt_phi", &missEt_phi, &b_missEt_phi);
  fChain->SetBranchAddress("nPF", &nPF, &b_nPF);
  fChain->SetBranchAddress("sumEnergyPF", &sumEnergyPF, &b_SumEnergyPF);
  fChain->SetBranchAddress("protonsArm", &protonsArm, &b_protonsArm);
  fChain->SetBranchAddress("protonsStation", &protonsStation, &b_protonsStation);
  fChain->SetBranchAddress("protonsRP", &protonsRP, &b_protonsRP);
  fChain->SetBranchAddress("protonsX", &protonsX, &b_protonsX);
  fChain->SetBranchAddress("protonsXUnc", &protonsXUnc, &b_protonsXUnc);
  fChain->SetBranchAddress("protonsY", &protonsY, &b_protonsY);
  fChain->SetBranchAddress("protonsYUnc", &protonsYUnc, &b_protonsYUnc);
  fChain->SetBranchAddress("protonsTime", &protonsTime, &b_protonsTime);
  fChain->SetBranchAddress("protonsTimeUnc", &protonsTimeUnc, &b_protonsTimeUnc);
  Notify();
}

Bool_t MissingMassNtupleAnalyzer::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void MissingMassNtupleAnalyzer::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t MissingMassNtupleAnalyzer::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef MissingMassNtupleAnalyzer_cxx
