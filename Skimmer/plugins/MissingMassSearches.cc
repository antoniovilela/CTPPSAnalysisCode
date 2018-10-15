// -*- C++ -*-
//
// Package:    CTPPSAnalysisCode/MissingMassSearches
// Class:      MissingMassSearches
// 
/**\class MissingMassSearches MissingMassSearches.cc CTPPSAnalysisCode/MissingMassSearches/plugins/MissingMassSearches.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Diego Figueiredo and Nicola Turini
//         Created:  Fri, 05 Oct 2018 09:08:27 GMT
//
//

// user include files
#include "MissingMassSearches.h"

//
// constructors and destructor
//
MissingMassSearches::MissingMassSearches(const edm::ParameterSet& iConfig):
  debug_               ( iConfig.getParameter<bool>( "debugging" ) ),
  mode_                ( iConfig.getParameter<std::string>( "mode" ) ),
  includeMuons_        ( iConfig.getParameter<bool>( "includeMuons" ) ),
  includeElectrons_    ( iConfig.getParameter<bool>( "includeElectrons" ) ),
  includeJets_         ( iConfig.getParameter<bool>( "includeJets" ) ),
  includePhotons_      ( iConfig.getParameter<bool>( "includePhotons" ) ),
  includeMET_          ( iConfig.getParameter<bool>( "includeMET" ) ),
  includeVertices_     ( iConfig.getParameter<bool>( "includeVertices" ) ),
  includePF_           ( iConfig.getParameter<bool>( "includeParticleFlow" ) ),
  includeProtons_      ( iConfig.getParameter<bool>( "includeProtons" ) ),
  triggersList_        ( iConfig.getParameter<std::vector<std::string>>              ( "triggersList" ) ),
  triggerResultsToken_ ( consumes<edm::TriggerResults>                               ( iConfig.getParameter<edm::InputTag>( "triggerResults" ) ) ),
  jetTokenA_           ( consumes<edm::View<pat::Jet>>                               ( iConfig.getParameter<edm::InputTag>( "JetAlgoA" ) ) ),
  jetTokenB_           ( consumes<edm::View<pat::Jet>>                               ( iConfig.getParameter<edm::InputTag>( "JetAlgoB" ) ) ),
  eleToken_            ( consumes<edm::View<pat::Electron>>                          ( iConfig.getParameter<edm::InputTag>( "electronTag" ) ) ),
  muonToken_           ( consumes<edm::View<pat::Muon>>                              ( iConfig.getParameter<edm::InputTag>( "muonTag" ) ) ),
  pfToken_             ( consumes<std::vector< reco::PFCandidate>>                   ( iConfig.getParameter<edm::InputTag>( "pfTag" ) ) ),
  packedToken_         ( consumes<std::vector< pat::PackedCandidate>>                ( iConfig.getParameter<edm::InputTag>( "packedTag" ) ) ),
  photonToken_         ( consumes<edm::View<pat::Photon>>                            ( iConfig.getParameter<edm::InputTag>( "photonTag" ) ) ),
  metToken_            ( consumes<edm::View<pat::MET>>                               ( iConfig.getParameter<edm::InputTag>( "metTag" ) ) ),
  vertexToken_         ( consumes<edm::View<reco::Vertex>>                           ( iConfig.getParameter<edm::InputTag>( "vertexTag" ) ) ),
  protonToken_         ( consumes<vector<CTPPSLocalTrackLite>>                       ( iConfig.getParameter<edm::InputTag>( "protonTag" ) ) ),
  pixelsppsToken_      ( consumes<edm::DetSetVector<CTPPSPixelLocalTrack>>           ( iConfig.getParameter<edm::InputTag>( "pixelsppsTag" ) ) ),
  timingppsToken_      ( consumes<edm::DetSetVector<CTPPSDiamondLocalTrack>>         ( iConfig.getParameter<edm::InputTag>( "timingppsTag" ) ) ),
  stripstotemToken_    ( consumes<edm::DetSetVector<TotemRPLocalTrack>>              ( iConfig.getParameter<edm::InputTag>( "stripstotemTag" ) ) )
{
  //now do what ever initialization is needed
  usesResource("TFileService");

  if ( mode_ == "Muon" || mode_ =="Muons" ) {
    includeElectrons_ = false;
    includeMuons_ = true;
  }
  else if ( mode_ == "Electron" || mode_=="Electrons" ) {
    includeElectrons_ = true;
    includeMuons_ = false;
  }
  else throw cms::Exception( "MissingMassSearches" ) << "'mode' parameter should either be:\n"
    << "   * 'Electron' or 'Muon' (for same-flavour leptons)";

  if(includeElectrons_){
    const edm::ParameterSet eleIdLabelSet = iConfig.getParameter<edm::ParameterSet>( "eleIdLabels" );
    eleMediumIdLabel_ = eleIdLabelSet.getParameter<edm::InputTag>( "mediumLabel" ).encode();
    eleTightIdLabel_ = eleIdLabelSet.getParameter<edm::InputTag>( "tightLabel" ).encode();
  }

}


MissingMassSearches::~MissingMassSearches()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//


// ------------ Cleaning  ------------
void MissingMassSearches::cleaning(){

  triggerVec.clear();
  jetsVecA.clear(); // AK4
  jetsVecB.clear(); // AK8
  electronsVec.clear();
  muonsVec.clear();
  packedVec.clear();
  pfVec.clear();
  photonsVec.clear();
  vtxVec.clear();
  protonsVec.clear();
  pixelsVec.clear();
  timingVec.clear();
  stripsVec.clear();

  jetsak4Cand.clear();
  jetsak8Cand.clear();

  MuonsJetsCand.clear();
  ElectronsJetsCand.clear();

}

// ------------ Cleaning  ------------
void MissingMassSearches::eventClear(){

  *run_=-999;
  *ev_=-999;
  *lumiblock_=-999;

  (*trigger_).clear();

  if(includeVertices_){
    (*vertex_x_).clear();
    (*vertex_y_).clear();
    (*vertex_z_).clear();

    (*vertex_ntrack_).clear();
    (*vertex_chi2_).clear();
    (*vertex_ndof_).clear();
  }

  if(includeJets_){
    (*jetsak4_pt_).clear();
    (*jetsak4_energy_).clear();
    (*jetsak4_phi_).clear();
    (*jetsak4_eta_).clear();
    (*jetsak4_vz_).clear();
    (*jetsak4_bdis_).clear();
    (*jetsak4_looseJetId_).clear();
    (*jetsak4_tightJetId_).clear();
    (*jetsak4_lepVetoJetId_).clear();
    (*jetsak8_pt_).clear();
    (*jetsak8_energy_).clear();
    (*jetsak8_phi_).clear();
    (*jetsak8_eta_).clear();
    (*jetsak8_vz_).clear();
    (*jetsak8_bdis_).clear();
    (*jetsak8_looseJetId_).clear();
    (*jetsak8_tightJetId_).clear();
    (*jetsak8_lepVetoJetId_).clear();
  }

  if(includeElectrons_ || includeMuons_){
    (*leptons_energy_).clear();
    (*leptons_pt_).clear();
    (*leptons_eta_).clear();
    (*leptons_phi_).clear();
    (*leptons_px_).clear();
    (*leptons_py_).clear();
    (*leptons_pz_).clear();
    (*leptons_charge_).clear();
    (*leptons_vx_).clear();
    (*leptons_vy_).clear();
    (*leptons_vz_).clear();
    (*leptons_looseId_).clear();
    (*leptons_mediumId_).clear();
    (*leptons_tightId_).clear();
    (*leptons_pfIsoMedium_).clear();
    (*leptons_miniIsoTight_).clear();
    (*leptons_pfIsoVeryTight_).clear();
    (*leptons_pfIso_).clear();
    (*leptons_tkIso_).clear();
  }

  if(includeProtons_){
    (*protonsArm_).clear();
    (*protonsStation_).clear();
    (*protonsRP_).clear();
    (*protonsX_).clear();
    (*protonsXUnc_).clear();
    (*protonsY_).clear();
    (*protonsYUnc_).clear();
    (*protonsTime_).clear();
    (*protonsTimeUnc_).clear();
  }

}

// ------------ running over muons  ------------
void MissingMassSearches::fetchMuons(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  MuonsEvent *mu;
  mu = new MuonsEvent(iEvent, iSetup, muonToken_);
  muonsVec = mu->GetMuons();

}

// ------------ running over electrons  ------------
void MissingMassSearches::fetchElectrons(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  ElectronsEvent *ele;
  ele = new ElectronsEvent(iEvent, iSetup, eleToken_);
  electronsVec = ele->GetElectrons();

}

// ------------ running over jets  ------------
void MissingMassSearches::fetchJets(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  JetsEvent *jetAlgoA;
  jetAlgoA = new JetsEvent(iEvent, iSetup, jetTokenA_);
  jetsVecA = jetAlgoA->GetJets();

  JetsEvent *jetAlgoB;
  jetAlgoB = new JetsEvent(iEvent, iSetup, jetTokenB_);
  jetsVecB = jetAlgoB->GetJets();

}

// ------------ running over photons  ------------
void MissingMassSearches::fetchPhotons(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  PhotonsEvent *photon;
  photon = new PhotonsEvent(iEvent, iSetup, photonToken_);
  photonsVec = photon->GetPhotons();

}

// ------------ running over MET  ------------
void MissingMassSearches::fetchMET(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  METEvent *missinget;
  missinget = new METEvent(iEvent, iSetup, metToken_);
  met = missinget->GetMET();

  *misset_ = met->et();
  *misset_phi_ = met->phi();

}

// ------------ running over vertices  ------------
void MissingMassSearches::fetchVertices(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  VerticesEvent *vtx;
  vtx = new VerticesEvent(iEvent, iSetup, vertexToken_);
  vtxVec = vtx->GetVertices();

  if(includeVertices_){
    for ( const auto& vtxEvt: vtxVec){
      (*vertex_x_).push_back(vtxEvt->x());
      (*vertex_y_).push_back(vtxEvt->y());
      (*vertex_z_).push_back(vtxEvt->z());
      (*vertex_ntrack_).push_back(vtxEvt->nTracks());
      (*vertex_chi2_).push_back(vtxEvt->chi2());
      (*vertex_ndof_).push_back(vtxEvt->ndof());
    }
  }

}

// ------------ running over PF  ------------
void MissingMassSearches::fetchPF(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  ParticleFlowEvent *packed;
  packed = new ParticleFlowEvent(iEvent, iSetup, packedToken_);
  packedVec = packed->GetPackedFlow();

  ParticleFlowEvent *pf;
  pf = new ParticleFlowEvent(iEvent, iSetup, pfToken_);
  pfVec = pf->GetParticleFlow();

  if(packedVec.empty()){
    float sumEnergyPF = 0;
    for ( const auto& pfEvt: pfVec){
      sumEnergyPF+=pfEvt->energy();
    }
    *nPF_=pfVec.size();
    *SumPF_energy_=sumEnergyPF;
  }else{
    float sumEnergyPF = 0;
    for ( const auto& packedEvt: packedVec){
      sumEnergyPF+=packedEvt->energy();
    }
    *nPF_=packedVec.size();
    *SumPF_energy_=sumEnergyPF;
  }

}

// ------------ running over protons  ------------
void MissingMassSearches::fetchProtons(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  ProtonsEvent *protons;
  protons = new ProtonsEvent(iEvent, iSetup, protonToken_);
  protonsVec = protons->GetProtonsLite();

  // Not in Miniaod, Proton Event take care of it!
  ProtonsEvent *pixels;
  pixels = new ProtonsEvent(iEvent, iSetup, pixelsppsToken_);
  pixelsVec = pixels->GetPixels();

  // Not in Miniaod, Proton Event take care of it!
  ProtonsEvent *timing;
  timing = new ProtonsEvent(iEvent, iSetup, timingppsToken_);
  timingVec = timing->GetTiming();

  // Not in Miniaod, Proton Event take care of it!
  ProtonsEvent *strips;
  strips = new ProtonsEvent(iEvent, iSetup, stripstotemToken_);
  stripsVec = strips->GetStrips();

  if(protonsVec.empty()){
    for ( const auto& pixelsEvt: pixelsVec ) {
      (*protonsArm_).push_back(pixelsEvt.second.arm());
      (*protonsStation_).push_back(pixelsEvt.second.station());
      (*protonsRP_).push_back(pixelsEvt.second.rp());
      (*protonsX_).push_back(pixelsEvt.first->getX0());
      (*protonsXUnc_).push_back(pixelsEvt.first->getX0Sigma());
      (*protonsY_).push_back(pixelsEvt.first->getY0());
      (*protonsYUnc_).push_back(pixelsEvt.first->getY0Sigma());
      (*protonsTime_).push_back(-1.);
      (*protonsTimeUnc_).push_back(-1.);
    }
    for ( const auto& stripsEvt: stripsVec ) {
      (*protonsArm_).push_back(stripsEvt.second.arm());
      (*protonsStation_).push_back(stripsEvt.second.station());
      (*protonsRP_).push_back(stripsEvt.second.rp());
      (*protonsX_).push_back(stripsEvt.first->getX0());
      (*protonsXUnc_).push_back(stripsEvt.first->getX0Sigma());
      (*protonsY_).push_back(stripsEvt.first->getY0());
      (*protonsYUnc_).push_back(stripsEvt.first->getY0Sigma());
      (*protonsTime_).push_back(-1.);
      (*protonsTimeUnc_).push_back(-1.);
    }
    for ( const auto& timingEvt: timingVec ) {
      (*protonsArm_).push_back(timingEvt.second.arm());
      (*protonsStation_).push_back(timingEvt.second.station());
      (*protonsRP_).push_back(timingEvt.second.rp());
      (*protonsX_).push_back(timingEvt.first->getX0());
      (*protonsXUnc_).push_back(timingEvt.first->getX0Sigma());
      (*protonsY_).push_back(timingEvt.first->getY0());
      (*protonsYUnc_).push_back(timingEvt.first->getY0Sigma());
      (*protonsTime_).push_back(timingEvt.first->getT());
      (*protonsTimeUnc_).push_back(timingEvt.first->getTSigma());
    }
  }else{
    for ( const auto& protonsEvt: protonsVec ) {
      const CTPPSDetId det_id( protonsEvt->getRPId() );
      (*protonsArm_).push_back(det_id.arm());
      (*protonsStation_).push_back(det_id.station());
      (*protonsRP_).push_back(det_id.rp());
      (*protonsX_).push_back(protonsEvt->getX());
      (*protonsXUnc_).push_back(protonsEvt->getXUnc());
      (*protonsY_).push_back(protonsEvt->getY());
      (*protonsYUnc_).push_back(protonsEvt->getYUnc());
      (*protonsTime_).push_back(protonsEvt->getTime());
      (*protonsTimeUnc_).push_back(protonsEvt->getTimeUnc());
    }
  }

}

// ------------ saving trigger  ------------
bool MissingMassSearches::fetchTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  bool trigger_fired = false;

  TriggerEvent *trigger;
  trigger = new TriggerEvent(iEvent, iSetup, triggerResultsToken_, triggersList_);
  triggerVec = trigger->GetTrigger();

  for ( const auto& triggerEvt: triggerVec) {
    if(triggerEvt==1) trigger_fired = true;
  }

  return trigger_fired;

}

// ------------ event tagger ------------
void MissingMassSearches::fetchEventTagger(const edm::Event& iEvent){

  // Force pT sorting... (CMSSW does the job)
  std::sort(muonsVec.begin(), muonsVec.end(), orderPt());
  std::sort(electronsVec.begin(), electronsVec.end(), orderPt());
  std::sort(jetsVecA.begin(), jetsVecA.end(), orderPt());
  std::sort(jetsVecB.begin(), jetsVecB.end(), orderPt());

  if(includeElectrons_){
    for ( const auto& leptonEvt: electronsVec) {
      (*leptons_energy_).push_back(leptonEvt->energy());
      (*leptons_pt_).push_back(leptonEvt->pt());
      (*leptons_eta_).push_back(leptonEvt->eta());
      (*leptons_phi_).push_back(leptonEvt->phi());
      (*leptons_px_).push_back(leptonEvt->px());
      (*leptons_py_).push_back(leptonEvt->py());
      (*leptons_pz_).push_back(leptonEvt->pz());
      (*leptons_charge_).push_back(leptonEvt->charge());
      (*leptons_vx_).push_back(leptonEvt->vertex().x());
      (*leptons_vy_).push_back(leptonEvt->vertex().y());
      (*leptons_vz_).push_back(leptonEvt->vertex().z());
      (*leptons_looseId_).push_back(leptonEvt->electronID("cutBasedElectronID-Fall17-94X-V1-loose"));
      (*leptons_mediumId_).push_back(leptonEvt->electronID("cutBasedElectronID-Fall17-94X-V1-medium"));
      (*leptons_tightId_).push_back(leptonEvt->electronID("cutBasedElectronID-Fall17-94X-V1-tight"));
      (*leptons_pfIsoMedium_).push_back(0);
      (*leptons_miniIsoTight_).push_back(0);
      (*leptons_pfIsoVeryTight_).push_back(0);
      (*leptons_pfIso_).push_back(0);
      (*leptons_tkIso_).push_back(0);
    }
    try{
      // storing jets AK4 candidates which are not matching each lepton candidate (highest and second highest pt)
      jetsak4Cand = Matching(electronsVec, jetsVecA, 0.5);
      // storing jets AK8 candidates which are not matching each lepton candidate (highest and second highest pt)
      jetsak8Cand = Matching(electronsVec, jetsVecB, 0.9);
    }catch(...){}
  }

  if(includeMuons_){
    for ( const auto& leptonEvt: muonsVec) {
      (*leptons_energy_).push_back(leptonEvt->energy());
      (*leptons_pt_).push_back(leptonEvt->pt());
      (*leptons_eta_).push_back(leptonEvt->eta());
      (*leptons_phi_).push_back(leptonEvt->phi());
      (*leptons_px_).push_back(leptonEvt->px());
      (*leptons_py_).push_back(leptonEvt->py());
      (*leptons_pz_).push_back(leptonEvt->pz());
      (*leptons_charge_).push_back(leptonEvt->charge());
      (*leptons_vx_).push_back(leptonEvt->vertex().x());
      (*leptons_vy_).push_back(leptonEvt->vertex().y());
      (*leptons_vz_).push_back(leptonEvt->vertex().z());
      (*leptons_looseId_).push_back(0);
      (*leptons_mediumId_).push_back(0);
      (*leptons_tightId_).push_back(0);

      double pfIso = (leptonEvt->pfIsolationR04().sumChargedHadronPt + max(0., leptonEvt->pfIsolationR04().sumNeutralHadronEt + leptonEvt->pfIsolationR04().sumPhotonEt - 0.5*leptonEvt->pfIsolationR04().sumPUPt))/leptonEvt->pt();
      double tkIso = (leptonEvt->isolationR03().sumPt)/(leptonEvt->pt());

      (*leptons_pfIsoMedium_).push_back(leptonEvt->passed(reco::Muon::PFIsoMedium));
      (*leptons_miniIsoTight_).push_back(leptonEvt->passed(reco::Muon::MiniIsoTight));
      (*leptons_pfIsoVeryTight_).push_back(leptonEvt->passed(reco::Muon::PFIsoVeryTight));
      (*leptons_pfIso_).push_back(pfIso);
      (*leptons_tkIso_).push_back(tkIso);
    }
    try{
      // storing jets AK4 candidates which are not matching each lepton candidate (highest and second highest pt)
      jetsak4Cand = Matching(muonsVec, jetsVecA, 0.5);
      // storing jets AK8 candidates which are not matching each lepton candidate (highest and second highest pt)
      jetsak8Cand = Matching(muonsVec, jetsVecB, 0.9);
    }catch(...){}
  }

  // Force pT sorting again...
  std::sort(jetsak4Cand.begin(), jetsak4Cand.end(), orderPt());
  std::sort(jetsak8Cand.begin(), jetsak8Cand.end(), orderPt());

  double NHFID, NEMFID, CHFID, MUFID, CEMFID, NumConstID, CHMID = -999;
  bool LooseJetID, TightJetID, tightLepVetoJetID = false;

  double vz_mean = -999.;
  double pt2_z = 0.;
  double pt2 = 0.;

  for ( const auto& jetsak4Evt: jetsak4Cand) {
    if(debug_) std::cout << "<Candidate> Jets (AK4) pT: " << jetsak4Evt->pt() << " [GeV], eta: " << jetsak4Evt->eta() << ", phi: "<< jetsak4Evt->phi() << std::endl;
    (*jetsak4_pt_).push_back(jetsak4Evt->pt());
    (*jetsak4_energy_).push_back(jetsak4Evt->energy());
    (*jetsak4_phi_).push_back(jetsak4Evt->phi());
    (*jetsak4_eta_).push_back(jetsak4Evt->eta());
    (*jetsak4_bdis_).push_back(jetsak4Evt->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));

    if(jetsak4Evt->isJPTJet() || jetsak4Evt->isPFJet()){
      NHFID  = jetsak4Evt->neutralHadronEnergyFraction();
      NEMFID = jetsak4Evt->neutralEmEnergyFraction();
      CHFID  = jetsak4Evt->chargedHadronEnergyFraction();
      MUFID  = jetsak4Evt->muonEnergyFraction();
      CEMFID = jetsak4Evt->chargedEmEnergyFraction();
      NumConstID = jetsak4Evt->chargedMultiplicity()+jetsak4Evt->neutralMultiplicity();
      CHMID      = jetsak4Evt->chargedMultiplicity();

      LooseJetID = ((NHFID<0.99 && NEMFID<0.99 && NumConstID>1 && CHFID>0 && CHMID>0 && CEMFID<0.99));
      TightJetID = ((NHFID<0.90 && NEMFID<0.90 && NumConstID>1 && CHFID>0 && CHMID>0 && CEMFID<0.99));
      tightLepVetoJetID=((NHFID<0.90 && NEMFID<0.90 && NumConstID>1 && MUFID<0.8 && CHFID>0 && CHMID>0 && CEMFID<0.90));

      (*jetsak4_looseJetId_).push_back(LooseJetID);
      (*jetsak4_tightJetId_).push_back(TightJetID);
      (*jetsak4_lepVetoJetId_).push_back(tightLepVetoJetID);
    }else{
      (*jetsak4_looseJetId_).push_back(false);
      (*jetsak4_tightJetId_).push_back(false);
      (*jetsak4_lepVetoJetId_).push_back(false);
    }

    reco::CompositePtrCandidate::daughters pfconst = jetsak4Evt->daughterPtrVector();
    for (reco::CompositePtrCandidate::daughters::const_iterator itpf = pfconst.begin(); itpf != pfconst.end(); ++itpf) {
      pt2 += (*itpf)->pt()*(*itpf)->pt();
      pt2_z += (*itpf)->pt()*(*itpf)->pt()*(*itpf)->vz();
    }

    if (pt2 > 0.) {
      vz_mean = pt2_z/pt2;
    }

    (*jetsak4_vz_).push_back(vz_mean);

  }

  NHFID = -999;
  NEMFID = -999;
  CHFID = -999;
  MUFID = -999;
  CEMFID = -999;
  NumConstID = -999;
  CHMID = -999;
  LooseJetID = false;
  TightJetID = false; 
  tightLepVetoJetID = false;

  vz_mean = -999.;
  pt2_z = 0.;
  pt2 = 0.;

  for ( const auto& jetsak8Evt: jetsak8Cand) {
    if(debug_) std::cout << "<Candidate> Jets (AK8) pT: " << jetsak8Evt->pt() << " [GeV], eta: " << jetsak8Evt->eta() << ", phi: "<< jetsak8Evt->phi() << std::endl;
    (*jetsak8_pt_).push_back(jetsak8Evt->pt());
    (*jetsak8_energy_).push_back(jetsak8Evt->energy());
    (*jetsak8_phi_).push_back(jetsak8Evt->phi());
    (*jetsak8_eta_).push_back(jetsak8Evt->eta());
    (*jetsak8_bdis_).push_back(jetsak8Evt->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));

    if(jetsak8Evt->isJPTJet() || jetsak8Evt->isPFJet()){
      NHFID  = jetsak8Evt->neutralHadronEnergyFraction();
      NEMFID = jetsak8Evt->neutralEmEnergyFraction();
      CHFID  = jetsak8Evt->chargedHadronEnergyFraction();
      MUFID  = jetsak8Evt->muonEnergyFraction();
      CEMFID = jetsak8Evt->chargedEmEnergyFraction();
      NumConstID = jetsak8Evt->chargedMultiplicity()+jetsak8Evt->neutralMultiplicity();
      CHMID      = jetsak8Evt->chargedMultiplicity();

      LooseJetID = ((NHFID<0.99 && NEMFID<0.99 && NumConstID>1 && CHFID>0 && CHMID>0 && CEMFID<0.99));
      TightJetID = ((NHFID<0.90 && NEMFID<0.90 && NumConstID>1 && CHFID>0 && CHMID>0 && CEMFID<0.99));
      tightLepVetoJetID=((NHFID<0.90 && NEMFID<0.90 && NumConstID>1 && MUFID<0.8 && CHFID>0 && CHMID>0 && CEMFID<0.90));

      (*jetsak8_looseJetId_).push_back(LooseJetID);
      (*jetsak8_tightJetId_).push_back(TightJetID);
      (*jetsak8_lepVetoJetId_).push_back(tightLepVetoJetID);
    }else{
      (*jetsak8_looseJetId_).push_back(false);
      (*jetsak8_tightJetId_).push_back(false);
      (*jetsak8_lepVetoJetId_).push_back(false);
    }

    reco::CompositePtrCandidate::daughters pfconst = jetsak8Evt->daughterPtrVector();
    for (reco::CompositePtrCandidate::daughters::const_iterator itpf = pfconst.begin(); itpf != pfconst.end(); ++itpf) {
      pt2 += (*itpf)->pt()*(*itpf)->pt();
      pt2_z += (*itpf)->pt()*(*itpf)->pt()*(*itpf)->vz();
    }

    if (pt2 > 0.) {
      vz_mean = pt2_z/pt2;
    }

    (*jetsak8_vz_).push_back(vz_mean);

  }

}

// ------------ debugging  ------------
void MissingMassSearches::Debug(){

  bool debug_deep = false;

  for ( const auto& triggerEvt: triggerVec ) {
    std::cout << "Trigger: " << triggerEvt << std::endl;
  }

  for ( const auto& jetsak4Evt: jetsVecA) {
    std::cout << "Jets (AK4) pT: " << jetsak4Evt->pt() << " [GeV], eta: " << jetsak4Evt->eta() << ", phi: "<< jetsak4Evt->phi() << std::endl;
  }

  for ( const auto& jetsak8Evt: jetsVecB) {
    std::cout << "Jets (AK8) pT: " << jetsak8Evt->pt() << " [GeV], eta: " << jetsak8Evt->eta() << ", phi: "<< jetsak8Evt->phi() << std::endl;
  }

  for ( const auto& electronsEvt: electronsVec ) {
    std::cout << "Electrons pT: " << electronsEvt->pt() << " [GeV], eta: " << electronsEvt->eta() << ", phi: "<< electronsEvt->phi() << std::endl;
  }

  for ( const auto& muonsEvt: muonsVec ) {
    std::cout << "Muons pT: " << muonsEvt->pt() << " [GeV], eta: " << muonsEvt->eta() << ", phi: "<< muonsEvt->phi() << std::endl;
  }

  for ( const auto& photonsEvt: photonsVec ) {
    std::cout << "Photons pT: " << photonsEvt->pt() << " [GeV], eta: " << photonsEvt->eta() << ", phi: "<< photonsEvt->phi() << std::endl;
  }

  std::cout << "MET eT: " << met->et() << " [GeV], phi: " << met->phi() << ", significance: " << met->significance() << std::endl; 

  if(debug_deep){
    for ( const auto& pfEvt: pfVec ) {
      std::cout << "PF pT: " << pfEvt->pt() << " [GeV], eta: " << pfEvt->eta() << ", phi: "<< pfEvt->phi() << std::endl;
    }

    for ( const auto& packedEvt: packedVec ) {
      std::cout << "Packed PF pT: " << packedEvt->pt() << " [GeV], eta: " << packedEvt->eta() << ", phi: "<< packedEvt->phi() << ", vertex (x,y,z) [cm]: " << packedEvt->vertex() << std::endl;
    }

    for ( const auto& vtxEvt: vtxVec ) {
      std::cout << "Vertices position (x,y,z) [cm]: " << vtxEvt->x() << ", " << vtxEvt->y() << ", " << vtxEvt->z() << std::endl;
    }
  }

  for ( const auto& protonsEvt: protonsVec ) {
    const CTPPSDetId det_id( protonsEvt->getRPId() );
    std::cout << "Protons (x,y)[mm]: (" << protonsEvt->getX() << ", " << protonsEvt->getY() << "), arm: " << det_id.arm() << ", station: " << det_id.station() << ", pot: " << det_id.rp() << std::endl;
  }

  for ( const auto& pixelsEvt: pixelsVec ) {
    std::cout << "Pixels (x,y)[mm]: (" << pixelsEvt.first->getX0() << ", " << pixelsEvt.first->getY0() << "), Arm: " << pixelsEvt.second.arm() << ", Pot: " << pixelsEvt.second.rp() << ", Station: " << pixelsEvt.second.station() << std::endl;
  }

  for ( const auto& timingEvt: timingVec ) {
    std::cout << "Timing (x,y)[mm]: (" << timingEvt.first->getX0() << ", " << timingEvt.first->getY0() << "), Arm: " << timingEvt.second.arm() << ", Pot: " << timingEvt.second.rp() << ", Station: " << timingEvt.second.station() << std::endl;
  }

  for ( const auto& stripsEvt: stripsVec ) {
    std::cout << "Strips (x,y)[mm]: (" << stripsEvt.first->getX0() << ", " << stripsEvt.first->getY0() << "), Arm: " << stripsEvt.second.arm() << ", Pot: " << stripsEvt.second.rp() << ", Station: " << stripsEvt.second.station() << std::endl;
  }

}

// ------------ method called for each event  ------------
  void
MissingMassSearches::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // Cleaning vectors per event
  cleaning();

  // Trigger fired!
  if(fetchTrigger(iEvent,iSetup)){;

    for (std::size_t i = 0; i != triggerVec.size(); ++i) {
      hltTriggerNamesHisto_->Fill( triggersList_[i].c_str(),triggerVec[i]);
      (*trigger_).push_back(triggerVec[i]);
    }

    *run_ = iEvent.id().run();
    *ev_ = iEvent.id().event();
    *lumiblock_ = iEvent.luminosityBlock();

    // Filling Vectors
    if(includeMuons_) fetchMuons(iEvent, iSetup);
    if(includeElectrons_) fetchElectrons(iEvent, iSetup);
    if(includeJets_) fetchJets(iEvent, iSetup);
    if(includePhotons_) fetchPhotons(iEvent, iSetup);
    if(includeMET_) fetchMET(iEvent, iSetup);
    if(includeVertices_) fetchVertices(iEvent, iSetup);
    if(includePF_) fetchPF(iEvent, iSetup);
    if(includeProtons_) fetchProtons(iEvent, iSetup);

    // Printout for Debugging
    if(debug_) Debug();

    // Event Algorith to tag two leading leptons and two highest pt jets not associated with the leptons!
    fetchEventTagger(iEvent);

    tree_->Fill();
  }

  eventClear();

}

// ------------ method called once each job just before starting event loop  ------------
  void 
MissingMassSearches::beginJob()
{

  edm::Service<TFileService> fs;
  tree_=fs->make<TTree>("analyzer","analyzer");

  // Control Histograms
  TFileDirectory triggerDir = fs->mkdir("TriggerInfo");
  hltTriggerNamesHisto_ = triggerDir.make<TH1F>("HLTTriggerNames","HLTTriggerNames",1,0,1);
  hltTriggerNamesHisto_->SetCanExtend(TH1::kXaxis);

  run_ = new int;
  ev_ = new long int;
  lumiblock_ = new int;

  trigger_ = new std::vector<int>;

  if(includeVertices_){
    vertex_x_ = new std::vector<float>;
    vertex_y_ = new std::vector<float>;
    vertex_z_ = new std::vector<float>;
    vertex_ntrack_ = new std::vector<float>;
    vertex_chi2_ = new std::vector<float>;
    vertex_ndof_ = new std::vector<float>;
  }

  if(includeJets_){
    jetsak4_pt_ = new std::vector<float>;
    jetsak4_energy_ = new std::vector<float>;
    jetsak4_phi_ = new std::vector<float>;
    jetsak4_eta_ = new std::vector<float>;
    jetsak4_vz_ = new std::vector<float>;
    jetsak4_bdis_ = new std::vector<float>;
    jetsak4_looseJetId_ = new std::vector<bool>;
    jetsak4_tightJetId_ = new std::vector<bool>;
    jetsak4_lepVetoJetId_ = new std::vector<bool>;
    jetsak8_pt_ = new std::vector<float>;
    jetsak8_energy_ = new std::vector<float>;
    jetsak8_phi_ = new std::vector<float>;
    jetsak8_eta_ = new std::vector<float>;
    jetsak8_vz_ = new std::vector<float>;
    jetsak8_bdis_ = new std::vector<float>;
    jetsak8_looseJetId_ = new std::vector<bool>;
    jetsak8_tightJetId_ = new std::vector<bool>;
    jetsak8_lepVetoJetId_ = new std::vector<bool>;
  }

  if(includeElectrons_ || includeMuons_){
    leptons_energy_ = new std::vector<float>;
    leptons_pt_ = new std::vector<float>;
    leptons_eta_ = new std::vector<float>;
    leptons_phi_ = new std::vector<float>;
    leptons_px_ = new std::vector<float>;
    leptons_py_ = new std::vector<float>;
    leptons_pz_ = new std::vector<float>;
    leptons_charge_ = new std::vector<float>;
    leptons_vx_ = new std::vector<float>;
    leptons_vy_ = new std::vector<float>;
    leptons_vz_ = new std::vector<float>;
    leptons_looseId_ = new std::vector<bool>;
    leptons_mediumId_ = new std::vector<bool>;
    leptons_tightId_ = new std::vector<bool>;
    leptons_pfIsoMedium_ = new std::vector<bool>;
    leptons_miniIsoTight_ = new std::vector<bool>;
    leptons_pfIsoVeryTight_ = new std::vector<bool>;
    leptons_pfIso_ = new std::vector<float>;
    leptons_tkIso_ = new std::vector<float>;
  }

  if(includeMET_){
    misset_ = new float;
    misset_phi_ = new float;
  }

  if(includePF_){
    nPF_ = new int;
    SumPF_energy_ = new float;
  }

  if(includeProtons_){
    protonsArm_ = new std::vector<int>;
    protonsStation_ = new std::vector<int>;
    protonsRP_ = new std::vector<int>;
    protonsX_ = new std::vector<float>;
    protonsXUnc_ = new std::vector<float>;
    protonsY_ = new std::vector<float>;
    protonsYUnc_ = new std::vector<float>;
    protonsTime_ = new std::vector<float>;
    protonsTimeUnc_ = new std::vector<float>;
  }

  tree_->Branch("run",run_,"run/I");
  tree_->Branch("event",ev_,"event/L");
  tree_->Branch("lumiblock",lumiblock_,"lumiblock/I");

  tree_->Branch("trigger",&trigger_);

  if(includeVertices_){
    tree_->Branch("vertex_x",&vertex_x_);
    tree_->Branch("vertex_y",&vertex_y_);
    tree_->Branch("vertex_z",&vertex_z_);
    tree_->Branch("vertex_ntrack",&vertex_ntrack_);
    tree_->Branch("vertex_chi2",&vertex_chi2_);
    tree_->Branch("vertex_ndof",&vertex_ndof_);
  }

  if(includeJets_){
    tree_->Branch("jetsak4_pt",&jetsak4_pt_);
    tree_->Branch("jetsak4_energy",&jetsak4_energy_);
    tree_->Branch("jetsak4_phi",&jetsak4_phi_);
    tree_->Branch("jetsak4_eta",&jetsak4_eta_);
    tree_->Branch("jetsak4_vz",&jetsak4_vz_);
    tree_->Branch("jetsak4_bdis",&jetsak4_bdis_);
    tree_->Branch("jetsak4_looseId",&jetsak4_looseJetId_);
    tree_->Branch("jetsak4_tightId",&jetsak4_tightJetId_);
    tree_->Branch("jetsak4_lepVeto",&jetsak4_lepVetoJetId_);
    tree_->Branch("jetsak8_pt",&jetsak8_pt_);
    tree_->Branch("jetsak8_energy",&jetsak8_energy_);
    tree_->Branch("jetsak8_phi",&jetsak8_phi_);
    tree_->Branch("jetsak8_eta",&jetsak8_eta_);
    tree_->Branch("jetsak8_vz",&jetsak8_vz_);
    tree_->Branch("jetsak8_bdis",&jetsak8_bdis_);
    tree_->Branch("jetsak8_looseId",&jetsak8_looseJetId_);
    tree_->Branch("jetsak8_tightId",&jetsak8_tightJetId_);
    tree_->Branch("jetsak8_lepVeto",&jetsak8_lepVetoJetId_);
  }

  if(includeElectrons_ || includeMuons_){
    tree_->Branch("leptons_energy",&leptons_energy_);
    tree_->Branch("leptons_pt",&leptons_pt_);
    tree_->Branch("leptons_eta",&leptons_eta_);
    tree_->Branch("leptons_phi",&leptons_phi_);
    tree_->Branch("leptons_px",&leptons_px_);
    tree_->Branch("leptons_py",&leptons_py_);
    tree_->Branch("leptons_pz",&leptons_pz_);
    tree_->Branch("leptons_charge",&leptons_charge_);
    tree_->Branch("leptons_vx",&leptons_vx_);
    tree_->Branch("leptons_vy",&leptons_vy_);
    tree_->Branch("leptons_vz",&leptons_vz_);
    tree_->Branch("leptons_looseId",&leptons_looseId_);
    tree_->Branch("leptons_mediumId",&leptons_mediumId_);
    tree_->Branch("leptons_tightId",&leptons_tightId_);
    tree_->Branch("leptons_pfIsoMedium_",&leptons_pfIsoMedium_);
    tree_->Branch("leptons_miniIsoTight_",&leptons_miniIsoTight_);
    tree_->Branch("leptons_pfIsoVeryTight_",&leptons_pfIsoVeryTight_);
    tree_->Branch("leptons_pfIso_",&leptons_pfIso_);
    tree_->Branch("leptons_tkIso_",&leptons_tkIso_);
  }

  if(includeMET_){
    tree_->Branch("missEt",misset_,"missEt/F");
    tree_->Branch("missEt_phi",misset_phi_,"missEt_phi/F");
  }

  if(includePF_){
    tree_->Branch("nPF",nPF_,"nPF/I");
    tree_->Branch("sumEnergyPF",SumPF_energy_,"SumEnergyPF/F");
  }

  if(includeProtons_){
    tree_->Branch("protonsArm",&protonsArm_);
    tree_->Branch("protonsStation",&protonsStation_);
    tree_->Branch("protonsRP",&protonsRP_);
    tree_->Branch("protonsX",&protonsX_);
    tree_->Branch("protonsXUnc",&protonsXUnc_);
    tree_->Branch("protonsY",&protonsY_);
    tree_->Branch("protonsYUnc",&protonsYUnc_);
    tree_->Branch("protonsTime",&protonsTime_);
    tree_->Branch("protonsTimeUnc",&protonsTimeUnc_);
  }

}

// ------------ method called once each job just after ending the event loop  ------------
  void 
MissingMassSearches::endJob() 
{

  delete run_;
  delete ev_;
  delete lumiblock_;

  delete trigger_;

  if(includeVertices_){
    delete vertex_x_;
    delete vertex_y_;
    delete vertex_z_;
    delete vertex_ntrack_;
    delete vertex_chi2_;
    delete vertex_ndof_;
  }

  if(includeJets_){
    delete jetsak4_pt_;
    delete jetsak4_energy_;
    delete jetsak4_phi_;
    delete jetsak4_eta_;
    delete jetsak4_vz_;
    delete jetsak4_bdis_;
    delete jetsak4_looseJetId_;
    delete jetsak4_tightJetId_;
    delete jetsak4_lepVetoJetId_;
    delete jetsak8_pt_;
    delete jetsak8_energy_;
    delete jetsak8_phi_;
    delete jetsak8_eta_;
    delete jetsak8_vz_;
    delete jetsak8_bdis_;
    delete jetsak8_looseJetId_;
    delete jetsak8_tightJetId_;
    delete jetsak8_lepVetoJetId_;
  }

  if(includeElectrons_ || includeMuons_){
    delete leptons_energy_;
    delete leptons_pt_;
    delete leptons_eta_;
    delete leptons_phi_;
    delete leptons_px_;
    delete leptons_py_;
    delete leptons_pz_;
    delete leptons_charge_;
    delete leptons_vx_;
    delete leptons_vy_;
    delete leptons_vz_;
    delete leptons_looseId_;
    delete leptons_mediumId_;
    delete leptons_tightId_;
    delete leptons_pfIsoMedium_;
    delete leptons_miniIsoTight_;
    delete leptons_pfIsoVeryTight_;
    delete leptons_pfIso_;
    delete leptons_tkIso_;
  }

  if(includeMET_){
    delete misset_;
    delete misset_phi_;
  }

  if(includePF_){
    delete nPF_;
    delete SumPF_energy_;
  }

  if(includeProtons_){
    delete protonsArm_;
    delete protonsStation_;
    delete protonsRP_;
    delete protonsX_;
    delete protonsXUnc_;
    delete protonsY_;
    delete protonsYUnc_;
    delete protonsTime_;
    delete protonsTimeUnc_;
  }

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MissingMassSearches::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


// ------------ templates ------------

// ------------ Invariant Mass
template <class T, class W>
math::XYZTLorentzVector MissingMassSearches::DiSystem(T obj1_, W obj2_){
  math::XYZTLorentzVector DiObj(0.,0.,0.,0.);
  DiObj += obj1_->p4();
  DiObj += obj2_->p4();
  return DiObj;
}

// ------------ Transverse Mass
template <class T, class W>
double MissingMassSearches::TransverseMass(T lepton_, W met_){

  double w_et = met_->et() + lepton_->pt();
  double w_px = met_->px() + lepton_->px();
  double w_py = met_->py() + lepton_->py();

  double tmass_ = w_et*w_et - w_px*w_px - w_py*w_py;
  tmass_ = (tmass_ > 0) ? sqrt(tmass_) : 0;

  return tmass_;

}

// ------------ Disentangle jets/leptons
template <class T, class W>
std::vector<const pat::Jet*> MissingMassSearches::Matching(T lepton_, W jet_, double radius){

  std::vector<const pat::Jet*> jetsCand_;
  jetsCand_.clear();

  // storing jets AK4 candidates which are not matching each lepton candidate (highest and second highest pt)
  try{
    for (std::size_t i = 0; i != jet_.size(); ++i) {
      if(!(lepton_.size()>1)) continue;
      if ((deltaR(lepton_[0]->eta(), lepton_[0]->phi(), jet_[i]->eta(), jet_[i]->phi()) > radius)&&(deltaR(lepton_[1]->eta(), lepton_[1]->phi(), jet_[i]->eta(), jet_[i]->phi())>radius)){
	jetsCand_.push_back(jet_[i]);
      }
    }
  }catch(...){}

  return jetsCand_;

}

//define this as a plug-in
DEFINE_FWK_MODULE(MissingMassSearches);
