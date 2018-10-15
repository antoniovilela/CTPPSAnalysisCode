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


// system include files
#include <memory>
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TAxis.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Event
#include "TriggerEvent.h"
#include "JetsEvent.h"
#include "ElectronsEvent.h"
#include "MuonsEvent.h"
#include "PhotonsEvent.h"
#include "METEvent.h"
#include "ParticleFlowEvent.h"
#include "VerticesEvent.h"
#include "ProtonsEvent.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MissingMassSearches : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit MissingMassSearches(const edm::ParameterSet&);
    ~MissingMassSearches();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() override;
    virtual void cleaning();
    virtual void eventClear();
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // Creating CMS Physics Object
    void fetchMuons(const edm::Event&, const edm::EventSetup&);
    void fetchElectrons(const edm::Event&, const edm::EventSetup&);
    void fetchJets(const edm::Event&, const edm::EventSetup&);
    void fetchPhotons(const edm::Event&, const edm::EventSetup&);
    void fetchMET(const edm::Event&, const edm::EventSetup&);
    void fetchVertices(const edm::Event&, const edm::EventSetup&);
    void fetchPF(const edm::Event&, const edm::EventSetup&);
    void fetchProtons(const edm::Event&, const edm::EventSetup&);
    bool fetchTrigger(const edm::Event&, const edm::EventSetup&);
    void fetchEventTagger(const edm::Event&);
    void Debug();

    TTree * tree_;

    template <class T, class W>
      math::XYZTLorentzVector DiSystem(T obj1, W obj2);

    template <class T, class W>
      double TransverseMass(T lepton1, W lepton2);

    template <class T, class W>
      std::vector<const pat::Jet*> Matching(T lepton, W jet, double radius);

    // ----------member data ---------------------------

    // Switches
    bool debug_;
    std::string mode_;
    bool includeMuons_, includeElectrons_;
    bool includeJets_;
    bool includePhotons_;
    bool includeMET_;
    bool includeVertices_;
    bool includePF_;     
    bool includeProtons_;

    // Trigger
    std::vector<std::string> triggersList_;
    edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;

    // Jets
    edm::EDGetTokenT<edm::View<pat::Jet> > jetTokenA_;
    edm::EDGetTokenT<edm::View<pat::Jet> > jetTokenB_;

    // Electrons
    edm::EDGetTokenT<edm::View<pat::Electron> > eleToken_;

    // Muons
    edm::EDGetTokenT<edm::View<pat::Muon> > muonToken_;

    // PF
    edm::EDGetTokenT<std::vector< reco::PFCandidate > > pfToken_;
    edm::EDGetTokenT<std::vector< pat::PackedCandidate > > packedToken_;

    // Photons
    edm::EDGetTokenT<edm::View<pat::Photon> > photonToken_;

    // MET
    edm::EDGetTokenT<edm::View<pat::MET> > metToken_;

    // Vertex
    edm::EDGetTokenT<edm::View<reco::Vertex> > vertexToken_;

    // Protons Lite
    edm::EDGetTokenT<vector<CTPPSLocalTrackLite> > protonToken_;

    // PPS/TOTEM Detectors
    edm::EDGetTokenT<edm::DetSetVector<CTPPSPixelLocalTrack> > pixelsppsToken_;
    edm::EDGetTokenT<edm::DetSetVector<CTPPSDiamondLocalTrack> > timingppsToken_;
    edm::EDGetTokenT<edm::DetSetVector<TotemRPLocalTrack> > stripstotemToken_;

    // Vectors with Collections
    std::vector<int> triggerVec;
    std::vector<const pat::Jet*> jetsVecA;
    std::vector<const pat::Jet*> jetsVecB;
    std::vector<const pat::Electron*> electronsVec;
    std::vector<const pat::Muon*> muonsVec;
    std::vector<const reco::PFCandidate*> pfVec;
    std::vector<const pat::PackedCandidate*> packedVec;
    std::vector<const pat::Photon*> photonsVec;
    std::vector<const reco::Vertex*> vtxVec;
    std::vector<const CTPPSLocalTrackLite*> protonsVec;
    std::vector<std::pair<const CTPPSPixelLocalTrack*, const CTPPSDetId>> pixelsVec;
    std::vector<std::pair<const CTPPSDiamondLocalTrack*, const CTPPSDetId>> timingVec;
    std::vector<std::pair<const TotemRPLocalTrack*, const CTPPSDetId>> stripsVec;

    std::vector<const pat::Jet*> jetsak4Cand;
    std::vector<const pat::Jet*> jetsak8Cand;

    edm::View<pat::MET>::const_iterator met;

    std::vector<std::pair<const pat::Muon*, const pat::Jet*>> MuonsJetsCand;
    std::vector<std::pair<const pat::Electron*, const pat::Jet*>> ElectronsJetsCand;

    // TTree Vectors

    long int *ev_;
    int *run_;
    int *lumiblock_;

    std::vector<int> *trigger_;

    std::vector<float> *vertex_x_;
    std::vector<float> *vertex_y_;
    std::vector<float> *vertex_z_;
    std::vector<float> *vertex_ntrack_;
    std::vector<float> *vertex_chi2_;
    std::vector<float> *vertex_ndof_;

    std::vector<float> *jetsak4_pt_;
    std::vector<float> *jetsak4_energy_;
    std::vector<float> *jetsak4_phi_;
    std::vector<float> *jetsak4_eta_;
    std::vector<float> *jetsak4_vz_;
    std::vector<float> *jetsak4_bdis_;
    std::vector<bool> *jetsak4_looseJetId_;
    std::vector<bool> *jetsak4_tightJetId_;
    std::vector<bool> *jetsak4_lepVetoJetId_;

    std::vector<float> *jetsak8_pt_;
    std::vector<float> *jetsak8_energy_;
    std::vector<float> *jetsak8_phi_;
    std::vector<float> *jetsak8_eta_;
    std::vector<float> *jetsak8_vz_;
    std::vector<float> *jetsak8_bdis_;
    std::vector<bool> *jetsak8_looseJetId_;
    std::vector<bool> *jetsak8_tightJetId_;
    std::vector<bool> *jetsak8_lepVetoJetId_;

    std::vector<float> *leptons_energy_;
    std::vector<float> *leptons_pt_;
    std::vector<float> *leptons_eta_;
    std::vector<float> *leptons_phi_;
    std::vector<float> *leptons_px_;
    std::vector<float> *leptons_py_;
    std::vector<float> *leptons_pz_;
    std::vector<float> *leptons_charge_;
    std::vector<float> *leptons_vx_;
    std::vector<float> *leptons_vy_;
    std::vector<float> *leptons_vz_;
    std::vector<bool> *leptons_looseId_;
    std::vector<bool> *leptons_mediumId_;
    std::vector<bool> *leptons_tightId_;
    std::vector<bool> *leptons_pfIsoMedium_;
    std::vector<bool> *leptons_miniIsoTight_;
    std::vector<bool> *leptons_pfIsoVeryTight_;
    std::vector<float> *leptons_pfIso_;
    std::vector<float> *leptons_tkIso_;

    float *misset_;
    float *misset_phi_;

    int *nPF_;
    float *SumPF_energy_;

    std::vector<int> *protonsArm_;
    std::vector<int> *protonsStation_;
    std::vector<int> *protonsRP_;
    std::vector<float> *protonsX_;
    std::vector<float> *protonsXUnc_;
    std::vector<float> *protonsY_;
    std::vector<float> *protonsYUnc_;
    std::vector<float> *protonsTime_;
    std::vector<float> *protonsTimeUnc_;

    // E/gamma identification
    edm::ParameterSet eleIdLabelSet_;
    std::string eleMediumIdLabel_, eleTightIdLabel_;

    // Cross-check histograms
    TH1F *hltTriggerNamesHisto_;

    // To be used with i.e. std::sort(v.begin(), v.end(), orderPT()), vector will be organized by pt.
    struct orderPt
    {
      template <class T, class W>
	inline bool operator() (T vec1, W vec2)
	{
	  return (vec1->pt() > vec2->pt());
	}
    };

    struct orderEta
    {
      template <class T, class W>
	inline bool operator() (T vec1, W vec2)
	{
	  return (vec1->eta() > vec2->eta());
	}
    };

    struct orderAbsolutPz
    {
      template <class T, class W>
	inline bool operator() (T vec1, W vec2)
	{
	  return (fabs(vec1->pz()) > fabs(vec2->pz()));
	}
    };

    struct orderVz
    {
      template <class T, class W>
	inline bool operator() (T vec1, W vec2)
	{
	  return (fabs(vec1->z()) > fabs(vec2->z()));
	}
    };

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

