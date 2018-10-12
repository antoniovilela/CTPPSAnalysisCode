// -*- C++ -*-
//
// Package:    Filter/ProtonTagFilter
// Class:      ProtonTagFilter
// 
/**\class ProtonTagFilter ProtonTagFilter.cc Filter/ProtonTagFilter/plugins/ProtonTagFilter.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Diego Figueiredo
//         Created:  Wed, 10 Oct 2018 00:10:06 GMT
//
//

#include "ProtonTagFilter.h"

ProtonTagFilter::ProtonTagFilter(const edm::ParameterSet& iConfig):
  debug_              ( iConfig.getParameter<bool>( "debugging" ) ),
  protonToken_        ( consumes<std::vector<CTPPSLocalTrackLite>>                  ( iConfig.getParameter<edm::InputTag>( "protonTag" ) ) ),
  pixelsppsToken_     ( consumes<edm::DetSetVector<CTPPSPixelLocalTrack>>           ( iConfig.getParameter<edm::InputTag>( "pixelsppsTag" ) ) ),
  timingppsToken_     ( consumes<edm::DetSetVector<CTPPSDiamondLocalTrack>>         ( iConfig.getParameter<edm::InputTag>( "timingppsTag" ) ) ),
  stripstotemToken_   ( consumes<edm::DetSetVector<TotemRPLocalTrack>>              ( iConfig.getParameter<edm::InputTag>( "stripstotemTag" ) ) )
{
  //now do what ever initialization is needed

}


ProtonTagFilter::~ProtonTagFilter()
{

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
  bool
ProtonTagFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;

  ProtonsEvent *protons;
  protons = new ProtonsEvent(iEvent, iSetup, protonToken_);
  protonsVec = protons->GetProtonsLite();

  if(protonsVec.empty()){
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
  }

  bool tracking_sec45 = false;
  bool tracking_sec56 = false;

  bool timing_sec45 = false;
  bool timing_sec56 = false;

  bool filter_event = false;

  /**
   * \brief Detectors Mapping
   *   => arm: 0 (sector 45), 1 (sector 56)
   *   => station: 0 (210m), 1 (cylyndrical pots), 2 (220m)
   *   => Roman Pot: 0 (near top), 1 (near bottom), 2 (near horizontal), 3 (far horizontal), 4 (far top), 5 (far bottom)
   *
   *   => pixels  arm 0/1, station 0/2 and pot 3
   *   => strips: arm 0/1, station 0/2 and pot 4/5
   *   => timing: arm 0/1, station 1 and pot 6.
   **/

  for ( const auto& protonsEvt: protonsVec ) {
    const CTPPSDetId det_id( protonsEvt->getRPId() );

    // Pixels at SEC45 
    if(det_id.arm()==0 && (det_id.station()==0 || det_id.station()==2) && det_id.rp()==3) tracking_sec45 = true;

    // Pixels at SEC56
    if(det_id.arm()==1 && (det_id.station()==0 || det_id.station()==2) && det_id.rp()==3) tracking_sec56 = true;

    // Strips at SEC45 
    if(det_id.arm()==0 && (det_id.station()==0 || det_id.station()==2) && (det_id.rp()==4 || det_id.rp()==5)) tracking_sec45 = true;

    // Strips at SEC56
    if(det_id.arm()==1 && (det_id.station()==0 || det_id.station()==2) && (det_id.rp()==4 || det_id.rp()==5)) tracking_sec56 = true;

    // Timing at SEC45 
    if(det_id.arm()==0 && det_id.station()==1 && det_id.rp()==6) timing_sec45 = true;

    // Timing at SEC56 
    if(det_id.arm()==1 && det_id.station()==1 && det_id.rp()==6) timing_sec56 = true;

    if(debug_){
      if(timing_sec45 || timing_sec56 || tracking_sec45 || tracking_sec56){
	std::cout << "Filter, proton x[mm]: "<< protonsEvt->getX() << ", arm: " << det_id.arm() << ", station: " << det_id.station() << ", pot: " << det_id.rp() << std::endl;
      }
    }

  }

  if(protonsVec.empty()){

    tracking_sec45 = false;
    tracking_sec56 = false;

    timing_sec45 = false;
    timing_sec56 = false;

    for ( const auto& pixelsEvt: pixelsVec ) {

      // Pixels at SEC45 
      if(pixelsEvt.second.arm()==0 && (pixelsEvt.second.station()==0 || pixelsEvt.second.station()==2) && pixelsEvt.second.rp()==3) tracking_sec45 = true;

      // Pixels at SEC56
      if(pixelsEvt.second.arm()==1 && (pixelsEvt.second.station()==0 || pixelsEvt.second.station()==2) && pixelsEvt.second.rp()==3) tracking_sec56 = true;

      if(debug_){
	if(tracking_sec45 || tracking_sec56){
	  std::cout << "Filter, Pixels x[mm]: " << pixelsEvt.first->getX0() << ", Arm: " << pixelsEvt.second.arm() << ", Pot: " << pixelsEvt.second.rp() << ", Station: " << pixelsEvt.second.station() << std::endl;
	}
      }

    }

    for ( const auto& timingEvt: timingVec ) {

      // Timing at SEC45 
      if(timingEvt.second.arm()==0 && timingEvt.second.station()==1 && timingEvt.second.rp()==6) timing_sec45 = true;

      // Timing at SEC56 
      if(timingEvt.second.arm()==1 && timingEvt.second.station()==1 && timingEvt.second.rp()==6) timing_sec56 = true;

      if(debug_){
	if(timing_sec45 || timing_sec56){
	  std::cout << "Filter, Timing x[mm]: " << timingEvt.first->getX0() << ", Arm: " << timingEvt.second.arm() << ", Pot: " << timingEvt.second.rp() << ", Station: " << timingEvt.second.station() << std::endl;
	}
      }

    }

    for ( const auto& stripsEvt: stripsVec ) {

      // Strips at SEC45 
      if(stripsEvt.second.arm()==0 && (stripsEvt.second.station()==0 || stripsEvt.second.station()==2) && (stripsEvt.second.rp()==4 || stripsEvt.second.rp()==5)) tracking_sec45 = true;

      // Strips at SEC56
      if(stripsEvt.second.arm()==1 && (stripsEvt.second.station()==0 || stripsEvt.second.station()==2) && (stripsEvt.second.rp()==4 || stripsEvt.second.rp()==5)) tracking_sec56 = true;

      if(debug_){
	if(tracking_sec45 || tracking_sec56){
	  std::cout << "Filter, Strips x[mm]: " << stripsEvt.first->getX0() << ", Arm: " << stripsEvt.second.arm() << ", Pot: " << stripsEvt.second.rp() << ", Station: " << stripsEvt.second.station() << std::endl;
	}
      }

    }
  }

  if(tracking_sec45 && tracking_sec56) filter_event=true;
  if(filter_event && debug_) std::cout << "\n\nEvent Selected!\n\n" << std::endl;

  return filter_event;

}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
  void
ProtonTagFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
ProtonTagFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
   void
   ProtonTagFilter::beginRun(edm::Run const&, edm::EventSetup const&)
   { 
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void
   ProtonTagFilter::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void
   ProtonTagFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void
   ProtonTagFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ProtonTagFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ProtonTagFilter);
