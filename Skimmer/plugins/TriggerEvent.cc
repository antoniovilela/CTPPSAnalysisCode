// -*- C++ -*-
//
// Package:    CTPPSAnalysisCode/TriggerEvent
// Class:      TriggerEvent
// 
/**\class TriggerEvent TriggerEvent.cc CTPPSAnalysisCode/TriggerEvent/plugins/TriggerEvent.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Diego Figueiredo
//         Created:  Fri, 05 Oct 2018 09:08:27 GMT
//
//

#include "TriggerEvent.h"

// Retrieve Trigger
std::vector<int> TriggerEvent::GetTrigger(){
  return triggerlist;
}

// Class Definition, loop over all trigger events and fill vector.
TriggerEvent::TriggerEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<edm::TriggerResults>& triggerResultsToken_, std::vector<std::string>& triggersList_){

  triggerlist.clear();

  try{
    edm::Handle<edm::TriggerResults> hltResults;
    iEvent.getByToken( triggerResultsToken_, hltResults);

    if( hltResults.isValid() ){

      unsigned int nSize = hltResults->size();
      const edm::TriggerNames& triggerNames = iEvent.triggerNames(*hltResults);

      size_t idxpath = 0;
      std::vector<std::string>::const_iterator hltpath = triggersList_.begin();
      std::vector<std::string>::const_iterator hltpaths_end = triggersList_.end();

      for(; hltpath != hltpaths_end; ++hltpath,++idxpath){

	std::string resolvedPathName;
	if( edm::is_glob( *hltpath ) ){
	  std::vector< std::vector<std::string>::const_iterator > matches = edm::regexMatch(triggerNames.triggerNames(), *hltpath);
	  if( matches.empty() ){
	    LogDebug("Configuration") << "Could not find trigger " << *hltpath << " in the path list.\n";
	  } 
	  else if( matches.size() > 1)
	    throw cms::Exception("Configuration") << "HLT path type " << *hltpath << " not unique\n";
	  else resolvedPathName = *(matches[0]);
	} else{
	  resolvedPathName = *hltpath;
	} 

	unsigned int idx_HLT = triggerNames.triggerIndex(resolvedPathName);

	if (idx_HLT < nSize){
	  int accept_HLT = ( hltResults->wasrun(idx_HLT) && hltResults->accept(idx_HLT) ) ? 1 : 0;
	  triggerlist.push_back(accept_HLT);
	}else{
	  triggerlist.push_back(-1);
	} 
      } 

    }else{
      edm::LogWarning("Trigger Event") << "trigger is not valid.\n";
    } 
  }catch(...){
    triggerlist.push_back(-1);
  }

}
