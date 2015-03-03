#ifndef TriggersNtuplizer_H
#define TriggersNtuplizer_H

#include "../interface/CandidateNtuplizer.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

class TriggersNtuplizer : public CandidateNtuplizer {

public:
   TriggersNtuplizer(edm::EDGetTokenT<edm::TriggerResults> tokens, edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> object, edm::EDGetTokenT<pat::PackedTriggerPrescales> prescale,  NtupleBranches* nBranches );
   ~TriggersNtuplizer( void );
   
  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );

private:
   edm::EDGetTokenT<edm::TriggerResults> 							HLTtriggersToken_;
   edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
   edm::EDGetTokenT<pat::PackedTriggerPrescales>            triggerPrescales_;
	
   edm::Handle< edm::TriggerResults> 								HLTtriggers_;
	edm::Handle<pat::TriggerObjectStandAloneCollection> 		triggerObjects;
	edm::Handle<pat::PackedTriggerPrescales> 						triggerPrescales;
      
};

#endif // TriggersNtuplizer_H



