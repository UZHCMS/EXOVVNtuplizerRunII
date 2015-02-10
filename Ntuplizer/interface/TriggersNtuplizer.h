#ifndef TriggersNtuplizer_H
#define TriggersNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class TriggersNtuplizer : public CandidateNtuplizer {

public:
   TriggersNtuplizer( std::vector<edm::EDGetTokenT<edm::TriggerResults>> tokens, NtupleBranches* nBranches );
   ~TriggersNtuplizer( void );
   
  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );

private:
   edm::EDGetTokenT<edm::TriggerResults	> HLTtriggersToken_;
   edm::Handle< edm::TriggerResults 	> HLTtriggers_;
      
};

#endif // TriggersNtuplizer_H



