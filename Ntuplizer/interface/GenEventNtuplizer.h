#ifndef GenEventNtuplizer_H
#define GenEventNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class GenEventNtuplizer : public CandidateNtuplizer {

public:
  GenEventNtuplizer( std::vector< edm::EDGetTokenT< GenEventInfoProduct > > tokens, NtupleBranches* nBranches );
  ~GenEventNtuplizer( void );
  
  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
  
private:
   edm::EDGetTokenT< GenEventInfoProduct > geneventToken_; 
     
   edm::Handle< GenEventInfoProduct >  geneventInfo_;
      
};

#endif // GenEventNtuplizer_H
