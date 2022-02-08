#ifndef CandidateNtuplizer_H
#define CandidateNtuplizer_H

#include "../interface/NtupleBranches.h"

class CandidateNtuplizer{

public:
   CandidateNtuplizer( NtupleBranches* nBranches );
   virtual ~CandidateNtuplizer( void );

   virtual bool fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){return true;};
   
protected:
   NtupleBranches* nBranches_;
      
};

#endif // CandidateNtuplizer_H
