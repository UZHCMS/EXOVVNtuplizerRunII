#ifndef METsNtuplizer_H
#define METsNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class TFormula;

class METsNtuplizer : public CandidateNtuplizer {

public:
   METsNtuplizer( std::vector<edm::EDGetTokenT<pat::METCollection>> tokens, NtupleBranches* nBranches );
   ~METsNtuplizer( void );
   
   void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
   
private:      
   //edm::InputTag METsRawLabel_;
   edm::EDGetTokenT<pat::METCollection> metInputToken_   ; 
     
   //edm::Handle< std::vector<pat::MET> > METsRaw_;
   edm::Handle<pat::METCollection> METs_	;
  
};

#endif // METsNtuplizer_H
