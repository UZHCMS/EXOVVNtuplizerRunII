#ifndef VerticesNtuplizer_H
#define VerticesNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class VerticesNtuplizer : public CandidateNtuplizer {

public:
  VerticesNtuplizer( std::vector<edm::EDGetTokenT<reco::VertexCollection>> tokens, NtupleBranches* nBranches );
  ~VerticesNtuplizer( void );
  
  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
  
private:
   edm::EDGetTokenT<reco::VertexCollection> vtxToken_   ;
   
   edm::Handle< reco::VertexCollection >  vertices_;
      
};

#endif // VerticesNtuplizer_H
