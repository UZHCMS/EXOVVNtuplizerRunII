#include "../interface/VerticesNtuplizer.h"

#include <cmath>

//===================================================================================================================
VerticesNtuplizer::VerticesNtuplizer( std::vector<edm::EDGetTokenT<reco::VertexCollection>> tokens, NtupleBranches* nBranches )
   : CandidateNtuplizer( nBranches )
   , vtxToken_( tokens[0] )
{

}

//===================================================================================================================
VerticesNtuplizer::~VerticesNtuplizer( void )
{

}

//===================================================================================================================
void VerticesNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  event.getByToken(vtxToken_, vertices_);
    
  nBranches_->nPVs = vertices_->size();
  
}
