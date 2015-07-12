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
    
  nBranches_->PV_N = vertices_->size();
  
  nBranches_->PV_filter = true;
  
  reco::VertexCollection::const_iterator firstGoodVertex = vertices_->end();
  int firstGoodVertexIdx = 0;
  for( reco::VertexCollection::const_iterator vtx = vertices_->begin(); vtx != vertices_->end(); ++vtx, ++firstGoodVertexIdx){
  
    nBranches_->PV_chi2.push_back(vtx->chi2());
    nBranches_->PV_ndof.push_back(vtx->ndof());
    nBranches_->PV_rho.push_back(vtx->position().Rho());
    nBranches_->PV_z.push_back(vtx->position().Z());
    
    bool isFake = (vtx->chi2()==0 && vtx->ndof()==0);
    if( !isFake && vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0) {
      firstGoodVertex = vtx;
      break;
    }
    
  }
  
  if ( firstGoodVertex==vertices_->end() ) nBranches_->PV_filter = false;
  
}
