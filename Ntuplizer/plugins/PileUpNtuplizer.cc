#include "../interface/PileUpNtuplizer.h"

//===================================================================================================================
PileUpNtuplizer::PileUpNtuplizer( std::vector< edm::EDGetTokenT< std::vector<PileupSummaryInfo> > > tokens, NtupleBranches* nBranches )
   : CandidateNtuplizer( nBranches )
   , pileUpToken_( tokens[0] )
{

}

//===================================================================================================================
PileUpNtuplizer::~PileUpNtuplizer( void )
{

}

//===================================================================================================================
void PileUpNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  event.getByToken(pileUpToken_, pileUpInfo_);  
  
  std::vector<PileupSummaryInfo>::const_iterator PVI;
  
  for( PVI = pileUpInfo_->begin(); PVI != pileUpInfo_->end(); ++PVI ) {
     nBranches_->nPuVtxTrue.push_back(PVI->getTrueNumInteractions());
     nBranches_->nPuVtx.push_back(PVI->getPU_NumInteractions());
     nBranches_->bX.push_back(PVI->getBunchCrossing());
  }

}
