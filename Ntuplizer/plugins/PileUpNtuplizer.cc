#include "../interface/PileUpNtuplizer.h"

//===================================================================================================================
PileUpNtuplizer::PileUpNtuplizer( std::vector< edm::EDGetTokenT< std::vector<PileupSummaryInfo> > > tokens, NtupleBranches* nBranches, std::map< std::string, bool >& runFlags )
   : CandidateNtuplizer( nBranches )
   , pileUpToken_( tokens[0] )
   , isJpsiMu_( runFlags["doJpsiMu"]  )
   , isJpsiEle_( runFlags["doJpsiEle"]  )
{

}

//===================================================================================================================
PileUpNtuplizer::~PileUpNtuplizer( void )
{

}

//===================================================================================================================
void PileUpNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
  //chunk to remove those events with no jspi if that analysis is chosen
  std::vector<int> doJpsi_;
  if(isJpsiEle_) {
    doJpsi_ = nBranches_->IsJpsiEle;
    //std::cout<<"im getting inside the electron part"<<std::endl;
  }else if(isJpsiMu_){
    doJpsi_ = nBranches_->IsJpsiMu;
    //std::cout<<"nbranch thing\t"<<size(isJpsi_)<<"; "<< isJpsi_[0]<<std::endl;
  }
  if(size(doJpsi_)>0) if(doJpsi_[0]==0) return;


  event.getByToken(pileUpToken_, pileUpInfo_);  



  std::vector<PileupSummaryInfo>::const_iterator PVI;
  
  for( PVI = pileUpInfo_->begin(); PVI != pileUpInfo_->end(); ++PVI ) {
     nBranches_->nPuVtxTrue.push_back(PVI->getTrueNumInteractions());
     nBranches_->nPuVtx.push_back(PVI->getPU_NumInteractions());
     nBranches_->bX.push_back(PVI->getBunchCrossing());
  }

}
