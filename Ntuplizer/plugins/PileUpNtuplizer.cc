#include "../interface/PileUpNtuplizer.h"

//===================================================================================================================
PileUpNtuplizer::PileUpNtuplizer( std::vector< edm::EDGetTokenT< std::vector<PileupSummaryInfo> > > tokens, NtupleBranches* nBranches, std::map< std::string, bool >& runFlags )
   : CandidateNtuplizer( nBranches )
   , pileUpToken_( tokens[0] )
   , isJpsiMu_( runFlags["doJpsiMu"]  )
   , isJpsiEle_( runFlags["doJpsiEle"]  )
   , isJpsiTau_( runFlags["doJpsiTau"]  )
{

}

//===================================================================================================================
PileUpNtuplizer::~PileUpNtuplizer( void )
{

}

//===================================================================================================================
bool PileUpNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
 
    //Skip events with no jspi if that analysis is chosen
   
//    bool rflag = false;
//    
//    if(isJpsiEle_ || isJpsiMu_){
//      if ( nBranches_->JpsiMu_B_pt.size() >=1) rflag = true;
//    }
//    
//    if(isJpsiTau_){
//      if ( nBranches_->JpsiTau_B_pt.size() >=1)  rflag = true;
//    }
//
//    if(rflag==false) return;


    event.getByToken(pileUpToken_, pileUpInfo_);  



    std::vector<PileupSummaryInfo>::const_iterator PVI;

    //    int ii=0;
    for( PVI = pileUpInfo_->begin(); PVI != pileUpInfo_->end(); ++PVI ) {
      //      std::cout << ii << " bx = " <<  PVI->getBunchCrossing() << ", nputrue = " << PVI->getTrueNumInteractions() << std::endl;
      if(PVI->getBunchCrossing() == 0) { 
	//	std::cout << "-----> enter!!!! : " << ii << std::endl;
        nBranches_->nPuVtxTrue.push_back(PVI->getTrueNumInteractions());
        nBranches_->nPuVtx.push_back(PVI->getPU_NumInteractions());
        nBranches_->bX.push_back(PVI->getBunchCrossing());
      }
      //      ii+=1;
    }

    return true;

}
