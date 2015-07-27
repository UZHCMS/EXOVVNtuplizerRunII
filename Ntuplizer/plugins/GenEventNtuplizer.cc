#include "../interface/GenEventNtuplizer.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

//===================================================================================================================
GenEventNtuplizer::GenEventNtuplizer( std::vector< edm::EDGetTokenT< GenEventInfoProduct > > tokens, NtupleBranches* nBranches )
   : CandidateNtuplizer( nBranches )
   , geneventToken_( tokens[0] )
{

}

//===================================================================================================================
GenEventNtuplizer::~GenEventNtuplizer( void )
{

}

//===================================================================================================================
void GenEventNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  event.getByToken(geneventToken_, geneventInfo_);  
  
  nBranches_->genWeight=geneventInfo_->weight();
  nBranches_->qScale=geneventInfo_->qScale();

}
