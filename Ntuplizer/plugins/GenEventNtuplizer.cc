#include "../interface/GenEventNtuplizer.h"

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
  nBranches_->PDF_x.push_back((geneventInfo_->pdf()->x).first);
  nBranches_->PDF_x.push_back((geneventInfo_->pdf()->x).second);
  nBranches_->PDF_xPDF.push_back((geneventInfo_->pdf()->xPDF).first);
  nBranches_->PDF_xPDF.push_back((geneventInfo_->pdf()->xPDF).second);
  nBranches_->PDF_id.push_back((geneventInfo_->pdf()->id).first);
  nBranches_->PDF_id.push_back((geneventInfo_->pdf()->id).second);
  
}
