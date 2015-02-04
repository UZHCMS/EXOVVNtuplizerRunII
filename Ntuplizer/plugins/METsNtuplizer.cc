#include "../interface/METsNtuplizer.h"

#include <TFormula.h>

//===================================================================================================================        
METsNtuplizer::METsNtuplizer( std::vector<edm::EDGetTokenT<pat::METCollection>> tokens, NtupleBranches* nBranches )
   : CandidateNtuplizer( nBranches ) 
   , metInputToken_( tokens[0] )
   //, METsRawLabel_( labels[0] )
{
}

//===================================================================================================================
METsNtuplizer::~METsNtuplizer( void )
{
}

//===================================================================================================================
void METsNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

    //event.getByLabel(METsRawLabel_, METsRaw_);  
    event.getByToken(metInputToken_, METs_ ); 

    /*for( unsigned int m = 0; m < METsRaw_->size(); ++m){   	   
      nBranches_->METraw_et .push_back((*METsRaw_)[m].et() );	    
      nBranches_->METraw_phi.push_back((*METsRaw_)[m].phi());
    }*/
    
    for( unsigned int m = 0; m < METs_->size(); ++m){
      nBranches_->MET_et .push_back((*METs_)[m].et());    
      nBranches_->MET_phi.push_back((*METs_)[m].phi());
    }    
   
}

