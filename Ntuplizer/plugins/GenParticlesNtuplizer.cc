#include "../interface/GenParticlesNtuplizer.h"
 
//===================================================================================================================        
GenParticlesNtuplizer::GenParticlesNtuplizer( std::vector<edm::EDGetTokenT<reco::GenParticleCollection>> tokens, NtupleBranches* nBranches ) 
   : CandidateNtuplizer( nBranches )
   , genParticlesToken_( tokens[0] )
{

}

//===================================================================================================================        
GenParticlesNtuplizer::~GenParticlesNtuplizer( void )
{
}

//===================================================================================================================        
  void GenParticlesNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
  
    event.getByToken(genParticlesToken_ , genParticles_); 

   /* here we want to save  gen particles info*/
   
    std::vector<int> vDau ;
    std::vector<int> vMoth;
    int nMoth = 0;
    int nDau  = 0;  
    nBranches_->genParticle_N = genParticles_->size();
    for( unsigned p=0; p<genParticles_->size(); ++p ){
      //if( (*genParticles_)[p].status() != 3 ) continue;
      vDau.clear(); vMoth.clear();
      nDau = 0; nMoth = 0;
      nBranches_->genParticle_pt    .push_back((*genParticles_)[p].pt()     );
      nBranches_->genParticle_px    .push_back((*genParticles_)[p].px()     );
      nBranches_->genParticle_py    .push_back((*genParticles_)[p].py()     );
      nBranches_->genParticle_pz    .push_back((*genParticles_)[p].pz()     );
      nBranches_->genParticle_eta   .push_back((*genParticles_)[p].eta()    );
      nBranches_->genParticle_mass  .push_back((*genParticles_)[p].mass()   );
      nBranches_->genParticle_phi   .push_back((*genParticles_)[p].phi()    );
      nBranches_->genParticle_e     .push_back((*genParticles_)[p].energy() );
      nBranches_->genParticle_status.push_back((*genParticles_)[p].status() );
      nBranches_->genParticle_pdgId .push_back((*genParticles_)[p].pdgId()  );


      // needed for the gen matching
      nBranches_->genParticle_isPrompt.push_back((*genParticles_)[p].statusFlags().isPrompt());
      nBranches_->genParticle_isDirectPromptTauDecayProduct.push_back((*genParticles_)[p].statusFlags().isDirectPromptTauDecayProduct());

      // needed for the MVA recoil correction
      //      nBranches_->genParticle_fromHardProcessFinalState.push_back((*genParticles_)[p].statusFlags().fromHardProcessFinalState());
      //      nBranches_->genParticle_isDirectHardProcessTauDecayProductFinalState.push_back((*genParticles_)[p].statusFlags().isDirectHardProcessTauDecayProductFinalState());
      nBranches_->genParticle_fromHardProcessFinalState.push_back((*genParticles_)[p].fromHardProcessFinalState());
      nBranches_->genParticle_isDirectHardProcessTauDecayProductFinalState.push_back((*genParticles_)[p].isDirectHardProcessTauDecayProductFinalState());


      for( unsigned int d=0; d<(*genParticles_)[p].numberOfDaughters(); ++d ){
        vDau.push_back( (*genParticles_)[p].daughter(d)->pdgId() );
	    nDau++;
      }
      for( unsigned int m=0; m<(*genParticles_)[p].numberOfMothers(); ++m ){
        vMoth.push_back( (*genParticles_)[p].mother(m)->pdgId() );
    	nMoth++;
      }
      nBranches_->genParticle_nDau  .push_back( nDau  );
      nBranches_->genParticle_nMoth .push_back( nMoth );      
      nBranches_->genParticle_mother.push_back( vMoth );
      nBranches_->genParticle_dau   .push_back( vDau  );      
    }

    
}

