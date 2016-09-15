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
    std::vector<float> vTauDau_pt ;
    std::vector<float> vTauDau_eta ;
    std::vector<float> vTauDau_phi ;
    std::vector<float> vTauDau_mass ;
    std::vector<int> vTauDau_pdgId ;
    std::vector<int> vMoth;
    int nMoth = 0;
    int nDau  = 0;  
    nBranches_->genParticle_N = genParticles_->size();
    for( unsigned p=0; p<genParticles_->size(); ++p ){
      //if( (*genParticles_)[p].status() != 3 ) continue;
      vDau.clear(); vMoth.clear();
      vTauDau_pt.clear();
      vTauDau_eta.clear();
      vTauDau_phi.clear();
      vTauDau_mass.clear();
      vTauDau_pdgId.clear();

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
      nBranches_->genParticle_fromHardProcessFinalState.push_back((*genParticles_)[p].fromHardProcessFinalState());
      nBranches_->genParticle_isDirectHardProcessTauDecayProductFinalState.push_back((*genParticles_)[p].isDirectHardProcessTauDecayProductFinalState());


      for( unsigned int d=0; d<(*genParticles_)[p].numberOfDaughters(); ++d ){
        vDau.push_back( (*genParticles_)[p].daughter(d)->pdgId() );
	nDau++;
      }


      if(abs((*genParticles_)[p].pdgId())==15 && (*genParticles_)[p].statusFlags().isPrompt()){
	for( unsigned int d=0; d<(*genParticles_)[p].numberOfDaughters(); ++d ){
	  vTauDau_pt.push_back( (*genParticles_)[p].daughter(d)->pt() );
	  vTauDau_eta.push_back( (*genParticles_)[p].daughter(d)->eta() );
	  vTauDau_phi.push_back( (*genParticles_)[p].daughter(d)->phi() );
	  vTauDau_mass.push_back( (*genParticles_)[p].daughter(d)->mass() );
	  vTauDau_pdgId.push_back( (*genParticles_)[p].daughter(d)->pdgId() );
	}
      }


      for( unsigned int m=0; m<(*genParticles_)[p].numberOfMothers(); ++m ){
        vMoth.push_back( (*genParticles_)[p].mother(m)->pdgId() );
    	nMoth++;
      }
      nBranches_->genParticle_nDau  .push_back( nDau  );
      nBranches_->genParticle_nMoth .push_back( nMoth );      
      nBranches_->genParticle_mother.push_back( vMoth );
      nBranches_->genParticle_dau   .push_back( vDau  );      

      // only for taus
      nBranches_->genParticle_taudau_pt   .push_back( vTauDau_pt  );      
      nBranches_->genParticle_taudau_eta   .push_back( vTauDau_eta  );      
      nBranches_->genParticle_taudau_phi   .push_back( vTauDau_phi  );      
      nBranches_->genParticle_taudau_mass   .push_back( vTauDau_mass  );      
      nBranches_->genParticle_taudau_pdgId   .push_back( vTauDau_pdgId  );      
    }

    
}

