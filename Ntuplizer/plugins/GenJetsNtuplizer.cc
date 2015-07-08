#include "../interface/GenJetsNtuplizer.h"

//===================================================================================================================        

GenJetsNtuplizer::GenJetsNtuplizer( edm::EDGetTokenT<reco::GenJetCollection> token, NtupleBranches* nBranches )

   : CandidateNtuplizer     ( nBranches )
   , genJetInputToken_	    ( token )

{
}

//===================================================================================================================
GenJetsNtuplizer::~GenJetsNtuplizer( void )
{
}

//===================================================================================================================
void GenJetsNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
  
  event.getByToken(genJetInputToken_   , genJets_    );
  
  nBranches_->genJetAK4_N = 0;
  
  for (const reco::GenJet &j : *genJets_) {
	  
      	 if (j.pt() < 20) continue;
      
	  nBranches_->genJetAK4_N++;
      
	  nBranches_->genJetAK4_pt     	    .push_back(j.pt());
	  nBranches_->genJetAK4_eta    	    .push_back(j.eta());
	  nBranches_->genJetAK4_mass   	    .push_back(j.mass());
	  nBranches_->genJetAK4_phi    	    .push_back(j.phi());   
	  nBranches_->genJetAK4_e      	    .push_back(j.energy());
          float visibleFraction = (j.energy()>0) ?(j.energy()-j.invisibleEnergy())/j.energy() : 0;
	  nBranches_->genJetNoNuAK4_pt     	    .push_back(j.pt()*visibleFraction);
	  nBranches_->genJetNoNuAK4_mass   	    .push_back(j.mass()*visibleFraction);
	  nBranches_->genJetNoNuAK4_e      	    .push_back(j.energy()*visibleFraction);

  }
 
}



