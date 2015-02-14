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
  
  nBranches_->ngenJetsAK4 = 0;
  
  for (const reco::GenJet &j : *genJets_) {
	  
      	 if (j.pt() < 20) continue;
      
	  nBranches_->ngenJetsAK4++;
      
	  nBranches_->genJetAK4_pt     	    .push_back(j.pt());
	  nBranches_->genJetAK4_eta    	    .push_back(j.eta());
	  nBranches_->genJetAK4_mass   	    .push_back(j.mass());
	  nBranches_->genJetAK4_phi    	    .push_back(j.phi());   
	  nBranches_->genJetAK4_e      	    .push_back(j.energy());
          float visibleFraction = 0;
          if(j.energy()>0) visibleFraction = (j.energy()-j.invisibleEnergy())/j.energy();
	  nBranches_->genJetNoNuAK4_pt     	    .push_back(j.pt()*visibleFraction);
	  nBranches_->genJetNoNuAK4_mass   	    .push_back(j.mass()*visibleFraction);
	  nBranches_->genJetNoNuAK4_e      	    .push_back(j.energy()*visibleFraction);

  }
 
}



