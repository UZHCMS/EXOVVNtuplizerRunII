#include "../interface/GenJetsNtuplizer.h"

//===================================================================================================================        

GenJetsNtuplizer::GenJetsNtuplizer( edm::EDGetTokenT<reco::GenJetCollection> token, edm::EDGetTokenT<pat::JetCollection> AK8token, NtupleBranches* nBranches, std::map< std::string, bool >& runFlags )

  : CandidateNtuplizer     ( nBranches )
  , genJetInputToken_	    ( token     )
  , genJetAK8InputToken_	  ( AK8token  )
  , isJpsiMu_( runFlags["doJpsiMu"]  )
  , isJpsiEle_( runFlags["doJpsiEle"]  )
{
}

//===================================================================================================================
GenJetsNtuplizer::~GenJetsNtuplizer( void )
{
}

//===================================================================================================================
void GenJetsNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
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

 
  bool doGenAK8  = event.getByToken(genJetAK8InputToken_, genJetsAK8_ );
  
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
  
  if( !doGenAK8) return;
  
  
  
  event.getByToken(genJetAK8InputToken_   , genJetsAK8_    );
  
  nBranches_->genJetAK8_N = 0;
  
  for (const pat::Jet &j : *genJetsAK8_) {
    // if (j.pt() < 20) continue;
    nBranches_->genJetAK8_N++;
    nBranches_->genJetAK8_pt     	    .push_back(j.pt());
    nBranches_->genJetAK8_eta    	    .push_back(j.eta());
    nBranches_->genJetAK8_mass   	    .push_back(j.mass());
    nBranches_->genJetAK8_phi    	    .push_back(j.phi());   
    nBranches_->genJetAK8_e      	    .push_back(j.energy());
    nBranches_->genJetAK8_prunedmass  .push_back(j.userFloat("ak8GenJetsPrunedMass"));
    nBranches_->genJetAK8_softdropmass.push_back(j.userFloat("ak8GenJetsSoftDropMass"));
  }
}



