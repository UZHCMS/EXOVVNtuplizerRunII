#include "../interface/GenParticlesNtuplizer.h"
 
//===================================================================================================================        
GenParticlesNtuplizer::GenParticlesNtuplizer( std::vector<edm::EDGetTokenT<reco::GenParticleCollection>> tokens, 
					      std::vector<edm::EDGetTokenT<bool>> tauSpinnerTokens_bool, 
					      std::vector<edm::EDGetTokenT<double>> tauSpinnerTokens_double, 
					      NtupleBranches* nBranches ) 

   : CandidateNtuplizer( nBranches )
   , genParticlesToken_( tokens[0] )
   , tauSpinnerWTisValidToken_ (tauSpinnerTokens_bool[0])
   , tauSpinnerWTToken_ (tauSpinnerTokens_double[0])
   , tauSpinnerWThminusToken_ (tauSpinnerTokens_double[1])
   , tauSpinnerWThplusToken_ (tauSpinnerTokens_double[2])
   , tauSpinnerTauPolFromZToken_ (tauSpinnerTokens_double[3])
   , tauSpinnerWRightToken_ (tauSpinnerTokens_double[4])
   , tauSpinnerWLeftToken_ (tauSpinnerTokens_double[5])
   , tauSpinnerIsRightLeftToken_ (tauSpinnerTokens_double[6])

{

}

//===================================================================================================================        
GenParticlesNtuplizer::~GenParticlesNtuplizer( void )
{
}

//===================================================================================================================        
  void GenParticlesNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
  
    event.getByToken(genParticlesToken_ , genParticles_); 


//    std::cout << "Retrieving TauSpinner information !!" << std::endl;
//
//    event.getByToken(tauSpinnerWTisValidToken_, tauSpinnerWTisValid_); 
//    event.getByToken(tauSpinnerWTToken_, tauSpinnerWT_); 
//    event.getByToken(tauSpinnerWThminusToken_ , tauSpinnerWThminus_);
//    event.getByToken(tauSpinnerWThplusToken_ , tauSpinnerWThplus_);
//    event.getByToken(tauSpinnerTauPolFromZToken_ , tauSpinnerTauPolFromZ_);
//    event.getByToken(tauSpinnerWRightToken_ , tauSpinnerWRight_);
//    event.getByToken(tauSpinnerWLeftToken_ , tauSpinnerWLeft_);
//    event.getByToken(tauSpinnerIsRightLeftToken_ , tauSpinnerIsRightLeft_);
//
//    
//    if(tauSpinnerWTisValid_.isValid()){
//      std::cout << "tauSpinnerWTisValid =" << *tauSpinnerWTisValid_  << std::endl;
//      std::cout << "tauSpinnerWT_ =" << *tauSpinnerWT_  << std::endl;
//      std::cout << "tauSpinnerWThminus_ =" << *tauSpinnerWThminus_  << std::endl;
//      std::cout << "tauSpinnerWThplus_ =" << *tauSpinnerWThplus_  << std::endl;
//      std::cout << "tauSpinnerTauPolFromZ_ =" << *tauSpinnerTauPolFromZ_  << std::endl;
//      std::cout << "tauSpinnerWRight_ =" << *tauSpinnerWRight_  << std::endl;
//      std::cout << "tauSpinnerWLeft_ =" << *tauSpinnerWLeft_  << std::endl;
//      std::cout << "tauSpinnerIsRightLeft_ =" << *tauSpinnerIsRightLeft_  << std::endl;
//      
//      nBranches_->TauSpinnerWTisValid  .push_back(*tauSpinnerWTisValid_);
//      nBranches_->TauSpinnerWT         .push_back(*tauSpinnerWT_);	
//      nBranches_->TauSpinnerWThminus   .push_back(*tauSpinnerWThminus_);	
//      nBranches_->TauSpinnerWThplus    .push_back(*tauSpinnerWThplus_);	
//      nBranches_->TauSpinnerTauPolFromZ.push_back(*tauSpinnerTauPolFromZ_);	
//      nBranches_->TauSpinnerWRight     .push_back(*tauSpinnerWRight_);	
//      nBranches_->TauSpinnerWLeft      .push_back(*tauSpinnerWLeft_);	
//      nBranches_->TauSpinnerIsRightLeft.push_back(*tauSpinnerIsRightLeft_);	
//
//    }




    //    std::cout << "tauSpinnerWTisValid valid =" << tauSpinnerWTisValid_.isValid()  << std::endl;
    //    std::cout << "tauSpinnerWTisValid =" << *tauSpinnerWTisValid_  << std::endl;
//    std::cout << "tauSpinnerWT_ =" << tauSpinnerWT_  << std::endl;
//    std::cout << "tauSpinnerWThminus_ =" << tauSpinnerWThminusToken_  << std::endl;
//    std::cout << "tauSpinnerWThplus_ =" << tauSpinnerWThplus_  << std::endl;
//    std::cout << "tauSpinnerTauPolFromZ_ =" << tauSpinnerTauPolFromZ_  << std::endl;
//    std::cout << "tauSpinnerWRight_ =" << tauSpinnerWRight_  << std::endl;
//    std::cout << "tauSpinnerWLeft_ =" << tauSpinnerWLeft_  << std::endl;
//    std::cout << "tauSpinnerIsRightLeft_ =" << tauSpinnerIsRightLeft_  << std::endl;


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
      nBranches_->genParticle_fromHardProcessFinalState.push_back((*genParticles_)[p].fromHardProcessFinalState());
      nBranches_->genParticle_isDirectHardProcessTauDecayProductFinalState.push_back((*genParticles_)[p].isDirectHardProcessTauDecayProductFinalState());


      for( unsigned int d=0; d<(*genParticles_)[p].numberOfDaughters(); ++d ){
        vDau.push_back( (*genParticles_)[p].daughter(d)->pdgId() );
	nDau++;
      }


      if(abs((*genParticles_)[p].pdgId())==15 && (*genParticles_)[p].statusFlags().isPrompt() && (*genParticles_)[p].status()==2){

	if(nDau>1){

	  bool flag_radioactive_gamma = false;
	  bool flag_radioactive_tau = false;
	  
	  for( unsigned int d=0; d<(*genParticles_)[p].numberOfDaughters(); ++d ){
	    Int_t taupdgId = abs((*genParticles_)[p].daughter(d)->pdgId());
	    
	    if(taupdgId==22) flag_radioactive_gamma = true;   
	    if(taupdgId==15) flag_radioactive_tau = true; 
	      
	  }

	  if(!(flag_radioactive_gamma && flag_radioactive_tau)){

	    TLorentzVector tau;
	    tau.SetPtEtaPhiM(0,0,0,0);
	    Int_t decaymode = -1;
	  
	    for( unsigned int d=0; d<(*genParticles_)[p].numberOfDaughters(); ++d ){
	      Float_t taupt = (*genParticles_)[p].daughter(d)->pt();
	      Float_t taueta = (*genParticles_)[p].daughter(d)->eta();
	      Float_t tauphi = (*genParticles_)[p].daughter(d)->phi();
	      Float_t taumass = (*genParticles_)[p].daughter(d)->mass();
	      Int_t taupdgId = abs((*genParticles_)[p].daughter(d)->pdgId());
	      
	      if(!(taupdgId >= 11 && taupdgId<=16)){
		TLorentzVector taudau;
		taudau.SetPtEtaPhiM(taupt, taueta, tauphi, taumass);
		tau += taudau;
		decaymode = 4;
	      }
	      if(taupdgId==11) decaymode = 2; // electron decay
	      if(taupdgId==13) decaymode = 3; // muon decay
	  
	    }

	    nBranches_->genParticle_tauvispt  .push_back( (float)tau.Pt()  );
	    nBranches_->genParticle_tauviseta  .push_back( (float)tau.Eta()  );
	    nBranches_->genParticle_tauvisphi  .push_back( (float)tau.Phi()  );
	    nBranches_->genParticle_tauvismass  .push_back( (float)tau.M()  );
	    nBranches_->genParticle_taudecay  .push_back( decaymode  );
	  }else{
	    nBranches_->genParticle_tauvispt  .push_back( -99.  );
	    nBranches_->genParticle_tauviseta  .push_back( -99.  );
	    nBranches_->genParticle_tauvisphi  .push_back( -99.  );
	    nBranches_->genParticle_tauvismass  .push_back( -99.  );
	    nBranches_->genParticle_taudecay  .push_back( 0  ); // self decay (tau -> tau)
	  }
	  
	}else{
	  nBranches_->genParticle_tauvispt  .push_back( -99.  );
	  nBranches_->genParticle_tauviseta  .push_back( -99.  );
	  nBranches_->genParticle_tauvisphi  .push_back( -99.  );
	  nBranches_->genParticle_tauvismass  .push_back( -99.  );
	  
	  nBranches_->genParticle_taudecay  .push_back( -1  ); // self decay (tau -> tau)
	  
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

    }

    
}

