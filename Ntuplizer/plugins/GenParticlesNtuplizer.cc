#include "../interface/GenParticlesNtuplizer.h"
 
//===================================================================================================================        
GenParticlesNtuplizer::GenParticlesNtuplizer( std::vector<edm::EDGetTokenT<reco::GenParticleCollection>> tokens, NtupleBranches* nBranches, std::map< std::string, bool >& runFlags//, TH1F* histGenWeights
 ) 

   : CandidateNtuplizer( nBranches )
   , genParticlesToken_( tokens[0] )
   , doGenHist_( runFlags["doGenHist"]  )
   , verbose_   (runFlags["verbose"])
   // , isBkgBSample_ (runFlags["isBkgBSample"])
//    , histGenWeights_ (histGenWeights)
{

}

//===================================================================================================================        
GenParticlesNtuplizer::~GenParticlesNtuplizer( void )
{
}




// Hints from https://github.com/cms-sw/cmssw/blob/00d3c78393965c7352859b8df14ed7d052dfff8c/PhysicsTools/HepMCCandAlgos/plugins/GenParticlePruner.cc

//void GenParticlesNtuplizer::recursiveDaughters(size_t index,
//					       int rank, 
//					       const reco::GenParticleCollection &src,
//					       std::vector<size_t> &allIndices,
//					       std::vector<int> &pdgs,
//					       std::vector<int> &layers,
//					       std::vector<float> &ppt,
//					       std::vector<float> &peta,
//					       std::vector<float> &pphi
//					       ) {
//
//  reco::GenParticleRefVector daughters = src[index].daughterRefVector();
//  // avoid infinite recursion if the daughters are set to "this" particle.
//  size_t cachedIndex = index;
//
//  for (reco::GenParticleRefVector::const_iterator i = daughters.begin(); i != daughters.end(); ++i) {
//
//    index = i->key();
//
//    // To also avoid infinite recursion if a "loop" is found in the daughter list,
//    // check to make sure the index hasn't already been added.
//    if (find(allIndices.begin(), allIndices.end(), index) == allIndices.end()) {
//
//      allIndices.push_back(index);
//      pdgs.push_back(src[index].pdgId());
//      layers.push_back(rank);
//      ppt.push_back(src[index].pt());
//      peta.push_back(src[index].eta());
//      pphi.push_back(src[index].phi());
//
////      for(int ii=0; ii<rank-1; ii++){
////	std::cout << "  "; 
////      }
////      std::cout << " +->" << src[index].pdgId() << std::endl;
//
//      if (cachedIndex != index) {
//        recursiveDaughters(index, rank+1, src, allIndices, pdgs, layers, ppt, peta, pphi);
//      }
//    }
//  }
//}







//void GenParticlesNtuplizer::allDaughters(reco::Candidate particle){
//
//  //	    rank += 1 ;
//  for(int i = 0; i < (int)particle.numberOfDaughters(); i++){
//    reco::GenParticle dau = GenParticle(particle.daughter(i));
//    std::cout << "-->" << dau.pdgId() << std::endl;
//    //	      dau.rank = rank
//    //		  daughters.append( dau )
//    allDaughters( dau );
//  }
//
//}



//===================================================================================================================        
bool GenParticlesNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
  

    event.getByToken(genParticlesToken_ , genParticles_); 

    if ( doGenHist_ ) {
        for( unsigned p=0; p<genParticles_->size(); ++p ){
            // Looking At B mesons who decay to Jpsi+X. Catalogue what else they decay to in addition to the Jpsi. Get the particle's pdgid, the pT, eta, and phi of it and the two muons from the jpsi, the jpsi's (aka dimuon) pt, eta, phi, and mass, and the B's visible pt, eta, phi, and mass
            if ( (  abs((*genParticles_)[p].pdgId()) >= 500
                    && abs((*genParticles_)[p].pdgId()) < 600 )
                 && (*genParticles_)[p].status() == 2 ) {
                TLorentzVector mu1, mu2, X;
          
                for (unsigned int d=0; d<(*genParticles_)[p].numberOfDaughters(); ++d ) {
                    if ( fabs((*genParticles_)[p].daughter(d)->pdgId()) == 12  || fabs((*genParticles_)[p].daughter(d)->pdgId()) == 14  || fabs((*genParticles_)[p].daughter(d)->pdgId()) == 16 ) {continue;} 
                    if ( (*genParticles_)[p].daughter(d)->pdgId() == 443 ) { 
                        // Loop over jpsi daughters & get the two mus. if there aren't two mus then skip it
                        for ( unsigned int jd=0; jd<(*genParticles_)[p].daughter(d)->numberOfDaughters(); ++jd ) {
                            if ( (*genParticles_)[p].daughter(d)->daughter(jd)->pdgId() == 13 ) {
                                mu1.SetPtEtaPhiM((*genParticles_)[p].daughter(d)->daughter(jd)->pt(),(*genParticles_)[p].daughter(d)->daughter(jd)->eta(),(*genParticles_)[p].daughter(d)->daughter(jd)->phi(), (*genParticles_)[p].daughter(d)->daughter(jd)->mass());
                            } 
                            else if ( (*genParticles_)[p].daughter(d)->daughter(jd)->pdgId() == -13 ) {
                                mu2.SetPtEtaPhiM((*genParticles_)[p].daughter(d)->daughter(jd)->pt(), (*genParticles_)[p].daughter(d)->daughter(jd)->eta(),(*genParticles_)[p].daughter(d)->daughter(jd)->phi(), (*genParticles_)[p].daughter(d)->daughter(jd)->mass());
                            }
                        }
                    } 
                }
                if (mu1.Pt() > 0 && mu2.Pt() > 0) {
                    float min_dau_pt =0;
                    int dau_index=0;
                    for (unsigned int d=0; d<(*genParticles_)[p].numberOfDaughters(); ++d ) {
                  
                        if ( fabs((*genParticles_)[p].daughter(d)->pdgId()) == 12   || fabs((*genParticles_)[p].daughter(d)->pdgId()) == 14  || fabs((*genParticles_)[p].daughter(d)->pdgId()) == 16 ) {continue;} 
                        if ( (*genParticles_)[p].daughter(d)->pdgId() == 443 ) {continue;}

                        if ((*genParticles_)[p].daughter(d)->pt() >min_dau_pt){
                            min_dau_pt= (*genParticles_)[p].daughter(d)->pt();
                            dau_index=d;
                        }
                    }

                  
                    // Histogram with bins labeled in Ntuplizer.cc
                    std::vector<int> bin={13,111,211,113,213,221,331,233,333,311,321,313,323,411,421,441,551,553};
                    std::vector<int>::iterator it =std::find(bin.begin(), bin.end(),fabs((*genParticles_)[p].daughter(dau_index)->pdgId())); 
                    if (it==bin.end() ) continue;
                    int index = std::distance(bin.begin(), it);
                    // mu+-, pi0,, pi+-,rho0,rho+-,eta, etaP,omega,phi,K0, K+, K*0,K*+, D+
                    //    D0, eta_c,eta_b, upsilon.
                  
                     
  
                     
                    // nBranches_->genParticle_Bdau_X_id->Fill(index); 
                    nBranches_->genParticle_Bdau_X_id->Fill(index); 
                  
                 
                    X.SetPtEtaPhiM((*genParticles_)[p].daughter(dau_index)->pt(),
                                   (*genParticles_)[p].daughter(dau_index)->eta(),
                                   (*genParticles_)[p].daughter(dau_index)->phi(),
                                   (*genParticles_)[p].daughter(dau_index)->mass());
               
                 
                }
                if (X.Pt() > 0) {
                    nBranches_->genParticle_Bdau_X_pt->Fill(X.Pt()); 
                    nBranches_->genParticle_Bdau_X_eta->Fill(X.Eta()); 
                    nBranches_->genParticle_Bdau_X_phi->Fill(X.Phi()); 
                    nBranches_->genParticle_Bdau_mu1_pt->Fill(mu1.Pt()); 
                    nBranches_->genParticle_Bdau_mu1_eta->Fill(mu1.Eta()); 
                    nBranches_->genParticle_Bdau_mu1_phi->Fill(mu1.Phi()); 
                    nBranches_->genParticle_Bdau_mu2_pt->Fill(mu2.Pt()); 
                    nBranches_->genParticle_Bdau_mu2_eta->Fill(mu2.Eta()); 
                    nBranches_->genParticle_Bdau_mu2_phi->Fill(mu2.Phi()); 
                    nBranches_->genParticle_Bdau_Jpsi_mass->Fill((mu1+mu2).M()); 
                    nBranches_->genParticle_Bdau_Jpsi_pt->Fill((mu1+mu2).Pt()); 
                    nBranches_->genParticle_Bdau_Jpsi_eta->Fill((mu1+mu2).Eta()); 
                    nBranches_->genParticle_Bdau_Jpsi_phi->Fill((mu1+mu2).Phi()); 
                    nBranches_->genParticle_Bvis_mass->Fill((mu1+mu2+X).M()); 
                    nBranches_->genParticle_Bvis_pt->Fill((mu1+mu2+X).Pt()); 
                    nBranches_->genParticle_Bvis_eta->Fill((mu1+mu2+X).Eta()); 
                    nBranches_->genParticle_Bvis_phi->Fill((mu1+mu2+X).Phi()); 
                }
            }
        }

    }


  
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
    
  
    /* here we want to save  gen particles info*/
   
    std::vector<int> vDau ;
    std::vector<int> vMoth;
    std::vector<float> ptMoth;
    int nMoth = 0;
    int nDau  = 0;  
    //nBranches_->genParticle_N = genParticles_->size(); // the genParticles are filtered below

    if(verbose_) std::cout << "------------------------------------------------------" << std::endl;


    for( unsigned p=0; p<genParticles_->size(); ++p ){
      
        vDau.clear(); vMoth.clear(); ptMoth.clear();
        nDau = 0; nMoth = 0;
      
        bool isPrompt( (*genParticles_)[p].statusFlags().isPrompt() );
        bool isDirectPromptTauDecayProduct( (*genParticles_)[p].statusFlags().isDirectPromptTauDecayProduct() );
        bool fromHardProcessFinalState( (*genParticles_)[p].fromHardProcessFinalState() );
        bool isDirectHardProcessTauDecayProductFinalState( (*genParticles_)[p].isDirectHardProcessTauDecayProductFinalState() );
        bool isLepton( abs((*genParticles_)[p].pdgId())>=11 && abs((*genParticles_)[p].pdgId())<=18 );
        bool isQuark( abs((*genParticles_)[p].pdgId())<=6 && abs((*genParticles_)[p].status())<=29 );
        bool isPhoton( abs((*genParticles_)[p].pdgId())==22 && (*genParticles_)[p].pt()>10. );
        bool isGluon( abs((*genParticles_)[p].pdgId())==22 && (*genParticles_)[p].pt()>10. );
        bool isWZH( abs((*genParticles_)[p].pdgId())>=23 && abs((*genParticles_)[p].pdgId())<=25 );
        bool isHeavyMeson( abs((*genParticles_)[p].pdgId())>0 && abs((*genParticles_)[p].pdgId())<=1000 );
        bool isHeavyBaryon( abs((*genParticles_)[p].pdgId())>=1000);
        bool isBSM( (abs((*genParticles_)[p].pdgId())>=30 && abs((*genParticles_)[p].pdgId())<=50) || abs((*genParticles_)[p].pdgId())>=1000000 );
        bool isB( (abs((*genParticles_)[p].pdgId())>=511 && abs((*genParticles_)[p].pdgId())<=545));
        bool isB2( abs((*genParticles_)[p].pdgId())==511 || abs((*genParticles_)[p].pdgId())==521 || abs((*genParticles_)[p].pdgId())==531 || abs((*genParticles_)[p].pdgId())==541 );
	//        bool isD( abs((*genParticles_)[p].pdgId())==411 || abs((*genParticles_)[p].pdgId())==421 || abs((*genParticles_)[p].pdgId())==431);
        bool isStatus2( (*genParticles_)[p].status()==2 );
        bool isStatus1( (*genParticles_)[p].status()==1 );

      //   // Implementation fo the weight for the B chain decay in the generic background B sample
      //   if (isBkgBSample_){        
	  // //            int motherID=0;
            
      //       if (abs((*genParticles_)[p].pdgId())==443 and abs((*genParticles_)[p].daughter(0)->pdgId())==13 ){
      //           int motherID = abs((GenParticlesNtuplizer::checkMom(&(*genParticles_)[p]))->pdgId());
      //           std:: cout<< " Jpsi status is "    <<  (*genParticles_)[p].status() << std::endl;
      //           std:: cout<< " Jpsi daugh is "    <<   (*genParticles_)[p].daughter(0)->pdgId() << std::endl;
      //           std:: cout<< " Jpsi mother is "    << motherID << std::endl;
                
      //           std::vector<int> B_hadron = {511,521,531,541,5112,5122,5132,5212,5232};   // at the beginning of the hist there is a bin for the all other possible decays   
                
      //           std::vector<int>::iterator it = std::find(B_hadron.begin(), B_hadron.end(), motherID);
      //           int index;
      //           if (it != B_hadron.end()) {
      //               index = std::distance(B_hadron.begin(), it);
      //               std:: cout<< "index is " << index <<" weight is " << histGenWeights_->GetBinContent(index+2)<< std::endl;
      //               nBranches_->genWeightBkgB = histGenWeights_->GetBinContent(index+2);
      //           } else {
      //               nBranches_->genWeightBkgB = histGenWeights_->GetBinContent(1); // in the first bin of the hist there is a generic 'other' for all the b decays not contained in the B_hadron vector;
      //           } 
      //       }
      //   } else { 
      //       nBranches_->genWeightBkgB = 1;
      //   }
    

//	std::cout << "[TESTGenParticlesNtuplizer] " << (*genParticles_)[p].pdgId() << " from ";
//	
//	for( unsigned int m=0; m<(*genParticles_)[p].numberOfMothers(); ++m ){
//	  if(verbose_) std::cout << (*genParticles_)[p].mother(m)->pdgId() << " ";
//	}
//	
//	if(verbose_) std::cout << "decays into ..." << std::endl;
//	
//	for( unsigned int d=0; d<(*genParticles_)[p].numberOfDaughters(); ++d ){
//	  if(verbose_) std::cout << "[TESTGenParticlesNtuplizer]   ---->  " << (*genParticles_)[p].daughter(d)->pdgId() << std::endl;
//	}
//	
//	std::cout << "[TESTGenParticlesNtuplizer] " << std::endl;


    
    
         if(!isLepton && !isQuark && !isPhoton && !isGluon && !isWZH && !isHeavyMeson && !isHeavyBaryon && !isBSM && !isDirectPromptTauDecayProduct && !fromHardProcessFinalState && !isDirectHardProcessTauDecayProductFinalState && !isB && !isStatus2 && !isStatus1) continue;
      
        //      nBranches_->genParticle_px    .push_back((*genParticles_)[p].px()     );
        //      nBranches_->genParticle_py    .push_back((*genParticles_)[p].py()     );
        //      nBranches_->genParticle_pz    .push_back((*genParticles_)[p].pz()     );
        //      nBranches_->genParticle_e     .push_back((*genParticles_)[p].energy() );
        nBranches_->genParticle_pt    .push_back((*genParticles_)[p].pt()     );
        nBranches_->genParticle_eta   .push_back((*genParticles_)[p].eta()    );
        nBranches_->genParticle_phi   .push_back((*genParticles_)[p].phi()    );
        nBranches_->genParticle_mass  .push_back((*genParticles_)[p].mass()   );
        nBranches_->genParticle_status.push_back((*genParticles_)[p].status() );
        nBranches_->genParticle_pdgId .push_back((*genParticles_)[p].pdgId()  );
	
	///////////////////////////////
	// Yuta 
	
	//	if(isB || isD){
	if(isB){
	  
	  if(verbose_) std::cout << "[GenParticlesNtuplizer] " << (*genParticles_)[p].pdgId() << " (" << (*genParticles_)[p].status() << ") from ";
	  
	  //	  Bool_t ismother_B = false;
	    
	  for( unsigned int m=0; m<(*genParticles_)[p].numberOfMothers(); ++m ){
	    if(verbose_) std::cout << (*genParticles_)[p].mother(m)->pdgId() << " (" << (*genParticles_)[p].mother(m)->status() <<") ";

	    //	    if( abs((*genParticles_)[p].mother(m)->pdgId())>=511 && abs((*genParticles_)[p].mother(m)->pdgId())<=545 ){
	    //	      ismother_B = true;
	    //	    }
	  }

	  //	  if(verbose_) std::cout << "[GenParticlesNtuplizer] ismother_B = " << ismother_B << " " << (*genParticles_)[p].numberOfDaughters() << " " << (*genParticles_)[p].daughter(0)->pdgId() << "  " << (*genParticles_)[p].pdgId()  << std::endl;
	  if(verbose_) std::cout << "decays into ..." << std::endl;

	  for( unsigned int d=0; d<(*genParticles_)[p].numberOfDaughters(); ++d ){
	    if(verbose_) std::cout << "[GenParticlesNtuplizer]   ---->  " << (*genParticles_)[p].daughter(d)->pdgId() << " (" << (*genParticles_)[p].daughter(d)->status()  << ")" << std::endl;
	  }
	}	  

//	  if( (*genParticles_)[p].numberOfDaughters()!=1 && 
//	      (*genParticles_)[p].numberOfMothers()==1 && 
//	      (*genParticles_)[p].mother(0)->pdgId()==(*genParticles_)[p].pdgId() 
//	      ){
//	    ismother_B = false; 
//	  }

	  

	if(isB2){
	  if(! ( (*genParticles_)[p].numberOfDaughters()==1 && (*genParticles_)[p].daughter(0)->pdgId()==(*genParticles_)[p].pdgId() )  && (*genParticles_)[p].status()==2){

	    int _rank = 1;
	    std::vector<size_t> allIndices;
	    std::vector<int> pdgs;
	    std::vector<int> layers;
	    std::vector<float> ppt;
	    std::vector<float> peta;
	    std::vector<float> pphi;
	    std::vector<int> isfinal;
	    
	    pdgs.push_back((*genParticles_)[p].pdgId());
	    layers.push_back(0);
	    ppt.push_back((*genParticles_)[p].pt());
	    peta.push_back((*genParticles_)[p].eta());
	    pphi.push_back((*genParticles_)[p].phi());
	    isfinal.push_back(0);

	    aux.recursiveDaughters(p, _rank, *genParticles_, allIndices, pdgs, layers, ppt, peta, pphi, isfinal, verbose_);


//	    std::cout << "******************* this is stored ******************" << std::endl;
//	    for(int ivec = 0; ivec < (int)pdgs.size(); ivec++){
//	      std::cout << ivec << " " << pdgs[ivec] << " " << layers[ivec] << " " << isfinal[ivec] << std::endl;
//	    }
//	    std::cout << "*****************************************************" << std::endl;

//	  for(int ivec = 0; ivec < (int)pdgs.size(); ivec++){
//	    nBranches_->genParticle_pmother.push_back( (*genParticles_)[p].pdgId() );
	    nBranches_->genParticle_pdgs.push_back( pdgs );
	    nBranches_->genParticle_layers.push_back( layers);
	    nBranches_->genParticle_ppt.push_back( ppt);
	    nBranches_->genParticle_peta.push_back( peta);
	    nBranches_->genParticle_pphi.push_back( pphi);
	    nBranches_->genParticle_isfinal.push_back( isfinal);
	    //	  }
	  }
	  

	  
	}
	

	///////////////////////////////



        // needed for the gen matching
        nBranches_->genParticle_isPrompt.push_back( isPrompt );
        nBranches_->genParticle_isDirectPromptTauDecayProduct.push_back( isDirectPromptTauDecayProduct );

        // needed for the MVA recoil correction
        nBranches_->genParticle_fromHardProcessFinalState.push_back( fromHardProcessFinalState );
        nBranches_->genParticle_isDirectHardProcessTauDecayProductFinalState.push_back( isDirectHardProcessTauDecayProductFinalState );


        for( unsigned int d=0; d<(*genParticles_)[p].numberOfDaughters(); ++d ){
            vDau.push_back( (*genParticles_)[p].daughter(d)->pdgId() );
            nDau++;
        }

        for( unsigned int m=0; m<(*genParticles_)[p].numberOfMothers(); ++m ){
            vMoth.push_back( (*genParticles_)[p].mother(m)->pdgId() );
            ptMoth.push_back( (*genParticles_)[p].mother(m)->pt() );
            nMoth++;
        }

        nBranches_->genParticle_nDau  .push_back( nDau  );
        nBranches_->genParticle_nMoth .push_back( nMoth );      
        nBranches_->genParticle_mother.push_back( vMoth );
        nBranches_->genParticle_mother_pt.push_back( ptMoth );
        nBranches_->genParticle_dau   .push_back( vDau  );      

    }

    nBranches_->genParticle_N = nBranches_->genParticle_pt.size(); // save number of save genParticles
    //if     nBranches_->genWeightBkgB = 1;    

    return true;
}
