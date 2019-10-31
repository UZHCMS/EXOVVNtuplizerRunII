#include "../interface/GenParticlesNtuplizer.h"
 
//===================================================================================================================        
GenParticlesNtuplizer::GenParticlesNtuplizer( std::vector<edm::EDGetTokenT<reco::GenParticleCollection>> tokens, NtupleBranches* nBranches, std::map< std::string, bool >& runFlags ) 

   : CandidateNtuplizer( nBranches )
   , genParticlesToken_( tokens[0] )
   , isJpsiMu_( runFlags["doJpsiMu"])
   , isJpsiEle_( runFlags["doJpsiEle"]  )
   , isJpsiTau_( runFlags["doJpsiTau"]  )
   , doGenHist_( runFlags["doGenHist"]  )
{

}

//===================================================================================================================        
GenParticlesNtuplizer::~GenParticlesNtuplizer( void )
{
}

//===================================================================================================================        
void GenParticlesNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
  

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
   
    bool rflag = false;
    
    if(isJpsiEle_ || isJpsiMu_){
      if ( nBranches_->JpsiMu_B_pt.size() >=1) rflag = true;
    }
    
    if(isJpsiTau_){
      if ( nBranches_->JpsiTau_B_pt.size() >=1)  rflag = true;
    }

    if(rflag==false) return;
    
  
    /* here we want to save  gen particles info*/
   
    std::vector<int> vDau ;
    std::vector<int> vMoth;
    int nMoth = 0;
    int nDau  = 0;  
    //nBranches_->genParticle_N = genParticles_->size(); // the genParticles are filtered below
    for( unsigned p=0; p<genParticles_->size(); ++p ){
      
        vDau.clear(); vMoth.clear();
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
        bool isStatus2( (*genParticles_)[p].status()==2 );
        bool isStatus1( (*genParticles_)[p].status()==1 );
      
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
            nMoth++;
        }

        nBranches_->genParticle_nDau  .push_back( nDau  );
        nBranches_->genParticle_nMoth .push_back( nMoth );      
        nBranches_->genParticle_mother.push_back( vMoth );
        nBranches_->genParticle_dau   .push_back( vDau  );      

    }

    nBranches_->genParticle_N = nBranches_->genParticle_pt.size(); // save number of save genParticles
    
}

