#include "../interface/JetsNtuplizer.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
//#include "AnalysisDataFormats/TopObjects/interface/CATopJetTagInfo.h"

//===================================================================================================================        
JetsNtuplizer::JetsNtuplizer( std::vector<edm::EDGetTokenT<pat::JetCollection>> tokens, edm::EDGetTokenT<reco::JetFlavourMatchingCollection> flavourToken, NtupleBranches* nBranches )
   : CandidateNtuplizer     ( nBranches )
   , jetInputToken_	    ( tokens[0] )
   , fatjetInputToken_	    ( tokens[1] )
   , prunedjetInputToken_   ( tokens[2] )
   , softdropjetInputToken_ ( tokens[3] )
   //, flavourToken_			( flavourToken 	) //For subjet flavour matching!! Not done yet.

{
   
}

//===================================================================================================================
JetsNtuplizer::~JetsNtuplizer( void )
{
}

//===================================================================================================================
void JetsNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
  
  event.getByToken(jetInputToken_   , jets_    );
  event.getByToken(fatjetInputToken_, fatjets_ );
 //event.getByToken(flavourToken_, jetMC );
  
  bool doPruning = event.getByToken(prunedjetInputToken_, prunedjets_ );
  bool doSoftDrop = event.getByToken(softdropjetInputToken_, softdropjets_ );
  
  /****************************************************************/
  nBranches_->njetsAK4 = 0;
  
  for (const pat::Jet &j : *jets_) {
	  
	  if (j.pt() < 20) continue;
	  
	  bool IDLoose = false;
	  
	  if( j.nConstituents() > 1 && 
              j.muonEnergyFraction() < 0.99 && 
              j.photonEnergyFraction() < 0.99 && 
              j.chargedEmEnergyFraction() < 0.99 &&
              j.neutralHadronEnergyFraction() < 0.99 && 
              j.chargedHadronEnergyFraction() > 0. ) IDLoose = true;
	  
	  if( !IDLoose ) continue;
	  nBranches_->njetsAK4++;
	  
	  nBranches_->jetAK4_pt     	    .push_back(j.pt());
	  nBranches_->jetAK4_eta    	    .push_back(j.eta());
	  nBranches_->jetAK4_mass   	    .push_back(j.mass());
	  nBranches_->jetAK4_phi    	    .push_back(j.phi());   
	  nBranches_->jetAK4_e      	    .push_back(j.energy());
	  nBranches_->jetAK4_IDLoose        .push_back(IDLoose);
	  nBranches_->jetAK4_cm     	    .push_back(j.chargedMultiplicity());
	  nBranches_->jetAK4_nm     	    .push_back(j.neutralMultiplicity());
	  nBranches_->jetAK4_che     	    .push_back(j.chargedHadronEnergy()+j.electronEnergy()+j.muonEnergy());
	  nBranches_->jetAK4_ne     	    .push_back(j.neutralHadronEnergy()+j.photonEnergy());
	  nBranches_->jetAK4_charge 	    .push_back(j.charge());
	  nBranches_->jetAK4_ssv    	    .push_back(j.bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));
	  nBranches_->jetAK4_cisv    	    .push_back(std::max(0.f,j.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags")));
	  nBranches_->jetAK4_tchp   	    .push_back(j.bDiscriminator("trackCountingHighPurBJetTags"));
	  nBranches_->jetAK4_tche   	    .push_back(j.bDiscriminator("trackCountingHighEffBJetTags"));
	  nBranches_->jetAK4_jp     	    .push_back(j.bDiscriminator("jetProbabilityBJetTags"));
	  nBranches_->jetAK4_jbp    	    .push_back(j.bDiscriminator("jetBProbabilityBJetTags"));
	  nBranches_->jetAK4_flavour	    .push_back(abs(j.partonFlavour()));
	  nBranches_->jetAK4_vtxMass	    .push_back(j.userFloat("vtxMass")); 
	  nBranches_->jetAK4_vtxNtracks	    .push_back(j.userFloat("vtxNtracks")); 
	  nBranches_->jetAK4_vtx3DVal	    .push_back(j.userFloat("vtx3DVal")); 
	  nBranches_->jetAK4_vtx3DSig	    .push_back(j.userFloat("vtx3DSig"));

  }
  
  /****************************************************************/
  nBranches_->njetsAK8 = 0;
  int nsubjets   = 0;
  //int nsoftdropsubjets = 0;
  
  std::vector<float> vPrunedSubjetpt     ;
  std::vector<float> vPrunedSubjeteta    ;
  std::vector<float> vPrunedSubjetmass   ;
  std::vector<float> vPrunedSubjetphi    ;
  std::vector<float> vPrunedSubjete	 ;
  std::vector<int  > vPrunedSubjetcharge ;
  std::vector<int  > vPrunedSubjetflavour;
  std::vector<float> vPrunedSubjetssv    ;
  std::vector<float> vPrunedSubjetcsv    ;
  std::vector<float> vPrunedSubjettchp   ;
  std::vector<float> vPrunedSubjettche   ; 
  std::vector<float> vPrunedSubjetjp     ;
  std::vector<float> vPrunedSubjetjbp    ; 
  
  //std::vector<float> vSoftDropSubjetpt     ;
  //std::vector<float> vSoftDropSubjeteta    ;
  //std::vector<float> vSoftDropSubjetmass   ;
  //std::vector<float> vSoftDropSubjetphi    ;
  //std::vector<float> vSoftDropSubjete      ;
  //std::vector<int  > vSoftDropSubjetcharge ;
  //std::vector<int  > vSoftDropSubjetflavour;
  //std::vector<float> vSoftDropSubjetssv    ;
  //std::vector<float> vSoftDropSubjetcsv    ;
  //std::vector<float> vSoftDropSubjettchp   ;
  //std::vector<float> vSoftDropSubjettche   ;
  //std::vector<float> vSoftDropSubjetjp     ;
  //std::vector<float> vSoftDropSubjetjbp    ;
  
  for (const pat::Jet &fj : *fatjets_) {
	  
	  if (fj.pt() < 80) continue;
	  
	  bool IDLoose = false;
	  
	  if( fj.nConstituents() > 1 && 
              fj.muonEnergyFraction() < 0.99 && 
              fj.photonEnergyFraction() < 0.99 && 
              fj.chargedEmEnergyFraction() < 0.99 &&
              fj.neutralHadronEnergyFraction() < 0.99 && 
              fj.chargedHadronEnergyFraction() > 0. ) IDLoose = true;

	  
	  if( !IDLoose ) continue;
	  
	  nBranches_->njetsAK8++;	       
                   
	  nBranches_->jetAK8_pt     	    .push_back(fj.pt());
	  nBranches_->jetAK8_eta    	    .push_back(fj.eta());
	  nBranches_->jetAK8_mass   	    .push_back(fj.mass());
	  nBranches_->jetAK8_phi    	    .push_back(fj.phi());
	  nBranches_->jetAK8_e      	    .push_back(fj.energy());
	  nBranches_->jetAK8_IDLoose        .push_back(IDLoose);
	  nBranches_->jetAK8_cm     	    .push_back(fj.chargedMultiplicity());
	  nBranches_->jetAK8_nm     	    .push_back(fj.neutralMultiplicity());     				       
	  nBranches_->jetAK8_che     	    .push_back(fj.chargedHadronEnergy()+fj.electronEnergy()+fj.muonEnergy());
	  nBranches_->jetAK8_ne     	    .push_back(fj.neutralHadronEnergy()+fj.photonEnergy());
	  nBranches_->jetAK8_charge 	    .push_back(fj.charge());					       
	  nBranches_->jetAK8_ssv    	    .push_back(fj.bDiscriminator("simpleSecondaryVertexHighPurBJetTags")); 
	  nBranches_->jetAK8_csv    	    .push_back(fj.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"));      
	  nBranches_->jetAK8_tchp   	    .push_back(fj.bDiscriminator("trackCountingHighPurBJetTags"));	       
	  nBranches_->jetAK8_tche   	    .push_back(fj.bDiscriminator("trackCountingHighEffBJetTags"));	       
	  nBranches_->jetAK8_jp     	    .push_back(fj.bDiscriminator("jetProbabilityBJetTags"));	       
	  nBranches_->jetAK8_jbp    	    .push_back(fj.bDiscriminator("jetBProbabilityBJetTags"));
	  //nBranches_->jetAK8_flavour		.push_back(abs(fj.partonFlavour()));
	  nBranches_->jetAK8_tau1	    .push_back(fj.userFloat("NjettinessAK8:tau1"));	   
	  nBranches_->jetAK8_tau2	    .push_back(fj.userFloat("NjettinessAK8:tau2"));
	  nBranches_->jetAK8_tau21	    .push_back((fj.userFloat("NjettinessAK8:tau2"))/(fj.userFloat("NjettinessAK8:tau1")));
	  nBranches_->jetAK8_tau3	    .push_back(fj.userFloat("NjettinessAK8:tau3")); 
	
	  nBranches_->jetAK8_prunedmass		.push_back(fj.userFloat("ak8PFJetsCHSPrunedLinks"));
	  nBranches_->jetAK8_softdropmass	.push_back(fj.userFloat("ak8PFJetsCHSSoftDropLinks"));
	  
	  nBranches_->jetAK8_trimmedmass	.push_back(fj.userFloat("ak8PFJetsCHSTrimmedLinks"));
	  nBranches_->jetAK8_filteredmass	.push_back(fj.userFloat("ak8PFJetsCHSFilteredLinks"));
	  
	  
	  nBranches_->njetsAK8pruned = 0;
	  
	  /****************************************************************/
	  if( doPruning ){
	  
	     TLorentzVector FatJet; FatJet.SetPtEtaPhiE( fj.pt(), fj.eta(), fj.phi(), fj.energy() );  
	  
	     //printf("FJ with pt %5.1f , softdrop mass %+4.2f\n", fj.pt(),  fj.userFloat("ak8PFJetsCHSSoftDropLinks"));

	     // Loop over pruned jets, store dR match. Used to obtain subjets.
	     float dRmin =  999. ; 
	     	   
	     pat::Jet prunedjet;
	  
	     for (const pat::Jet &pj : *prunedjets_) {
	     	     TLorentzVector jetPruned; jetPruned.SetPtEtaPhiE( pj.pt(), pj.eta(), pj.phi(), pj.energy() );   

	     	     float dRtmp   = FatJet.DeltaR(jetPruned);
	     	     if( dRtmp < dRmin && dRtmp < 0.8 ){
	     		     dRmin     = dRtmp;
	     		     prunedjet = pj;
	     	     }
	     }
	  
             //printf("PJ with pt %5.1f , mass %+4.2f\n", prunedjet.pt(),  prunedjet.mass());

	     nBranches_->jetAK8pruned_pt     .push_back(prunedjet.pt());
	     nBranches_->jetAK8pruned_eta    .push_back(prunedjet.eta());
	     nBranches_->jetAK8pruned_mass   .push_back(prunedjet.mass());
	     nBranches_->jetAK8pruned_phi    .push_back(prunedjet.phi());
	     nBranches_->jetAK8pruned_e      .push_back(prunedjet.energy());
	     //nBranches_->jetAK8pruned_flavour.push_back(abs(prunedjet.partonFlavour()));
	     nBranches_->jetAK8pruned_charge .push_back(prunedjet.charge());
	     nBranches_->jetAK8pruned_ssv    .push_back(prunedjet.bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));
	     nBranches_->jetAK8pruned_csv    .push_back(prunedjet.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags"));  //Obs! b-tag for fatjets not implemented yet. Recipe from Dinko when ready
	     nBranches_->jetAK8pruned_tchp   .push_back(prunedjet.bDiscriminator("trackCountingHighPurBJetTags"));
	     nBranches_->jetAK8pruned_tche   .push_back(prunedjet.bDiscriminator("trackCountingHighEffBJetTags"));
	     nBranches_->jetAK8pruned_jp     .push_back(prunedjet.bDiscriminator("jetProbabilityBJetTags"));
	     nBranches_->jetAK8pruned_jbp    .push_back(prunedjet.bDiscriminator("jetBProbabilityBJetTags"));
	  
	     nBranches_->njetsAK8pruned++;
	  
	     vPrunedSubjetpt	 .clear();
	     vPrunedSubjeteta	 .clear();    
	     vPrunedSubjetmass   .clear();   
	     vPrunedSubjetphi	 .clear();    
	     vPrunedSubjete	 .clear();
	     vPrunedSubjetcharge .clear();
	     vPrunedSubjetflavour.clear();
	     vPrunedSubjetssv	 .clear();    
	     vPrunedSubjetcsv	 .clear();    
	     vPrunedSubjettchp   .clear();   
	     vPrunedSubjettche   .clear();
	     vPrunedSubjetjp	 .clear();
	     vPrunedSubjetjbp	 .clear();	
	  
	     //Get subjets from matched pruned jet
	     std::vector<reco::CandidatePtr> pruneddaus(prunedjet.daughterPtrVector());
  
	     nsubjets = 0;

	     for (unsigned int k = 0; k < pruneddaus.size(); ++k) {
	     	     const reco::PFJet &prunedsubjet = dynamic_cast<const reco::PFJet &>(*pruneddaus[k]);
	     	     //std::cout<<"prunedsubjet.pt()  "<<prunedsubjet.pt() <<std::endl;
	      
	     	     if( prunedsubjet.pt() < 0.01 ) continue;
	     	     
	     	     nsubjets++;
      
	     	     //vPrunedSubjetflavour.push_back(flavor);
	     	     vPrunedSubjetpt	 .push_back(prunedsubjet.pt()); 
	     	     vPrunedSubjeteta	 .push_back(prunedsubjet.eta()); 
	     	     vPrunedSubjetmass   .push_back(prunedsubjet.mass()); 
	     	     vPrunedSubjetphi	 .push_back(prunedsubjet.phi()); 
	     	     vPrunedSubjete	 .push_back(prunedsubjet.energy());
	     	     vPrunedSubjetcharge .push_back(prunedsubjet.charge()); 
	     	     // vPrunedSubjetssv    .push_back(prunedsubjet.bDiscriminator("simpleSecondaryVertexHighPurBJetTags") ); //Obs!! b-tag for subjets not implemented yet
	     	     // vPrunedSubjetcsv    .push_back(prunedsubjet.bDiscriminator("combinedSecondaryVertexBJetTags"	 ) );
	     	     // vPrunedSubjettchp   .push_back(prunedsubjet.bDiscriminator("trackCountingHighPurBJetTags"	)  );
	     	     // vPrunedSubjettche   .push_back(prunedsubjet.bDiscriminator("trackCountingHighEffBJetTags"	)  );
	     	     // vPrunedSubjetjp     .push_back(prunedsubjet.bDiscriminator("jetProbabilityBJetTags"		    )	   );
	     	     // vPrunedSubjetjbp    .push_back(prunedsubjet.bDiscriminator("jetBProbabilityBJetTags"		    )	   );
	     }

	     nBranches_->nsubjets		.push_back(nsubjets		   );
	     nBranches_->subjetAK8pruned_pt	.push_back(vPrunedSubjetpt     );  
	     nBranches_->subjetAK8pruned_eta	.push_back(vPrunedSubjeteta    );  
	     nBranches_->subjetAK8pruned_mass	.push_back(vPrunedSubjetmass   );  
	     nBranches_->subjetAK8pruned_phi	.push_back(vPrunedSubjetphi    );  
	     nBranches_->subjetAK8pruned_e	.push_back(vPrunedSubjete      );
	     nBranches_->subjetAK8pruned_charge .push_back(vPrunedSubjetcharge ); 
	     //nBranches_->subjetAK8pruned_flavour.push_back(vPrunedSubjetflavour);  
	     nBranches_->subjetAK8pruned_ssv	.push_back(vPrunedSubjetssv    );  
	     nBranches_->subjetAK8pruned_csv	.push_back(vPrunedSubjetcsv    );  
	     nBranches_->subjetAK8pruned_tchp	.push_back(vPrunedSubjettchp   );  
	     nBranches_->subjetAK8pruned_tche	.push_back(vPrunedSubjettche   );
	     nBranches_->subjetAK8pruned_jp	.push_back(vPrunedSubjetjp     );
	     nBranches_->subjetAK8pruned_jbp	.push_back(vPrunedSubjetjbp    ); 
	  
	  }
	  
	  /****************************************************************/
	  if( doSoftDrop ){
	  
	     // Loop over softdrop jets, store match. Used to compare ak8PFJetCHS softdrop mass (ak8PFJetsCHSSoftDropLinks) as this is with JEC
	     float dRmin =  999. ; 
	     pat::Jet softdropjet;
	
	     TLorentzVector FatJet; FatJet.SetPtEtaPhiE( fj.pt(), fj.eta(), fj.phi(), fj.energy() );
	     
	     for (const pat::Jet &sdj : *softdropjets_) {
	     	     TLorentzVector jetSoftDrop; jetSoftDrop.SetPtEtaPhiE( sdj.pt(), sdj.eta(), sdj.phi(), sdj.energy() );   
	     	     float dRtmp   = FatJet.DeltaR(jetSoftDrop);
	     	     if( dRtmp < dRmin && dRtmp < 0.8 ){
	     		     dRmin  = dRtmp;
	     		     softdropjet = sdj;
	     	     }
	     }
	     nBranches_->jetAK8softdrop_mass   .push_back(softdropjet.mass());
	     //printf("SJ with pt %5.1f , mass %+4.2f\n", softdropjet.pt(),  softdropjet.mass());
  	  
	    //
	    // vSoftDropSubjetpt      .clear();
	    // vSoftDropSubjeteta    .clear();
	    // vSoftDropSubjetmass   .clear();
	    // vSoftDropSubjetphi    .clear();
	    // vSoftDropSubjete       .clear();
	    // vSoftDropSubjetcharge .clear();
	    // vSoftDropSubjetflavour.clear();
	    // vSoftDropSubjetssv    .clear();
	    // vSoftDropSubjetcsv    .clear();
	    // vSoftDropSubjettchp   .clear();
	    // vSoftDropSubjettche   .clear();
	    // vSoftDropSubjetjp      .clear();
	    // vSoftDropSubjetjbp .clear();
	    //
 	    // std::vector<reco::CandidatePtr> softdropdaus(softdropjet.daughterPtrVector());
 	    	 
 	    // nsoftdropsubjets = 0;
 	    //
 	    //   for (unsigned int k = 0; k < softdropdaus.size(); ++k) {
 	    //  	 std::cout<<" softdropdaus.size()  "<< softdropdaus.size() <<std::endl;
 	    //   }
 	    //
 	    //  	 if( softdropsubjet.pt() < 0.01 ) continue;
 	    //
 	    //  	 nsoftdropsubjets++;
 	    //
 	    //
 	    //  	 //vPrunedSubjetflavour.push_back(flavor);
 	    //  	 vSoftDropSubjetpt     .push_back(softdropsubjet.pt());
 	    //  	 vSoftDropSubjeteta    .push_back(softdropsubjet.eta());
 	    //  	 vSoftDropSubjetmass   .push_back(softdropsubjet.mass());
 	    //  	 vSoftDropSubjetphi    .push_back(softdropsubjet.phi());
 	    //  	 vSoftDropSubjete      .push_back(softdropsubjet.energy());
 	    //  	 vSoftDropSubjetcharge .push_back(softdropsubjet.charge());
 	    //  	 // vSoftDropSubjetssv    .push_back(softdropsubjet.bDiscriminator("simpleSecondaryVertexHighPurBJetTags")     );
 	    //  	 // vSoftDropSubjetcsv    .push_back(softdropsubjet.bDiscriminator("combinedSecondaryVertexBJetTags"	 )     );
 	    //  	 // vSoftDropSubjettchp   .push_back(softdropsubjet.bDiscriminator("trackCountingHighPurBJetTags"	    )  );
 	    //  	 // vSoftDropSubjettche   .push_back(softdropsubjet.bDiscriminator("trackCountingHighEffBJetTags"	    )  );
 	    //  	 // vSoftDropSubjetjp	  .push_back(softdropsubjet.bDiscriminator("jetProbabilityBJetTags"		)      );
 	    //  	 // vSoftDropSubjetjbp    .push_back(softdropsubjet.bDiscriminator("jetBProbabilityBJetTags"		)      );
 	    //   }
 	    //
 	    //   nBranches_->nsoftdropsubjets	    .push_back(nsoftdropsubjets      );
 	    //   nBranches_->subjetAK8softdrop_pt     .push_back(vSoftDropSubjetpt     );
 	    //   nBranches_->subjetAK8softdrop_eta    .push_back(vSoftDropSubjeteta    );
 	    //   nBranches_->subjetAK8softdrop_mass   .push_back(vSoftDropSubjetmass   );
 	    //   nBranches_->subjetAK8softdrop_phi    .push_back(vSoftDropSubjetphi    );
 	    //   nBranches_->subjetAK8softdrop_e      .push_back(vSoftDropSubjete      );
 	    //   nBranches_->subjetAK8softdrop_charge .push_back(vSoftDropSubjetcharge );
 	    //   //nBranches_->subjetAK8softdrop_flavour.push_back(vSoftDropSubjetflavour);
 	    //   nBranches_->subjetAK8softdrop_ssv    .push_back(vSoftDropSubjetssv    );
 	    //   nBranches_->subjetAK8softdrop_csv    .push_back(vSoftDropSubjetcsv    );
 	    //   nBranches_->subjetAK8softdrop_tchp   .push_back(vSoftDropSubjettchp   );
 	    //   nBranches_->subjetAK8softdrop_tche   .push_back(vSoftDropSubjettche   );
 	    //   nBranches_->subjetAK8softdrop_jp     .push_back(vSoftDropSubjetjp     );
 	    //   nBranches_->subjetAK8softdrop_jbp    .push_back(vSoftDropSubjetjbp    );
	    
	  }   
	  	   
  }
  
}



