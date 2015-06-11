#include "../interface/JetsNtuplizer.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
//#include "AnalysisDataFormats/TopObjects/interface/CATopJetTagInfo.h"

//===================================================================================================================        

JetsNtuplizer::JetsNtuplizer( std::vector<edm::EDGetTokenT<pat::JetCollection>> tokens, std::vector<std::string> jecAK4Labels, std::vector<std::string> jecAK8Labels, edm::EDGetTokenT<reco::JetFlavourMatchingCollection> flavourToken, edm::EDGetTokenT<double> rhoToken, edm::EDGetTokenT<reco::VertexCollection> verticeToken, NtupleBranches* nBranches )

  : CandidateNtuplizer     ( nBranches )

  , jetInputToken_	    ( tokens[0] )
  , fatjetInputToken_	    ( tokens[1] )
  , prunedjetInputToken_    ( tokens[2] )
  , softdropjetInputToken_  ( tokens[3] )
  , rhoToken_	       	    ( rhoToken  )
  , verticeToken_     	    ( verticeToken  )	    
   //, flavourToken_			( flavourToken 	) //For subjet flavour matching!! Not done yet.
    
{
	
  jecAK4PayloadNames_ = jecAK4Labels;
  jecAK4PayloadNames_.pop_back();
	
  jecAK8PayloadNames_ = jecAK8Labels;
  jecAK8PayloadNames_.pop_back();

      
  initJetCorrFactors();
   
}

//===================================================================================================================
JetsNtuplizer::~JetsNtuplizer( void )
{
}
//===================================================================================================================
bool JetsNtuplizer::looseJetID( const pat::Jet& j ) {
  // In sync with: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_8_TeV_data_a
  return ( j.nConstituents() > 1 && 
    j.photonEnergyFraction() < 0.99 && 
      j.neutralHadronEnergyFraction() < 0.99 && 
        j.muonEnergyFraction() < 0.8 && 
          j.electronEnergyFraction() < 0.9 && 
            (j.chargedHadronMultiplicity() > 0 || fabs(j.eta())>2.4 ) && 
              (j.chargedEmEnergyFraction() < 0.99 || fabs(j.eta())>2.4 ) &&
                (j.chargedHadronEnergyFraction() > 0. || fabs(j.eta())>2.4 ) );
}
//===================================================================================================================
void JetsNtuplizer::initJetCorrFactors( void ){

  std::vector<JetCorrectorParameters> vPar;
  for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8PayloadNames_.begin(), payloadEnd = jecAK8PayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
    JetCorrectorParameters pars(*ipayload);
    vPar.push_back(pars);
  }
  
  // Make the FactorizedJetCorrector
  jecAK8_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
 
  vPar.clear();
  for ( std::vector<std::string>::const_iterator payloadBegin = jecAK4PayloadNames_.begin(), payloadEnd = jecAK4PayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
    JetCorrectorParameters pars(*ipayload);
    vPar.push_back(pars);
  }
  
  // Make the FactorizedJetCorrector
  jecAK4_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
}

//===================================================================================================================
void JetsNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
  
  event.getByToken(jetInputToken_   , jets_    );
  event.getByToken(fatjetInputToken_, fatjets_ );
  event.getByToken(rhoToken_	 	, rho_	   );
  event.getByToken(verticeToken_ 	, vertices_ );
  //event.getByToken(flavourToken_, jetMC );
  
  nBranches_->rho = *(rho_.product());

  bool doPruning  = event.getByToken(prunedjetInputToken_, prunedjets_ );
  bool doSoftDrop = event.getByToken(softdropjetInputToken_, softdropjets_ );
  
  /****************************************************************/
  nBranches_->njetsAK4 = 0;
  
  for (const pat::Jet &j : *jets_) {
	  
    if (j.pt() < 20) continue;
      
    bool IDLoose = looseJetID(j);
	     
    //if( !IDLoose ) continue;

    reco::Candidate::LorentzVector uncorrJet = j.correctedP4(0);

    jecAK4_->setJetEta( uncorrJet.eta()	);
    jecAK4_->setJetPt ( uncorrJet.pt()	);
    jecAK4_->setJetE  ( uncorrJet.energy() );
    jecAK4_->setRho   ( nBranches_->rho  );
    jecAK4_->setNPV   ( vertices_->size()  );
    jecAK4_->setJetA  ( j.jetArea()	);
    double corr = jecAK4_->getCorrection();
    
    //jecAK4Unc_->setJetEta( uncorrJet.eta() );
    //jecAK4Unc_->setJetPt( corr * uncorrJet.pt() );
    //double corrUp = corr * (1 + fabs(jecAK4Unc_->getUncertainty(1)));
    //jecAK4Unc_->setJetEta( uncorrJet.eta() );
    //jecAK4Unc_->setJetPt( corr * uncorrJet.pt() );
    //double corrDown = corr * ( 1 - fabs(jecAK4Unc_->getUncertainty(-1)) );

    nBranches_->njetsAK4++;

    nBranches_->jetAK4_pt     	    .push_back(corr*uncorrJet.pt());      
    nBranches_->jetAK4_eta    	    .push_back(j.eta());
    nBranches_->jetAK4_mass   	    .push_back(corr*uncorrJet.mass());
    nBranches_->jetAK4_phi    	    .push_back(j.phi());   
    nBranches_->jetAK4_e      	    .push_back(corr*uncorrJet.energy());
    nBranches_->jetAK4_jec    	    .push_back(corr);
    //nBranches_->jetAK4_jecUp        	.push_back(corrUp);
    //nBranches_->jetAK4_jecDown      	.push_back(corrDown);
    nBranches_->jetAK4_IDLoose      .push_back(IDLoose);
    nBranches_->jetAK4_cm     	    .push_back(j.chargedMultiplicity());
    nBranches_->jetAK4_nm     	    .push_back(j.neutralMultiplicity());
    nBranches_->jetAK4_muf     	    .push_back(j.muonEnergyFraction());
    nBranches_->jetAK4_phf     	    .push_back(j.photonEnergyFraction());
    nBranches_->jetAK4_emf     	    .push_back(j.chargedEmEnergyFraction());
    nBranches_->jetAK4_nhf     	    .push_back(j.neutralHadronEnergyFraction());
    nBranches_->jetAK4_chf     	    .push_back(j.chargedHadronEnergyFraction());
    nBranches_->jetAK4_area         .push_back(j.jetArea());	  
    nBranches_->jetAK4_che     	    .push_back(j.chargedHadronEnergy()+j.electronEnergy()+j.muonEnergy());
    nBranches_->jetAK4_ne     	    .push_back(j.neutralHadronEnergy()+j.photonEnergy());
    nBranches_->jetAK4_charge 	    .push_back(j.charge());
    nBranches_->jetAK4_ssv    	    .push_back(j.bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags"));
    nBranches_->jetAK4_cisv    	    .push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    nBranches_->jetAK4_tchp   	    .push_back(j.bDiscriminator("pfTrackCountingHighPurBJetTags"));
    nBranches_->jetAK4_tche   	    .push_back(j.bDiscriminator("pfTrackCountingHighEffBJetTags"));
    nBranches_->jetAK4_jp     	    .push_back(j.bDiscriminator("pfJetProbabilityBJetTags"));
    nBranches_->jetAK4_jbp    	    .push_back(j.bDiscriminator("pfJetBProbabilityBJetTags"));
    nBranches_->jetAK4_flavour	    .push_back(abs(j.partonFlavour()));
    nBranches_->jetAK4_vtxMass	    .push_back(j.userFloat("vtxMass")); 
    nBranches_->jetAK4_vtxNtracks   .push_back(j.userFloat("vtxNtracks")); 
    nBranches_->jetAK4_vtx3DVal	    .push_back(j.userFloat("vtx3DVal")); 
    nBranches_->jetAK4_vtx3DSig	    .push_back(j.userFloat("vtx3DSig"));
	
  }
  
  /****************************************************************/
  nBranches_->njetsAK8 = 0;
  int nsubjets   = 0;
  // int nsoftdropsubjets = 0;
  
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
  
  std::vector<float> vSoftDropSubjetpt     ;
  std::vector<float> vSoftDropSubjeteta    ;
  std::vector<float> vSoftDropSubjetmass   ;
  std::vector<float> vSoftDropSubjetphi    ;
  std::vector<float> vSoftDropSubjete      ;
  std::vector<int  > vSoftDropSubjetcharge ;
  std::vector<int  > vSoftDropSubjetflavour;
  std::vector<float> vSoftDropSubjetssv    ;
  std::vector<float> vSoftDropSubjetcsv    ;
  std::vector<float> vSoftDropSubjettchp   ;
  std::vector<float> vSoftDropSubjettche   ;
  std::vector<float> vSoftDropSubjetjp     ;
  std::vector<float> vSoftDropSubjetjbp    ;
     
  for (const pat::Jet &fj : *fatjets_) {
	  
    reco::Candidate::LorentzVector uncorrJet = fj.correctedP4(0);

    jecAK8_->setJetEta( uncorrJet.eta()	       );
    jecAK8_->setJetPt ( uncorrJet.pt()	       );
    jecAK8_->setJetE  ( uncorrJet.energy()       );
    jecAK8_->setJetA  ( fj.jetArea() );
    jecAK8_->setRho   ( nBranches_->rho	     );
    jecAK8_->setNPV   ( vertices_->size()        );
    double corr = jecAK8_->getCorrection();

    //jecAK8Unc_->setJetEta( uncorrJet.eta() );
    //jecAK8Unc_->setJetPt( corr * uncorrJet.pt() );
    //double corrUp = corr * (1 + fabs(jecAK8Unc_->getUncertainty(1)));
    //jecAK8Unc_->setJetEta( uncorrJet.eta() );
    //jecAK8Unc_->setJetPt( corr * uncorrJet.pt() );
    //double corrDown = corr * ( 1 - fabs(jecAK8Unc_->getUncertainty(-1)) );
    
    if( corr*uncorrJet.pt() < 80. ) continue;
    
    if (!doSoftDrop){
      
      vSoftDropSubjetpt.clear();
      vSoftDropSubjeteta.clear();
      vSoftDropSubjetmass.clear();
      vSoftDropSubjetphi.clear();
      vSoftDropSubjete.clear();
      vSoftDropSubjetcharge.clear();
      vSoftDropSubjetflavour.clear();
      vSoftDropSubjetssv.clear();
      vSoftDropSubjetcsv.clear();
      vSoftDropSubjettchp.clear();
      vSoftDropSubjettche.clear();
      vSoftDropSubjetjp.clear();
      vSoftDropSubjetjbp.clear();
      
      nsubjets = 0;
      
      const std::vector<edm::Ptr<pat::Jet> > &wSubjets = fj.subjets("SoftDrop");
    
    	for ( const pat::Jet & softdropsubjet : wSubjets ) {
         if( softdropsubjet.pt() < 0.01 ) continue;
         
         nsubjets++;

         vSoftDropSubjetflavour.push_back(abs(softdropsubjet.partonFlavour()));
         vSoftDropSubjetpt.push_back(softdropsubjet.pt());
         vSoftDropSubjeteta.push_back(softdropsubjet.eta());
         vSoftDropSubjetmass.push_back(softdropsubjet.mass());
         vSoftDropSubjetphi.push_back(softdropsubjet.phi());
         vSoftDropSubjete.push_back(softdropsubjet.energy());
         vSoftDropSubjetflavour.push_back(abs(softdropsubjet.partonFlavour()));
         vSoftDropSubjetcharge.push_back(softdropsubjet.charge());
         vSoftDropSubjetssv.push_back(softdropsubjet.bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags")); 
         vSoftDropSubjetcsv.push_back(softdropsubjet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
         vSoftDropSubjettchp.push_back(softdropsubjet.bDiscriminator("pfTrackCountingHighPurBJetTags") );
         vSoftDropSubjettche.push_back(softdropsubjet.bDiscriminator("pfTrackCountingHighEffBJetTags") );
         vSoftDropSubjetjp.push_back(softdropsubjet.bDiscriminator("pfJetProbabilityBJetTags") );
         vSoftDropSubjetjbp.push_back(softdropsubjet.bDiscriminator("pfJetBProbabilityBJetTags") );
         
       } 

      nBranches_->nsoftdropsubjets.push_back(nsubjets);
      nBranches_->subjetAK8softdrop_pt.push_back(vSoftDropSubjetpt);
      nBranches_->subjetAK8softdrop_eta.push_back(vSoftDropSubjeteta);
      nBranches_->subjetAK8softdrop_mass.push_back(vSoftDropSubjetmass);
      nBranches_->subjetAK8softdrop_phi.push_back(vSoftDropSubjetphi);
      nBranches_->subjetAK8softdrop_e.push_back(vSoftDropSubjete);
      nBranches_->subjetAK8softdrop_charge.push_back(vSoftDropSubjetcharge);
      nBranches_->subjetAK8softdrop_flavour.push_back(vSoftDropSubjetflavour);
      nBranches_->subjetAK8softdrop_ssv.push_back(vSoftDropSubjetssv);
      nBranches_->subjetAK8softdrop_csv.push_back(vSoftDropSubjetcsv);
      nBranches_->subjetAK8softdrop_tchp.push_back(vSoftDropSubjettchp);
      nBranches_->subjetAK8softdrop_tche.push_back(vSoftDropSubjettche);
      nBranches_->subjetAK8softdrop_jp.push_back(vSoftDropSubjetjp);
      nBranches_->subjetAK8softdrop_jbp.push_back(vSoftDropSubjetjbp);
    
    
    }
    	                    
    bool IDLoose = looseJetID(fj);
    //if( !IDLoose ) continue;

    nBranches_->njetsAK8++;	       
    nBranches_->jetAK8_pt     	    .push_back(corr*uncorrJet.pt());                   
    nBranches_->jetAK8_eta    	    .push_back(fj.eta());
    nBranches_->jetAK8_mass   	    .push_back(corr*uncorrJet.mass());
    nBranches_->jetAK8_phi    	    .push_back(fj.phi());
    nBranches_->jetAK8_e      	    .push_back(corr*uncorrJet.energy());
    nBranches_->jetAK8_jec    	    .push_back(corr);
    nBranches_->jetAK8_IDLoose      .push_back(IDLoose);
    nBranches_->jetAK8_muf     	    .push_back(fj.muonEnergyFraction());
    nBranches_->jetAK8_phf     	    .push_back(fj.photonEnergyFraction());
    nBranches_->jetAK8_emf     	    .push_back(fj.chargedEmEnergyFraction());
    nBranches_->jetAK8_nhf     	    .push_back(fj.neutralHadronEnergyFraction());
    nBranches_->jetAK8_chf     	    .push_back(fj.chargedHadronEnergyFraction());
    nBranches_->jetAK8_area         .push_back(fj.jetArea());
    nBranches_->jetAK8_cm     	    .push_back(fj.chargedMultiplicity());
    nBranches_->jetAK8_nm     	    .push_back(fj.neutralMultiplicity());     				       
    nBranches_->jetAK8_che     	    .push_back(fj.chargedHadronEnergy()+fj.electronEnergy()+fj.muonEnergy());
    nBranches_->jetAK8_ne     	    .push_back(fj.neutralHadronEnergy()+fj.photonEnergy());
    nBranches_->jetAK8_charge 	    .push_back(fj.charge());					 
    nBranches_->jetAK8_Hbbtag    	  .push_back(fj.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags"));       
    nBranches_->jetAK8_ssv    	    .push_back(fj.bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags")); 
    nBranches_->jetAK8_csv    	    .push_back(fj.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));      
    nBranches_->jetAK8_tchp   	    .push_back(fj.bDiscriminator("pfTrackCountingHighPurBJetTags"));	       
    nBranches_->jetAK8_tche   	    .push_back(fj.bDiscriminator("pfTrackCountingHighEffBJetTags"));	       
    nBranches_->jetAK8_jp     	    .push_back(fj.bDiscriminator("pfJetProbabilityBJetTags"));	       
    nBranches_->jetAK8_jbp    	    .push_back(fj.bDiscriminator("pfJetBProbabilityBJetTags"));
    nBranches_->jetAK8_flavour      .push_back(abs(fj.partonFlavour()));
    nBranches_->jetAK8_tau1         .push_back(fj.userFloat("NjettinessAK8:tau1"));	   
    nBranches_->jetAK8_tau2         .push_back(fj.userFloat("NjettinessAK8:tau2"));
    nBranches_->jetAK8_tau3         .push_back(fj.userFloat("NjettinessAK8:tau3")); 
    nBranches_->jetAK8_prunedmass   .push_back(fj.userFloat("ak8PFJetsCHSPrunedMass"));
    nBranches_->jetAK8_softdropmass .push_back(fj.userFloat("ak8PFJetsCHSSoftDropMass"));
	  	 
    TLorentzVector FatJet; FatJet.SetPtEtaPhiE( fj.pt(), fj.eta(), fj.phi(), fj.energy() );  

    if( doPruning ){

      // nBranches_->njetsAK8pruned = 0;

      // Loop over pruned jets, store dR match. Used to obtain subjets.
      float dRmin =  999. ;

      pat::Jet prunedjet;

      for (const pat::Jet &pj : *prunedjets_) {

        TLorentzVector jetPruned; jetPruned.SetPtEtaPhiE( pj.pt(), pj.eta(), pj.phi(), pj.energy() );

        float dRtmp   = FatJet.DeltaR(jetPruned);
        if( dRtmp < dRmin && dRtmp < 0.8 ){ dRmin = dRtmp; prunedjet = pj;}
        else continue;

      }

      //Compute JEC for pruned mass
      reco::Candidate::LorentzVector uncorrPrunedJet = prunedjet.correctedP4(0);

      jecAK8_->setJetEta( uncorrPrunedJet.eta()     );
      jecAK8_->setJetPt ( uncorrPrunedJet.pt()     );
      jecAK8_->setJetE  ( uncorrPrunedJet.energy()     );
      jecAK8_->setJetA  ( prunedjet.jetArea() );
      jecAK8_->setRho   ( nBranches_->rho       );
      jecAK8_->setNPV   ( vertices_->size()       );
      double prunedcorr = jecAK8_->getCorrection();

      nBranches_->jetAK8_prunedmassCorr.push_back(prunedcorr*fj.userFloat("ak8PFJetsCHSPrunedMass"));
      nBranches_->jetAK8pruned_jec   .push_back(prunedcorr);

      /****************************************************************/

      // nBranches_->jetAK8pruned_pt     .push_back(prunedjet.pt());
      // nBranches_->jetAK8pruned_eta    .push_back(prunedjet.eta());
      // nBranches_->jetAK8pruned_mass   .push_back(prunedjet.mass());
      // nBranches_->jetAK8pruned_phi    .push_back(prunedjet.phi());
      // nBranches_->jetAK8pruned_e      .push_back(prunedjet.energy());
      // nBranches_->jetAK8pruned_jec    .push_back(prunedjet.mass());
      // nBranches_->jetAK8pruned_flavour.push_back(abs(prunedjet.partonFlavour()));
      // nBranches_->jetAK8pruned_charge .push_back(prunedjet.charge());
      // nBranches_->jetAK8pruned_ssv    .push_back(prunedjet.bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags"));
      // nBranches_->jetAK8pruned_csv    .push_back(prunedjet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
      // nBranches_->jetAK8pruned_tchp   .push_back(prunedjet.bDiscriminator("pfTrackCountingHighPurBJetTags"));
      // nBranches_->jetAK8pruned_tche   .push_back(prunedjet.bDiscriminator("pfTrackCountingHighEffBJetTags"));
      // nBranches_->jetAK8pruned_jp    .push_back(prunedjet.bDiscriminator("pfJetProbabilityBJetTags"));
      // nBranches_->jetAK8pruned_jbp    .push_back(prunedjet.bDiscriminator("pfJetBProbabilityBJetTags"));
      // nBranches_->njetsAK8pruned++;

      vPrunedSubjetpt     .clear();
      vPrunedSubjeteta    .clear();
      vPrunedSubjetmass   .clear();
      vPrunedSubjetphi    .clear();
      vPrunedSubjete      .clear();
      vPrunedSubjetcharge .clear();
      vPrunedSubjetflavour.clear();
      vPrunedSubjetssv    .clear();
      vPrunedSubjetcsv    .clear();
      vPrunedSubjettchp   .clear();
      vPrunedSubjettche   .clear();
      vPrunedSubjetjp     .clear();
      vPrunedSubjetjbp    .clear();

      //Get subjets from matched pruned jet
      std::vector<reco::CandidatePtr> pruneddaus(prunedjet.daughterPtrVector());

      nsubjets = 0;

      for (unsigned int k = 0; k < pruneddaus.size(); ++k) {

        // const reco::PFJet &prunedsubjet = dynamic_cast<const reco::PFJet &>(*pruneddaus[k]);
        
        const pat::Jet &prunedsubjet = dynamic_cast<const pat::Jet &>(*pruneddaus[k]);   

        if( prunedsubjet.pt() < 0.01 ) continue;

        nsubjets++;

        // vPrunedSubjetflavour.push_back(flavor);
        vPrunedSubjetpt.push_back(prunedsubjet.pt());
        vPrunedSubjeteta.push_back(prunedsubjet.eta());
        vPrunedSubjetmass.push_back(prunedsubjet.mass());
        vPrunedSubjetphi.push_back(prunedsubjet.phi());
        vPrunedSubjete.push_back(prunedsubjet.energy());
        vPrunedSubjetflavour.push_back(abs(prunedsubjet.partonFlavour()));
        vPrunedSubjetcharge.push_back(prunedsubjet.charge());
        vPrunedSubjetssv.push_back(prunedsubjet.bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags")); 
        vPrunedSubjetcsv.push_back(prunedsubjet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
        vPrunedSubjettchp.push_back(prunedsubjet.bDiscriminator("pfTrackCountingHighPurBJetTags") );
        vPrunedSubjettche.push_back(prunedsubjet.bDiscriminator("pfTrackCountingHighEffBJetTags") );
        vPrunedSubjetjp.push_back(prunedsubjet.bDiscriminator("pfJetProbabilityBJetTags") );
        vPrunedSubjetjbp.push_back(prunedsubjet.bDiscriminator("pfJetBProbabilityBJetTags") );     
        
      }

      nBranches_->nprunedsubjets.push_back(nsubjets         );
      nBranches_->subjetAK8pruned_pt.push_back(vPrunedSubjetpt     );
      nBranches_->subjetAK8pruned_eta.push_back(vPrunedSubjeteta    );
      nBranches_->subjetAK8pruned_mass.push_back(vPrunedSubjetmass   );
      nBranches_->subjetAK8pruned_phi.push_back(vPrunedSubjetphi    );
      nBranches_->subjetAK8pruned_e.push_back(vPrunedSubjete      );
      nBranches_->subjetAK8pruned_charge.push_back(vPrunedSubjetcharge );
      nBranches_->subjetAK8pruned_flavour.push_back(vPrunedSubjetflavour);
      nBranches_->subjetAK8pruned_ssv.push_back(vPrunedSubjetssv    );
      nBranches_->subjetAK8pruned_csv.push_back(vPrunedSubjetcsv    );
      nBranches_->subjetAK8pruned_tchp.push_back(vPrunedSubjettchp   );
      nBranches_->subjetAK8pruned_tche.push_back(vPrunedSubjettche   );
      nBranches_->subjetAK8pruned_jp.push_back(vPrunedSubjetjp     );
      nBranches_->subjetAK8pruned_jbp.push_back(vPrunedSubjetjbp    );
    }
	  
    /****************************************************************/
    if( doSoftDrop ){

      // Loop over softdrop jets
      float dRmin =  999. ;
      pat::Jet softdropjet;

      for (const pat::Jet &sdj : *softdropjets_) {
        TLorentzVector jetSoftDrop; jetSoftDrop.SetPtEtaPhiE( sdj.pt(), sdj.eta(), sdj.phi(), sdj.energy() );

        float dRtmp   = FatJet.DeltaR(jetSoftDrop);
        if( dRtmp < dRmin && dRtmp < 0.8 ){ dRmin  = dRtmp; softdropjet = sdj;}
        else continue;
      }

      reco::Candidate::LorentzVector uncorrSoftDropJet = softdropjet.correctedP4(0);

      jecAK8_->setJetEta( uncorrSoftDropJet.eta()   );
      jecAK8_->setJetPt ( uncorrSoftDropJet.pt()    );
      jecAK8_->setJetE  ( uncorrSoftDropJet.energy());
      jecAK8_->setJetA  ( softdropjet.jetArea()              );
      jecAK8_->setRho   ( nBranches_->rho     );
      jecAK8_->setNPV   ( vertices_->size()     );
      double softdropcorr = jecAK8_->getCorrection();

      nBranches_->jetAK8_softdropmassCorr.push_back(softdropcorr*fj.userFloat("ak8PFJetsCHSSoftDropMass"));
      nBranches_->jetAK8softdrop_jec  .push_back(softdropcorr);
      
      vSoftDropSubjetpt.clear();
      vSoftDropSubjeteta.clear();
      vSoftDropSubjetmass.clear();
      vSoftDropSubjetphi.clear();
      vSoftDropSubjete.clear();
      vSoftDropSubjetcharge.clear();
      vSoftDropSubjetflavour.clear();
      vSoftDropSubjetssv.clear();
      vSoftDropSubjetcsv.clear();
      vSoftDropSubjettchp.clear();
      vSoftDropSubjettche.clear();
      vSoftDropSubjetjp.clear();
      vSoftDropSubjetjbp.clear();

      //Get subjets from matched pruned jet
      std::vector<reco::CandidatePtr> softdropdaus(softdropjet.daughterPtrVector());

      nsubjets = 0;

      for (unsigned int k = 0; k < softdropdaus.size(); ++k) {

        // const reco::PFJet &softdropsubjet = dynamic_cast<const reco::PFJet &>(*softdropdaus[k]);
        
        const pat::Jet &softdropsubjet = dynamic_cast<const pat::Jet &>(*softdropdaus[k]);   

        if( softdropsubjet.pt() < 0.01 ) continue;

        nsubjets++;

        // vSoftdropSubjetflavour.push_back(flavor);
        vSoftDropSubjetpt.push_back(softdropsubjet.pt());
        vSoftDropSubjeteta.push_back(softdropsubjet.eta());
        vSoftDropSubjetmass.push_back(softdropsubjet.mass());
        vSoftDropSubjetphi.push_back(softdropsubjet.phi());
        vSoftDropSubjete.push_back(softdropsubjet.energy());
        vSoftDropSubjetflavour.push_back(abs(softdropsubjet.partonFlavour()));
        vSoftDropSubjetcharge.push_back(softdropsubjet.charge());
        vSoftDropSubjetssv.push_back(softdropsubjet.bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags")); 
        vSoftDropSubjetcsv.push_back(softdropsubjet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
        vSoftDropSubjettchp.push_back(softdropsubjet.bDiscriminator("pfTrackCountingHighPurBJetTags") );
        vSoftDropSubjettche.push_back(softdropsubjet.bDiscriminator("pfTrackCountingHighEffBJetTags") );
        vSoftDropSubjetjp.push_back(softdropsubjet.bDiscriminator("pfJetProbabilityBJetTags") );
        vSoftDropSubjetjbp.push_back(softdropsubjet.bDiscriminator("pfJetBProbabilityBJetTags") );       
      }

      nBranches_->nsoftdropsubjets.push_back(nsubjets);
      nBranches_->subjetAK8softdrop_pt.push_back(vSoftDropSubjetpt);
      nBranches_->subjetAK8softdrop_eta.push_back(vSoftDropSubjeteta);
      nBranches_->subjetAK8softdrop_mass.push_back(vSoftDropSubjetmass);
      nBranches_->subjetAK8softdrop_phi.push_back(vSoftDropSubjetphi);
      nBranches_->subjetAK8softdrop_e.push_back(vSoftDropSubjete);
      nBranches_->subjetAK8softdrop_charge.push_back(vSoftDropSubjetcharge);
      nBranches_->subjetAK8softdrop_flavour.push_back(vSoftDropSubjetflavour);
      nBranches_->subjetAK8softdrop_ssv.push_back(vSoftDropSubjetssv);
      nBranches_->subjetAK8softdrop_csv.push_back(vSoftDropSubjetcsv);
      nBranches_->subjetAK8softdrop_tchp.push_back(vSoftDropSubjettchp);
      nBranches_->subjetAK8softdrop_tche.push_back(vSoftDropSubjettche);
      nBranches_->subjetAK8softdrop_jp.push_back(vSoftDropSubjetjp);
      nBranches_->subjetAK8softdrop_jbp.push_back(vSoftDropSubjetjbp);
    } 
  }
}



