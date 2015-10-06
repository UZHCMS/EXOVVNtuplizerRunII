#include "../interface/JetsNtuplizer.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"


//===================================================================================================================        

JetsNtuplizer::JetsNtuplizer( std::vector<edm::EDGetTokenT<pat::JetCollection>> tokens, std::vector<std::string> jecAK4Labels, std::vector<std::string> jecAK8Labels, std::vector<std::string> jecAK8GroomedLabels, std::vector<std::string> jecAK8PuppiLabels, edm::EDGetTokenT<reco::JetFlavourMatchingCollection> flavourToken, edm::EDGetTokenT<double> rhoToken, edm::EDGetTokenT<reco::VertexCollection> verticeToken, NtupleBranches* nBranches, std::map< std::string, bool >& runFlags )

  : CandidateNtuplizer     ( nBranches )

  , jetInputToken_	    ( tokens[0] )
  , fatjetInputToken_	    ( tokens[1] )
  , prunedjetInputToken_    ( tokens[2] )
  , softdropjetInputToken_  ( tokens[3] )
  , trimmedjetInputToken_   ( tokens[4] )
  , puppijetInputToken_     ( tokens[5] )
  , rhoToken_	       	    ( rhoToken  )
  , verticeToken_     	    ( verticeToken  )
  , doAK4Jets_ (runFlags["doAK4Jets"])
  , doAK8Jets_ (runFlags["doAK8Jets"])
  , runOnMC_ (runFlags["runOnMC"])
   //, flavourToken_			( flavourToken 	) //For subjet flavour matching!! Not done yet.
    
{
	
  doCorrOnTheFly_ = false;	
  if( jecAK4Labels.size() != 0 && jecAK8Labels.size() != 0 && jecAK8GroomedLabels.size() != 0 && jecAK8PuppiLabels.size() != 0 ){	
  
     jecAK4PayloadNames_ = jecAK4Labels;
     jecAK4PayloadNames_.pop_back();
	
     jecAK8PayloadNames_ = jecAK8Labels;
     jecAK8PayloadNames_.pop_back();

     jecAK8GroomedPayloadNames_ = jecAK8GroomedLabels;
     jecAK8GroomedPayloadNames_.pop_back();

     jecAK8PuppiPayloadNames_ = jecAK8PuppiLabels;
     jecAK8PuppiPayloadNames_.pop_back();

      
     initJetCorrFactors();
     
     doCorrOnTheFly_ = true;	
     
  }
   
}

//===================================================================================================================
JetsNtuplizer::~JetsNtuplizer( void )
{
}

//===================================================================================================================
bool JetsNtuplizer::looseJetID( const pat::Jet& j ) {

  /*return ( j.nConstituents() > 1 && 
    j.photonEnergyFraction() < 0.99 && 
      j.neutralHadronEnergyFraction() < 0.99 && 
        j.muonEnergyFraction() < 0.8 && 
          j.electronEnergyFraction() < 0.9 && 
            (j.chargedHadronMultiplicity() > 0 || fabs(j.eta())>2.4 ) && 
              (j.chargedEmEnergyFraction() < 0.99 || fabs(j.eta())>2.4 ) &&
                (j.chargedHadronEnergyFraction() > 0. || fabs(j.eta())>2.4 ) );*/

  //In sync with: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_8_TeV_data_a
  //and with dijet group https://github.com/CMSDIJET/DijetRootTreeMaker/blob/e650ba19e9e9bc676754a948298bb5cf850f4ecc/plugins/DijetTreeProducer.cc#L869

  double eta = j.eta();		
  double chf = j.chargedHadronEnergyFraction();
  double nhf = j.neutralHadronEnergyFraction(); // + j.HFHadronEnergyFraction();
  double muf = j.muonEnergy()/(j.jecFactor(0) * j.energy());  
  double nemf = j.neutralEmEnergyFraction();
  double cemf = j.chargedEmEnergyFraction();
  int chMult = j.chargedMultiplicity();
  int neMult = j.neutralMultiplicity();
  int npr    = chMult + neMult;
  int NumConst = npr;

  return (nhf<0.99 && nemf<0.99 && NumConst>1 && muf < 0.8) && ((fabs(eta) <= 2.4 && chf>0 && chMult>0 && cemf<0.99) || fabs(eta)>2.4);        
      		
}
//===================================================================================================================
bool JetsNtuplizer::tightJetID( const pat::Jet& j ) {



  //In sync with: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_8_TeV_data_a
  double eta = j.eta();		
  double chf = j.chargedHadronEnergyFraction();
  double nhf = j.neutralHadronEnergyFraction(); // + j.HFHadronEnergyFraction();
  double muf = j.muonEnergy()/(j.jecFactor(0) * j.energy());  
  double nemf = j.neutralEmEnergyFraction();
  double cemf = j.chargedEmEnergyFraction();
  int chMult = j.chargedMultiplicity();
  int neMult = j.neutralMultiplicity();
  int npr    = chMult + neMult;
  int NumConst = npr;
     
  return (nhf<0.90 && nemf<0.90 && NumConst>1 && muf<0.8) && ((fabs(eta)<=2.4 && chf>0 && chMult>0 && cemf<0.90) || fabs(eta)>2.4);  		
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
  for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8GroomedPayloadNames_.begin(), payloadEnd = jecAK8GroomedPayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
    JetCorrectorParameters pars(*ipayload);
    vPar.push_back(pars);
  }
  
  // Make the FactorizedJetCorrector
  jecAK8Groomed_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );

  vPar.clear();
  for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8PuppiPayloadNames_.begin(), payloadEnd = jecAK8PuppiPayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
    JetCorrectorParameters pars(*ipayload);
    vPar.push_back(pars);
  }
  
  // Make the FactorizedJetCorrector
  jecAK8Puppi_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );

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
  
  event.getByToken(jetInputToken_   , jets_     );
  event.getByToken(fatjetInputToken_, fatjets_  );
  event.getByToken(rhoToken_	    , rho_	);
  event.getByToken(verticeToken_    , vertices_ );
  //event.getByToken(flavourToken_, jetMC );
  
  nBranches_->rho = *(rho_.product());

  bool doPruning  = event.getByToken(prunedjetInputToken_, prunedjets_ );
  bool doSoftDrop = event.getByToken(softdropjetInputToken_, softdropjets_ );
  bool doTrimming  = event.getByToken(trimmedjetInputToken_, trimmedjets_ );
  bool doPuppi  = event.getByToken(puppijetInputToken_, puppijets_ );
  
  bool isMC = runOnMC_;

  /****************************************************************/
  if (doAK4Jets_) {
    nBranches_->jetAK4_N = 0;
  
    for (const pat::Jet &j : *jets_) {
	        
    bool IDLoose = looseJetID(j);
    bool IDTight = tightJetID(j);
    //if( !IDLoose ) continue;


      reco::Candidate::LorentzVector uncorrJet;
      double corr = 1;
    
      if( doCorrOnTheFly_ ){
    
         uncorrJet = j.correctedP4(0);
         jecAK4_->setJetEta( uncorrJet.eta()	);
         jecAK4_->setJetPt ( uncorrJet.pt()	);
         jecAK4_->setJetE  ( uncorrJet.energy() );
         jecAK4_->setRho   ( nBranches_->rho  );
         jecAK4_->setNPV   ( vertices_->size()  );
         jecAK4_->setJetA  ( j.jetArea()	);
         corr = jecAK4_->getCorrection();
       
      }
      else{
         uncorrJet = j.p4();
      }
      //jecAK4Unc_->setJetEta( uncorrJet.eta() );
      //jecAK4Unc_->setJetPt( corr * uncorrJet.pt() );
      //double corrUp = corr * (1 + fabs(jecAK4Unc_->getUncertainty(1)));
      //jecAK4Unc_->setJetEta( uncorrJet.eta() );
      //jecAK4Unc_->setJetPt( corr * uncorrJet.pt() );
      //double corrDown = corr * ( 1 - fabs(jecAK4Unc_->getUncertainty(-1)) );

      if (corr*uncorrJet.pt() < 20) continue;
      nBranches_->jetAK4_N++;

      nBranches_->jetAK4_pt     	    .push_back(corr*uncorrJet.pt());      
      nBranches_->jetAK4_eta    	    .push_back(j.eta());
      nBranches_->jetAK4_mass   	    .push_back(corr*uncorrJet.mass());
      nBranches_->jetAK4_phi    	    .push_back(j.phi());   
      nBranches_->jetAK4_e      	    .push_back(corr*uncorrJet.energy());
      nBranches_->jetAK4_jec    	    .push_back(corr);
      //nBranches_->jetAK4_jecUp        	.push_back(corrUp);
      //nBranches_->jetAK4_jecDown      	.push_back(corrDown);
      nBranches_->jetAK4_IDLoose      .push_back(IDLoose);
      nBranches_->jetAK4_IDTight      .push_back(IDTight);
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
      
      nBranches_->jetAK4_hf_hf  	    .push_back(j.HFHadronEnergyFraction());
      nBranches_->jetAK4_hf_emf  	    .push_back(j.HFEMEnergyFraction());
      nBranches_->jetAK4_hof     	    .push_back(j.hoEnergyFraction());
      
      nBranches_->jetAK4_chm     	    .push_back(j.chargedHadronMultiplicity());
      nBranches_->jetAK4_neHadMult    .push_back(j.neutralHadronMultiplicity());
      nBranches_->jetAK4_phoMult      .push_back(j.photonMultiplicity());
      
      nBranches_->jetAK4_nemf    	    .push_back(j.neutralEmEnergyFraction());
      nBranches_->jetAK4_cemf    	    .push_back(j.chargedEmEnergyFraction());
      
      nBranches_->jetAK4_charge 	    .push_back(j.charge());
      // nBranches_->jetAK4_ssv          .push_back(j.bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags"));
      nBranches_->jetAK4_cisv    	    .push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
      // nBranches_->jetAK4_tchp         .push_back(j.bDiscriminator("pfTrackCountingHighPurBJetTags"));
      // nBranches_->jetAK4_tche         .push_back(j.bDiscriminator("pfTrackCountingHighEffBJetTags"));
      // nBranches_->jetAK4_jp           .push_back(j.bDiscriminator("pfJetProbabilityBJetTags"));
      // nBranches_->jetAK4_jbp          .push_back(j.bDiscriminator("pfJetBProbabilityBJetTags"));
      nBranches_->jetAK4_vtxMass	    .push_back(j.userFloat("vtxMass")); 
      nBranches_->jetAK4_vtxNtracks   .push_back(j.userFloat("vtxNtracks")); 
      nBranches_->jetAK4_vtx3DVal	    .push_back(j.userFloat("vtx3DVal")); 
      nBranches_->jetAK4_vtx3DSig	    .push_back(j.userFloat("vtx3DSig"));
      if(isMC){
         int genP_pdgId = j.genParton() ? j.genParton()->pdgId() : -99; //det default to -99 when no genParton is found!!!
          nBranches_->jetAK4_partonFlavour.push_back(j.partonFlavour()); // b,c can be identified by either j.hadronFlavour() (hadron definition) or j.partonFlavour() (physics definition)
          nBranches_->jetAK4_hadronFlavour.push_back(j.hadronFlavour()); // g,u,d,s can be identified by j.partonFlavour() (physics definition)
          nBranches_->jetAK4_genParton_pdgID  .push_back(genP_pdgId);
          nBranches_->jetAK4_nbHadrons         .push_back((j.jetFlavourInfo().getbHadrons()).size()); //gluon splitting to bb can be identified roughly by (j.partonFlavour()==21 &% len(j.jetFlavourInfo().getbHadrons())==2 )
          nBranches_->jetAK4_ncHadrons         .push_back((j.jetFlavourInfo().getcHadrons()).size());
        }
    }
  } //doAK4Jets_
  
  /****************************************************************/
  if (doAK8Jets_) {
    nBranches_->jetAK8_N = 0;
    int nsubjets   = 0;
    // int nsoftdropsubjets = 0;
  
    std::vector<float> vPrunedSubjetpt     ;
    std::vector<float> vPrunedSubjeteta    ;
    std::vector<float> vPrunedSubjetmass   ;
    std::vector<float> vPrunedSubjetphi    ;
    std::vector<float> vPrunedSubjete	   ;
    std::vector<int  > vPrunedSubjetcharge ;
    std::vector<int  > vPrunedSubjetPartonFlavour;
    std::vector<int  > vPrunedSubjetHadronFlavour;
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
    std::vector<int  > vSoftDropSubjetPartonFlavour;
    std::vector<int  > vSoftDropSubjetHadronFlavour;
    std::vector<float> vSoftDropSubjetssv    ;
    std::vector<float> vSoftDropSubjetcsv    ;
    std::vector<float> vSoftDropSubjettchp   ;
    std::vector<float> vSoftDropSubjettche   ;
    std::vector<float> vSoftDropSubjetjp     ;
    std::vector<float> vSoftDropSubjetjbp    ;
     
    for (const pat::Jet &fj : *fatjets_) {
	  

      reco::Candidate::LorentzVector uncorrJet;
      double corr = 1;

      if( doCorrOnTheFly_ ){
    
         uncorrJet = fj.correctedP4(0);
         jecAK8_->setJetEta( uncorrJet.eta() );
         jecAK8_->setJetPt ( uncorrJet.pt() );
         jecAK8_->setJetE  ( uncorrJet.energy() );
         jecAK8_->setJetA  ( fj.jetArea() );
         jecAK8_->setRho   ( nBranches_->rho );
         jecAK8_->setNPV   ( vertices_->size() );
         corr = jecAK8_->getCorrection();
    
      }
      else{
         uncorrJet = fj.p4();
      }

      //jecAK8Unc_->setJetEta( uncorrJet.eta() );
      //jecAK8Unc_->setJetPt( corr * uncorrJet.pt() );
      //double corrUp = corr * (1 + fabs(jecAK8Unc_->getUncertainty(1)));
      //jecAK8Unc_->setJetEta( uncorrJet.eta() );
      //jecAK8Unc_->setJetPt( corr * uncorrJet.pt() );
      //double corrDown = corr * ( 1 - fabs(jecAK8Unc_->getUncertainty(-1)) );
    
      if( corr*uncorrJet.pt() < 100. ) continue;
            
      if (!doSoftDrop){
      
        vSoftDropSubjetpt.clear();
        vSoftDropSubjeteta.clear();
        vSoftDropSubjetmass.clear();
        vSoftDropSubjetphi.clear();
        vSoftDropSubjete.clear();
        vSoftDropSubjetcharge.clear();
        vSoftDropSubjetPartonFlavour.clear();
        vSoftDropSubjetHadronFlavour.clear();
        // vSoftDropSubjetssv.clear();
        vSoftDropSubjetcsv.clear();
        // vSoftDropSubjettchp.clear();
        // vSoftDropSubjettche.clear();
        // vSoftDropSubjetjp.clear();
        // vSoftDropSubjetjbp.clear();
      
        nsubjets = 0;
      
        const std::vector<edm::Ptr<pat::Jet> > &wSubjets = fj.subjets("SoftDrop");
    
      	for ( const pat::Jet & softdropsubjet : wSubjets ) {
           if( softdropsubjet.pt() < 0.01 ) continue;
         
           nsubjets++;

           vSoftDropSubjetPartonFlavour.push_back(abs(softdropsubjet.partonFlavour()));
           vSoftDropSubjetpt.push_back(softdropsubjet.pt());
           vSoftDropSubjeteta.push_back(softdropsubjet.eta());
           vSoftDropSubjetmass.push_back(softdropsubjet.mass());
           vSoftDropSubjetphi.push_back(softdropsubjet.phi());
           vSoftDropSubjete.push_back(softdropsubjet.energy());
           vSoftDropSubjetPartonFlavour.push_back(abs(softdropsubjet.partonFlavour()));
           vSoftDropSubjetHadronFlavour.push_back(abs(softdropsubjet.hadronFlavour()));
           vSoftDropSubjetcharge.push_back(softdropsubjet.charge());
           // vSoftDropSubjetssv.push_back(softdropsubjet.bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags"));
           vSoftDropSubjetcsv.push_back(softdropsubjet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
           // vSoftDropSubjettchp.push_back(softdropsubjet.bDiscriminator("pfTrackCountingHighPurBJetTags") );
           // vSoftDropSubjettche.push_back(softdropsubjet.bDiscriminator("pfTrackCountingHighEffBJetTags") );
           // vSoftDropSubjetjp.push_back(softdropsubjet.bDiscriminator("pfJetProbabilityBJetTags") );
           // vSoftDropSubjetjbp.push_back(softdropsubjet.bDiscriminator("pfJetBProbabilityBJetTags") );
         
         } 

        nBranches_->subjetAK8_softdrop_N.push_back(nsubjets);
        nBranches_->subjetAK8_softdrop_pt.push_back(vSoftDropSubjetpt);
        nBranches_->subjetAK8_softdrop_eta.push_back(vSoftDropSubjeteta);
        nBranches_->subjetAK8_softdrop_mass.push_back(vSoftDropSubjetmass);
        nBranches_->subjetAK8_softdrop_phi.push_back(vSoftDropSubjetphi);
        nBranches_->subjetAK8_softdrop_e.push_back(vSoftDropSubjete);
        nBranches_->subjetAK8_softdrop_charge.push_back(vSoftDropSubjetcharge);
        // nBranches_->subjetAK8_softdrop_ssv.push_back(vSoftDropSubjetssv);
        nBranches_->subjetAK8_softdrop_csv.push_back(vSoftDropSubjetcsv);
        // nBranches_->subjetAK8_softdrop_tchp.push_back(vSoftDropSubjettchp);
        // nBranches_->subjetAK8_softdrop_tche.push_back(vSoftDropSubjettche);
        // nBranches_->subjetAK8_softdrop_jp.push_back(vSoftDropSubjetjp);
        // nBranches_->subjetAK8_softdrop_jbp.push_back(vSoftDropSubjetjbp);
        if(isMC){
          nBranches_->subjetAK8_softdrop_partonFlavour.push_back(vSoftDropSubjetPartonFlavour);
          nBranches_->subjetAK8_softdrop_hadronFlavour.push_back(vSoftDropSubjetHadronFlavour);
        }
    
    
      }
    	                    
      bool IDLoose = looseJetID(fj);
      bool IDTight = tightJetID(fj);
      //if( !IDLoose ) continue;
      
      nBranches_->jetAK8_N++;	       
      nBranches_->jetAK8_pt     	    .push_back(corr*uncorrJet.pt());                   
      nBranches_->jetAK8_eta    	    .push_back(fj.eta());
      nBranches_->jetAK8_mass   	    .push_back(corr*uncorrJet.mass());
      nBranches_->jetAK8_phi    	    .push_back(fj.phi());
      nBranches_->jetAK8_e      	    .push_back(corr*uncorrJet.energy());
      nBranches_->jetAK8_jec    	    .push_back(corr);
      nBranches_->jetAK8_IDLoose      	    .push_back(IDLoose);
      nBranches_->jetAK8_IDTight      	    .push_back(IDTight);
      nBranches_->jetAK8_muf     	    .push_back(fj.muonEnergyFraction());
      nBranches_->jetAK8_phf     	    .push_back(fj.photonEnergyFraction());
      nBranches_->jetAK8_emf     	    .push_back(fj.chargedEmEnergyFraction());
      nBranches_->jetAK8_nhf     	    .push_back(fj.neutralHadronEnergyFraction());
      nBranches_->jetAK8_chf     	    .push_back(fj.chargedHadronEnergyFraction());
      nBranches_->jetAK8_area               .push_back(fj.jetArea());
      nBranches_->jetAK8_cm     	    .push_back(fj.chargedMultiplicity());
      nBranches_->jetAK8_nm     	    .push_back(fj.neutralMultiplicity());     				       
      nBranches_->jetAK8_che     	    .push_back(fj.chargedHadronEnergy()+fj.electronEnergy()+fj.muonEnergy());
      nBranches_->jetAK8_ne     	    .push_back(fj.neutralHadronEnergy()+fj.photonEnergy());
      
      nBranches_->jetAK8_hf_hf  	    .push_back(fj.HFHadronEnergyFraction());
      nBranches_->jetAK8_hf_emf  	    .push_back(fj.HFEMEnergyFraction());
      nBranches_->jetAK8_hof     	    .push_back(fj.hoEnergyFraction());
      
      nBranches_->jetAK8_chm     	    .push_back(fj.chargedHadronMultiplicity());
      nBranches_->jetAK8_neHadMult          .push_back(fj.neutralHadronMultiplicity());
      nBranches_->jetAK8_phoMult            .push_back(fj.photonMultiplicity());
      
      nBranches_->jetAK8_nemf    	    .push_back(fj.neutralEmEnergyFraction());
      nBranches_->jetAK8_cemf    	    .push_back(fj.chargedEmEnergyFraction());
      
      nBranches_->jetAK8_charge 	    .push_back(fj.charge());					 
      nBranches_->jetAK8_Hbbtag             .push_back(fj.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags"));       
      // nBranches_->jetAK8_ssv          .push_back(fj.bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags"));
      nBranches_->jetAK8_csv                .push_back(fj.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
      // nBranches_->jetAK8_tchp         .push_back(fj.bDiscriminator("pfTrackCountingHighPurBJetTags"));
      // nBranches_->jetAK8_tche         .push_back(fj.bDiscriminator("pfTrackCountingHighEffBJetTags"));
      // nBranches_->jetAK8_jp           .push_back(fj.bDiscriminator("pfJetProbabilityBJetTags"));
      // nBranches_->jetAK8_jbp          .push_back(fj.bDiscriminator("pfJetBProbabilityBJetTags"));
      nBranches_->jetAK8_tau1            .push_back(fj.userFloat("NjettinessAK8:tau1"));	   
      nBranches_->jetAK8_tau2            .push_back(fj.userFloat("NjettinessAK8:tau2"));
      nBranches_->jetAK8_tau3            .push_back(fj.userFloat("NjettinessAK8:tau3")); 
      nBranches_->jetAK8_pruned_mass     .push_back(fj.userFloat("ak8PFJetsCHSPrunedMass"));
      nBranches_->jetAK8_softdrop_mass   .push_back(fj.userFloat("ak8PFJetsCHSSoftDropMass"));
      if(isMC){
        int genP_pdgId_fj = fj.genParton() ? fj.genParton()->pdgId() : -99; //det default to -99 when no genParton is found!!!
        nBranches_->jetAK8_partonFlavour.push_back(fj.partonFlavour()); // b,c can be identified by either j.hadronFlavour() (hadron definition) or j.partonFlavour() (physics definition)
        nBranches_->jetAK8_hadronFlavour.push_back(fj.hadronFlavour()); // g,u,d,s can be identified by j.partonFlavour() (physics definition)
        nBranches_->jetAK8_genParton_pdgID  .push_back(genP_pdgId_fj);
        nBranches_->jetAK8_nbHadrons         .push_back((fj.jetFlavourInfo().getbHadrons()).size()); //gluon splitting to bb can be identified roughly by (j.partonFlavour()==21 &% len(j.jetFlavourInfo().getbHadrons())==2 )
        nBranches_->jetAK8_ncHadrons         .push_back((fj.jetFlavourInfo().getcHadrons()).size());
      }
     
      TLorentzVector FatJet; FatJet.SetPtEtaPhiE( fj.pt(), fj.eta(), fj.phi(), fj.energy() );  

      if( doPruning ){

        // nBranches_->jetAK8_pruned_N = 0;

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
        reco::Candidate::LorentzVector uncorrPrunedJet;
        double prunedcorr = 1;
      
        if( doCorrOnTheFly_ ){
      
           uncorrPrunedJet = prunedjet.correctedP4(0);

           jecAK8Groomed_->setJetEta( uncorrPrunedJet.eta()    );
           jecAK8Groomed_->setJetPt ( uncorrPrunedJet.pt()     );
           jecAK8Groomed_->setJetE  ( uncorrPrunedJet.energy() );
           jecAK8Groomed_->setJetA  ( prunedjet.jetArea()      );
           jecAK8Groomed_->setRho   ( nBranches_->rho          );
           jecAK8Groomed_->setNPV   ( vertices_->size()        );
           prunedcorr = jecAK8Groomed_->getCorrection();
           nBranches_->jetAK8_pruned_massCorr.push_back(prunedcorr*fj.userFloat("ak8PFJetsCHSPrunedMass"));
           nBranches_->jetAK8_pruned_jec.push_back(prunedcorr);
         
        }
        else{
           nBranches_->jetAK8_pruned_massCorr.push_back(fj.userFloat("ak8PFJetsCHSPrunedMassCorrected"));
           nBranches_->jetAK8_pruned_jec.push_back(fj.userFloat("ak8PFJetsCHSPrunedMassCorrected")/fj.userFloat("ak8PFJetsCHSPrunedMass"));
        }

        /****************************************************************/

        // nBranches_->jetAK8_pruned_pt     .push_back(prunedjet.pt());
        // nBranches_->jetAK8_pruned_eta    .push_back(prunedjet.eta());
        // nBranches_->jetAK8_pruned_mass   .push_back(prunedjet.mass());
        // nBranches_->jetAK8_pruned_phi    .push_back(prunedjet.phi());
        // nBranches_->jetAK8_pruned_e      .push_back(prunedjet.energy());
        // nBranches_->jetAK8_pruned_jec    .push_back(prunedjet.mass());
        // nBranches_->jetAK8_pruned_flavour.push_back(abs(prunedjet.partonFlavour()));
        // nBranches_->jetAK8_pruned_charge .push_back(prunedjet.charge());
        // nBranches_->jetAK8_pruned_ssv    .push_back(prunedjet.bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags"));
        // nBranches_->jetAK8_pruned_csv    .push_back(prunedjet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
        // nBranches_->jetAK8_pruned_tchp   .push_back(prunedjet.bDiscriminator("pfTrackCountingHighPurBJetTags"));
        // nBranches_->jetAK8_pruned_tche   .push_back(prunedjet.bDiscriminator("pfTrackCountingHighEffBJetTags"));
        // nBranches_->jetAK8_pruned_jp    .push_back(prunedjet.bDiscriminator("pfJetProbabilityBJetTags"));
        // nBranches_->jetAK8_pruned_jbp    .push_back(prunedjet.bDiscriminator("pfJetBProbabilityBJetTags"));
        // nBranches_->jetAK8_pruned_N++;

        vPrunedSubjetpt     .clear();
        vPrunedSubjeteta    .clear();
        vPrunedSubjetmass   .clear();
        vPrunedSubjetphi    .clear();
        vPrunedSubjete      .clear();
        vPrunedSubjetcharge .clear();
        vPrunedSubjetPartonFlavour.clear();
        vPrunedSubjetHadronFlavour.clear();
        // vPrunedSubjetssv    .clear();
        vPrunedSubjetcsv    .clear();
        // vPrunedSubjettchp   .clear();
        // vPrunedSubjettche   .clear();
        // vPrunedSubjetjp     .clear();
        // vPrunedSubjetjbp    .clear();

        //Get subjets from matched pruned jet
        std::vector<reco::CandidatePtr> pruneddaus(prunedjet.daughterPtrVector());

        nsubjets = 0;

        for (unsigned int k = 0; k < pruneddaus.size(); ++k) {

          // const reco::PFJet &prunedsubjet = dynamic_cast<const reco::PFJet &>(*pruneddaus[k]);
        
          const pat::Jet &prunedsubjet = dynamic_cast<const pat::Jet &>(*pruneddaus[k]);   

          if( prunedsubjet.pt() < 0.01 ) continue;

          nsubjets++;

          // vPrunedSubjetPartonFlavour.push_back(flavor);
          vPrunedSubjetpt.push_back(prunedsubjet.pt());
          vPrunedSubjeteta.push_back(prunedsubjet.eta());
          vPrunedSubjetmass.push_back(prunedsubjet.mass());
          vPrunedSubjetphi.push_back(prunedsubjet.phi());
          vPrunedSubjete.push_back(prunedsubjet.energy());
          vPrunedSubjetPartonFlavour.push_back(abs(prunedsubjet.partonFlavour()));
          vPrunedSubjetHadronFlavour.push_back(abs(prunedsubjet.hadronFlavour()));
          vPrunedSubjetcharge.push_back(prunedsubjet.charge());
          // vPrunedSubjetssv.push_back(prunedsubjet.bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags"));
          vPrunedSubjetcsv.push_back(prunedsubjet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
          // vPrunedSubjettchp.push_back(prunedsubjet.bDiscriminator("pfTrackCountingHighPurBJetTags") );
          // vPrunedSubjettche.push_back(prunedsubjet.bDiscriminator("pfTrackCountingHighEffBJetTags") );
          // vPrunedSubjetjp.push_back(prunedsubjet.bDiscriminator("pfJetProbabilityBJetTags") );
          // vPrunedSubjetjbp.push_back(prunedsubjet.bDiscriminator("pfJetBProbabilityBJetTags") );   
        
        }

        nBranches_->subjetAK8_pruned_N.push_back(nsubjets                     );
        nBranches_->subjetAK8_pruned_pt.push_back(vPrunedSubjetpt          );
        nBranches_->subjetAK8_pruned_eta.push_back(vPrunedSubjeteta        );
        nBranches_->subjetAK8_pruned_mass.push_back(vPrunedSubjetmass      );
        nBranches_->subjetAK8_pruned_phi.push_back(vPrunedSubjetphi        );
        nBranches_->subjetAK8_pruned_e.push_back(vPrunedSubjete            );
        nBranches_->subjetAK8_pruned_charge.push_back(vPrunedSubjetcharge  );
        // nBranches_->subjetAK8_pruned_ssv.push_back(vPrunedSubjetssv        );
        nBranches_->subjetAK8_pruned_csv.push_back(vPrunedSubjetcsv        );
        // nBranches_->subjetAK8_pruned_tchp.push_back(vPrunedSubjettchp      );
        // nBranches_->subjetAK8_pruned_tche.push_back(vPrunedSubjettche      );
        // nBranches_->subjetAK8_pruned_jp.push_back(vPrunedSubjetjp          );
        // nBranches_->subjetAK8_pruned_jbp.push_back(vPrunedSubjetjbp        );
        if(isMC){
          nBranches_->subjetAK8_pruned_partonFlavour.push_back(vPrunedSubjetPartonFlavour);
          nBranches_->subjetAK8_pruned_hadronFlavour.push_back(vPrunedSubjetHadronFlavour);          
        }
        
      }
      else{
        nBranches_->jetAK8_pruned_massCorr.push_back(-99);
        nBranches_->jetAK8_pruned_jec.push_back(-99);
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

        reco::Candidate::LorentzVector uncorrSoftDropJet;
        double softdropcorr = 1;
      
        if( doCorrOnTheFly_ ){
      
           uncorrSoftDropJet = softdropjet.correctedP4(0);

           jecAK8Groomed_->setJetEta( uncorrSoftDropJet.eta()    );
           jecAK8Groomed_->setJetPt ( uncorrSoftDropJet.pt()     );
           jecAK8Groomed_->setJetE  ( uncorrSoftDropJet.energy() );
           jecAK8Groomed_->setJetA  ( softdropjet.jetArea()      );
           jecAK8Groomed_->setRho   ( nBranches_->rho            );
           jecAK8Groomed_->setNPV   ( vertices_->size()          );
           softdropcorr = jecAK8Groomed_->getCorrection();
           nBranches_->jetAK8_softdrop_massCorr.push_back(softdropcorr*fj.userFloat("ak8PFJetsCHSSoftDropMass"));
           nBranches_->jetAK8_softdrop_jec  .push_back(softdropcorr);
        }
        else{
           nBranches_->jetAK8_softdrop_massCorr.push_back(fj.userFloat("ak8PFJetsCHSSoftDropMassCorrected"));
           nBranches_->jetAK8_softdrop_jec  .push_back(fj.userFloat("ak8PFJetsCHSSoftDropMassCorrected")/fj.userFloat("ak8PFJetsCHSSoftDropMass"));
        }
      
        vSoftDropSubjetpt.clear();
        vSoftDropSubjeteta.clear();
        vSoftDropSubjetmass.clear();
        vSoftDropSubjetphi.clear();
        vSoftDropSubjete.clear();
        vSoftDropSubjetcharge.clear();
        vSoftDropSubjetPartonFlavour.clear();
        vSoftDropSubjetHadronFlavour.clear();
        // vSoftDropSubjetssv.clear();
        vSoftDropSubjetcsv.clear();
        // vSoftDropSubjettchp.clear();
        // vSoftDropSubjettche.clear();
        // vSoftDropSubjetjp.clear();
        // vSoftDropSubjetjbp.clear();

        //Get subjets from matched pruned jet
        std::vector<reco::CandidatePtr> softdropdaus(softdropjet.daughterPtrVector());

        nsubjets = 0;

        for (unsigned int k = 0; k < softdropdaus.size(); ++k) {

          // const reco::PFJet &softdropsubjet = dynamic_cast<const reco::PFJet &>(*softdropdaus[k]);
        
          const pat::Jet &softdropsubjet = dynamic_cast<const pat::Jet &>(*softdropdaus[k]);   

          if( softdropsubjet.pt() < 0.01 ) continue;

          nsubjets++;

          // vSoftDropSubjetPartonFlavour.push_back(flavor);
          vSoftDropSubjetpt.push_back(softdropsubjet.pt());
          vSoftDropSubjeteta.push_back(softdropsubjet.eta());
          vSoftDropSubjetmass.push_back(softdropsubjet.mass());
          vSoftDropSubjetphi.push_back(softdropsubjet.phi());
          vSoftDropSubjete.push_back(softdropsubjet.energy());
          vSoftDropSubjetPartonFlavour.push_back(abs(softdropsubjet.partonFlavour()));
          vSoftDropSubjetHadronFlavour.push_back(abs(softdropsubjet.hadronFlavour()));
          vSoftDropSubjetcharge.push_back(softdropsubjet.charge());
          // vSoftDropSubjetssv.push_back(softdropsubjet.bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags"));
          vSoftDropSubjetcsv.push_back(softdropsubjet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
          // vSoftDropSubjettchp.push_back(softdropsubjet.bDiscriminator("pfTrackCountingHighPurBJetTags") );
          // vSoftDropSubjettche.push_back(softdropsubjet.bDiscriminator("pfTrackCountingHighEffBJetTags") );
          // vSoftDropSubjetjp.push_back(softdropsubjet.bDiscriminator("pfJetProbabilityBJetTags") );
          // vSoftDropSubjetjbp.push_back(softdropsubjet.bDiscriminator("pfJetBProbabilityBJetTags") );
        }

        nBranches_->subjetAK8_softdrop_N.push_back(nsubjets);
        nBranches_->subjetAK8_softdrop_pt.push_back(vSoftDropSubjetpt);
        nBranches_->subjetAK8_softdrop_eta.push_back(vSoftDropSubjeteta);
        nBranches_->subjetAK8_softdrop_mass.push_back(vSoftDropSubjetmass);
        nBranches_->subjetAK8_softdrop_phi.push_back(vSoftDropSubjetphi);
        nBranches_->subjetAK8_softdrop_e.push_back(vSoftDropSubjete);
        nBranches_->subjetAK8_softdrop_charge.push_back(vSoftDropSubjetcharge);
        // nBranches_->subjetAK8_softdrop_ssv.push_back(vSoftDropSubjetssv);
        nBranches_->subjetAK8_softdrop_csv.push_back(vSoftDropSubjetcsv);
        // nBranches_->subjetAK8_softdrop_tchp.push_back(vSoftDropSubjettchp);
        // nBranches_->subjetAK8_softdrop_tche.push_back(vSoftDropSubjettche);
        // nBranches_->subjetAK8_softdrop_jp.push_back(vSoftDropSubjetjp);
        // nBranches_->subjetAK8_softdrop_jbp.push_back(vSoftDropSubjetjbp);
        if(isMC){
          nBranches_->subjetAK8_softdrop_partonFlavour.push_back(vSoftDropSubjetPartonFlavour);
          nBranches_->subjetAK8_softdrop_hadronFlavour.push_back(vSoftDropSubjetHadronFlavour);          
        }
      
      }
      else{
        nBranches_->jetAK8_softdrop_massCorr.push_back(-99);
        nBranches_->jetAK8_softdrop_jec  .push_back(-99);
      } 

      /****************************************************************/
      if( doPuppi ){

        // nBranches_->jetAK8_puppi__N = 0;

        // Loop over  jets, store dR match. Used to obtain subjets.
        float dRmin =  999. ;

        pat::Jet puppijet;

        for (const pat::Jet &pj : *puppijets_) {

          TLorentzVector jet; jet.SetPtEtaPhiE( pj.pt(), pj.eta(), pj.phi(), pj.energy() );

          float dRtmp   = FatJet.DeltaR(jet);
          if( dRtmp < dRmin && dRtmp < 0.8 ){ dRmin = dRtmp; puppijet = pj;}
          else continue;

        }

        nBranches_->jetAK8_puppi_tau1	 .push_back(puppijet.userFloat("NjettinessAK8Puppi:tau1"));	 
        nBranches_->jetAK8_puppi_tau2	 .push_back(puppijet.userFloat("NjettinessAK8Puppi:tau2"));
        nBranches_->jetAK8_puppi_tau3	 .push_back(puppijet.userFloat("NjettinessAK8Puppi:tau3")); 
        nBranches_->jetAK8_puppi_pruned_mass.push_back(puppijet.userFloat("ak8PFJetsPuppiPrunedMass"));
        nBranches_->jetAK8_puppi_softdrop_mass.push_back(puppijet.userFloat("ak8PFJetsPuppiSoftDropMass"));

        //Compute JEC for pruned mass
        reco::Candidate::LorentzVector uncorrJet;
        double corr = 1;
      
        if( doCorrOnTheFly_ ){
           if(puppijet.jecSetsAvailable())
             uncorrJet = puppijet.correctedP4(0);
	   else
             uncorrJet = puppijet.p4();

           jecAK8Puppi_->setJetEta( uncorrJet.eta()    );
           jecAK8Puppi_->setJetPt ( uncorrJet.pt()     );
           jecAK8Puppi_->setJetE  ( uncorrJet.energy() );
           jecAK8Puppi_->setJetA  ( puppijet.jetArea()      );
           jecAK8Puppi_->setRho   ( nBranches_->rho          );
           jecAK8Puppi_->setNPV   ( vertices_->size()        );
           corr = jecAK8Puppi_->getCorrection();
           nBranches_->jetAK8_puppi_pruned_massCorr.push_back(corr*puppijet.userFloat("ak8PFJetsPuppiPrunedMass"));
           nBranches_->jetAK8_puppi_pruned_jec.push_back(corr);
           nBranches_->jetAK8_puppi_softdrop_massCorr.push_back(corr*puppijet.userFloat("ak8PFJetsPuppiSoftDropMass"));
           nBranches_->jetAK8_puppi_softdrop_jec.push_back(corr);
        }
        else{
           nBranches_->jetAK8_puppi_pruned_massCorr.push_back(puppijet.userFloat("ak8PFJetsPuppiPrunedMassCorrected"));
           nBranches_->jetAK8_puppi_pruned_jec.push_back(puppijet.userFloat("ak8PFJetsPuppiPrunedMassCorrected")/puppijet.userFloat("ak8PFJetsPuppiPrunedMass"));
           nBranches_->jetAK8_puppi_softdrop_massCorr.push_back(puppijet.userFloat("ak8PFJetsPuppiSoftDropMassCorrected"));
           nBranches_->jetAK8_puppi_softdrop_jec.push_back(puppijet.userFloat("ak8PFJetsPuppiSoftDropMassCorrected")/puppijet.userFloat("ak8PFJetsPuppiSoftDropMass"));
        }

      }

    } // ak8 jet loop

  } //doAK8Jets

  if( doTrimming ){

    for (const pat::Jet &tj : *trimmedjets_) {

        reco::Candidate::LorentzVector uncorrTrimmedJet = tj.correctedP4(0);

        nBranches_->jetAK10_trimmed_mass .push_back(uncorrTrimmedJet.mass());
        nBranches_->jetAK10_ecf1	 .push_back(tj.userFloat("ECFAK10:ecf1"));	 
        nBranches_->jetAK10_ecf2	 .push_back(tj.userFloat("ECFAK10:ecf2"));
        nBranches_->jetAK10_ecf3	 .push_back(tj.userFloat("ECFAK10:ecf3")); 

        //Compute JEC for trimmed mass
        double trimmedcorr = 1;
      
        if( doCorrOnTheFly_ ){
      
           jecAK8Groomed_->setJetEta( uncorrTrimmedJet.eta()    );
           jecAK8Groomed_->setJetPt ( uncorrTrimmedJet.pt()     );
           jecAK8Groomed_->setJetE  ( uncorrTrimmedJet.energy() );
           jecAK8Groomed_->setJetA  ( tj.jetArea()      );
           jecAK8Groomed_->setRho   ( nBranches_->rho          );
           jecAK8Groomed_->setNPV   ( vertices_->size()        );
           trimmedcorr = jecAK8Groomed_->getCorrection();
           nBranches_->jetAK10_trimmed_massCorr.push_back(trimmedcorr*uncorrTrimmedJet.mass());
           nBranches_->jetAK10_trimmed_jec.push_back(trimmedcorr);
         
        }
        else{
           nBranches_->jetAK10_trimmed_massCorr.push_back(tj.mass());
	   if(uncorrTrimmedJet.mass()>0)
           nBranches_->jetAK10_trimmed_jec.push_back(tj.mass()/uncorrTrimmedJet.mass());
	   else
           nBranches_->jetAK10_trimmed_jec.push_back(1.);
        }
    }
  }

}



