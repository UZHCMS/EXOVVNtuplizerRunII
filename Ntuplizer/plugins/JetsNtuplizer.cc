#include "../interface/JetsNtuplizer.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include <CondFormats/JetMETObjects/interface/JetResolutionObject.h>
//===================================================================================================================        

JetsNtuplizer::JetsNtuplizer( std::vector<edm::EDGetTokenT<pat::JetCollection>> tokens, std::vector<std::string> jecAK4Labels, std::vector<std::string> jecAK8Labels, std::vector<std::string> jecAK8GroomedLabels, std::vector<std::string> jecAK8PuppiLabels, edm::EDGetTokenT<reco::JetFlavourMatchingCollection> flavourToken, edm::EDGetTokenT<double> rhoToken, edm::EDGetTokenT<reco::VertexCollection> verticeToken, NtupleBranches* nBranches, std::map< std::string, bool >& runFlags,  std::vector<std::string>   jerAK8chsFileLabel, std::vector<std::string>   jerAK4chsFileLabel, std::vector<std::string>   jerAK8PuppiFileLabel,  std::vector<std::string>   jerAK4PuppiFileLabel )

  : CandidateNtuplizer      ( nBranches )
  , jetInputToken_	    ( tokens[0] )
  , fatjetInputToken_	    ( tokens[1] )
  , prunedjetInputToken_    ( tokens[2] )
  , softdropjetInputToken_  ( tokens[3] )
  , trimmedjetInputToken_   ( tokens[4] )
  , puppijetInputToken_     ( tokens[5] )
  , rhoToken_	       	    ( rhoToken  )
  , verticeToken_     	    ( verticeToken )
  , doAK4Jets_ (runFlags["doAK4Jets"])
  , doAK8Jets_ (runFlags["doAK8Jets"])
  , doPuppiRecluster_ (runFlags["doPuppiRecluster"])
  , runOnMC_   (runFlags["runOnMC"])
   //, flavourToken_			( flavourToken 	) //For subjet flavour matching!! Not done yet.
    
{
	
  doCorrOnTheFly_ = false;	
  if( jecAK4Labels.size() > 1 && jecAK8Labels.size() > 1 && jecAK8GroomedLabels.size() != 0 && jecAK8PuppiLabels.size() != 0 ){	
  
     jecAK4PayloadNames_ = jecAK4Labels;
     jecAK4PayloadNames_.pop_back();
	
     jecAK8PayloadNames_ = jecAK8Labels;
     jecAK8PayloadNames_.pop_back();

     jecAK8GroomedPayloadNames_ = jecAK8GroomedLabels;
     //jecAK8GroomedPayloadNames_.pop_back();

     jecAK8PuppiPayloadNames_ = jecAK8PuppiLabels;
     //jecAK8PuppiPayloadNames_.pop_back();
      
     initJetCorrFactors();
     
     doCorrOnTheFly_ = true;	
     
  }

  jecAK4UncName_ = jecAK4Labels.back();
  jecAK8UncName_ = jecAK8Labels.back();

  jerAK4chsName_res_ = jerAK4chsFileLabel.front();
  jerAK8chsName_res_ = jerAK8chsFileLabel.front();
  jerAK4chsName_sf_ = jerAK4chsFileLabel.back();
  jerAK8chsName_sf_ = jerAK8chsFileLabel.back();
  
  jerAK8PuppiName_res_ = jerAK8PuppiFileLabel.front();
  jerAK8PuppiName_sf_ = jerAK8PuppiFileLabel.back();
  initJetCorrUncertainty();
     
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

//  double eta = j.eta();		
//  double chf = j.chargedHadronEnergyFraction();
//  double nhf = j.neutralHadronEnergyFraction(); // + j.HFHadronEnergyFraction();
//  double muf = j.muonEnergy()/(j.jecFactor(0) * j.energy());  
//  double nemf = j.neutralEmEnergyFraction();
//  double cemf = j.chargedEmEnergyFraction();
//  int chMult = j.chargedMultiplicity();
//  int neMult = j.neutralMultiplicity();
//  int npr    = chMult + neMult;
//  int NumConst = npr;
//
//  return (nhf<0.99 && nemf<0.99 && NumConst>1 && muf < 0.8) && ((fabs(eta) <= 2.4 && chf>0 && chMult>0 && cemf<0.99) || fabs(eta)>2.4);      

  // Change to 13 TeV definition
  // https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data

  double eta = j.eta();		
  double chf = j.chargedHadronEnergyFraction();
  double nhf = j.neutralHadronEnergyFraction(); // + j.HFHadronEnergyFraction();
  double nemf = j.neutralEmEnergyFraction();
  double cemf = j.chargedEmEnergyFraction();
  int chMult = j.chargedMultiplicity();
  int neMult = j.neutralMultiplicity();
  int npr    = chMult + neMult;
  int NumConst = npr;

  if(abs(eta) <= 2.7){
    return (nhf<0.99 && nemf<0.99 && NumConst>1) && ((abs(eta)<=2.4 && chf>0 && chMult>0 && cemf<0.99) || abs(eta)>2.4);
  }else if(abs(eta) <= 3.0){
    return (nemf<0.90 && neMult>2);
  }else{
    return (nemf<0.90 && neMult>10);
  }

      		
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
     
  if(abs(eta) <= 2.7){
    return (nhf<0.90 && nemf<0.90 && NumConst>1 && muf<0.8) && ((fabs(eta)<=2.4 && chf>0 && chMult>0 && cemf<0.90) || fabs(eta)>2.4);  		
  }else if(abs(eta) <= 3.0){
    return (nemf<0.90 && neMult>2);
  }else{
    return (nemf<0.90 && neMult>10);
  }

}



bool JetsNtuplizer::tightJetIDWithoutLepVeto( const pat::Jet& j ) {

  // Change to 13 TeV definition
  // https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data

  double eta = j.eta();		
  double chf = j.chargedHadronEnergyFraction();
  double nhf = j.neutralHadronEnergyFraction(); // + j.HFHadronEnergyFraction();
  double nemf = j.neutralEmEnergyFraction();
  double cemf = j.chargedEmEnergyFraction();
  int chMult = j.chargedMultiplicity();
  int neMult = j.neutralMultiplicity();
  int npr    = chMult + neMult;
  int NumConst = npr;

  if(abs(eta) <= 2.7){
    return (nhf<0.90 && nemf<0.90 && NumConst>1) && ((abs(eta)<=2.4 && chf>0 && chMult>0 && cemf<0.99) || abs(eta)>2.4);
  }else if(abs(eta) <= 3.0){
    return (nemf<0.90 && neMult>2);
  }else{
    return (nemf<0.90 && neMult>10);
  }


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
void JetsNtuplizer::initJetCorrUncertainty( void ){

  jecAK8Unc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecAK8UncName_) );
  jecAK4Unc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecAK4UncName_) );

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
  bool doPuppi  = doPuppiRecluster_;
  if (doPuppi) doPuppi=event.getByToken(puppijetInputToken_, puppijets_ );
  bool isMC = runOnMC_;

  /****************************************************************/
  if (doAK4Jets_) {
  

    JME::JetResolution resolution = JME::JetResolution(jerAK4chsName_res_);
    JME::JetResolutionScaleFactor resolution_sf = JME::JetResolutionScaleFactor(jerAK4chsName_sf_);


    nBranches_->jetAK4_N = 0;
  
    for (const pat::Jet &j : *jets_) {
	        
     bool IDLoose = looseJetID(j);
     bool IDTight = tightJetID(j);
     bool IDTightWithoutLepVeto = tightJetIDWithoutLepVeto(j);

     reco::Candidate::LorentzVector uncorrJet;
     double corr = 1;
     double corrUp = 1;
     double corrDown = 1;   
     if( doCorrOnTheFly_ ){
   
       uncorrJet = j.correctedP4(0);
       jecAK4_->setJetEta( uncorrJet.eta()    );
       jecAK4_->setJetPt ( uncorrJet.pt()     );
       jecAK4_->setJetE  ( uncorrJet.energy() );
       jecAK4_->setRho	( nBranches_->rho  );
       jecAK4_->setNPV	( vertices_->size()  );
       jecAK4_->setJetA  ( j.jetArea()	     );
       corr = jecAK4_->getCorrection();
    
     }
     else{
       uncorrJet = j.p4();
     }

     jecAK4Unc_->setJetEta( j.correctedP4(0).eta() );
     jecAK4Unc_->setJetPt( corr * j.correctedP4(0).pt() );
     corrUp = corr * (1 + fabs(jecAK4Unc_->getUncertainty(1)));
     jecAK4Unc_->setJetEta( j.correctedP4(0).eta() );
     jecAK4Unc_->setJetPt( corr * j.correctedP4(0).pt() );
     corrDown = corr * ( 1 - fabs(jecAK4Unc_->getUncertainty(-1)) );

     if (corr*uncorrJet.pt() < 20) continue;    

     nBranches_->jetAK4_N++;

     nBranches_->jetAK4_pt	  .push_back(corr*uncorrJet.pt());	
     nBranches_->jetAK4_eta	  .push_back(j.eta());
     nBranches_->jetAK4_mass	  .push_back(corr*uncorrJet.mass());
     nBranches_->jetAK4_phi	  .push_back(j.phi());   
     nBranches_->jetAK4_e	  .push_back(corr*uncorrJet.energy());
     nBranches_->jetAK4_jec	  .push_back(corr);
     nBranches_->jetAK4_jecUp 	  .push_back(corrUp);
     nBranches_->jetAK4_jecDown	  .push_back(corrDown);
     nBranches_->jetAK4_IDLoose   .push_back(IDLoose);
     nBranches_->jetAK4_IDTight   .push_back(IDTight);
     nBranches_->jetAK4_IDTightWithoutLepVeto   .push_back(IDTightWithoutLepVeto);
     nBranches_->jetAK4_PUIDdiscriminat.push_back(j.userFloat("pileupJetIdUpdated:fullDiscriminant")); 
     int fullId=j.userInt("pileupJetIdUpdated:fullId");
     nBranches_->jetAK4_PUIDloose .push_back(fullId & (1 << 2)); 
     nBranches_->jetAK4_PUIDmedium.push_back(fullId & (1 << 1)); 
     nBranches_->jetAK4_PUIDtight .push_back(fullId & (1 << 0)); 
     nBranches_->jetAK4_cm	  .push_back(j.chargedMultiplicity());
     nBranches_->jetAK4_nm	  .push_back(j.neutralMultiplicity());
     nBranches_->jetAK4_muf	  .push_back(j.muonEnergyFraction());
     nBranches_->jetAK4_phf	  .push_back(j.photonEnergyFraction());
     nBranches_->jetAK4_emf	  .push_back(j.chargedEmEnergyFraction());
     nBranches_->jetAK4_nhf	  .push_back(j.neutralHadronEnergyFraction());
     nBranches_->jetAK4_chf	  .push_back(j.chargedHadronEnergyFraction());
     nBranches_->jetAK4_area	  .push_back(j.jetArea());	
     nBranches_->jetAK4_che	  .push_back(j.chargedHadronEnergy()+j.electronEnergy()+j.muonEnergy());
     nBranches_->jetAK4_ne	  .push_back(j.neutralHadronEnergy()+j.photonEnergy());   
     nBranches_->jetAK4_hf_hf	  .push_back(j.HFHadronEnergyFraction());
     nBranches_->jetAK4_hf_emf    .push_back(j.HFEMEnergyFraction());
     nBranches_->jetAK4_hof	  .push_back(j.hoEnergyFraction());   
     nBranches_->jetAK4_chm	  .push_back(j.chargedHadronMultiplicity());
     nBranches_->jetAK4_neHadMult .push_back(j.neutralHadronMultiplicity());
     nBranches_->jetAK4_phoMult   .push_back(j.photonMultiplicity());	
     nBranches_->jetAK4_nemf	  .push_back(j.neutralEmEnergyFraction());
     nBranches_->jetAK4_cemf	  .push_back(j.chargedEmEnergyFraction());   
     nBranches_->jetAK4_charge    .push_back(j.charge());
     nBranches_->jetAK4_csv	  .push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
     nBranches_->jetAK4_vtxMass   .push_back(j.userFloat("vtxMass")); 
     nBranches_->jetAK4_vtxNtracks.push_back(j.userFloat("vtxNtracks")); 
     nBranches_->jetAK4_vtx3DVal  .push_back(j.userFloat("vtx3DVal")); 
     nBranches_->jetAK4_vtx3DSig  .push_back(j.userFloat("vtx3DSig"));
     
     if(isMC){
     
       int genP_pdgId = j.genParton() ? j.genParton()->pdgId() : -99; // default to -99 when no genParton is found!!!
       nBranches_->jetAK4_partonFlavour  .push_back(j.partonFlavour()); //Algorithmic definition! g,u,d,s can be identified by j.partonFlavour() 
       nBranches_->jetAK4_hadronFlavour  .push_back(j.hadronFlavour()); //Hadron flavour. b,c can be identified by either j.hadronFlavour() (hadron definition) or j.partonFlavour() (parton definition)
       nBranches_->jetAK4_genParton_pdgID.push_back(genP_pdgId); //Physics definition! see https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#MC_Truth
       nBranches_->jetAK4_nbHadrons	 .push_back((j.jetFlavourInfo().getbHadrons()).size()); 
       nBranches_->jetAK4_ncHadrons	 .push_back((j.jetFlavourInfo().getcHadrons()).size());
       

       JME::JetParameters parameters;
       parameters.setJetPt(corr*uncorrJet.pt());
       parameters.setJetEta(j.eta());
       parameters.setRho(nBranches_->rho);
       float jer_res = resolution.getResolution(parameters);
       float jer_sf = resolution_sf.getScaleFactor(parameters);
       float jer_sf_up = resolution_sf.getScaleFactor(parameters, Variation::UP);
       float jer_sf_down = resolution_sf.getScaleFactor(parameters, Variation::DOWN);
       nBranches_->jetAK4_jer_sf.push_back(jer_sf);
       nBranches_->jetAK4_jer_sf_up.push_back(jer_sf_up);
       nBranches_->jetAK4_jer_sf_down.push_back(jer_sf_down);
       nBranches_->jetAK4_jer_sigma_pt.push_back(jer_res);
     }
     
    }//close ak4 jets loop
    
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
    std::vector<int  > vPrunedSubjetGenPartonPdgId;
    std::vector<int  > vPrunedSubjetnbHadrons;
    std::vector<int  > vPrunedSubjetncHadrons;
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
    std::vector<int  > vSoftDropSubjetGenPartonPdgId;
    std::vector<int  > vSoftDropSubjetnbHadrons;
    std::vector<int  > vSoftDropSubjetncHadrons;
    std::vector<int  > vSoftDropSubjetPartonFlavour;
    std::vector<int  > vSoftDropSubjetHadronFlavour;
    std::vector<float> vSoftDropSubjetssv    ;
    std::vector<float> vSoftDropSubjetcsv    ;
    std::vector<float> vSoftDropSubjettchp   ;
    std::vector<float> vSoftDropSubjettche   ;
    std::vector<float> vSoftDropSubjetjp     ;
    std::vector<float> vSoftDropSubjetjbp    ;
     
    std::vector<float> vPuppiSoftDropSubjetpt     ;
    std::vector<float> vPuppiSoftDropSubjeteta    ;
    std::vector<float> vPuppiSoftDropSubjetmass   ;
    std::vector<float> vPuppiSoftDropSubjetphi    ;
    std::vector<float> vPuppiSoftDropSubjete      ;
    std::vector<int  > vPuppiSoftDropSubjetcharge ;
    std::vector<int  > vPuppiSoftDropSubjetGenPartonPdgId;
    std::vector<int  > vPuppiSoftDropSubjetnbHadrons;
    std::vector<int  > vPuppiSoftDropSubjetncHadrons;
    std::vector<int  > vPuppiSoftDropSubjetPartonFlavour;
    std::vector<int  > vPuppiSoftDropSubjetHadronFlavour;
    std::vector<float> vPuppiSoftDropSubjetssv    ;
    std::vector<float> vPuppiSoftDropSubjetcsv    ;
    std::vector<float> vPuppiSoftDropSubjettchp   ;
    std::vector<float> vPuppiSoftDropSubjettche   ;
    std::vector<float> vPuppiSoftDropSubjetjp     ;
    std::vector<float> vPuppiSoftDropSubjetjbp    ;
     

    JME::JetResolution resolution_ak8 = JME::JetResolution(jerAK8chsName_res_);
    JME::JetResolutionScaleFactor resolution_ak8_sf = JME::JetResolutionScaleFactor(jerAK8chsName_sf_);
    

    for (const pat::Jet &fj : *fatjets_) {
	  
      reco::Candidate::LorentzVector uncorrJet;
      double corr = 1;
      double corrUp = 1;
      double corrDown = 1;
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

      jecAK8Unc_->setJetEta( fj.correctedP4(0).eta() );
      jecAK8Unc_->setJetPt( corr * fj.correctedP4(0).pt() );
      corrUp = corr * (1 + fabs(jecAK8Unc_->getUncertainty(1)));
      jecAK8Unc_->setJetEta( fj.correctedP4(0).eta() );
      jecAK8Unc_->setJetPt( corr * fj.correctedP4(0).pt() );
      corrDown = corr * ( 1 - fabs(jecAK8Unc_->getUncertainty(-1)) );   
      
      if( corr*uncorrJet.pt() < 100. ) continue;

      bool IDLoose = looseJetID(fj);
      bool IDTight = tightJetID(fj);

      nBranches_->jetAK8_N++;	       
      nBranches_->jetAK8_pt     	    .push_back(corr*uncorrJet.pt());                   
      nBranches_->jetAK8_eta    	    .push_back(fj.eta());
      nBranches_->jetAK8_mass   	    .push_back(corr*uncorrJet.mass());
      nBranches_->jetAK8_phi    	    .push_back(fj.phi());
      nBranches_->jetAK8_e      	    .push_back(corr*uncorrJet.energy());
      nBranches_->jetAK8_jec    	    .push_back(corr);
      nBranches_->jetAK8_jecUp    	    .push_back(corrUp);
      nBranches_->jetAK8_jecDown    	    .push_back(corrDown);
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
      nBranches_->jetAK8_csv                .push_back(fj.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
      nBranches_->jetAK8_tau1               .push_back(fj.userFloat("NjettinessAK8:tau1"));	  
      nBranches_->jetAK8_tau2               .push_back(fj.userFloat("NjettinessAK8:tau2"));
      nBranches_->jetAK8_tau3               .push_back(fj.userFloat("NjettinessAK8:tau3")); 
      nBranches_->jetAK8_pruned_mass        .push_back(fj.userFloat("ak8PFJetsCHSPrunedMass"));
      nBranches_->jetAK8_softdrop_mass      .push_back(fj.userFloat("ak8PFJetsCHSSoftDropMass"));
      
      if(isMC){
      
        int genP_pdgId_fj = fj.genParton() ? fj.genParton()->pdgId() : -99; // default to -99 when no genParton is found!!!
        nBranches_->jetAK8_partonFlavour  .push_back(fj.partonFlavour()); //Algorithmic definition! g,u,d,s can be identified by j.partonFlavour() 
        nBranches_->jetAK8_hadronFlavour  .push_back(fj.hadronFlavour()); //Hadron flavour. b,c can be identified by either j.hadronFlavour() (hadron definition) or j.partonFlavour() (parton definition)
        nBranches_->jetAK8_genParton_pdgID.push_back(genP_pdgId_fj); //Physics definition! see https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#MC_Truth
        nBranches_->jetAK8_nbHadrons      .push_back((fj.jetFlavourInfo().getbHadrons()).size());
        nBranches_->jetAK8_ncHadrons      .push_back((fj.jetFlavourInfo().getcHadrons()).size());

	JME::JetParameters parameters_ak8;
	parameters_ak8.setJetPt(corr*uncorrJet.pt());
	parameters_ak8.setJetEta(fj.eta());
	parameters_ak8.setRho(nBranches_->rho);
      
	float jer_res= resolution_ak8.getResolution(parameters_ak8);
	float jer_sf = resolution_ak8_sf.getScaleFactor(parameters_ak8);
	float jer_sf_up = resolution_ak8_sf.getScaleFactor(parameters_ak8, Variation::UP);
	float jer_sf_down = resolution_ak8_sf.getScaleFactor(parameters_ak8, Variation::DOWN);
	nBranches_->jetAK8_jer_sf.push_back(jer_sf);
	nBranches_->jetAK8_jer_sf_up.push_back(jer_sf_up);
	nBranches_->jetAK8_jer_sf_down.push_back(jer_sf_down);
	nBranches_->jetAK8_jer_sigma_pt.push_back(jer_res);
	
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

	if(prunedjet.pt()>0.1) {

        //Compute JEC for pruned mass
        reco::Candidate::LorentzVector uncorrPrunedJet;
        double prunedcorr = 1;
	double prunedcorrUp = 1;
	double prunedcorrDown = 1;
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

      	jecAK8Unc_->setJetEta( prunedjet.correctedP4(0).eta() );
      	jecAK8Unc_->setJetPt( prunedcorr * prunedjet.correctedP4(0).pt() );
      	prunedcorrUp = prunedcorr * (1 + fabs(jecAK8Unc_->getUncertainty(1)));
      	jecAK8Unc_->setJetEta( prunedjet.correctedP4(0).eta() );
      	jecAK8Unc_->setJetPt( prunedcorr * prunedjet.correctedP4(0).pt() );
      	prunedcorrDown = prunedcorr * ( 1 - fabs(jecAK8Unc_->getUncertainty(-1)) );
      
        nBranches_->jetAK8_pruned_jecUp.push_back(prunedcorrUp);
	nBranches_->jetAK8_pruned_jecDown.push_back(prunedcorrDown);
	
        /****************************************************************/

        vPrunedSubjetpt     .clear();
        vPrunedSubjeteta    .clear();
        vPrunedSubjetmass   .clear();
        vPrunedSubjetphi    .clear();
        vPrunedSubjete      .clear();
        vPrunedSubjetcharge .clear();
        vPrunedSubjetGenPartonPdgId.clear();
        vPrunedSubjetnbHadrons.clear();
        vPrunedSubjetncHadrons.clear();
        vPrunedSubjetPartonFlavour.clear();
        vPrunedSubjetHadronFlavour.clear();
        vPrunedSubjetcsv    .clear();

        //Get subjets from matched pruned jet
        std::vector<reco::CandidatePtr> pruneddaus(prunedjet.daughterPtrVector());

        nsubjets = 0;

        for (unsigned int k = 0; k < pruneddaus.size(); ++k) {

          // const reco::PFJet &prunedsubjet = dynamic_cast<const reco::PFJet &>(*pruneddaus[k]);
        
          const pat::Jet &prunedsubjet = dynamic_cast<const pat::Jet &>(*pruneddaus[k]);   

          if( prunedsubjet.pt() < 0.01 ) continue;

          nsubjets++;

          vPrunedSubjetpt.push_back(prunedsubjet.pt());
          vPrunedSubjeteta.push_back(prunedsubjet.eta());
          vPrunedSubjetmass.push_back(prunedsubjet.mass());
          vPrunedSubjetphi.push_back(prunedsubjet.phi());
          vPrunedSubjete.push_back(prunedsubjet.energy());
          
          if(isMC){
            int genP_pdgId_prunedsubjet = prunedsubjet.genParton() ? prunedsubjet.genParton()->pdgId() : -99; // default to -99 when no genParton is found!!!
            vPrunedSubjetGenPartonPdgId.push_back(genP_pdgId_prunedsubjet);
            vPrunedSubjetnbHadrons.push_back(prunedsubjet.jetFlavourInfo().getbHadrons().size());
            vPrunedSubjetncHadrons.push_back(prunedsubjet.jetFlavourInfo().getcHadrons().size());
            vPrunedSubjetPartonFlavour.push_back(abs(prunedsubjet.partonFlavour()));
            vPrunedSubjetHadronFlavour.push_back(abs(prunedsubjet.hadronFlavour()));
          }
          
          vPrunedSubjetcharge.push_back(prunedsubjet.charge());
          vPrunedSubjetcsv.push_back(prunedsubjet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );  
        
        }

        nBranches_->jetAK8_subjet_pruned_N     .push_back(nsubjets);
        nBranches_->jetAK8_subjet_pruned_pt    .push_back(vPrunedSubjetpt);
        nBranches_->jetAK8_subjet_pruned_eta   .push_back(vPrunedSubjeteta);
        nBranches_->jetAK8_subjet_pruned_mass  .push_back(vPrunedSubjetmass);
        nBranches_->jetAK8_subjet_pruned_phi   .push_back(vPrunedSubjetphi);
        nBranches_->jetAK8_subjet_pruned_e     .push_back(vPrunedSubjete);
        nBranches_->jetAK8_subjet_pruned_charge.push_back(vPrunedSubjetcharge);
        nBranches_->jetAK8_subjet_pruned_csv   .push_back(vPrunedSubjetcsv);

        if(isMC){
	
          nBranches_->jetAK8_subjet_pruned_genParton_pdgID.push_back(vPrunedSubjetGenPartonPdgId);
          nBranches_->jetAK8_subjet_pruned_nbHadrons	 .push_back(vPrunedSubjetnbHadrons); 
          nBranches_->jetAK8_subjet_pruned_ncHadrons	 .push_back(vPrunedSubjetncHadrons);
          nBranches_->jetAK8_subjet_pruned_partonFlavour.push_back(vPrunedSubjetPartonFlavour);
          nBranches_->jetAK8_subjet_pruned_hadronFlavour.push_back(vPrunedSubjetHadronFlavour); 
	           
        }
	
	}
	
      }//doPruning
      else{

        double prunedcorr = 1;
	double prunedcorrUp = 1;
	double prunedcorrDown = 1;
        if( doCorrOnTheFly_ ){
      
           reco::Candidate::LorentzVector uncorrPrunedJet = fj.correctedP4(0);

           jecAK8Groomed_->setJetEta( uncorrPrunedJet.eta()    );
           jecAK8Groomed_->setJetPt ( uncorrPrunedJet.pt()     );
           jecAK8Groomed_->setJetE  ( uncorrPrunedJet.energy() );
           jecAK8Groomed_->setJetA  ( fj.jetArea()             );
           jecAK8Groomed_->setRho   ( nBranches_->rho          );
           jecAK8Groomed_->setNPV   ( vertices_->size()        );
           prunedcorr = jecAK8Groomed_->getCorrection();
           nBranches_->jetAK8_pruned_massCorr.push_back(prunedcorr*fj.userFloat("ak8PFJetsCHSPrunedMass"));
           nBranches_->jetAK8_pruned_jec.push_back(prunedcorr);
	          
        }
	else{
	      
           nBranches_->jetAK8_pruned_massCorr.push_back(-99);
           nBranches_->jetAK8_pruned_jec.push_back(-99);
	   
	}

        jecAK8Unc_->setJetEta( fj.correctedP4(0).eta() );
        jecAK8Unc_->setJetPt( prunedcorr * fj.correctedP4(0).pt() );
        prunedcorrUp = prunedcorr * (1 + fabs(jecAK8Unc_->getUncertainty(1)));
        jecAK8Unc_->setJetEta( fj.correctedP4(0).eta() );
        jecAK8Unc_->setJetPt( prunedcorr * fj.correctedP4(0).pt() );
        prunedcorrDown = prunedcorr * ( 1 - fabs(jecAK8Unc_->getUncertainty(-1)) ); 
      
        nBranches_->jetAK8_pruned_jecUp.push_back(prunedcorrUp);
	nBranches_->jetAK8_pruned_jecDown.push_back(prunedcorrDown);
			
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

	if(softdropjet.pt()>0.1) {

        reco::Candidate::LorentzVector uncorrSoftDropJet;
        double softdropcorr = 1;
	double softdropcorrUp = 1;
	double softdropcorrDown = 1;
      
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
           nBranches_->jetAK8_softdrop_jec.push_back(fj.userFloat("ak8PFJetsCHSSoftDropMassCorrected")/fj.userFloat("ak8PFJetsCHSSoftDropMass"));
        }

      	jecAK8Unc_->setJetEta( softdropjet.correctedP4(0).eta() );
      	jecAK8Unc_->setJetPt( softdropcorr * softdropjet.correctedP4(0).pt() );
      	softdropcorrUp = softdropcorr * (1 + fabs(jecAK8Unc_->getUncertainty(1)));
      	jecAK8Unc_->setJetEta( softdropjet.correctedP4(0).eta() );
      	jecAK8Unc_->setJetPt( softdropcorr * softdropjet.correctedP4(0).pt() );
      	softdropcorrDown = softdropcorr * ( 1 - fabs(jecAK8Unc_->getUncertainty(-1)) );
      
        nBranches_->jetAK8_softdrop_jecUp.push_back(softdropcorrUp);
	nBranches_->jetAK8_softdrop_jecDown.push_back(softdropcorrDown);
      
      	jecAK8Unc_->setJetEta( softdropjet.correctedP4(0).eta() );
      	jecAK8Unc_->setJetPt( softdropcorr * softdropjet.correctedP4(0).pt() );
      	softdropcorrUp = softdropcorr * (1 + fabs(jecAK8Unc_->getUncertainty(1)));
      	jecAK8Unc_->setJetEta( softdropjet.correctedP4(0).eta() );
      	jecAK8Unc_->setJetPt( softdropcorr * softdropjet.correctedP4(0).pt() );
      	softdropcorrDown = softdropcorr * ( 1 - fabs(jecAK8Unc_->getUncertainty(-1)) );
      
        nBranches_->jetAK8_softdrop_jecUp.push_back(softdropcorrUp);
	nBranches_->jetAK8_softdrop_jecDown.push_back(softdropcorrDown);
	
        vSoftDropSubjetpt.clear();
        vSoftDropSubjeteta.clear();
        vSoftDropSubjetmass.clear();
        vSoftDropSubjetphi.clear();
        vSoftDropSubjete.clear();
        vSoftDropSubjetcharge.clear();
        vSoftDropSubjetGenPartonPdgId.clear();
        vSoftDropSubjetnbHadrons.clear();
        vSoftDropSubjetncHadrons.clear();
        vSoftDropSubjetPartonFlavour.clear();
        vSoftDropSubjetHadronFlavour.clear();
        vSoftDropSubjetcsv.clear();

        //Get subjets from matched pruned jet
        std::vector<reco::CandidatePtr> softdropdaus(softdropjet.daughterPtrVector());

        nsubjets = 0;

        for (unsigned int k = 0; k < softdropdaus.size(); ++k) {

          // const reco::PFJet &softdropsubjet = dynamic_cast<const reco::PFJet &>(*softdropdaus[k]);
        
          const pat::Jet &softdropsubjet = dynamic_cast<const pat::Jet &>(*softdropdaus[k]);   

          if( softdropsubjet.pt() < 0.01 ) continue;

          nsubjets++;

          vSoftDropSubjetpt.push_back(softdropsubjet.pt());
          vSoftDropSubjeteta.push_back(softdropsubjet.eta());
          vSoftDropSubjetmass.push_back(softdropsubjet.mass());
          vSoftDropSubjetphi.push_back(softdropsubjet.phi());
          vSoftDropSubjete.push_back(softdropsubjet.energy());
          
          if(isMC){
            int genP_pdgId_softdropsubjet = softdropsubjet.genParton() ? softdropsubjet.genParton()->pdgId() : -99; // default to -99 when no genParton is found!!!
            vSoftDropSubjetGenPartonPdgId.push_back(genP_pdgId_softdropsubjet);
            vSoftDropSubjetnbHadrons.push_back(softdropsubjet.jetFlavourInfo().getbHadrons().size());
            vSoftDropSubjetncHadrons.push_back(softdropsubjet.jetFlavourInfo().getcHadrons().size());
            vSoftDropSubjetPartonFlavour.push_back(abs(softdropsubjet.partonFlavour()));
            vSoftDropSubjetHadronFlavour.push_back(abs(softdropsubjet.hadronFlavour()));
          }
          
          
          vSoftDropSubjetPartonFlavour.push_back(abs(softdropsubjet.partonFlavour()));
          vSoftDropSubjetHadronFlavour.push_back(abs(softdropsubjet.hadronFlavour()));
          vSoftDropSubjetcharge.push_back(softdropsubjet.charge());
          vSoftDropSubjetcsv.push_back(softdropsubjet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
	  
        }

        nBranches_->jetAK8_subjet_softdrop_N.push_back(nsubjets);
        nBranches_->jetAK8_subjet_softdrop_pt.push_back(vSoftDropSubjetpt);
        nBranches_->jetAK8_subjet_softdrop_eta.push_back(vSoftDropSubjeteta);
        nBranches_->jetAK8_subjet_softdrop_mass.push_back(vSoftDropSubjetmass);
        nBranches_->jetAK8_subjet_softdrop_phi.push_back(vSoftDropSubjetphi);
        nBranches_->jetAK8_subjet_softdrop_e.push_back(vSoftDropSubjete);
        nBranches_->jetAK8_subjet_softdrop_charge.push_back(vSoftDropSubjetcharge);
        nBranches_->jetAK8_subjet_softdrop_csv.push_back(vSoftDropSubjetcsv);

        if(isMC){
          nBranches_->jetAK8_subjet_softdrop_genParton_pdgID.push_back(vSoftDropSubjetGenPartonPdgId);
          nBranches_->jetAK8_subjet_softdrop_nbHadrons	 .push_back(vSoftDropSubjetnbHadrons); 
          nBranches_->jetAK8_subjet_softdrop_ncHadrons	 .push_back(vSoftDropSubjetncHadrons);
          nBranches_->jetAK8_subjet_softdrop_partonFlavour.push_back(vSoftDropSubjetPartonFlavour);
          nBranches_->jetAK8_subjet_softdrop_hadronFlavour.push_back(vSoftDropSubjetHadronFlavour);          
        }
	
	}
      
      }
      else{
        double softdropcorr = 1;
        if( doCorrOnTheFly_ ){
      
           reco::Candidate::LorentzVector uncorrsoftdropJet = fj.correctedP4(0);

           jecAK8Groomed_->setJetEta( uncorrsoftdropJet.eta()    );
           jecAK8Groomed_->setJetPt ( uncorrsoftdropJet.pt()     );
           jecAK8Groomed_->setJetE  ( uncorrsoftdropJet.energy() );
           jecAK8Groomed_->setJetA  ( fj.jetArea()             );
           jecAK8Groomed_->setRho   ( nBranches_->rho          );
           jecAK8Groomed_->setNPV   ( vertices_->size()        );
           softdropcorr = jecAK8Groomed_->getCorrection();
           nBranches_->jetAK8_softdrop_massCorr.push_back(softdropcorr*fj.userFloat("ak8PFJetsCHSSoftDropMass"));
           nBranches_->jetAK8_softdrop_jec.push_back(softdropcorr);
	          
        }
	else{
	      
           nBranches_->jetAK8_softdrop_massCorr.push_back(-99);
           nBranches_->jetAK8_softdrop_jec.push_back(-99);
	   
	}

        nBranches_->jetAK8_softdrop_jecUp.push_back(corrUp);
	nBranches_->jetAK8_softdrop_jecDown.push_back(corrDown);

        vSoftDropSubjetpt.clear();
        vSoftDropSubjeteta.clear();
        vSoftDropSubjetmass.clear();
        vSoftDropSubjetphi.clear();
        vSoftDropSubjete.clear();
        vSoftDropSubjetcharge.clear();
        vSoftDropSubjetPartonFlavour.clear();
        vSoftDropSubjetHadronFlavour.clear();
        vSoftDropSubjetcsv.clear();
      
        nsubjets = 0;
      
        const std::vector<edm::Ptr<pat::Jet> > &wSubjets = fj.subjets("SoftDrop");
    
      	for ( const pat::Jet & softdropsubjet : wSubjets ) {
	
           if( softdropsubjet.pt() < 0.01 ) continue;
         
           nsubjets++;

           vSoftDropSubjetpt.push_back(softdropsubjet.pt());
           vSoftDropSubjeteta.push_back(softdropsubjet.eta());
           vSoftDropSubjetmass.push_back(softdropsubjet.mass());
           vSoftDropSubjetphi.push_back(softdropsubjet.phi());
           vSoftDropSubjete.push_back(softdropsubjet.energy());
           vSoftDropSubjetPartonFlavour.push_back(abs(softdropsubjet.partonFlavour()));
           vSoftDropSubjetHadronFlavour.push_back(abs(softdropsubjet.hadronFlavour()));
           vSoftDropSubjetcharge.push_back(softdropsubjet.charge());
           vSoftDropSubjetcsv.push_back(softdropsubjet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
         
         } 

        nBranches_->jetAK8_subjet_softdrop_N.push_back(nsubjets);
        nBranches_->jetAK8_subjet_softdrop_pt.push_back(vSoftDropSubjetpt);
        nBranches_->jetAK8_subjet_softdrop_eta.push_back(vSoftDropSubjeteta);
        nBranches_->jetAK8_subjet_softdrop_mass.push_back(vSoftDropSubjetmass);
        nBranches_->jetAK8_subjet_softdrop_phi.push_back(vSoftDropSubjetphi);
        nBranches_->jetAK8_subjet_softdrop_e.push_back(vSoftDropSubjete);
        nBranches_->jetAK8_subjet_softdrop_charge.push_back(vSoftDropSubjetcharge);
        nBranches_->jetAK8_subjet_softdrop_csv.push_back(vSoftDropSubjetcsv);

        if(isMC){
	
          nBranches_->jetAK8_subjet_softdrop_partonFlavour.push_back(vSoftDropSubjetPartonFlavour);
          nBranches_->jetAK8_subjet_softdrop_hadronFlavour.push_back(vSoftDropSubjetHadronFlavour);
	  
        }
			
      }

      /****************************************************************/
      if( doPuppi ){
	JME::JetResolution resolution_ak8Puppi = JME::JetResolution(jerAK8PuppiName_res_);
	JME::JetResolutionScaleFactor resolution_ak8Puppi_sf = JME::JetResolutionScaleFactor(jerAK8PuppiName_sf_);

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
	
	if(puppijet.pt()>0.1) {
	
	  nBranches_->jetAK8_puppi_tau1	 .push_back(puppijet.userFloat("NjettinessAK8Puppi:tau1"));	 
	  nBranches_->jetAK8_puppi_tau2	 .push_back(puppijet.userFloat("NjettinessAK8Puppi:tau2"));
	  nBranches_->jetAK8_puppi_tau3	 .push_back(puppijet.userFloat("NjettinessAK8Puppi:tau3")); 
	  //nBranches_->jetAK8_puppi_pruned_mass.push_back(puppijet.userFloat("ak8PFJetsPuppiPrunedMass"));
	  nBranches_->jetAK8_puppi_softdrop_mass.push_back(puppijet.userFloat("ak8PFJetsPuppiSoftDropMass"));

	  //Compute JEC for pruned mass
	  reco::Candidate::LorentzVector uncorrJet=puppijet.p4();
	  double corr = 1;
      
	  if( doCorrOnTheFly_ ){
	    if(puppijet.jecSetsAvailable())
	      uncorrJet = puppijet.correctedP4(0);

	    jecAK8Puppi_->setJetEta( uncorrJet.eta()    );
	    jecAK8Puppi_->setJetPt ( uncorrJet.pt()     );
	    jecAK8Puppi_->setJetE  ( uncorrJet.energy() );
	    jecAK8Puppi_->setJetA  ( puppijet.jetArea()      );
	    jecAK8Puppi_->setRho   ( nBranches_->rho          );
	    jecAK8Puppi_->setNPV   ( vertices_->size()        );
	    corr = jecAK8Puppi_->getCorrection();
	    //nBranches_->jetAK8_puppi_pruned_massCorr.push_back(corr*puppijet.userFloat("ak8PFJetsPuppiPrunedMass"));
	    //nBranches_->jetAK8_puppi_pruned_jec.push_back(corr);
	    nBranches_->jetAK8_puppi_softdrop_massCorr.push_back(corr*puppijet.userFloat("ak8PFJetsPuppiSoftDropMass"));
	    nBranches_->jetAK8_puppi_softdrop_jec.push_back(corr);
	  }
	  else{
	    //nBranches_->jetAK8_puppi_pruned_massCorr.push_back(puppijet.userFloat("ak8PFJetsPuppiPrunedMassCorrected"));
	    //nBranches_->jetAK8_puppi_pruned_jec.push_back(puppijet.userFloat("ak8PFJetsPuppiPrunedMassCorrected")/puppijet.userFloat("ak8PFJetsPuppiPrunedMass"));
	    nBranches_->jetAK8_puppi_softdrop_massCorr.push_back(puppijet.userFloat("ak8PFJetsPuppiSoftDropMassCorrected"));
	    nBranches_->jetAK8_puppi_softdrop_jec.push_back(puppijet.userFloat("ak8PFJetsPuppiSoftDropMassCorrected")/puppijet.userFloat("ak8PFJetsPuppiSoftDropMass"));
	  }

	  nBranches_->jetAK8_puppi_pt     	    .push_back(corr*uncorrJet.pt());                   
	  nBranches_->jetAK8_puppi_eta    	    .push_back(puppijet.eta());
	  nBranches_->jetAK8_puppi_mass   	    .push_back(corr*uncorrJet.mass());
	  nBranches_->jetAK8_puppi_phi    	    .push_back(puppijet.phi());
	  nBranches_->jetAK8_puppi_e      	    .push_back(corr*uncorrJet.energy());

	  if (isMC){
	    JME::JetParameters parameters_ak8Puppi;
	    parameters_ak8Puppi.setJetPt(corr*uncorrJet.pt());
	    parameters_ak8Puppi.setJetEta(puppijet.eta());
	    parameters_ak8Puppi.setRho(nBranches_->rho);
	  
	    float jer_res = resolution_ak8Puppi.getResolution(parameters_ak8Puppi);
	    float jer_sf = resolution_ak8Puppi_sf.getScaleFactor(parameters_ak8Puppi);
	    float jer_sf_up = resolution_ak8Puppi_sf.getScaleFactor(parameters_ak8Puppi, Variation::UP);
	    float jer_sf_down = resolution_ak8Puppi_sf.getScaleFactor(parameters_ak8Puppi, Variation::DOWN);
	    nBranches_->jetAK8Puppi_jer_sf.push_back(jer_sf);
	    nBranches_->jetAK8Puppi_jer_sf_up.push_back(jer_sf_up);
	    nBranches_->jetAK8Puppi_jer_sf_down.push_back(jer_sf_down);
	    nBranches_->jetAK8Puppi_jer_sigma_pt.push_back(jer_res);
	  }
        } else {
        nBranches_->jetAK8_puppi_tau1	 .push_back(-99);	 
        nBranches_->jetAK8_puppi_tau2	 .push_back(-99);
        nBranches_->jetAK8_puppi_tau3	 .push_back(-99); 
        //nBranches_->jetAK8_puppi_pruned_mass.push_back(-99);
        nBranches_->jetAK8_puppi_softdrop_mass.push_back(-99);
        //nBranches_->jetAK8_puppi_pruned_massCorr.push_back(-99);
        //nBranches_->jetAK8_puppi_pruned_jec.push_back(-99);
        nBranches_->jetAK8_puppi_softdrop_massCorr.push_back(-99);
        nBranches_->jetAK8_puppi_softdrop_jec.push_back(-99);
        nBranches_->jetAK8_puppi_pt     	    .push_back(-99);                   
        nBranches_->jetAK8_puppi_eta    	    .push_back(-99);
        nBranches_->jetAK8_puppi_mass   	    .push_back(-99);
        nBranches_->jetAK8_puppi_phi    	    .push_back(-99);
        nBranches_->jetAK8_puppi_e      	    .push_back(-99);
	}

        vPuppiSoftDropSubjetpt.clear();
        vPuppiSoftDropSubjeteta.clear();
        vPuppiSoftDropSubjetmass.clear();
        vPuppiSoftDropSubjetphi.clear();
        vPuppiSoftDropSubjete.clear();
        vPuppiSoftDropSubjetcharge.clear();
        vPuppiSoftDropSubjetPartonFlavour.clear();
        vPuppiSoftDropSubjetHadronFlavour.clear();
        vPuppiSoftDropSubjetcsv.clear();
      
        nsubjets = 0;
      
        const std::vector<edm::Ptr<pat::Jet> > &wSubjets = puppijet.subjets();
    
      	for ( const pat::Jet & puppi_softdropsubjet : wSubjets ) {
	
           if( puppi_softdropsubjet.pt() < 0.01 ) continue;
         
           nsubjets++;

           vPuppiSoftDropSubjetpt.push_back(puppi_softdropsubjet.pt());
           vPuppiSoftDropSubjeteta.push_back(puppi_softdropsubjet.eta());
           vPuppiSoftDropSubjetmass.push_back(puppi_softdropsubjet.mass());
           vPuppiSoftDropSubjetphi.push_back(puppi_softdropsubjet.phi());
           vPuppiSoftDropSubjete.push_back(puppi_softdropsubjet.energy());
           vPuppiSoftDropSubjetPartonFlavour.push_back(abs(puppi_softdropsubjet.partonFlavour()));
           vPuppiSoftDropSubjetHadronFlavour.push_back(abs(puppi_softdropsubjet.hadronFlavour()));
           vPuppiSoftDropSubjetcharge.push_back(puppi_softdropsubjet.charge());
           vPuppiSoftDropSubjetcsv.push_back(puppi_softdropsubjet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
         
         } 

        nBranches_->jetAK8_subjet_puppi_softdrop_N.push_back(nsubjets);
        nBranches_->jetAK8_subjet_puppi_softdrop_pt.push_back(vPuppiSoftDropSubjetpt);
        nBranches_->jetAK8_subjet_puppi_softdrop_eta.push_back(vPuppiSoftDropSubjeteta);
        nBranches_->jetAK8_subjet_puppi_softdrop_mass.push_back(vPuppiSoftDropSubjetmass);
        nBranches_->jetAK8_subjet_puppi_softdrop_phi.push_back(vPuppiSoftDropSubjetphi);
        nBranches_->jetAK8_subjet_puppi_softdrop_e.push_back(vPuppiSoftDropSubjete);
        nBranches_->jetAK8_subjet_puppi_softdrop_charge.push_back(vPuppiSoftDropSubjetcharge);
        nBranches_->jetAK8_subjet_puppi_softdrop_csv.push_back(vPuppiSoftDropSubjetcsv);

        if(isMC){
	
          nBranches_->jetAK8_subjet_puppi_softdrop_partonFlavour.push_back(vPuppiSoftDropSubjetPartonFlavour);
          nBranches_->jetAK8_subjet_puppi_softdrop_hadronFlavour.push_back(vPuppiSoftDropSubjetHadronFlavour);
	 

        }


      } //doPuppi
      else{
	JME::JetResolution resolution_ak8Puppi = JME::JetResolution(jerAK8PuppiName_res_);
	JME::JetResolutionScaleFactor resolution_ak8Puppi_sf = JME::JetResolutionScaleFactor(jerAK8PuppiName_sf_);



        double puppi_pt  = fj.userFloat("ak8PFJetsPuppiValueMap:pt");
        double puppi_mass  = fj.userFloat("ak8PFJetsPuppiValueMap:mass");
        double puppi_eta  = fj.userFloat("ak8PFJetsPuppiValueMap:eta");
        double puppi_phi  = fj.userFloat("ak8PFJetsPuppiValueMap:phi");
        nBranches_->jetAK8_puppi_pt     	    .push_back(puppi_pt);                   
        nBranches_->jetAK8_puppi_eta    	    .push_back(puppi_eta);
        nBranches_->jetAK8_puppi_mass   	    .push_back(puppi_mass);
        nBranches_->jetAK8_puppi_phi    	    .push_back(puppi_phi);
        TLorentzVector puppi;
        puppi.SetPtEtaPhiM(puppi_pt,puppi_eta,puppi_phi,puppi_mass);
        nBranches_->jetAK8_puppi_e      	    .push_back(puppi.E());
        nBranches_->jetAK8_puppi_tau1	 .push_back(fj.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1"));	 
        nBranches_->jetAK8_puppi_tau2	 .push_back(fj.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2"));
        nBranches_->jetAK8_puppi_tau3	 .push_back(fj.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3")); 

	if(isMC){
	  JME::JetParameters parameters_ak8Puppi;
	  parameters_ak8Puppi.setJetPt(puppi_pt);
	  parameters_ak8Puppi.setJetEta(puppi_eta);
	  parameters_ak8Puppi.setRho(nBranches_->rho);
	  
	  float jer_res = resolution_ak8Puppi.getResolution(parameters_ak8Puppi);
	  float jer_sf = resolution_ak8Puppi_sf.getScaleFactor(parameters_ak8Puppi);
	  float jer_sf_up = resolution_ak8Puppi_sf.getScaleFactor(parameters_ak8Puppi, Variation::UP);
	  float jer_sf_down = resolution_ak8Puppi_sf.getScaleFactor(parameters_ak8Puppi, Variation::DOWN);
	  nBranches_->jetAK8Puppi_jer_sf.push_back(jer_sf);
	  nBranches_->jetAK8Puppi_jer_sf_up.push_back(jer_sf_up);
	  nBranches_->jetAK8Puppi_jer_sf_down.push_back(jer_sf_down);
	  nBranches_->jetAK8Puppi_jer_sigma_pt.push_back(jer_res);
	}

        TLorentzVector puppi_softdrop, puppi_softdrop_subjet;
        auto const & sdSubjetsPuppi = fj.subjets("SoftDropPuppi");
        for ( auto const & it : sdSubjetsPuppi ) {
          puppi_softdrop_subjet.SetPtEtaPhiM(it->correctedP4(0).pt(),it->correctedP4(0).eta(),it->correctedP4(0).phi(),it->correctedP4(0).mass());
          puppi_softdrop+=puppi_softdrop_subjet;
        }

        nBranches_->jetAK8_puppi_softdrop_mass.push_back(puppi_softdrop.M());

        double puppi_softdropcorr = 1;
        if( doCorrOnTheFly_ ){
           // Using puppi corrections for softdrop puppi jets. Approximation!
           jecAK8Puppi_->setJetEta( puppi_softdrop.Eta()    );
           jecAK8Puppi_->setJetPt ( puppi_softdrop.Pt()     );
           jecAK8Puppi_->setJetE  ( puppi_softdrop.E() );
           jecAK8Puppi_->setJetA  ( fj.jetArea()             );
           jecAK8Puppi_->setRho   ( nBranches_->rho          );
           jecAK8Puppi_->setNPV   ( vertices_->size()        );
           puppi_softdropcorr = jecAK8Puppi_->getCorrection();
           nBranches_->jetAK8_puppi_softdrop_massCorr.push_back(puppi_softdropcorr*puppi_softdrop.M());
           nBranches_->jetAK8_puppi_softdrop_jec.push_back(puppi_softdropcorr);
	          
        }
	else{
	      
           nBranches_->jetAK8_puppi_softdrop_massCorr.push_back(-99);
           nBranches_->jetAK8_puppi_softdrop_jec.push_back(-99);
	   
	}

        vPuppiSoftDropSubjetpt.clear();
        vPuppiSoftDropSubjeteta.clear();
        vPuppiSoftDropSubjetmass.clear();
        vPuppiSoftDropSubjetphi.clear();
        vPuppiSoftDropSubjete.clear();
        vPuppiSoftDropSubjetcharge.clear();
        vPuppiSoftDropSubjetPartonFlavour.clear();
        vPuppiSoftDropSubjetHadronFlavour.clear();
        vPuppiSoftDropSubjetcsv.clear();
      
        nsubjets = 0;
      
        const std::vector<edm::Ptr<pat::Jet> > &wSubjets = fj.subjets("SoftDropPuppi");
    
      	for ( const pat::Jet & puppi_softdropsubjet : wSubjets ) {
	
           if( puppi_softdropsubjet.pt() < 0.01 ) continue;
         
           nsubjets++;

           vPuppiSoftDropSubjetpt.push_back(puppi_softdropsubjet.pt());
           vPuppiSoftDropSubjeteta.push_back(puppi_softdropsubjet.eta());
           vPuppiSoftDropSubjetmass.push_back(puppi_softdropsubjet.mass());
           vPuppiSoftDropSubjetphi.push_back(puppi_softdropsubjet.phi());
           vPuppiSoftDropSubjete.push_back(puppi_softdropsubjet.energy());
           vPuppiSoftDropSubjetPartonFlavour.push_back(abs(puppi_softdropsubjet.partonFlavour()));
           vPuppiSoftDropSubjetHadronFlavour.push_back(abs(puppi_softdropsubjet.hadronFlavour()));
           vPuppiSoftDropSubjetcharge.push_back(puppi_softdropsubjet.charge());
           vPuppiSoftDropSubjetcsv.push_back(puppi_softdropsubjet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
         
         } 

        nBranches_->jetAK8_subjet_puppi_softdrop_N.push_back(nsubjets);
        nBranches_->jetAK8_subjet_puppi_softdrop_pt.push_back(vPuppiSoftDropSubjetpt);
        nBranches_->jetAK8_subjet_puppi_softdrop_eta.push_back(vPuppiSoftDropSubjeteta);
        nBranches_->jetAK8_subjet_puppi_softdrop_mass.push_back(vPuppiSoftDropSubjetmass);
        nBranches_->jetAK8_subjet_puppi_softdrop_phi.push_back(vPuppiSoftDropSubjetphi);
        nBranches_->jetAK8_subjet_puppi_softdrop_e.push_back(vPuppiSoftDropSubjete);
        nBranches_->jetAK8_subjet_puppi_softdrop_charge.push_back(vPuppiSoftDropSubjetcharge);
        nBranches_->jetAK8_subjet_puppi_softdrop_csv.push_back(vPuppiSoftDropSubjetcsv);

        if(isMC){
	
          nBranches_->jetAK8_subjet_puppi_softdrop_partonFlavour.push_back(vPuppiSoftDropSubjetPartonFlavour);
          nBranches_->jetAK8_subjet_puppi_softdrop_hadronFlavour.push_back(vPuppiSoftDropSubjetHadronFlavour);
	 

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



