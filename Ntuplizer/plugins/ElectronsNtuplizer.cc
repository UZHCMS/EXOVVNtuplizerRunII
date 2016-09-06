 #include "../interface/ElectronsNtuplizer.h"
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "Math/VectorUtil.h"

//===================================================================================================================
ElectronsNtuplizer::ElectronsNtuplizer( edm::EDGetTokenT<edm::View<pat::Electron>  >          electronToken, 
                                        edm::EDGetTokenT<reco::VertexCollection>              verticeToken , 
					edm::EDGetTokenT<double>                              rhoToken     ,
                                        std::vector< edm::EDGetTokenT<edm::ValueMap<bool> > > eleIDtokens  ,
					edm::EDGetTokenT<pat::TauCollection>                  boostedtauToken  ,
					NtupleBranches*                                       nBranches  ,
					std::map< std::string, bool >&                        runFlags 
					)
	: CandidateNtuplizer       ( nBranches      )
	, electronToken_           ( electronToken  )
	, verticeToken_	           ( verticeToken   )	
	, rhoToken_	           ( rhoToken	    )
	, electronVetoIdMapToken_  ( eleIDtokens[0] )
	, electronLooseIdMapToken_ ( eleIDtokens[1] )
	, electronMediumIdMapToken_( eleIDtokens[2] )
	, electronTightIdMapToken_ ( eleIDtokens[3] )
	, electronHEEPIdMapToken_  ( eleIDtokens[4] )
        , electronHEEPId51MapToken_( eleIDtokens[5] )
	, boostedtauToken_		   ( boostedtauToken    )
	, doBoostedTaus_   	   ( runFlags["doBoostedTaus"]  )
{

}

//===================================================================================================================
ElectronsNtuplizer::~ElectronsNtuplizer( void )
{

}

//===================================================================================================================
float ElectronsNtuplizer::dEtaInSeed( const pat::Electron &ele ){

	return ele.superCluster().isNonnull() && ele.superCluster()->seed().isNonnull() ? ele.deltaEtaSuperClusterTrackAtVtx() - ele.superCluster()->eta() + ele.superCluster()->seed()->eta() : std::numeric_limits<float>::max();
}


//===================================================================================================================

// YT added : 17 Aug 2016. This is for Non-Triggering MVA, used for tautau analysis, in general
bool isNonTrigElectronID(pat::Electron ele)
{
  Float_t eta = fabs(ele.superCluster()->eta());
  Float_t pt = ele.pt();

  Float_t mva = 0;
  if(ele.hasUserFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values")){
    mva = ele.userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values");
  }else{
    std::cout << "Not available nonTrig" << std::endl;
    return 0;
  }

  //  std::cout << "pT, eta, mva = " << pt << " " << eta << " " << mva << std::endl;

  if(pt > 10.){
    if(eta < 0.8) return mva > 0.967083;
    else if(eta < 1.479) return mva > 0.929117;
    else{
      return mva > 0.726311;
    }
  }else if(pt <= 10){
    if(eta < 0.8) return mva > 0.287435;
    else if(eta < 1.479) return mva > 0.221846;
    else{
      return mva > -0.303263;
    }
  }else{
    std::cout << "Not happens" << std::endl;
    return 0;
  }
  
}


float isNonTrigElectron(pat::Electron ele)
{
  
  float mva = 0;
  if(ele.hasUserFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values")){
    mva = ele.userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values");
    return mva;
  }else{
    std::cout << "Not available nonTrig" << std::endl;
    return -99;
  }

  return -99;
  
}



float ElectronCorrPFIso(pat::Electron ele, double Aeff03, float rho, edm::Handle<pat::TauCollection> 	     taus_  ){
  double TauSumChargedHadronPt = 0.;
  double TauSumNeutralHadronEt = 0.;
  double TauSumPhotonEt        = 0.;
  double dRmin = 0.4;
  pat::TauRef matchedTau;
  
  size_t numTaus = taus_->size();
  for(size_t tauIndex = 0; tauIndex < numTaus; ++tauIndex){
    pat::TauRef tau(taus_, tauIndex);
    double dR = reco::deltaR(ele.eta(), ele.phi(), tau->eta(), tau->phi());
				if ( dR < dRmin &&
				     tau->pt()>20 && 
				     fabs(tau->eta())<2.4 && 
				     tau->tauID("decayModeFindingNewDMs")>0.5  && 
				     // tau->tauID("againstMuonLoose")>0.5
				     // && 
				     // tau->tauID("againstElectronLoose")>0.5 && 
				     tau->tauID("byVLooseIsolationMVArun2v1PWnewDMwLT")>0.5
				     ) {
				  matchedTau = tau;
				  dRmin = dR;
				}
  }
  if(matchedTau.isNonnull()){
    
    for(size_t Ind1=0; Ind1<matchedTau->signalChargedHadrCands().size(); Ind1++){
      double dRConst = reco::deltaR(ele.eta(), ele.phi(), matchedTau->signalChargedHadrCands()[Ind1]->eta(), matchedTau->signalChargedHadrCands()[Ind1]->phi());
      if (dRConst <0.3)	TauSumChargedHadronPt = TauSumChargedHadronPt + matchedTau->signalChargedHadrCands()[Ind1]->pt();
    }
    for(size_t Ind2=0; Ind2<matchedTau->signalNeutrHadrCands().size(); Ind2++){
      double dRConst = reco::deltaR(ele.eta(), ele.phi(), matchedTau->signalNeutrHadrCands()[Ind2]->eta(), matchedTau->signalNeutrHadrCands()[Ind2]->phi()); 
      if (dRConst <0.3)	TauSumNeutralHadronEt = TauSumNeutralHadronEt + matchedTau->signalNeutrHadrCands()[Ind2]->pt();
    }
    for(size_t Ind3=0; Ind3<matchedTau->signalGammaCands().size(); Ind3++){
      double dRConst = reco::deltaR(ele.eta(), ele.phi(), matchedTau->signalGammaCands()[Ind3]->eta(), matchedTau->signalGammaCands()[Ind3]->phi()); 
      if (dRConst <0.3)	TauSumPhotonEt = TauSumPhotonEt + matchedTau->signalGammaCands()[Ind3]->pt();
    }
  }
  float sumChargedHadronPt = std::max(0., ele.chargedHadronIso()-TauSumChargedHadronPt);
  float sumNeutralEt       = std::max(0., ele.neutralHadronIso()-TauSumNeutralHadronEt+ele.photonIso()-TauSumPhotonEt);		
  float RhoCorrectedIso03Boost = sumChargedHadronPt + std::max(sumNeutralEt - rho * Aeff03, 0.);
  return 	RhoCorrectedIso03Boost;	
}


//===================================================================================================================


void ElectronsNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
  // bool doTausBoosted_ = event.getByToken( tauInputToken_ , taus_ ); 
 
    
   event.getByToken(electronToken_ , electrons_ ); 
   event.getByToken(verticeToken_  , vertices_  );
   event.getByToken(rhoToken_	   , rho_       );
   event.getByToken(boostedtauToken_   , taus_      );
 
   event.getByToken(electronHEEPIdMapToken_  , heep_id_decisions   );
   event.getByToken(electronHEEPId51MapToken_, heep_id51_decisions );
   event.getByToken(electronVetoIdMapToken_  , veto_id_decisions   );
   event.getByToken(electronLooseIdMapToken_ , loose_id_decisions  );
   event.getByToken(electronMediumIdMapToken_, medium_id_decisions );
   event.getByToken(electronTightIdMapToken_ , tight_id_decisions  );
  
   // Find the first vertex in the collection that passes good quality criteria
   // reco::VertexCollection::const_iterator firstGoodVertex = vertices_->end();
   reco::VertexCollection::const_iterator firstVertex = vertices_->begin();
   reco::VertexCollection::const_iterator firstGoodVertex = vertices_->begin();

   int firstGoodVertexIdx = 0;
   for( reco::VertexCollection::const_iterator vtx = vertices_->begin(); vtx != vertices_->end(); ++vtx, ++firstGoodVertexIdx){
     bool isFake = (vtx->chi2()==0 && vtx->ndof()==0);
     // Check the goodness
     if( !isFake && vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0) {
       firstGoodVertex = vtx;
       break;
     }
     
   }

   int nele = 0;

   for (const pat::Electron &ele : *electrons_) {
   	                 
    nBranches_->el_pdgId	   .push_back(ele.pdgId());
    nBranches_->el_charge	   .push_back(ele.charge());
    nBranches_->el_e		   .push_back(ele.energy());
    nBranches_->el_eta 	           .push_back(ele.eta());
    nBranches_->el_superCluster_eta.push_back(ele.superCluster()->eta());
    nBranches_->el_mass	   	   .push_back(ele.mass());
    nBranches_->el_pt  	   	   .push_back(ele.pt());  
    nBranches_->el_phi 	   	   .push_back(ele.phi()); 
    
    /*======= ISO ==========*/
    double rho = *(rho_.product()); 
    
    double Aeff04 = 0.5;
    if ( fabs(ele.superCluster()->eta()) < 1.0 ) Aeff04 = 0.208;
    if ( fabs(ele.superCluster()->eta()) > 1.0 && fabs(ele.superCluster()->eta()) < 1.479 ) Aeff04 = 0.209;
    if ( fabs(ele.superCluster()->eta()) > 1.479 && fabs(ele.superCluster()->eta()) < 2.0 ) Aeff04 = 0.115;
    if ( fabs(ele.superCluster()->eta()) > 2.0 && fabs(ele.superCluster()->eta()) < 2.2 ) Aeff04 = 0.143;
    if ( fabs(ele.superCluster()->eta()) > 2.2 && fabs(ele.superCluster()->eta()) < 2.3 ) Aeff04 = 0.183;
    if ( fabs(ele.superCluster()->eta()) > 2.3 && fabs(ele.superCluster()->eta()) < 2.4 ) Aeff04 = 0.194;
    if ( fabs(ele.superCluster()->eta()) > 2.4 ) Aeff04 = 0.261;

    double Aeff03 = 0.5;
    if ( fabs(ele.superCluster()->eta()) < 1.0 ) Aeff03 = 0.13;
    if ( fabs(ele.superCluster()->eta()) > 1.0 && fabs(ele.superCluster()->eta()) < 1.479 ) Aeff03 = 0.14;
    if ( fabs(ele.superCluster()->eta()) > 1.479 && fabs(ele.superCluster()->eta()) < 2.0 ) Aeff03 = 0.07;
    if ( fabs(ele.superCluster()->eta()) > 2.0 && fabs(ele.superCluster()->eta()) < 2.2 ) Aeff03 = 0.09;
    if ( fabs(ele.superCluster()->eta()) > 2.2 && fabs(ele.superCluster()->eta()) < 2.3 ) Aeff03 = 0.11;
    if ( fabs(ele.superCluster()->eta()) > 2.3 && fabs(ele.superCluster()->eta()) < 2.4 ) Aeff03 = 0.11;
    if ( fabs(ele.superCluster()->eta()) > 2.4 ) Aeff03 = 0.14;

    float  DeltaCorrectedIso = (ele.chargedHadronIso() + std::max(0., ele.neutralHadronIso() + ele.photonIso() - 0.5*ele.puChargedHadronIso()))/ele.pt();
    float  RhoCorrectedIso04 = ele.chargedHadronIso() + std::max(ele.neutralHadronIso() + ele.photonIso() - rho*Aeff04, 0.);
    float  RhoCorrectedIso03 = ele.chargedHadronIso() + std::max(ele.neutralHadronIso() + ele.photonIso() - rho*Aeff03, 0.);

    nBranches_->el_pfDeltaCorrRelIso	    .push_back(DeltaCorrectedIso);
    nBranches_->el_pfRhoCorrRelIso04	    .push_back(RhoCorrectedIso04);
    nBranches_->el_pfRhoCorrRelIso03	    .push_back(RhoCorrectedIso03);
    nBranches_->el_pfRelIso		    .push_back((ele.chargedHadronIso() + ele.neutralHadronIso()+ ele.photonIso())/ele.pt());
    nBranches_->el_photonIso		    .push_back(ele.photonIso());
    nBranches_->el_neutralHadIso	    .push_back(ele.neutralHadronIso());
    nBranches_->el_chargedHadIso	    .push_back(ele.chargedHadronIso());
    nBranches_->el_trackIso		    .push_back(ele.trackIso());
    
    float  DeltaCorrectedIsoBoost = (ele.userIsolation(pat::PfChargedHadronIso) + std::max(0., ele.userIsolation(pat::PfNeutralHadronIso) + ele.userIsolation(pat::PfGammaIso) - 0.5*ele.userIsolation(pat::PfPUChargedHadronIso)))/ele.pt();
    float  RhoCorrectedIso04Boost = ele.userIsolation(pat::PfChargedHadronIso) + std::max(ele.userIsolation(pat::PfNeutralHadronIso) + ele.userIsolation(pat::PfGammaIso) - rho*Aeff04, 0.);
    float  RhoCorrectedIso03Boost = ele.userIsolation(pat::PfChargedHadronIso) + std::max(ele.userIsolation(pat::PfNeutralHadronIso) + ele.userIsolation(pat::PfGammaIso) - rho*Aeff03, 0.);   
    
    nBranches_->el_pfDeltaCorrRelIsoBoost.push_back(DeltaCorrectedIsoBoost);
    nBranches_->el_pfRhoCorrRelIso04Boost.push_back(RhoCorrectedIso04Boost);
    nBranches_->el_pfRhoCorrRelIso03Boost.push_back(RhoCorrectedIso03Boost);
    nBranches_->el_pfRelIsoBoost	  .push_back((ele.userIsolation(pat::PfChargedHadronIso) + ele.userIsolation(pat::PfNeutralHadronIso) + ele.userIsolation(pat::PfGammaIso))/ele.pt());   
    nBranches_->el_photonIsoBoost        .push_back(ele.userIsolation(pat::PfGammaIso));
    nBranches_->el_neutralHadIsoBoost    .push_back(ele.userIsolation(pat::PfNeutralHadronIso));
    nBranches_->el_chargedHadIsoBoost    .push_back(ele.userIsolation(pat::PfChargedHadronIso));
    nBranches_->el_SemileptonicPFIso	  .push_back(RhoCorrectedIso03);// /ele.pt()							     );
    if ( doBoostedTaus_ )  nBranches_->el_SemileptonicCorrPFIso .push_back( ElectronCorrPFIso(ele, Aeff03, rho , taus_ ));// /ele.pt()
     
    /*======= IDs ==========*/   	    
    float et = ele.energy()!=0. ? ele.et()/ele.energy()*ele.caloEnergy() : 0.;
    nBranches_->el_et                      .push_back(et);
    nBranches_->el_passConversionVeto      .push_back(ele.passConversionVeto());
    nBranches_->el_full5x5_sigmaIetaIeta   .push_back(ele.full5x5_sigmaIetaIeta());
    nBranches_->el_dEtaIn		   .push_back(ele.deltaEtaSuperClusterTrackAtVtx());
    nBranches_->el_dPhiIn		   .push_back(ele.deltaPhiSuperClusterTrackAtVtx());
    nBranches_->el_hOverE		   .push_back(ele.hcalOverEcal());
    nBranches_->el_dz  		           .push_back(ele.gsfTrack()->dz((*firstGoodVertex).position()));
    nBranches_->el_d0  	                   .push_back((-1)*ele.gsfTrack()->dxy((*firstGoodVertex).position()));
    nBranches_->el_dz_allvertices          .push_back(ele.gsfTrack()->dz((*firstVertex).position()));
    nBranches_->el_d0_allvertices          .push_back((-1)*ele.gsfTrack()->dxy((*firstVertex).position()));
    nBranches_->el_expectedMissingInnerHits.push_back(ele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS));
    
    nBranches_->el_dr03EcalRecHitSumEt.push_back(ele.dr03EcalRecHitSumEt());
    nBranches_->el_dr03HcalDepth1TowerSumEt.push_back(ele.dr03HcalDepth1TowerSumEt());
    nBranches_->el_rho.push_back(*(rho_.product()));
    nBranches_->el_ecalDriven.push_back(ele.ecalDriven());
    nBranches_->el_dEtaInSeed.push_back(dEtaInSeed( ele ));
    nBranches_->el_full5x5_e2x5Max.push_back(ele.full5x5_e2x5Max());
    nBranches_->el_full5x5_e5x5.push_back(ele.full5x5_e5x5());
    nBranches_->el_full5x5_e1x5.push_back(ele.full5x5_e1x5());
    nBranches_->el_dr03TkSumPt.push_back(ele.dr03TkSumPt());
    nBranches_->el_superCluster_e.push_back(ele.superCluster()->energy());
    nBranches_->el_hadronicOverEm.push_back(ele.hadronicOverEm());

    reco::GsfElectron::PflowIsolationVariables pfIso = ele.pfIsolationVariables();
    nBranches_->el_relIsoWithDBeta	    .push_back((pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt ))/ele.pt());
   
    if( ele.ecalEnergy() == 0 ){
    	    nBranches_->el_ooEmooP.push_back(1e30);
    }
    else if( !std::isfinite(ele.ecalEnergy())){
    	    nBranches_->el_ooEmooP.push_back(1e30);
    }
    else{
    	    nBranches_->el_ooEmooP.push_back(fabs(1.0/ele.ecalEnergy() - ele.eSuperClusterOverP()/ele.ecalEnergy() ));
    }

    const auto el = electrons_->ptrAt(nele);

    bool isVetoElectron   = (*veto_id_decisions)[el]  ;
    bool isLooseElectron  = (*loose_id_decisions)[el] ;
    bool isMediumElectron = (*medium_id_decisions)[el];
    bool isTightElectron  = (*tight_id_decisions)[el] ;
    bool isHeepElectron   = (*heep_id_decisions)[el]  ;
    bool isHeep51Electron = (*heep_id51_decisions)[el];
  
    nBranches_->el_isVetoElectron  .push_back(isVetoElectron);
    nBranches_->el_isLooseElectron .push_back(isLooseElectron);
    nBranches_->el_isMediumElectron.push_back(isMediumElectron);
    nBranches_->el_isTightElectron .push_back(isTightElectron);
    nBranches_->el_nonTrigMVAID	   .push_back(isNonTrigElectronID(ele)); 
    nBranches_->el_nonTrigMVA	   .push_back(isNonTrigElectron(ele)); 
    nBranches_->el_isHeepElectron  .push_back(isHeepElectron);  
    nBranches_->el_isHeep51Electron.push_back(isHeep51Electron);  

    if ( doBoostedTaus_ ) {
      bool isVetoElectronBoosted   = eleIDpassedBoosted("Veto"  ,ele) ;
      bool isLooseElectronBoosted  = eleIDpassedBoosted("Loose" ,ele) ;
      bool isMediumElectronBoosted = eleIDpassedBoosted("Medium",ele) ;
      bool isTightElectronBoosted  = eleIDpassedBoosted("Tight" ,ele) ;
      bool isHeepElectronBoosted   = eleIDpassedBoosted("Heep"  ,ele) ;
      bool isHeep51ElectronBoosted = eleIDpassedBoosted("Heep51",ele) ;
      
      
      
      nBranches_->el_isVetoElectronBoosted  .push_back(isVetoElectronBoosted);
      nBranches_->el_isLooseElectronBoosted .push_back(isLooseElectronBoosted);
      nBranches_->el_isMediumElectronBoosted.push_back(isMediumElectronBoosted);
      nBranches_->el_isTightElectronBoosted .push_back(isTightElectronBoosted);
      nBranches_->el_isHeepElectronBoosted  .push_back(isHeepElectronBoosted);  
      nBranches_->el_isHeep51ElectronBoosted.push_back(isHeep51ElectronBoosted);  
    }
    



    /***************************************************************************/

    ++nele;
        
   }

   nBranches_->el_N = nele;

}


bool ElectronsNtuplizer::eleIDpassedBoosted(std::string id, const pat::Electron &ele ){
  ///NOT REQUIRING ISOLATION FOR TAU ANALYSIS
 
  // Find the first vertex in the collection that passes good quality criteria
  // reco::VertexCollection::const_iterator firstGoodVertex = vertices_->end();
  reco::VertexCollection::const_iterator firstGoodVertex = vertices_->begin();
  int firstGoodVertexIdx = 0;
  for( reco::VertexCollection::const_iterator vtx = vertices_->begin(); vtx != vertices_->end(); ++vtx, ++firstGoodVertexIdx){
    bool isFake = (vtx->chi2()==0 && vtx->ndof()==0);
    // Check the goodness
    if( !isFake && vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0) {
      firstGoodVertex = vtx;
      break;
    }
    
  }
   
  // Ele veto, medium and tight ID starts here!!! 
  // All ID variables
  Float_t dEtaIn_;
  Float_t dPhiIn_;
  Float_t hOverE_;
  Float_t full5x5_sigmaIetaIeta_;
  // Float_t relIsoWithDBeta_;
  Float_t ooEmooP_;
  Float_t d0_;
  Float_t dz_;
  Float_t eta;
  Int_t   expectedMissingInnerHits_;
  Int_t   passConversionVeto_;
  bool    isVetoElectron   = false;  
  bool    isLooseElectron   = false;
  bool    isMediumElectron = false;
  bool    isTightElectron  = false;

  eta = ele.superCluster()->eta();  
  dEtaIn_ = ele.deltaEtaSuperClusterTrackAtVtx();
  dPhiIn_ = ele.deltaPhiSuperClusterTrackAtVtx();
  hOverE_ = ele.hcalOverEcal();
  full5x5_sigmaIetaIeta_ = ele.full5x5_sigmaIetaIeta();
  // |1/E-1/p| = |1/E - EoverPinner/E| is computed below
  // The if protects against ecalEnergy == inf or zero (always
  // the case for electrons below 5 GeV in miniAOD)
  if( ele.ecalEnergy() == 0 ){
	  ooEmooP_ = 1e30;
  }
  else if( !std::isfinite(ele.ecalEnergy())){
	  ooEmooP_ = 1e30;
  }
  else{
	  ooEmooP_ = fabs(1.0/ele.ecalEnergy() - ele.eSuperClusterOverP()/ele.ecalEnergy() );
  }
  // Isolation
  //reco::GsfElectron::PflowIsolationVariables pfIso = ele.pfIsolationVariables();
  // Compute isolation with delta beta correction for PU
  // float absiso = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
  //  relIsoWithDBeta_ = absiso/ele.pt();
  // Impact parameter
  d0_ = (-1) * ele.gsfTrack()->dxy((*firstGoodVertex).position() );
  dz_ = ele.gsfTrack()->dz( (*firstGoodVertex).position() );

  // Conversion rejection
  expectedMissingInnerHits_ = ele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
  passConversionVeto_ = ele.passConversionVeto(); 
 
  //Barrel cuts 
  if(fabs(eta) <= 1.479){
  
	  if( passConversionVeto_	           &&
	      full5x5_sigmaIetaIeta_    < 0.0114 &&
	      fabs(dEtaIn_)	        < 0.0152 &&
	      fabs(dPhiIn_)	        < 0.216 &&
	      hOverE_		        < 0.181 &&
	      // relIsoWithDBeta_          < 0.126 &&
	      ooEmooP_  	        < 0.207 &&
	      fabs(d0_) 	        < 0.0564 &&
	      fabs(dz_) 	        < 0.472 &&			      
	      expectedMissingInnerHits_ <= 2.000000 &&
	      passConversionVeto_ 
	      ) isVetoElectron = true;

	   if( passConversionVeto_	           &&
	      full5x5_sigmaIetaIeta_    < 0.0103   &&
	      fabs(dEtaIn_)	        < 0.0105 &&
	      fabs(dPhiIn_)	        < 0.115 &&
	      hOverE_		        < 0.104 &&
	      // relIsoWithDBeta_          < 0.0893 &&
	      ooEmooP_  	        < 0.102 &&
	      fabs(d0_) 	        < 0.0261 &&
	      fabs(dz_) 	        < 0.41 &&			      
	      expectedMissingInnerHits_ <= 2.000000 &&
	      passConversionVeto_ 
	      ) isLooseElectron = true;
	  
	  if( passConversionVeto_		   &&
	      full5x5_sigmaIetaIeta_    < 0.0101 &&
	      fabs(dEtaIn_)		< 0.0103 &&
	      fabs(dPhiIn_)		< 0.0336 &&
	      hOverE_			< 0.0876 &&
	      //relIsoWithDBeta_  	< 0.0766 &&
	      ooEmooP_  		< 0.0174 &&
	      fabs(d0_) 		< 0.011811 &&
	      fabs(dz_) 		< 0.0373 &&				   
	      expectedMissingInnerHits_ <= 2.000000 &&
	      passConversionVeto_ 
	      ) isMediumElectron = true;
	  
	  if( passConversionVeto_		   &&
	      full5x5_sigmaIetaIeta_    < 0.0101 &&
	      fabs(dEtaIn_)		< 0.00926 &&
	      fabs(dPhiIn_)		< 0.0336 &&
	      hOverE_			< 0.0597 &&
	      // relIsoWithDBeta_  	< 0.0354 &&
	      ooEmooP_  		< 0.012 &&
	      fabs(d0_) 		< 0.0111 &&
	      fabs(dz_) 		< 0.0466 &&				   
	      expectedMissingInnerHits_ <= 2.000000 &&
	      passConversionVeto_ 
	      ) isTightElectron = true;	  
	      
  } 
  //Endcap cut
  else if(fabs(eta) > 1.479 && fabs(eta) < 2.5){
		  
	  if( passConversionVeto_		   &&
	      full5x5_sigmaIetaIeta_	< 0.0352 &&
	      fabs(dEtaIn_)		< 0.0113 &&
	      fabs(dPhiIn_)		< 0.237 &&
	      hOverE_			< 0.1116 &&
	      //relIsoWithDBeta_  	< 0.144 &&
	      ooEmooP_  		< 0.174 &&
	      fabs(d0_) 		< 0.222 &&
	      fabs(dz_) 		< 0.921 &&				   
	      expectedMissingInnerHits_ <= 3.000000 
	      ) isVetoElectron = true;
	
	  if( passConversionVeto_		   &&
	      full5x5_sigmaIetaIeta_    < 0.0301   &&
	      fabs(dEtaIn_)		< 0.00814   &&
	      fabs(dPhiIn_)		< 0.182   &&
	      hOverE_			< 0.0897    &&
	      //relIsoWithDBeta_  	< 0.121    &&
	      ooEmooP_  		< 0.126   &&
	      fabs(d0_) 		< 0.118   &&
	      fabs(dz_) 		< 0.822   &&
	      expectedMissingInnerHits_ < 1.00000 
	      ) isLooseElectron = true;  
	
	  if( passConversionVeto_		   &&
	      full5x5_sigmaIetaIeta_	< 0.0283 &&
	      fabs(dEtaIn_)		< 0.00733 &&
	      fabs(dPhiIn_)		< 0.114 &&
	      hOverE_			< 0.0678 &&
	      //relIsoWithDBeta_  	< 0.0678 &&
	      ooEmooP_  		< 0.0898 &&
	      fabs(d0_) 		< 0.0739 &&
	      fabs(dz_) 		< 0.602 &&				   
	      expectedMissingInnerHits_ <= 1.000000
	      ) isMediumElectron = true;  

	  if( passConversionVeto_		   &&
	      full5x5_sigmaIetaIeta_	< 0.0279 &&
	      fabs(dEtaIn_)		< 0.00724 &&
	      fabs(dPhiIn_)		< 0.0918 &&
	      hOverE_			< 0.0615 &&
	      //relIsoWithDBeta_  	< 0.0646 &&
	      ooEmooP_  		< 0.00999 &&
	      fabs(d0_) 		< 0.0351 &&
	      fabs(dz_) 		< 0.417 &&				   
	      expectedMissingInnerHits_ <= 1.000000
	      ) isTightElectron = true;  
	      
  }

  bool isHeepElectron = false;
  bool isHeep51Electron = false;

  float et = ele.energy()!=0. ? ele.et()/ele.energy()*ele.caloEnergy() : 0.;
  // double iso;
  // double isoCut;
  // double rho = *(rho_.product()); 
  double dxy = ( vertices_->size() ? ele.gsfTrack()->dxy((*firstGoodVertex).position()) :  ele.gsfTrack()->dxy() );
  
  if (ele.gsfTrack().isNonnull()){
  
    if( et > 35. ) {
    
       //barrel electrons
       if( fabs(eta) < 1.4442 ){
       
    	  // iso = ele.dr03EcalRecHitSumEt() + ele.dr03HcalDepth1TowerSumEt();
    	  // isoCut = 2 + 0.03*et + 0.28*rho;	    
    	  if(  ele.ecalDriven() == 1 && dEtaInSeed( ele ) < 0.004 && 
	       ele.deltaPhiSuperClusterTrackAtVtx() < 0.06 && 
    	      (ele.full5x5_e2x5Max()/ele.full5x5_e5x5() > 0.94 || ele.full5x5_e1x5()/ele.full5x5_e5x5() > 0.83) &&
    	       ele.dr03TkSumPt() < 5. && 
	       ele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 &&
    	       // iso < isoCut &&
	       fabs(dxy) < 0.02 )
    	  {
    	     if (ele.hadronicOverEm() < (2./ele.superCluster()->energy()+0.05)) isHeep51Electron = true;
    	     if (ele.hadronicOverEm() < (1./ele.superCluster()->energy()+0.05)) isHeepElectron = true;
    	  }
	  
       }
       //endcap electrons
       if( fabs(eta) > 1.566 && fabs(eta) < 2.5 ){
       
    	  // iso = ele.dr03EcalRecHitSumEt() + ele.dr03HcalDepth1TowerSumEt();
    	  // if( et <= 50 ) isoCut = 2.5 + 0.28*rho;
    	  // else isoCut = 2.5+0.03*(et-50.) + 0.28*rho;	 
    	  if( ele.ecalDriven() == 1 && 
	      dEtaInSeed( ele ) < 0.006 && 
	      ele.deltaPhiSuperClusterTrackAtVtx() < 0.06 && 
    	      ele.full5x5_sigmaIetaIeta() < 0.03 && 
    	      ele.dr03TkSumPt() < 5. && 
	      ele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 &&
    	      // iso < isoCut &&
	      fabs(dxy) < 0.05 )
    	  {
    	     if (ele.hadronicOverEm() < (12.5/ele.superCluster()->energy()+0.05)) isHeep51Electron = true;
    	     if (ele.hadronicOverEm() < (5./ele.superCluster()->energy()+0.05)) isHeepElectron = true;
    	  }
	  
       }
       	 
    }
    	   
  }
	     
  if( id == "Veto" ) return isVetoElectron;
  else if( id == "Loose" ) return isLooseElectron;
  else if( id == "Medium" ) return isMediumElectron;
  else if( id == "Tight" ) return isTightElectron;
  else if( id == "Heep51" ) return isHeep51Electron;
  else if( id == "Heep" ) return isHeepElectron;
  else return false;
   
}

