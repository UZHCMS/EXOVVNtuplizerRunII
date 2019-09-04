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
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"


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
					edm::EDGetTokenT<edm::ValueMap<float> >               mvaValuesMapToken,
					edm::EDGetTokenT<edm::ValueMap<int> >                 mvaCategoriesMapToken,
					edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> ebRecHitsToken, 	
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
        , electronHLTIdMapToken_   ( eleIDtokens[4] )
	, electronHEEPIdMapToken_  ( eleIDtokens[5] )
	, electronMVAMediumIdMapToken_( eleIDtokens[6] )
	, electronMVATightIdMapToken_ ( eleIDtokens[7] )
	, mvaValuesMapToken_( mvaValuesMapToken )
	, mvaCategoriesMapToken_( mvaCategoriesMapToken )
	, ebRecHitsToken_ ( ebRecHitsToken )
	, isJpsiMu_( runFlags["doJpsiMu"]  )
	, isJpsiEle_( runFlags["doJpsiEle"]  )
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


float ElectronsNtuplizer::EleEInverseMinusPInverse( const pat::Electron &ele ){

  const float ecal_energy_inverse = 1.0/ele.ecalEnergy();
  const float eSCoverP = ele.eSuperClusterOverP();

  return std::abs(1.0 - eSCoverP)*ecal_energy_inverse;
}


//===================================================================================================================


void ElectronsNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
  
 //chunk to remove those events with no jspi if that analysis is chosen
  std::vector<int> doJpsi_;
  if(isJpsiEle_) {
    doJpsi_ = nBranches_->IsJpsiEle;
    //std::cout<<"im getting inside the electron part"<<std::endl;
  }else if(isJpsiMu_){
    doJpsi_ = nBranches_->IsJpsiMu;
    //std::cout<<"nbranch thing\t"<<size(isJpsi_)<<"; "<< isJpsi_[0]<<std::endl;
  }

  event.getByToken(electronToken_ , electrons_ ); 
  event.getByToken(verticeToken_  , vertices_  );
  event.getByToken(rhoToken_	   , rho_       ); 
  event.getByToken(electronVetoIdMapToken_  , veto_id_decisions   );
  event.getByToken(electronLooseIdMapToken_ , loose_id_decisions  );
  event.getByToken(electronMediumIdMapToken_, medium_id_decisions );
  event.getByToken(electronTightIdMapToken_ , tight_id_decisions  );
  event.getByToken(electronHLTIdMapToken_, hlt_id_decisions );
  event.getByToken(electronHEEPIdMapToken_  , heep_id_decisions   );
  event.getByToken(electronMVAMediumIdMapToken_, mva_medium_id_decisions );
  event.getByToken(electronMVATightIdMapToken_ , mva_tight_id_decisions  );
  event.getByToken(mvaValuesMapToken_ , mva_value);
  event.getByToken(mvaCategoriesMapToken_ , mva_categories);
  event.getByToken(ebRecHitsToken_, _ebRecHits);
  
   
  
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
    if ( fabs(ele.superCluster()->eta()) < 1.0 ) Aeff03 = 0.1566;
    if ( fabs(ele.superCluster()->eta()) > 1.0 && fabs(ele.superCluster()->eta()) < 1.479 ) Aeff03 = 0.1628;
    if ( fabs(ele.superCluster()->eta()) > 1.479 && fabs(ele.superCluster()->eta()) < 2.0 ) Aeff03 = 0.1073;
    if ( fabs(ele.superCluster()->eta()) > 2.0 && fabs(ele.superCluster()->eta()) < 2.2 ) Aeff03 = 0.0854;
    if ( fabs(ele.superCluster()->eta()) > 2.2 && fabs(ele.superCluster()->eta()) < 2.3 ) Aeff03 = 0.1051;
    if ( fabs(ele.superCluster()->eta()) > 2.3 && fabs(ele.superCluster()->eta()) < 2.4 ) Aeff03 = 0.1204;
    if ( fabs(ele.superCluster()->eta()) > 2.4 ) Aeff03 = 0.1524;

    float  DeltaCorrectedIso = (ele.chargedHadronIso() + std::max(0., ele.neutralHadronIso() + ele.photonIso() - 0.5*ele.puChargedHadronIso()))/ele.pt();
    float  RhoCorrectedIso04 = ele.chargedHadronIso() + std::max(ele.neutralHadronIso() + ele.photonIso() - rho*Aeff04, 0.);
    float  RhoCorrectedIso03 = ele.chargedHadronIso() + std::max(ele.neutralHadronIso() + ele.photonIso() - rho*Aeff03, 0.);

    nBranches_->el_pfDeltaCorrRelIso	    .push_back(DeltaCorrectedIso);
    nBranches_->el_pfRhoCorrRelIso04	    .push_back(RhoCorrectedIso04);
    nBranches_->el_pfRhoCorrRelIso03	    .push_back(RhoCorrectedIso03);
    // nBranches_->el_pfRelIso		    .push_back((ele.chargedHadronIso() + ele.neutralHadronIso()+ ele.photonIso())/ele.pt());
    nBranches_->el_photonIso		    .push_back(ele.photonIso());
    nBranches_->el_neutralHadIso	    .push_back(ele.neutralHadronIso());
    nBranches_->el_chargedHadIso	    .push_back(ele.chargedHadronIso());
    //nBranches_->el_trackIso		    .push_back(ele.trackIso());
    // float  DeltaCorrectedIsoBoost = (ele.userIsolation(pat::PfChargedHadronIso) + std::max(0., ele.userIsolation(pat::PfNeutralHadronIso) + ele.userIsolation(pat::PfGammaIso) - 0.5*ele.userIsolation(pat::PfPUChargedHadronIso)))/ele.pt();
    // float  RhoCorrectedIso04Boost = ele.userIsolation(pat::PfChargedHadronIso) + std::max(ele.userIsolation(pat::PfNeutralHadronIso) + ele.userIsolation(pat::PfGammaIso) - rho*Aeff04, 0.);
    // float  RhoCorrectedIso03Boost = ele.userIsolation(pat::PfChargedHadronIso) + std::max(ele.userIsolation(pat::PfNeutralHadronIso) + ele.userIsolation(pat::PfGammaIso) - rho*Aeff03, 0.);   
    
    // nBranches_->el_pfDeltaCorrRelIsoBoost.push_back(DeltaCorrectedIsoBoost);
    // nBranches_->el_pfRhoCorrRelIso04Boost.push_back(RhoCorrectedIso04Boost);
    // nBranches_->el_pfRhoCorrRelIso03Boost.push_back(RhoCorrectedIso03Boost);
    // nBranches_->el_pfRelIsoBoost	  .push_back((ele.userIsolation(pat::PfChargedHadronIso) + ele.userIsolation(pat::PfNeutralHadronIso) + ele.userIsolation(pat::PfGammaIso))/ele.pt());   
    // nBranches_->el_photonIsoBoost        .push_back(ele.userIsolation(pat::PfGammaIso));
    // nBranches_->el_neutralHadIsoBoost    .push_back(ele.userIsolation(pat::PfNeutralHadronIso));
    // nBranches_->el_chargedHadIsoBoost    .push_back(ele.userIsolation(pat::PfChargedHadronIso));
    nBranches_->el_SemileptonicPFIso	  .push_back(RhoCorrectedIso03);// /ele.pt()			  
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
    nBranches_->el_expectedMissingInnerHits.push_back(ele.gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS));
    //numberOfHits gives compilation error in 94X->changed in numberOfAllHits
    
    nBranches_->el_dr03EcalRecHitSumEt.push_back(ele.dr03EcalRecHitSumEt());
    nBranches_->el_dr03HcalDepth1TowerSumEt.push_back(ele.dr03HcalDepth1TowerSumEt());
    nBranches_->el_rho.push_back(*(rho_.product()));
    nBranches_->el_ecalDriven.push_back(ele.ecalDriven());
    nBranches_->el_dEtaInSeed.push_back(dEtaInSeed( ele ));
    nBranches_->el_full5x5_e2x5Max.push_back(ele.full5x5_e2x5Max());
    nBranches_->el_full5x5_e5x5.push_back(ele.full5x5_e5x5());
    nBranches_->el_full5x5_e1x5.push_back(ele.full5x5_e1x5());
    nBranches_->el_full5x5_r9.push_back(ele.full5x5_r9());
    nBranches_->el_dr03TkSumPt.push_back(ele.dr03TkSumPt());
    nBranches_->el_superCluster_e.push_back(ele.superCluster()->energy());
    nBranches_->el_hadronicOverEm.push_back(ele.hadronicOverEm());
    
    // Seed energy for slew rate corrections
    DetId detid = ele.superCluster()->seed()->seed();
    const EcalRecHit * rh = NULL;
    double seedE(0.);
    if (detid.subdetId() == EcalBarrel) {
        auto rh_i =  _ebRecHits->find(detid);
        if( rh_i != _ebRecHits->end()) rh =  &(*rh_i);
        else rh = NULL;
    }
    if(rh==NULL) seedE = -1.;
    else {seedE = rh->energy();}
    nBranches_->el_seedEnergy.push_back(seedE);
    

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
    bool isHltElectron = (*hlt_id_decisions)[el];
    bool isMVAMediumElectron = (*mva_medium_id_decisions)[el];
    bool isMVATightElectron  = (*mva_tight_id_decisions)[el] ;
  
    nBranches_->el_isVetoElectron     .push_back(isVetoElectron);
    nBranches_->el_isLooseElectron    .push_back(isLooseElectron);
    nBranches_->el_isMediumElectron   .push_back(isMediumElectron);
    nBranches_->el_isTightElectron    .push_back(isTightElectron);
    nBranches_->el_isHeepElectron     .push_back(isHeepElectron);  
    nBranches_->el_isHltElectron      .push_back(isHltElectron);  
    nBranches_->el_isMVAMediumElectron.push_back(isMVAMediumElectron);
    nBranches_->el_isMVATightElectron .push_back(isMVATightElectron);
    nBranches_->el_MVAscore           .push_back((*mva_value)[el]);
    nBranches_->el_MVAcategory        .push_back((*mva_categories)[el]);

    bool isVetoElectronWithoutIPandIsolation   = eleIDpassedWithoutIPandIsolation("Veto"  ,ele,rho) ;
    bool isLooseElectronWithoutIPandIsolation  = eleIDpassedWithoutIPandIsolation("Loose" ,ele,rho) ;
    bool isMediumElectronWithoutIPandIsolation = eleIDpassedWithoutIPandIsolation("Medium",ele,rho) ;
    bool isTightElectronWithoutIPandIsolation  = eleIDpassedWithoutIPandIsolation("Tight" ,ele,rho) ;
    bool isHeepElectronWithoutIPandIsolation   = eleIDpassedWithoutIPandIsolation("Heep"  ,ele,rho) ;
    
    nBranches_->el_isVetoElectronWithoutIPandIsolation  .push_back(isVetoElectronWithoutIPandIsolation);
    nBranches_->el_isLooseElectronWithoutIPandIsolation .push_back(isLooseElectronWithoutIPandIsolation);
    nBranches_->el_isMediumElectronWithoutIPandIsolation.push_back(isMediumElectronWithoutIPandIsolation);
    nBranches_->el_isTightElectronWithoutIPandIsolation .push_back(isTightElectronWithoutIPandIsolation);
    nBranches_->el_isHeepElectronWithoutIPandIsolation  .push_back(isHeepElectronWithoutIPandIsolation);  
    


    ++nele;
        
   }

   nBranches_->el_N = nele;

}






bool ElectronsNtuplizer::eleIDpassedWithoutIPandIsolation(std::string id, const pat::Electron &ele, float rho){

  // No requirements on isolation AND d0/dz cuts for tau analysis
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Working_points_for_2016_data_for

  bool    isVetoElectron   = false;  
  bool    isLooseElectron   = false;
  bool    isMediumElectron = false;
  bool    isTightElectron  = false;

  Float_t full5x5_sigmaIetaIeta_ = ele.full5x5_sigmaIetaIeta(); 
  Float_t eta = ele.superCluster()->eta();  
  Float_t dEtaInSeed_ = dEtaInSeed( ele );
  Float_t dPhiIn_ = std::abs(ele.deltaPhiSuperClusterTrackAtVtx());
  Float_t hOverE_ = ele.hadronicOverEm();
  Float_t E_ =ele.correctedEcalEnergy();/// to be checked
  Float_t ooEmooP_ = EleEInverseMinusPInverse(ele);
  Int_t expectedMissingInnerHits_ = ele.gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);
  Int_t passConversionVeto_ = ele.passConversionVeto(); 
 
  bool pass_C_0_requirement=false;
  Float_t C_0=0.;
  
  //Barrel cuts 
  if(fabs(eta) <= 1.479){
  
    if(
       passConversionVeto_ && 
       full5x5_sigmaIetaIeta_    < 0.0128 &&
       fabs(dEtaInSeed_)        < 0.00523 &&
       fabs(dPhiIn_)	        < 0.159 &&
       //hOverE_		        < 0.05 &&
       ooEmooP_  	        < 0.193 &&
       expectedMissingInnerHits_ <= 2.
       ){ isVetoElectron = true; C_0=0.05;}

    if(
       passConversionVeto_ && 
       full5x5_sigmaIetaIeta_    < 0.0105 &&
       fabs(dEtaInSeed_)        < 0.00387 &&
       fabs(dPhiIn_)	        < 0.0716 &&
       //hOverE_		        < 0.05 &&
       ooEmooP_  	        < 0.129 &&
       expectedMissingInnerHits_ <= 1.
       ){ isLooseElectron = true;C_0=0.05;}

    if(
       passConversionVeto_ && 
       full5x5_sigmaIetaIeta_    < 0.0105 &&
       fabs(dEtaInSeed_)        < 0.00365 &&
       fabs(dPhiIn_)	        < 0.0588 &&
       // hOverE_		        < 0.026 &&
       ooEmooP_  	        < 0.0327 &&
       expectedMissingInnerHits_ <= 1.
       ){ isMediumElectron = true;C_0=0.026;}

    if(
       passConversionVeto_ && 
       full5x5_sigmaIetaIeta_    < 0.0104 &&
       fabs(dEtaInSeed_)        < 0.00353 &&
       fabs(dPhiIn_)	        < 0.0499 &&
       //hOverE_		        < 0.026 &&
       ooEmooP_  	        < 0.0278 &&
       expectedMissingInnerHits_ <= 1.
       ){ isTightElectron = true;C_0=0.026;}
  
    if ( hOverE_< C_0 +1.12 *E_+ 0.0368*rho/E_) pass_C_0_requirement =true; 
  	      
  } 
  //Endcap cut
  else if(fabs(eta) > 1.479 && fabs(eta) < 2.5){


    if(
       passConversionVeto_ && 
       full5x5_sigmaIetaIeta_    < 0.0445 &&
       fabs(dEtaInSeed_)        < 0.00984 &&
       fabs(dPhiIn_)	        < 0.157 &&
       // hOverE_		        < 0.0962 &&
       ooEmooP_  	        < 0.0962 &&
       expectedMissingInnerHits_ <= 3.
       ){ isVetoElectron = true;C_0=0.05;}

    if(
       passConversionVeto_ && 
       full5x5_sigmaIetaIeta_    < 0.0356 &&
       fabs(dEtaInSeed_)        < 0.0072 &&
       fabs(dPhiIn_)	        < 0.147 &&
       //hOverE_		        < 0.0414 &&
       ooEmooP_  	        < 0.0875 &&
       expectedMissingInnerHits_ <= 1.
       ) {isLooseElectron = true;C_0=0.0414;}

    if(
       passConversionVeto_ && 
       full5x5_sigmaIetaIeta_    < 0.0309 &&
       fabs(dEtaInSeed_)        < 0.00625 &&
       fabs(dPhiIn_)	        < 0.0355 &&
       //hOverE_		        < 0.026 &&
       ooEmooP_  	        < 0.0335 &&
       expectedMissingInnerHits_ <= 1.
       ) {isMediumElectron = true;C_0=0.026;}

    if(
       passConversionVeto_ && 
       full5x5_sigmaIetaIeta_    < 0.0305 &&
       fabs(dEtaInSeed_)        < 0.00567 &&
       fabs(dPhiIn_)	        < 0.0165 &&
       //hOverE_		        < 0.026 &&
       ooEmooP_  	        < 0.0158 &&
       expectedMissingInnerHits_ <= 1.
       ) {isTightElectron = true;C_0=0.026;}
    
    if ( hOverE_< C_0 +1.12 *E_+ 0.201*rho/E_) pass_C_0_requirement =true; 
    
  }


  // HEEP electron IDs
  // https://twiki.cern.ch/twiki/bin/view/CMS/HEEPElectronIdentificationRun2
  bool isHeepElectron = false;

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
	       ele.gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 &&
    	       // iso < isoCut &&
	       fabs(dxy) < 0.02 )
    	  {
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
	      ele.gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 &&
    	      // iso < isoCut &&
	      fabs(dxy) < 0.05 )
    	  {
    	     if (ele.hadronicOverEm() < (5./ele.superCluster()->energy()+0.05)) isHeepElectron = true;
    	  }
	  
       }
       	 
    }
    	   
  }
	     
  if( id == "Veto" ) return (isVetoElectron && pass_C_0_requirement);
  else if( id == "Loose" ) return (isLooseElectron && pass_C_0_requirement);
  else if( id == "Medium" ) return (isMediumElectron && pass_C_0_requirement);
  else if( id == "Tight" ) return (isTightElectron && pass_C_0_requirement);
  else if( id == "Heep" ) return isHeepElectron;
  else return false;
   
}
