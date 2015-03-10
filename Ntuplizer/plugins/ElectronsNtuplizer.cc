#include "../interface/ElectronsNtuplizer.h"

// #include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"
// #include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
// #include <cmath>

//===================================================================================================================
ElectronsNtuplizer::ElectronsNtuplizer( edm::EDGetTokenT<pat::ElectronCollection> electronToken, edm::EDGetTokenT<reco::VertexCollection> verticeToken, edm::EDGetTokenT<double> rhoToken, NtupleBranches* nBranches )
   : CandidateNtuplizer( nBranches     )
   , electronToken_    ( electronToken )
   , verticeToken_     ( verticeToken  )	
   , rhoToken_	       ( rhoToken      )
{

}

//===================================================================================================================
ElectronsNtuplizer::~ElectronsNtuplizer( void )
{

}

//===================================================================================================================
float ElectronsNtuplizer::dEtaInSeed( const pat::Electron &ele ){
  return ele.superCluster().isNonnull() && ele.superCluster()->seed().isNonnull() ? 
    ele.deltaEtaSuperClusterTrackAtVtx() - ele.superCluster()->eta() + ele.superCluster()->seed()->eta() : std::numeric_limits<float>::max();
}

//===================================================================================================================
void ElectronsNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  event.getByToken(electronToken_, electrons_); 
  event.getByToken(verticeToken_ , vertices_ );
  event.getByToken(rhoToken_	 , rho_	     );

  int nele = 0;
  
  for (const pat::Electron &ele : *electrons_) {
    
    bool isHEEP = false;
    bool isHEEPv50 = false;

    float et = ele.energy()!=0. ? ele.et()/ele.energy()*ele.caloEnergy() : 0.;
    
    float eta = ele.superCluster()->eta();
    double iso;
    double isoCut;
    double rho = *(rho_.product()); 
    double dxy = ( vertices_->size() ? ele.gsfTrack()->dxy(vertices_->at(0).position()) :  ele.gsfTrack()->dxy() );
        
    if (ele.gsfTrack().isNonnull()){
        
       if( et > 35. ) {
          if( fabs(eta) < 1.4442 ){
             iso = ele.dr03EcalRecHitSumEt() + ele.dr03HcalDepth1TowerSumEt();
             isoCut = 2 + 0.03*et + 0.28*rho;	  
             if( ele.ecalDriven() == 1 && dEtaInSeed( ele ) < 0.004 && ele.deltaPhiSuperClusterTrackAtVtx() < 0.06 && 
	         ele.hadronicOverEm() < (2./ele.superCluster()->energy()+0.05) && 
		 (ele.full5x5_e2x5Max()/ele.full5x5_e5x5() > 0.94 || ele.full5x5_e1x5()/ele.full5x5_e5x5() > 0.83) &&
		 ele.dr03TkSumPt() < 5. && ele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 &&
		 iso < isoCut && fabs(dxy) < 0.02 ) isHEEP = true;
          }
          if( fabs(eta) > 1.566 && fabs(eta) < 2.5 ){
             iso = ele.dr03EcalRecHitSumEt() + ele.dr03HcalDepth1TowerSumEt();
             if( et <= 50 )
             	isoCut = 2.5 + 0.28*rho;
             else
             	isoCut = 2.5+0.03*(et-50.) + 0.28*rho;       
             if( ele.ecalDriven() == 1 && dEtaInSeed( ele ) < 0.006 && ele.deltaPhiSuperClusterTrackAtVtx() < 0.06 && 
	         ele.hadronicOverEm() < (12.5/ele.superCluster()->energy()+0.05) && ele.full5x5_sigmaIetaIeta() < 0.03 && 
		 ele.dr03TkSumPt() < 5. && ele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 &&
		 iso < isoCut && fabs(dxy) < 0.05 ) isHEEP = true;
          }  
	  
       }
       
    }

    if (ele.gsfTrack().isNonnull()){
        
       if( et > 35. ) {
          if( fabs(eta) < 1.4442 ){
             iso = ele.dr03EcalRecHitSumEt() + ele.dr03HcalDepth1TowerSumEt();
             isoCut = 2 + 0.03*et + 0.28*rho;	  
             if( ele.ecalDriven() == 1 && ele.deltaEtaSuperClusterTrackAtVtx() < std::max(0.016-1E-4*et,0.004) && ele.deltaPhiSuperClusterTrackAtVtx() < 0.06 && 
	         ele.hadronicOverEm() < (2./ele.superCluster()->energy()+0.05) && 
		 (ele.full5x5_e2x5Max()/ele.full5x5_e5x5() > 0.94 || ele.full5x5_e1x5()/ele.full5x5_e5x5() > 0.83) &&
		 ele.dr03TkSumPt() < 5. && ele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 &&
		 iso < isoCut && fabs(dxy) < 0.02 ) isHEEPv50 = true;
          }
          if( fabs(eta) > 1.566 && fabs(eta) < 2.5 ){
             iso = ele.dr03EcalRecHitSumEt() + ele.dr03HcalDepth1TowerSumEt();
             if( et <= 50 )
             	isoCut = 2.5 + 0.28*rho;
             else
             	isoCut = 2.5+0.03*(et-50.) + 0.28*rho;       
             if( ele.ecalDriven() == 1 && ele.deltaEtaSuperClusterTrackAtVtx() < std::max(0.015-8.5E-5*et,0.006) && ele.deltaPhiSuperClusterTrackAtVtx() < 0.06 && 
	         ele.hadronicOverEm() < (12.5/ele.superCluster()->energy()+0.05) && ele.full5x5_sigmaIetaIeta() < 0.03 && 
		 ele.dr03TkSumPt() < 5. && ele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 &&
		 iso < isoCut && fabs(dxy) < 0.05 ) isHEEPv50 = true;
          }  
	  
       }
       
    }
     
    //if( !isHEEP ) continue;
    nele++;
          

    nBranches_->lep_isHEEP			.push_back( isHEEP );
    nBranches_->lep_isHEEPv50		.push_back( isHEEPv50 );
    nBranches_->lep_isHighPtMuon	.push_back(-99);       
	 nBranches_->lep_isTightMuon	.push_back(-99); 
	 nBranches_->lep_isLooseMuon	.push_back(-99);       
    nBranches_->lep_type			.push_back(ele.pdgId());
    nBranches_->lep_charge			.push_back(ele.charge());
    nBranches_->lep_e				.push_back(ele.energy());
    nBranches_->lep_eta				.push_back(ele.superCluster()->eta());
    nBranches_->lep_mass			.push_back(ele.mass());
    nBranches_->lep_pt				.push_back(ele.pt());  
    nBranches_->lep_phi				.push_back(ele.phi());    

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

    nBranches_->lep_pfDeltaCorrRelIso.push_back(DeltaCorrectedIso);
    nBranches_->lep_pfRhoCorrRelIso04.push_back(RhoCorrectedIso04);
    nBranches_->lep_pfRhoCorrRelIso03.push_back(RhoCorrectedIso03);
    nBranches_->lep_pfRelIso	     .push_back((ele.chargedHadronIso() + ele.neutralHadronIso()+ ele.photonIso())/ele.pt());
    nBranches_->lep_photonIso	     .push_back(ele.photonIso());
    nBranches_->lep_neutralHadIso    .push_back(ele.neutralHadronIso());
    nBranches_->lep_chargedHadIso    .push_back(ele.chargedHadronIso());
    nBranches_->lep_trackIso	     .push_back(ele.trackIso());
 
  }
  
  nBranches_->nlep += nele;

}
