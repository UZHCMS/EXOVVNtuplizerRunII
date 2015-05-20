#include "../interface/ElectronsNtuplizer.h"
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// #include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
// #include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "Math/VectorUtil.h"

// #include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"
// #include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
// #include <cmath>

//===================================================================================================================
ElectronsNtuplizer::ElectronsNtuplizer( edm::EDGetTokenT<pat::ElectronCollection> electronToken, edm::EDGetTokenT<reco::VertexCollection> verticeToken, edm::EDGetTokenT<double> rhoToken, NtupleBranches* nBranches )
	: CandidateNtuplizer( nBranches     )
		, electronToken_					( electronToken )
		, verticeToken_						( verticeToken  )	
		, rhoToken_								( rhoToken      )
					// , electronVetoIdMapToken_	( eleIDtokens[0])
					// , electronLooseIdMapToken_	( eleIDtokens[1])
					// , electronMediumIdMapToken_( eleIDtokens[2])
					// , electronTightIdMapToken_	( eleIDtokens[3])

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
	
	// Get stuff for conversions
	 

	event.getByToken(electronToken_, electrons_); 
	event.getByToken(verticeToken_ , vertices_ );
	event.getByToken(rhoToken_	 , rho_	     );
	
	// event.getByToken(electronVetoIdMapToken_		,veto_id_decisions		);
	// 	event.getByToken(electronLooseIdMapToken_		,loose_id_decisions		);
	// 	event.getByToken(electronMediumIdMapToken_	,medium_id_decisions	);
	//   event.getByToken(electronTightIdMapToken_		,tight_id_decisions		);

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
          

		nBranches_->lep_isHEEP				.push_back( isHEEP );
		nBranches_->lep_isHEEPv50			.push_back( isHEEPv50 );
		nBranches_->lep_isHighPtMuon	.push_back(-99);       
		nBranches_->lep_isTightMuon		.push_back(-99); 
		nBranches_->lep_isLooseMuon		.push_back(-99);       
		nBranches_->lep_type			.push_back(ele.pdgId());
		nBranches_->lep_charge		.push_back(ele.charge());
		nBranches_->lep_e					.push_back(ele.energy());
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

		nBranches_->lep_pfDeltaCorrRelIso	.push_back(DeltaCorrectedIso);
		nBranches_->lep_pfRhoCorrRelIso04	.push_back(RhoCorrectedIso04);
		nBranches_->lep_pfRhoCorrRelIso03	.push_back(RhoCorrectedIso03);
		nBranches_->lep_pfRelIso					.push_back((ele.chargedHadronIso() + ele.neutralHadronIso()+ ele.photonIso())/ele.pt());
		nBranches_->lep_photonIso					.push_back(ele.photonIso());
		nBranches_->lep_neutralHadIso			.push_back(ele.neutralHadronIso());
		nBranches_->lep_chargedHadIso			.push_back(ele.chargedHadronIso());
		nBranches_->lep_trackIso					.push_back(ele.trackIso());
		
		
		// Ele veto, medium and tight ID starts here!!!	
		// All ID variables
		Float_t dEtaIn_;
		Float_t dPhiIn_;
		Float_t hOverE_;
		Float_t full5x5_sigmaIetaIeta_;
		Float_t relIsoWithDBeta_;
		Float_t ooEmooP_;
		Float_t d0_;
		Float_t dz_;
		Int_t   expectedMissingInnerHits_;
		Int_t   passConversionVeto_;
		bool		isVetoElectron		= false;
		bool		isMediumElectron	= false;
		bool		isTightElectron		= false;

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
		reco::GsfElectron::PflowIsolationVariables pfIso = ele.pfIsolationVariables();
		// Compute isolation with delta beta correction for PU
		float absiso = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
		relIsoWithDBeta_ = absiso/ele.pt();
		// Impact parameter
		d0_ = (-1) * ele.gsfTrack()->dxy(vertices_->at(0).position() );
		dz_ = ele.gsfTrack()->dz( vertices_->at(0).position() );
		// Conversion rejection
		expectedMissingInnerHits_ = ele.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
		passConversionVeto_ = ele.passConversionVeto();	
		
		//Barrel cuts
		
		if(fabs(eta) <= 1.479){
			if( passConversionVeto_									&&
				full5x5_sigmaIetaIeta_		< 0.011100	&&
				fabs(dEtaIn_)							< 0.016315	&&
				fabs(dPhiIn_)							< 0.252044	&&
				hOverE_										< 0.345843	&&
				relIsoWithDBeta_					< 0.164369	&&
				ooEmooP_									< 0.248070	&&
				fabs(d0_)									< 0.060279	&&
				fabs(dz_)									< 0.800538	&&				
				expectedMissingInnerHits_	<= 2.000000
					) isVetoElectron = true;
			
			if( passConversionVeto_									&&
				full5x5_sigmaIetaIeta_		< 0.010399	&&
				fabs(dEtaIn_)							< 0.007641	&&
				fabs(dPhiIn_)							< 0.032643	&&
				hOverE_										< 0.060662	&&
				relIsoWithDBeta_					< 0.097213	&&
				ooEmooP_									< 0.153897	&&
				fabs(d0_)									< 0.011811	&&
				fabs(dz_)									< 0.070775	&&				
				expectedMissingInnerHits_	<= 1.000000
					) isMediumElectron = true;
			
			if( passConversionVeto_									&&
				full5x5_sigmaIetaIeta_		< 0.010181	&&
				fabs(dEtaIn_)							< 0.006574	&&
				fabs(dPhiIn_)							< 0.022868	&&
				hOverE_										< 0.037553	&&
				relIsoWithDBeta_					< 0.074355	&&
				ooEmooP_									< 0.131191	&&
				fabs(d0_)									< 0.009924	&&
				fabs(dz_)									< 0.015310	&&				
				expectedMissingInnerHits_	<= 1.000000
					) isTightElectron = true;	
		}
		
		//Endcap cuts
		else
			if(fabs(eta) > 1.479 && fabs(eta) < 2.5){
				
				if( passConversionVeto_									&&
					full5x5_sigmaIetaIeta_		< 0.033987	&&
					fabs(dEtaIn_)							< 0.010671	&&
					fabs(dPhiIn_)							< 0.245263	&&
					hOverE_										< 0.134691	&&
					relIsoWithDBeta_					< 0.212604	&&
					ooEmooP_									< 0.157160	&&
					fabs(d0_)									< 0.273097	&&
					fabs(dz_)									< 0.885860	&&				
					expectedMissingInnerHits_	<= 3.000000
						) isVetoElectron = true;
				
				if( passConversionVeto_								&&
					full5x5_sigmaIetaIeta_		< 0.0318	&&
					fabs(dEtaIn_)							< 0.0108	&&
					fabs(dPhiIn_)							< 0.0455	&&
					hOverE_										< 0.097		&&
					relIsoWithDBeta_					< 0.254		&&
					ooEmooP_									< 0.1201	&&
					fabs(d0_)									< 0.0845	&&
					fabs(dz_)									< 0.7523	&&
					expectedMissingInnerHits_	< 1.020000 
						) isMediumElectron = true;  
				
				if( passConversionVeto_									&&
					full5x5_sigmaIetaIeta_		< 0.028766	&&
					fabs(dEtaIn_)							< 0.005681	&&
					fabs(dPhiIn_)							< 0.032046	&&
					hOverE_										< 0.081902	&&
					relIsoWithDBeta_					< 0.090185	&&
					ooEmooP_									< 0.106055	&&
					fabs(d0_)									< 0.027261	&&
					fabs(dz_)									< 0.147154	&&				
					expectedMissingInnerHits_	<= 1.000000
						) isTightElectron = true;  
			}
			
			nBranches_->lep_passConversionVeto				.push_back(passConversionVeto_);
			nBranches_->lep_full5x5_sigmaIetaIeta			.push_back(full5x5_sigmaIetaIeta_);
			nBranches_->lep_dEtaIn										.push_back(dEtaIn_);
			nBranches_->lep_dPhiIn										.push_back(dPhiIn_);
			nBranches_->lep_hOverE										.push_back(hOverE_);
			nBranches_->lep_relIsoWithDBeta						.push_back(relIsoWithDBeta_);
			nBranches_->lep_ooEmooP										.push_back(ooEmooP_);
			nBranches_->lep_d0												.push_back(d0_);
			nBranches_->lep_dz												.push_back(dz_);
			nBranches_->lep_expectedMissingInnerHits	.push_back(expectedMissingInnerHits_);
			nBranches_->lep_isVetoElectron						.push_back(isVetoElectron);
			nBranches_->lep_isMediumElectron					.push_back(isMediumElectron);
			nBranches_->lep_isTightElectron						.push_back(isTightElectron);
		
	}
  
	nBranches_->nlep += nele;

}
