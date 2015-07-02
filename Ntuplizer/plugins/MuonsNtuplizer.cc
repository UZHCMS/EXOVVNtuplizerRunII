#include "../interface/MuonsNtuplizer.h"

#include <cmath>

//===================================================================================================================        
MuonsNtuplizer::MuonsNtuplizer( edm::EDGetTokenT<pat::MuonCollection>    muonToken   , 
                                edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
				edm::EDGetTokenT<double>                 rhoToken    , 
				NtupleBranches* nBranches )
	: CandidateNtuplizer( nBranches    )
	, muonToken_	    ( muonToken    )
	, verticeToken_     ( verticeToken )
	, rhoToken_	    ( rhoToken     )
	   
{
}

//===================================================================================================================
MuonsNtuplizer::~MuonsNtuplizer( void )
{

}

//===================================================================================================================
float MuonPFIso(pat::Muon muon, bool highpt){

  float sumChargedHadronPt = muon.pfIsolationR04().sumChargedHadronPt;
  float sumNeutralHadronEt = muon.pfIsolationR04().sumNeutralHadronEt;
  float sumPhotonEt = muon.pfIsolationR04().sumPhotonEt;
  float sumPUPt = muon.pfIsolationR04().sumPUPt;
  float iso = (sumChargedHadronPt + std::max( 0. ,sumNeutralHadronEt + sumPhotonEt - 0.5 * sumPUPt));// / muon.pt()
  //  if(highpt){
  //     reco::TrackRef cktTrack = (muon::tevOptimized(muon, 200, 17., 40., 0.25)).first;
  //     iso = (sumChargedHadronPt+ std::max ( 0. , sumNeutralHadronEt + sumPhotonEt - 0.5* sumPUPt ) )/cktTrack->pt();
  //   }
  return iso;
}

//===================================================================================================================
float MuonCorrPFIso(pat::Muon muon, bool highpt){

  float sumChargedHadronPt = muon.userIsolation(pat::PfChargedHadronIso);//pfIsolationR04().sumChargedHadronPt;
  float sumNeutralHadronEt = muon.userIsolation(pat::PfNeutralHadronIso);//pfIsolationR04().sumNeutralHadronEt;
  float sumPhotonEt = muon.userIsolation(pat::PfGammaIso);//pfIsolationR04().sumPhotonEt;
  float sumPUPt = muon.userIsolation(pat::User2Iso);//pfIsolationR04().sumPUPt;
  float iso = (sumChargedHadronPt+ std::max( 0., sumNeutralHadronEt + sumPhotonEt - 0.5 * sumPUPt ));// /muon.pt()
  // if(highpt){ 
  //     reco::TrackRef cktTrack = (muon::tevOptimized(muon, 200, 17., 40., 0.25)).first;
  //     iso = (sumChargedHadronPt+ std::max(0., sumNeutralHadronEt + sumPhotonEt - 0.5 * sumPUPt))/cktTrack->pt();
  //   }
  return iso;
}
  
//===================================================================================================================
void MuonsNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  event.getByToken(muonToken_	, muons_    ); 
  event.getByToken(verticeToken_, vertices_ ); 
  event.getByToken(rhoToken_	, rho_      );

  int nmus = 0;

  for (const pat::Muon &mu : *muons_) {
      
    nBranches_->lep_type   	     	    .push_back(mu.pdgId() );
    nBranches_->lep_charge 	     	    .push_back(mu.charge());
    nBranches_->lep_e	   	     	    .push_back(mu.energy());
    nBranches_->lep_eta    	     	    .push_back(mu.eta()   );
    nBranches_->lep_etaSC	            .push_back(mu.eta()   ); 
    nBranches_->lep_mass   	     	    .push_back(mu.mass()  );
    nBranches_->lep_pt     	     	    .push_back(mu.pt()    );
    nBranches_->lep_phi    		    .push_back(mu.phi()   );
    nBranches_->lep_isHeepElectron	    .push_back(-99);
    nBranches_->lep_isHeep51Electron	    .push_back(-99);
    nBranches_->lep_isLooseElectron         .push_back(-99);
    nBranches_->lep_passConversionVeto      .push_back(-99);
    nBranches_->lep_full5x5_sigmaIetaIeta   .push_back(-99);
    nBranches_->lep_dEtaIn		    .push_back(-99);
    nBranches_->lep_dPhiIn		    .push_back(-99);
    nBranches_->lep_hOverE		    .push_back(-99);
    nBranches_->lep_relIsoWithDBeta	    .push_back(-99);
    nBranches_->lep_ooEmooP		    .push_back(-99);
    nBranches_->lep_expectedMissingInnerHits.push_back(-99);
    nBranches_->lep_isVetoElectron	    .push_back(-99);
    nBranches_->lep_isMediumElectron	    .push_back(-99);
    nBranches_->lep_isTightElectron	    .push_back(-99);
    nBranches_->lep_TauType                 .push_back(0);      

    /*========== IDs ==============*/    
    nBranches_->lep_isHighPtMuon.push_back(mu.isHighPtMuon(vertices_->at(0)));
    nBranches_->lep_isTightMuon .push_back(mu.isTightMuon(vertices_->at(0)));
    nBranches_->lep_isLooseMuon .push_back(mu.isLooseMuon());
    nBranches_->lep_isPFMuon    .push_back(mu.isPFMuon());   

    double rho = *(rho_.product());     
    float deltaR = 0.3;
    double energy = TMath::Pi()*deltaR*deltaR*rho;
    float dxy = fabs(mu.muonBestTrack()->dxy(vertices_->at(0).position()));

    nBranches_->lep_d0  	.push_back(dxy);
    nBranches_->lep_dz  	.push_back(mu.muonBestTrack()->dz(vertices_->at(0).position()));
    nBranches_->lep_isGlobalMuon.push_back(mu.isGlobalMuon());  
    nBranches_->lep_isSoftMuon  .push_back(mu.isSoftMuon(vertices_->at(0)));  
      
    double normChi2	   = -99;
    int    trackerHits     = -99;
    int    pixelHits	   = -99;
    int    globalMuonHits  = -99;
  
    if( mu.isGlobalMuon() ) 
      normChi2=mu.normChi2();
  
    if( !mu.track().isNull() )
      trackerHits = (mu.track())->hitPattern().trackerLayersWithMeasurement();
  
    if( !mu.innerTrack().isNull() )
      pixelHits = (mu.innerTrack())->hitPattern().numberOfValidPixelHits();
  
    if( !mu.globalTrack().isNull() )
      globalMuonHits = (mu.globalTrack())->hitPattern().numberOfValidMuonHits();
  
    nBranches_->lep_normChi2	   .push_back(normChi2);
    nBranches_->lep_trackerHits    .push_back(trackerHits);
    nBranches_->lep_matchedStations.push_back(mu.numberOfMatchedStations());
    nBranches_->lep_pixelHits	   .push_back(pixelHits);
    nBranches_->lep_globalHits     .push_back(globalMuonHits);
        
    /*===== ISO ====*/
    deltaR = 0.3;
    energy = TMath::Pi()*deltaR*deltaR*rho;    
    nBranches_->lep_pfRhoCorrRelIso03.push_back((mu.chargedHadronIso() + std::max(0., mu.neutralHadronIso() + mu.photonIso() - energy))/mu.pt());
    nBranches_->lep_pfRhoCorrRelIso03Boost.push_back((mu.userIsolation(pat::PfChargedHadronIso) + std::max(0., mu.userIsolation(pat::PfNeutralHadronIso) + mu.userIsolation(pat::PfGammaIso) - energy))/mu.pt());
    
    deltaR = 0.4;
    energy = TMath::Pi()*deltaR*deltaR*rho;    
    nBranches_->lep_pfRhoCorrRelIso04.push_back((mu.chargedHadronIso() + std::max(0., mu.neutralHadronIso() + mu.photonIso() - energy))/mu.pt());
    nBranches_->lep_pfRhoCorrRelIso04Boost.push_back((mu.userIsolation(pat::PfChargedHadronIso) + std::max(0., mu.userIsolation(pat::PfNeutralHadronIso) + mu.userIsolation(pat::PfGammaIso) - energy))/mu.pt());
    
    nBranches_->lep_pfDeltaCorrRelIso     .push_back((mu.chargedHadronIso() + std::max(0., mu.neutralHadronIso() + mu.photonIso() - 0.5*mu.puChargedHadronIso()))/mu.pt());
    nBranches_->lep_pfRelIso	          .push_back((mu.chargedHadronIso() + mu.neutralHadronIso()+ mu.photonIso())/mu.pt()) ; 
    nBranches_->lep_photonIso	          .push_back(mu.photonIso());
    nBranches_->lep_neutralHadIso         .push_back(mu.neutralHadronIso());
    nBranches_->lep_chargedHadIso         .push_back(mu.chargedHadronIso());
    nBranches_->lep_trackIso	          .push_back(mu.trackIso());
    nBranches_->lep_pfDeltaCorrRelIsoBoost.push_back((mu.userIsolation(pat::PfChargedHadronIso) + std::max(0., mu.userIsolation(pat::PfNeutralHadronIso) + mu.userIsolation(pat::PfGammaIso) - 0.5*mu.userIsolation(pat::PfPUChargedHadronIso)))/mu.pt());
    nBranches_->lep_pfRelIsoBoost	  .push_back((mu.userIsolation(pat::PfChargedHadronIso) + mu.userIsolation(pat::PfNeutralHadronIso)+ mu.userIsolation(pat::PfGammaIso))/mu.pt()) ; 
    nBranches_->lep_photonIsoBoost	  .push_back(mu.userIsolation(pat::PfGammaIso));
    nBranches_->lep_neutralHadIsoBoost    .push_back(mu.userIsolation(pat::PfNeutralHadronIso));
    nBranches_->lep_chargedHadIsoBoost    .push_back(mu.userIsolation(pat::PfChargedHadronIso));
    nBranches_->lep_SemileptonicPFIso	  .push_back(MuonPFIso(mu,true));  
    nBranches_->lep_SemileptonicCorrPFIso .push_back(MuonCorrPFIso(mu,true));

    /*=======================*/

    nmus++;

    /*=======================*/  
    nBranches_->decayModeFindingNewDMs  		   .push_back(-99);
    nBranches_->decayModeFinding			   .push_back(-99);
    nBranches_->byLooseCombinedIsolationDeltaBetaCorr3Hits .push_back(-99);
    nBranches_->byMediumCombinedIsolationDeltaBetaCorr3Hits.push_back(-99);
    nBranches_->byTightCombinedIsolationDeltaBetaCorr3Hits .push_back(-99);
    nBranches_->byCombinedIsolationDeltaBetaCorrRaw3Hits   .push_back(-99);
    nBranches_->chargedIsoPtSum 			   .push_back(-99);
    nBranches_->neutralIsoPtSum 			   .push_back(-99);
    nBranches_->puCorrPtSum				   .push_back(-99);
    nBranches_->byIsolationMVA3oldDMwoLTraw		   .push_back(-99); 
    nBranches_->byVLooseIsolationMVA3oldDMwoLT  	   .push_back(-99);
    nBranches_->byLooseIsolationMVA3oldDMwoLT		   .push_back(-99);
    nBranches_->byMediumIsolationMVA3oldDMwoLT  	   .push_back(-99);
    nBranches_->byTightIsolationMVA3oldDMwoLT		   .push_back(-99);
    nBranches_->byVTightIsolationMVA3oldDMwoLT  	   .push_back(-99);
    nBranches_->byVVTightIsolationMVA3oldDMwoLT 	   .push_back(-99);
    nBranches_->byIsolationMVA3oldDMwLTraw		   .push_back(-99);
    nBranches_->byVLooseIsolationMVA3oldDMwLT		   .push_back(-99);
    nBranches_->byLooseIsolationMVA3oldDMwLT		   .push_back(-99);
    nBranches_->byMediumIsolationMVA3oldDMwLT		   .push_back(-99);
    nBranches_->byTightIsolationMVA3oldDMwLT		   .push_back(-99);
    nBranches_->byVTightIsolationMVA3oldDMwLT		   .push_back(-99);
    nBranches_->byVVTightIsolationMVA3oldDMwLT  	   .push_back(-99);
    nBranches_->byIsolationMVA3newDMwoLTraw		   .push_back(-99);
    nBranches_->byVLooseIsolationMVA3newDMwoLT  	   .push_back(-99);
    nBranches_->byLooseIsolationMVA3newDMwoLT		   .push_back(-99);
    nBranches_->byMediumIsolationMVA3newDMwoLT  	   .push_back(-99);
    nBranches_->byTightIsolationMVA3newDMwoLT		   .push_back(-99);
    nBranches_->byVTightIsolationMVA3newDMwoLT  	   .push_back(-99);
    nBranches_->byVVTightIsolationMVA3newDMwoLT 	   .push_back(-99);
    nBranches_->byIsolationMVA3newDMwLTraw		   .push_back(-99);
    nBranches_->byVLooseIsolationMVA3newDMwLT		   .push_back(-99);
    nBranches_->byLooseIsolationMVA3newDMwLT		   .push_back(-99);
    nBranches_->byMediumIsolationMVA3newDMwLT		   .push_back(-99);
    nBranches_->byTightIsolationMVA3newDMwLT		   .push_back(-99);
    nBranches_->byVTightIsolationMVA3newDMwLT		   .push_back(-99);
    nBranches_->byVVTightIsolationMVA3newDMwLT  	   .push_back(-99);
    nBranches_->againstElectronLoose			   .push_back(-99);
    nBranches_->againstElectronMedium			   .push_back(-99);
    nBranches_->againstElectronTight			   .push_back(-99);
    nBranches_->againstElectronMVA5raw  		   .push_back(-99);
    nBranches_->againstElectronMVA5category		   .push_back(-99);
    nBranches_->againstElectronVLooseMVA5		   .push_back(-99);
    nBranches_->againstElectronLooseMVA5		   .push_back(-99);
    nBranches_->againstElectronMediumMVA5		   .push_back(-99);
    nBranches_->againstElectronTightMVA5		   .push_back(-99);
    nBranches_->againstElectronVTightMVA5		   .push_back(-99);
    nBranches_->againstMuonLoose			   .push_back(-99);
    nBranches_->againstMuonMedium			   .push_back(-99);
    nBranches_->againstMuonTight			   .push_back(-99);
    nBranches_->againstMuonLoose2			   .push_back(-99);
    nBranches_->againstMuonMedium2			   .push_back(-99);
    nBranches_->againstMuonTight2			   .push_back(-99);
    nBranches_->againstMuonLoose3			   .push_back(-99);
    nBranches_->againstMuonTight3			   .push_back(-99);
    nBranches_->againstMuonMVAraw			   .push_back(-99);
    nBranches_->againstMuonLooseMVA			   .push_back(-99);
    nBranches_->againstMuonMediumMVA			   .push_back(-99);
    nBranches_->againstMuonTightMVA			   .push_back(-99);  
    /*=======================*/

  } 

  nBranches_->nlep +=  nmus;
    
}
