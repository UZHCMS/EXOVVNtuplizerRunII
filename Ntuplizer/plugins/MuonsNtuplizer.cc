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
      
    nBranches_->mu_pdgId   	     	    .push_back(mu.pdgId() );
    nBranches_->mu_charge 	     	    .push_back(mu.charge());
    nBranches_->mu_e	   	     	    .push_back(mu.energy());
    nBranches_->mu_eta    	     	    .push_back(mu.eta()   );
    nBranches_->mu_mass   	     	    .push_back(mu.mass()  );
    nBranches_->mu_pt     	     	    .push_back(mu.pt()    );
    nBranches_->mu_phi    		    .push_back(mu.phi()   );

    /*========== IDs ==============*/    
    nBranches_->mu_isHighPtMuon.push_back(mu.isHighPtMuon(vertices_->at(0)));
    nBranches_->mu_isTightMuon .push_back(mu.isTightMuon(vertices_->at(0)));
    nBranches_->mu_isLooseMuon .push_back(mu.isLooseMuon());
    nBranches_->mu_isPFMuon    .push_back(mu.isPFMuon());   

    double rho = *(rho_.product());     
    float deltaR = 0.3;
    double energy = TMath::Pi()*deltaR*deltaR*rho;
    float dxy = fabs(mu.muonBestTrack()->dxy(vertices_->at(0).position()));

    nBranches_->mu_d0  	.push_back(dxy);
    nBranches_->mu_dz  	.push_back(mu.muonBestTrack()->dz(vertices_->at(0).position()));
    nBranches_->mu_isGlobalMuon.push_back(mu.isGlobalMuon());  
    nBranches_->mu_isSoftMuon  .push_back(mu.isSoftMuon(vertices_->at(0)));  
      
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
  
    nBranches_->mu_normChi2	   .push_back(normChi2);
    nBranches_->mu_trackerHits    .push_back(trackerHits);
    nBranches_->mu_matchedStations.push_back(mu.numberOfMatchedStations());
    nBranches_->mu_pixelHits	   .push_back(pixelHits);
    nBranches_->mu_globalHits     .push_back(globalMuonHits);
        
    /*===== ISO ====*/
    deltaR = 0.3;
    energy = TMath::Pi()*deltaR*deltaR*rho;    
    nBranches_->mu_pfRhoCorrRelIso03.push_back((mu.chargedHadronIso() + std::max(0., mu.neutralHadronIso() + mu.photonIso() - energy))/mu.pt());
    nBranches_->mu_pfRhoCorrRelIso03Boost.push_back((mu.userIsolation(pat::PfChargedHadronIso) + std::max(0., mu.userIsolation(pat::PfNeutralHadronIso) + mu.userIsolation(pat::PfGammaIso) - energy))/mu.pt());
    
    deltaR = 0.4;
    energy = TMath::Pi()*deltaR*deltaR*rho;    
    nBranches_->mu_pfRhoCorrRelIso04.push_back((mu.chargedHadronIso() + std::max(0., mu.neutralHadronIso() + mu.photonIso() - energy))/mu.pt());
    nBranches_->mu_pfRhoCorrRelIso04Boost.push_back((mu.userIsolation(pat::PfChargedHadronIso) + std::max(0., mu.userIsolation(pat::PfNeutralHadronIso) + mu.userIsolation(pat::PfGammaIso) - energy))/mu.pt());
    
    nBranches_->mu_pfDeltaCorrRelIso     .push_back((mu.chargedHadronIso() + std::max(0., mu.neutralHadronIso() + mu.photonIso() - 0.5*mu.puChargedHadronIso()))/mu.pt());
    nBranches_->mu_pfRelIso	          .push_back((mu.chargedHadronIso() + mu.neutralHadronIso()+ mu.photonIso())/mu.pt()) ; 
    nBranches_->mu_photonIso	          .push_back(mu.photonIso());
    nBranches_->mu_neutralHadIso         .push_back(mu.neutralHadronIso());
    nBranches_->mu_chargedHadIso         .push_back(mu.chargedHadronIso());
    nBranches_->mu_trackIso	          .push_back(mu.trackIso());
    nBranches_->mu_pfDeltaCorrRelIsoBoost.push_back((mu.userIsolation(pat::PfChargedHadronIso) + std::max(0., mu.userIsolation(pat::PfNeutralHadronIso) + mu.userIsolation(pat::PfGammaIso) - 0.5*mu.userIsolation(pat::PfPUChargedHadronIso)))/mu.pt());
    nBranches_->mu_pfRelIsoBoost	  .push_back((mu.userIsolation(pat::PfChargedHadronIso) + mu.userIsolation(pat::PfNeutralHadronIso)+ mu.userIsolation(pat::PfGammaIso))/mu.pt()) ; 
    nBranches_->mu_photonIsoBoost	  .push_back(mu.userIsolation(pat::PfGammaIso));
    nBranches_->mu_neutralHadIsoBoost    .push_back(mu.userIsolation(pat::PfNeutralHadronIso));
    nBranches_->mu_chargedHadIsoBoost    .push_back(mu.userIsolation(pat::PfChargedHadronIso));
    nBranches_->mu_SemileptonicPFIso	  .push_back(MuonPFIso(mu,true));  
    nBranches_->mu_SemileptonicCorrPFIso .push_back(MuonCorrPFIso(mu,true));

    /*=======================*/

    ++nmus;



  } 

  nBranches_->mu_N =  nmus;
    
}
