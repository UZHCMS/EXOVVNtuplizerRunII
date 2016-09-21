#include "../interface/TausNtuplizer.h"

#include <cmath>

#include "TVector3.h"
#include "TLorentzVector.h"

//===================================================================================================================        
TausNtuplizer::TausNtuplizer( edm::EDGetTokenT<pat::TauCollection> tauToken,edm::EDGetTokenT<pat::TauCollection> tauBoostedTauToken,
			      edm::EDGetTokenT<double> rhoToken, 
			      edm::EDGetTokenT<reco::VertexCollection> verticeToken,
			      NtupleBranches* nBranches,
			      std::map< std::string, bool >& runFlags )
  : CandidateNtuplizer( nBranches )
  , tauInputToken_ (tauToken)
  , tauBoostedTauInputToken_ ( tauBoostedTauToken)
  , rhoToken_	       	    ( rhoToken  )
  , verticeToken_     	    ( verticeToken  )	 
  , doBoostedTaus_     	    ( runFlags["doBoostedTaus"]  )	 
 
{
  
}

//===================================================================================================================
TausNtuplizer::~TausNtuplizer( void )
{

}

//===================================================================================================================
void TausNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
 
 

  if ( !doBoostedTaus_ ) return;
  event.getByToken( tauBoostedTauInputToken_ , boostedTaus_  ); 
  event.getByToken( rhoToken_	 	   , rho_      );
    event.getByToken( verticeToken_ 	   , vertices_ );
    event.getByToken(  tauInputToken_ , taus_ ); 
    /********************************************************************/    
    for( size_t t = 0; t < taus_->size(); ++t ){
      
      pat::Tau tau = (*taus_)[t];
    
      nBranches_->tau_pdgId  	     	      .push_back(tau.pdgId());
      nBranches_->tau_charge 	     	      .push_back(tau.charge());
      nBranches_->tau_e      	     	      .push_back(tau.energy());
      nBranches_->tau_eta    	     	      .push_back(tau.eta());
      nBranches_->tau_mass   	     	      .push_back(tau.mass());
      nBranches_->tau_pt     	     	      .push_back(tau.pt());
      nBranches_->tau_phi    	     	      .push_back(tau.phi());
      nBranches_->tau_d0	      	      .push_back(tau.dxy());
      nBranches_->tau_TauType	     	      .push_back(1);
      nBranches_->tau_decayMode	     	      .push_back(tau.decayMode());

      // YT added : 17 Aug 2016
      pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
      nBranches_->tau_dz	      	      .push_back(fabs(packedLeadTauCand->dz()));

      // For the tau polarization measurement
      if(tau.decayMode() >= 1 && tau.decayMode() <= 4){
	Float_t pt_neutral = 0;
	Float_t pt_charged = 0;
	
	// This crashes ... don't know why ... 
	//	for(auto it = tau.signalCands().begin(), end = tau.signalCands().end(); it!=end; ++it){
	//	  std::cout << (*it)->pdgId() << std::endl;
	//	}

	for(int ii = 0; ii < (int)tau.signalCands().size(); ii++){
	  Int_t pdg = abs(tau.signalCands()[ii]->pdgId());
	  Float_t candpt = tau.signalCands()[ii]->pt();
	  
	  if(pdg == 11 || pdg == 22 || pdg == 130) pt_neutral += candpt;
	  if(pdg == 211){
	    if(candpt > pt_charged){
	      pt_charged = candpt;
	    }
	  }
	}
	
	if(pt_charged > 0){
	  nBranches_->tau_chargedPionPt.push_back(pt_charged);
	  nBranches_->tau_neutralPionPt.push_back(pt_neutral);
	}else{
	  nBranches_->tau_chargedPionPt.push_back(-1);
	  nBranches_->tau_neutralPionPt.push_back(-1);
	}
      }else{
	nBranches_->tau_chargedPionPt.push_back(-99);
	nBranches_->tau_neutralPionPt.push_back(-99);
      }

         
      /*====================== ISO ========================*/	         
      double rho = *(rho_.product());
      float  deltaR = 0.3;
      double energy = TMath::Pi()*deltaR*deltaR*rho;      
      nBranches_->tau_pfRhoCorrRelIso03.push_back((tau.chargedHadronIso() + std::max(0.,tau.neutralHadronIso() + tau.photonIso() - energy))/tau.pt());      
      nBranches_->tau_pfRhoCorrRelIso03Boost.push_back((tau.userIsolation(pat::PfChargedHadronIso) + std::max(0.,tau.userIsolation(pat::PfNeutralHadronIso) + tau.userIsolation(pat::PfGammaIso) - energy))/tau.pt()); 

      deltaR = 0.4;
      energy = TMath::Pi()*deltaR*deltaR*rho;      
      nBranches_->tau_pfRhoCorrRelIso04.push_back((tau.chargedHadronIso() + std::max(0.,tau.neutralHadronIso() + tau.photonIso() - energy))/tau.pt());
      nBranches_->tau_pfRhoCorrRelIso04Boost.push_back((tau.userIsolation(pat::PfChargedHadronIso) + std::max(0.,tau.userIsolation(pat::PfNeutralHadronIso) + tau.userIsolation(pat::PfGammaIso) - energy))/tau.pt()); 
            
      nBranches_->tau_pfDeltaCorrRelIso     .push_back((tau.chargedHadronIso() + std::max(0., tau.neutralHadronIso() + tau.photonIso() - 0.5*tau.puChargedHadronIso()))/tau.pt());
      nBranches_->tau_pfRelIso              .push_back((tau.chargedHadronIso() + tau.neutralHadronIso()+ tau.photonIso())/tau.pt()) ; 
      nBranches_->tau_photonIso             .push_back(tau.photonIso());
      nBranches_->tau_neutralHadIso         .push_back(tau.neutralHadronIso());
      nBranches_->tau_chargedHadIso         .push_back(tau.chargedHadronIso());
      nBranches_->tau_trackIso              .push_back(tau.trackIso());      
      nBranches_->tau_pfDeltaCorrRelIsoBoost.push_back((tau.userIsolation(pat::PfChargedHadronIso) + std::max(0., tau.userIsolation(pat::PfNeutralHadronIso) + tau.userIsolation(pat::PfGammaIso) - 0.5*tau.userIsolation(pat::PfPUChargedHadronIso)))/tau.pt());
      nBranches_->tau_pfRelIsoBoost         .push_back((tau.userIsolation(pat::PfChargedHadronIso) + tau.userIsolation(pat::PfNeutralHadronIso)+ tau.userIsolation(pat::PfGammaIso))/tau.pt()) ; 
      nBranches_->tau_photonIsoBoost        .push_back(tau.userIsolation(pat::PfGammaIso));
      nBranches_->tau_neutralHadIsoBoost    .push_back(tau.userIsolation(pat::PfNeutralHadronIso));
      nBranches_->tau_chargedHadIsoBoost    .push_back(tau.userIsolation(pat::PfChargedHadronIso));

      /*====================== IDs ========================*/	                         
      nBranches_->tau_decayModeFindingNewDMs                     .push_back(tau.tauID("decayModeFindingNewDMs"  		    ));
      nBranches_->tau_decayModeFinding	                         .push_back(tau.tauID("decayModeFinding"			    ));
      nBranches_->tau_byLooseCombinedIsolationDeltaBetaCorr3Hits .push_back(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"  ));
      nBranches_->tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.push_back(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits" ));
      nBranches_->tau_byTightCombinedIsolationDeltaBetaCorr3Hits .push_back(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"  ));
      nBranches_->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits   .push_back(tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"    ));
      nBranches_->tau_chargedIsoPtSum                            .push_back(tau.tauID("chargedIsoPtSum" 			    ));
      nBranches_->tau_neutralIsoPtSum                            .push_back(tau.tauID("neutralIsoPtSum" 			    ));
      nBranches_->tau_puCorrPtSum                                .push_back(tau.tauID("puCorrPtSum"				    ));
      
      nBranches_->tau_chargedIsoPtSumdR03                        .push_back(tau.tauID("chargedIsoPtSumdR03" 			    ));
      nBranches_->tau_footprintCorrectiondR03                    .push_back(tau.tauID("footprintCorrectiondR03"                         ));
      nBranches_->tau_neutralIsoPtSumdR03                        .push_back(tau.tauID("neutralIsoPtSumdR03" 			    ));
      nBranches_->tau_neutralIsoPtSumWeight                      .push_back(tau.tauID("neutralIsoPtSumWeight" 			    ));
      nBranches_->tau_neutralIsoPtSumWeightdR03                  .push_back(tau.tauID("neutralIsoPtSumWeightdR03" 			    ));
      nBranches_->tau_photonPtSumOutsideSignalConedR03           .push_back(tau.tauID("photonPtSumOutsideSignalConedR03"               ));

      nBranches_->tau_byIsolationMVArun2v1DBdR03oldDMwLTraw      .push_back(tau.tauID("byIsolationMVArun2v1DBdR03oldDMwLTraw"));
      nBranches_->tau_byIsolationMVArun2v1DBnewDMwLTraw          .push_back(tau.tauID("byIsolationMVArun2v1DBnewDMwLTraw"));
      nBranches_->tau_byIsolationMVArun2v1DBoldDMwLTraw          .push_back(tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw"				    ));
      nBranches_->tau_byIsolationMVArun2v1PWdR03oldDMwLTraw      .push_back(tau.tauID("byIsolationMVArun2v1PWdR03oldDMwLTraw"				    ));
      nBranches_->tau_byIsolationMVArun2v1PWnewDMwLTraw          .push_back(tau.tauID("byIsolationMVArun2v1PWnewDMwLTraw"				    ));
      nBranches_->tau_byIsolationMVArun2v1PWoldDMwLTraw          .push_back(tau.tauID("byIsolationMVArun2v1PWoldDMwLTraw"				    ));
      nBranches_->tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT    .push_back(tau.tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT"				    ));
      nBranches_->tau_byLooseIsolationMVArun2v1DBnewDMwLT        .push_back(tau.tauID("byLooseIsolationMVArun2v1DBnewDMwLT"				    ));
      nBranches_->tau_byLooseIsolationMVArun2v1DBoldDMwLT        .push_back(tau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT"			    ));
      nBranches_->tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT    .push_back(tau.tauID("byLooseIsolationMVArun2v1PWdR03oldDMwLT"				    ));
      nBranches_->tau_byLooseIsolationMVArun2v1PWnewDMwLT        .push_back(	tau.tauID("byLooseIsolationMVArun2v1PWnewDMwLT"			    ));
      nBranches_->tau_byLooseIsolationMVArun2v1PWoldDMwLT        .push_back(	tau.tauID("byLooseIsolationMVArun2v1PWoldDMwLT"			    ));
      
      nBranches_->tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT    .push_back(tau.tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT"				    ));
      nBranches_->tau_byMediumIsolationMVArun2v1DBnewDMwLT        .push_back(	tau.tauID("byMediumIsolationMVArun2v1DBnewDMwLT"			    ));
      nBranches_->tau_byMediumIsolationMVArun2v1DBoldDMwLT        .push_back(tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT"				    ));
      nBranches_->tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT    .push_back(tau.tauID("byMediumIsolationMVArun2v1PWdR03oldDMwLT"				    ));
      nBranches_->tau_byMediumIsolationMVArun2v1PWnewDMwLT        .push_back(tau.tauID("byMediumIsolationMVArun2v1PWnewDMwLT"				    ));
      nBranches_->tau_byMediumIsolationMVArun2v1PWoldDMwLT        .push_back(tau.tauID("byMediumIsolationMVArun2v1PWoldDMwLT"				    ));

    
      nBranches_->tau_byTightIsolationMVArun2v1DBdR03oldDMwLT     .push_back(tau.tauID("byTightIsolationMVArun2v1DBdR03oldDMwLT"				    ));
      nBranches_->tau_byTightIsolationMVArun2v1DBnewDMwLT         .push_back(tau.tauID("byTightIsolationMVArun2v1DBnewDMwLT"			    ));
      nBranches_->tau_byTightIsolationMVArun2v1DBoldDMwLT         .push_back(tau.tauID("byTightIsolationMVArun2v1DBoldDMwLT"				    ));
      nBranches_->tau_byTightIsolationMVArun2v1PWdR03oldDMwLT     .push_back(tau.tauID("byTightIsolationMVArun2v1PWdR03oldDMwLT"				    ));
      nBranches_->tau_byTightIsolationMVArun2v1PWnewDMwLT         .push_back(tau.tauID("byTightIsolationMVArun2v1PWnewDMwLT"			    ));
      nBranches_->tau_byTightIsolationMVArun2v1PWoldDMwLT         .push_back(tau.tauID("byTightIsolationMVArun2v1PWoldDMwLT"			    ));
      nBranches_->tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT    .push_back(tau.tauID("byVLooseIsolationMVArun2v1DBdR03oldDMwLT"				    ));

      nBranches_->tau_byVLooseIsolationMVArun2v1DBnewDMwLT        .push_back(	tau.tauID("byVLooseIsolationMVArun2v1DBnewDMwLT"			    ));
      nBranches_->tau_byVLooseIsolationMVArun2v1DBoldDMwLT        .push_back(tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"			    ));
      nBranches_->tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT    .push_back(tau.tauID("byVLooseIsolationMVArun2v1PWdR03oldDMwLT"			    ));
      nBranches_->tau_byVLooseIsolationMVArun2v1PWnewDMwLT        .push_back(tau.tauID("byVLooseIsolationMVArun2v1PWnewDMwLT"			    ));
      nBranches_->tau_byVLooseIsolationMVArun2v1PWoldDMwLT        .push_back(tau.tauID("byVLooseIsolationMVArun2v1PWoldDMwLT"			    ));
      nBranches_->tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT    .push_back(tau.tauID("byVTightIsolationMVArun2v1DBdR03oldDMwLT"			    ));
      nBranches_->tau_byVTightIsolationMVArun2v1DBnewDMwLT        .push_back(tau.tauID("byVTightIsolationMVArun2v1DBnewDMwLT"			    ));
      nBranches_->tau_byVTightIsolationMVArun2v1DBoldDMwLT        .push_back(tau.tauID("byVTightIsolationMVArun2v1DBoldDMwLT"			    ));

      nBranches_->tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT     .push_back(tau.tauID("byVTightIsolationMVArun2v1PWdR03oldDMwLT"			    ));
      nBranches_->tau_byVTightIsolationMVArun2v1PWnewDMwLT         .push_back(tau.tauID("byVTightIsolationMVArun2v1PWnewDMwLT"			    ));
      nBranches_->tau_byVTightIsolationMVArun2v1PWoldDMwLT         .push_back(tau.tauID("byVTightIsolationMVArun2v1PWoldDMwLT"			    ));
      nBranches_->tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT    .push_back(tau.tauID("byVVTightIsolationMVArun2v1DBdR03oldDMwLT"			    ));
      nBranches_->tau_byVVTightIsolationMVArun2v1DBnewDMwLT        .push_back(tau.tauID("byVVTightIsolationMVArun2v1DBnewDMwLT"			    ));
      nBranches_->tau_byVVTightIsolationMVArun2v1DBoldDMwLT        .push_back(tau.tauID("byVVTightIsolationMVArun2v1DBoldDMwLT"			    ));
      nBranches_->tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT    .push_back(tau.tauID("byVVTightIsolationMVArun2v1PWdR03oldDMwLT"			    ));
      nBranches_->tau_byVVTightIsolationMVArun2v1PWnewDMwLT        .push_back(tau.tauID("byVVTightIsolationMVArun2v1PWnewDMwLT"			    ));

      nBranches_->tau_byVVTightIsolationMVArun2v1PWoldDMwLT        .push_back(tau.tauID("byVVTightIsolationMVArun2v1PWoldDMwLT"			    ));




      nBranches_->tau_againstElectronMVA6raw                     .push_back(tau.tauID("againstElectronMVA6Raw"  		    ));
      nBranches_->tau_againstElectronMVA6category                .push_back(tau.tauID("againstElectronMVA6category"		    ));
      nBranches_->tau_againstElectronVLooseMVA6                  .push_back(tau.tauID("againstElectronVLooseMVA6"		    ));
      nBranches_->tau_againstElectronLooseMVA6                   .push_back(tau.tauID("againstElectronLooseMVA6"		    ));
      nBranches_->tau_againstElectronMediumMVA6                  .push_back(tau.tauID("againstElectronMediumMVA6"		    ));
      nBranches_->tau_againstElectronTightMVA6                   .push_back(tau.tauID("againstElectronTightMVA6"		    ));
      nBranches_->tau_againstElectronVTightMVA6                  .push_back(tau.tauID("againstElectronVTightMVA6"		    ));
      

      nBranches_->tau_againstMuonLoose3                          .push_back(tau.tauID("againstMuonLoose3"			    ));
      nBranches_->tau_againstMuonTight3                          .push_back(tau.tauID("againstMuonTight3"			    ));
      
      // nBranches_->tau_byPileupWeightedIsolationRaw3Hits          .push_back(tau.tauID("byPileupWeightedIsolationRaw3Hits"           ));
      // nBranches_->tau_byLoosePileupWeightedIsolation3Hits        .push_back(tau.tauID("byLoosePileupWeightedIsolation3Hits"         ));
      // nBranches_->tau_byMediumPileupWeightedIsolation3Hits       .push_back(tau.tauID("byMediumPileupWeightedIsolation3Hits"        ));
      // nBranches_->tau_byTightPileupWeightedIsolation3Hits        .push_back(tau.tauID("byTightPileupWeightedIsolation3Hits"         ));
      nBranches_->tau_byPhotonPtSumOutsideSignalCone             .push_back(tau.tauID("byPhotonPtSumOutsideSignalCone"              ));
      nBranches_->tau_footprintCorrection                        .push_back(tau.tauID("footprintCorrection"                         ));



      /*======================================================*/            

  }

  /********************************************************************/      
  for( size_t t = 0; t < boostedTaus_->size(); ++t ){
  
      pat::Tau boostedTau = (*boostedTaus_)[t];
              
      nBranches_->tau_pdgId  		      .push_back(boostedTau.pdgId());
      nBranches_->tau_charge 		      .push_back(boostedTau.charge());
      nBranches_->tau_e      		      .push_back(boostedTau.energy());
      nBranches_->tau_eta    		      .push_back(boostedTau.eta()); 
      nBranches_->tau_mass   		      .push_back(boostedTau.mass());
      nBranches_->tau_pt     		      .push_back(boostedTau.pt());
      nBranches_->tau_phi    		      .push_back(boostedTau.phi());
      nBranches_->tau_TauType		      .push_back(2);
      nBranches_->tau_d0	              .push_back(boostedTau.dxy());
      nBranches_->tau_decayMode		      .push_back(boostedTau.decayMode());

      // YT added : 17 Aug 2016
      pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(boostedTau.leadChargedHadrCand().get());
      nBranches_->tau_dz	      	      .push_back(fabs(packedLeadTauCand->dz()));

      if(boostedTau.decayMode() >= 1 && boostedTau.decayMode() <= 4){
	Float_t pt_neutral = 0;
	Float_t pt_charged = 0;

	for(int ii = 0; ii < (int)boostedTau.signalCands().size(); ii++){
	  Int_t pdg = abs(boostedTau.signalCands()[ii]->pdgId());
	  Float_t candpt = boostedTau.signalCands()[ii]->pt();
	  
	  if(pdg == 11 || pdg == 22 || pdg == 130) pt_neutral += candpt;
	  if(pdg == 211){
	    if(candpt > pt_charged){
	      pt_charged = candpt;
	    }
	  }
	}
	
	if(pt_charged > 0){
	  nBranches_->tau_chargedPionPt.push_back(pt_charged);
	  nBranches_->tau_neutralPionPt.push_back(pt_neutral);
	}else{
	  nBranches_->tau_chargedPionPt.push_back(-1);
	  nBranches_->tau_neutralPionPt.push_back(-1);
	}
      }else{
	nBranches_->tau_chargedPionPt.push_back(-99);
	nBranches_->tau_neutralPionPt.push_back(-99);
      }



      /*====================== ISO ========================*/	          
      double rho = *(rho_.product());
      float  deltaR = 0.3;
      double energy = TMath::Pi()*deltaR*deltaR*rho;      
      nBranches_->tau_pfRhoCorrRelIso03.push_back((boostedTau.chargedHadronIso() + std::max(0.,boostedTau.neutralHadronIso() + boostedTau.photonIso() - energy))/boostedTau.pt());      
      nBranches_->tau_pfRhoCorrRelIso03Boost.push_back((boostedTau.userIsolation(pat::PfChargedHadronIso) + std::max(0.,boostedTau.userIsolation(pat::PfNeutralHadronIso) + boostedTau.userIsolation(pat::PfGammaIso) - energy))/boostedTau.pt()); 

      deltaR = 0.4;
      energy = TMath::Pi()*deltaR*deltaR*rho;      
      nBranches_->tau_pfRhoCorrRelIso04.push_back((boostedTau.chargedHadronIso() + std::max(0.,boostedTau.neutralHadronIso() + boostedTau.photonIso() - energy))/boostedTau.pt());
      nBranches_->tau_pfRhoCorrRelIso04Boost.push_back((boostedTau.userIsolation(pat::PfChargedHadronIso) + std::max(0.,boostedTau.userIsolation(pat::PfNeutralHadronIso) + boostedTau.userIsolation(pat::PfGammaIso) - energy))/boostedTau.pt()); 
            
      nBranches_->tau_pfDeltaCorrRelIso     .push_back((boostedTau.chargedHadronIso() + std::max(0., boostedTau.neutralHadronIso() + boostedTau.photonIso() - 0.5*boostedTau.puChargedHadronIso()))/boostedTau.pt());
      nBranches_->tau_pfRelIso              .push_back((boostedTau.chargedHadronIso() + boostedTau.neutralHadronIso()+ boostedTau.photonIso())/boostedTau.pt()) ; 
      nBranches_->tau_photonIso             .push_back(boostedTau.photonIso());
      nBranches_->tau_neutralHadIso         .push_back(boostedTau.neutralHadronIso());
      nBranches_->tau_chargedHadIso         .push_back(boostedTau.chargedHadronIso());
      nBranches_->tau_trackIso              .push_back(boostedTau.trackIso());            
      nBranches_->tau_pfDeltaCorrRelIsoBoost.push_back((boostedTau.userIsolation(pat::PfChargedHadronIso) + std::max(0., boostedTau.userIsolation(pat::PfNeutralHadronIso) + boostedTau.userIsolation(pat::PfGammaIso) - 0.5*boostedTau.userIsolation(pat::PfPUChargedHadronIso)))/boostedTau.pt());
      nBranches_->tau_pfRelIsoBoost         .push_back((boostedTau.userIsolation(pat::PfChargedHadronIso) + boostedTau.userIsolation(pat::PfNeutralHadronIso)+ boostedTau.userIsolation(pat::PfGammaIso))/boostedTau.pt()) ; 
      nBranches_->tau_photonIsoBoost        .push_back(boostedTau.userIsolation(pat::PfGammaIso));
      nBranches_->tau_neutralHadIsoBoost    .push_back(boostedTau.userIsolation(pat::PfNeutralHadronIso));
      nBranches_->tau_chargedHadIsoBoost    .push_back(boostedTau.userIsolation(pat::PfChargedHadronIso));

      /*====================== IDs ========================*/	                                     
      
      nBranches_->tau_decayModeFindingNewDMs                     .push_back(boostedTau.tauID("decayModeFindingNewDMs"  		    ));
      nBranches_->tau_decayModeFinding	                         .push_back(boostedTau.tauID("decayModeFinding"			    ));
      nBranches_->tau_byLooseCombinedIsolationDeltaBetaCorr3Hits .push_back(boostedTau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"  ));
      nBranches_->tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.push_back(boostedTau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits" ));
      nBranches_->tau_byTightCombinedIsolationDeltaBetaCorr3Hits .push_back(boostedTau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"  ));
      nBranches_->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits   .push_back(boostedTau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"    ));
      nBranches_->tau_chargedIsoPtSum                            .push_back(boostedTau.tauID("chargedIsoPtSum" 			    ));
      nBranches_->tau_neutralIsoPtSum                            .push_back(boostedTau.tauID("neutralIsoPtSum" 			    ));
      nBranches_->tau_puCorrPtSum                                .push_back(boostedTau.tauID("puCorrPtSum"				    ));
       nBranches_->tau_chargedIsoPtSumdR03                       .push_back(boostedTau.tauID("chargedIsoPtSumdR03" 			    ));
      nBranches_->tau_footprintCorrectiondR03                    .push_back(boostedTau.tauID("footprintCorrectiondR03"                         ));
      nBranches_->tau_neutralIsoPtSumdR03                        .push_back(boostedTau.tauID("neutralIsoPtSumdR03" 			    ));
      nBranches_->tau_neutralIsoPtSumWeight                      .push_back(boostedTau.tauID("neutralIsoPtSumWeight" 			    ));
      nBranches_->tau_neutralIsoPtSumWeightdR03                  .push_back(boostedTau.tauID("neutralIsoPtSumWeightdR03" 			    ));
      nBranches_->tau_photonPtSumOutsideSignalConedR03           .push_back(boostedTau.tauID("photonPtSumOutsideSignalConedR03"               ));

      nBranches_->tau_byIsolationMVArun2v1DBdR03oldDMwLTraw      .push_back(boostedTau.tauID("byIsolationMVArun2v1DBdR03oldDMwLTraw"));
      nBranches_->tau_byIsolationMVArun2v1DBnewDMwLTraw          .push_back(boostedTau.tauID("byIsolationMVArun2v1DBnewDMwLTraw"));
      nBranches_->tau_byIsolationMVArun2v1DBoldDMwLTraw          .push_back(boostedTau.tauID("byIsolationMVArun2v1DBoldDMwLTraw"				    ));
      nBranches_->tau_byIsolationMVArun2v1PWdR03oldDMwLTraw      .push_back(boostedTau.tauID("byIsolationMVArun2v1PWdR03oldDMwLTraw"				    ));
      nBranches_->tau_byIsolationMVArun2v1PWnewDMwLTraw          .push_back(boostedTau.tauID("byIsolationMVArun2v1PWnewDMwLTraw"				    ));
      nBranches_->tau_byIsolationMVArun2v1PWoldDMwLTraw          .push_back(boostedTau.tauID("byIsolationMVArun2v1PWoldDMwLTraw"				    ));
      nBranches_->tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT    .push_back(boostedTau.tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT"				    ));
      nBranches_->tau_byLooseIsolationMVArun2v1DBnewDMwLT        .push_back(boostedTau.tauID("byLooseIsolationMVArun2v1DBnewDMwLT"				    ));
      nBranches_->tau_byLooseIsolationMVArun2v1DBoldDMwLT        .push_back(boostedTau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT"			    ));
      nBranches_->tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT    .push_back(boostedTau.tauID("byLooseIsolationMVArun2v1PWdR03oldDMwLT"				    ));
      nBranches_->tau_byLooseIsolationMVArun2v1PWnewDMwLT        .push_back(	boostedTau.tauID("byLooseIsolationMVArun2v1PWnewDMwLT"			    ));
      nBranches_->tau_byLooseIsolationMVArun2v1PWoldDMwLT        .push_back(	boostedTau.tauID("byLooseIsolationMVArun2v1PWoldDMwLT"			    ));
       nBranches_->tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT    .push_back(boostedTau.tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT"				    ));
      nBranches_->tau_byMediumIsolationMVArun2v1DBnewDMwLT        .push_back(	boostedTau.tauID("byMediumIsolationMVArun2v1DBnewDMwLT"			    ));
      nBranches_->tau_byMediumIsolationMVArun2v1DBoldDMwLT        .push_back(boostedTau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT"				    ));
      nBranches_->tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT    .push_back(boostedTau.tauID("byMediumIsolationMVArun2v1PWdR03oldDMwLT"				    ));
      nBranches_->tau_byMediumIsolationMVArun2v1PWnewDMwLT        .push_back(boostedTau.tauID("byMediumIsolationMVArun2v1PWnewDMwLT"				    ));
      nBranches_->tau_byMediumIsolationMVArun2v1PWoldDMwLT        .push_back(boostedTau.tauID("byMediumIsolationMVArun2v1PWoldDMwLT"				    ));

     
      nBranches_->tau_byTightIsolationMVArun2v1DBdR03oldDMwLT     .push_back(boostedTau.tauID("byTightIsolationMVArun2v1DBdR03oldDMwLT"				    ));
      nBranches_->tau_byTightIsolationMVArun2v1DBnewDMwLT         .push_back(boostedTau.tauID("byTightIsolationMVArun2v1DBnewDMwLT"			    ));
      nBranches_->tau_byTightIsolationMVArun2v1DBoldDMwLT         .push_back(boostedTau.tauID("byTightIsolationMVArun2v1DBoldDMwLT"				    ));
      nBranches_->tau_byTightIsolationMVArun2v1PWdR03oldDMwLT     .push_back(boostedTau.tauID("byTightIsolationMVArun2v1PWdR03oldDMwLT"				    ));
      nBranches_->tau_byTightIsolationMVArun2v1PWnewDMwLT         .push_back(boostedTau.tauID("byTightIsolationMVArun2v1PWnewDMwLT"			    ));
      nBranches_->tau_byTightIsolationMVArun2v1PWoldDMwLT         .push_back(boostedTau.tauID("byTightIsolationMVArun2v1PWoldDMwLT"			    ));
      nBranches_->tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT    .push_back(boostedTau.tauID("byVLooseIsolationMVArun2v1DBdR03oldDMwLT"				    ));

      nBranches_->tau_byVLooseIsolationMVArun2v1DBnewDMwLT        .push_back(	boostedTau.tauID("byVLooseIsolationMVArun2v1DBnewDMwLT"			    ));
      nBranches_->tau_byVLooseIsolationMVArun2v1DBoldDMwLT        .push_back(boostedTau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"			    ));
      nBranches_->tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT    .push_back(boostedTau.tauID("byVLooseIsolationMVArun2v1PWdR03oldDMwLT"			    ));
      nBranches_->tau_byVLooseIsolationMVArun2v1PWnewDMwLT        .push_back(boostedTau.tauID("byVLooseIsolationMVArun2v1PWnewDMwLT"			    ));
      nBranches_->tau_byVLooseIsolationMVArun2v1PWoldDMwLT        .push_back(boostedTau.tauID("byVLooseIsolationMVArun2v1PWoldDMwLT"			    ));
      nBranches_->tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT    .push_back(boostedTau.tauID("byVTightIsolationMVArun2v1DBdR03oldDMwLT"			    ));
      nBranches_->tau_byVTightIsolationMVArun2v1DBnewDMwLT        .push_back(boostedTau.tauID("byVTightIsolationMVArun2v1DBnewDMwLT"			    ));
      nBranches_->tau_byVTightIsolationMVArun2v1DBoldDMwLT        .push_back(boostedTau.tauID("byVTightIsolationMVArun2v1DBoldDMwLT"			    ));

      nBranches_->tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT     .push_back(boostedTau.tauID("byVTightIsolationMVArun2v1PWdR03oldDMwLT"			    ));
      nBranches_->tau_byVTightIsolationMVArun2v1PWnewDMwLT         .push_back(boostedTau.tauID("byVTightIsolationMVArun2v1PWnewDMwLT"			    ));
      nBranches_->tau_byVTightIsolationMVArun2v1PWoldDMwLT         .push_back(boostedTau.tauID("byVTightIsolationMVArun2v1PWoldDMwLT"			    ));
      nBranches_->tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT    .push_back(boostedTau.tauID("byVVTightIsolationMVArun2v1DBdR03oldDMwLT"			    ));
      nBranches_->tau_byVVTightIsolationMVArun2v1DBnewDMwLT        .push_back(boostedTau.tauID("byVVTightIsolationMVArun2v1DBnewDMwLT"			    ));
      nBranches_->tau_byVVTightIsolationMVArun2v1DBoldDMwLT        .push_back(boostedTau.tauID("byVVTightIsolationMVArun2v1DBoldDMwLT"			    ));
      nBranches_->tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT    .push_back(boostedTau.tauID("byVVTightIsolationMVArun2v1PWdR03oldDMwLT"			    ));
      nBranches_->tau_byVVTightIsolationMVArun2v1PWnewDMwLT        .push_back(boostedTau.tauID("byVVTightIsolationMVArun2v1PWnewDMwLT"			    ));

      nBranches_->tau_byVVTightIsolationMVArun2v1PWoldDMwLT        .push_back(boostedTau.tauID("byVVTightIsolationMVArun2v1PWoldDMwLT"			    ));




      nBranches_->tau_againstElectronMVA6raw                     .push_back(boostedTau.tauID("againstElectronMVA6Raw"  		    ));
      nBranches_->tau_againstElectronMVA6category                .push_back(boostedTau.tauID("againstElectronMVA6category"		    ));
      nBranches_->tau_againstElectronVLooseMVA6                  .push_back(boostedTau.tauID("againstElectronVLooseMVA6"		    ));
      nBranches_->tau_againstElectronLooseMVA6                   .push_back(boostedTau.tauID("againstElectronLooseMVA6"		    ));
      nBranches_->tau_againstElectronMediumMVA6                  .push_back(boostedTau.tauID("againstElectronMediumMVA6"		    ));
      nBranches_->tau_againstElectronTightMVA6                   .push_back(boostedTau.tauID("againstElectronTightMVA6"		    ));
      nBranches_->tau_againstElectronVTightMVA6                  .push_back(boostedTau.tauID("againstElectronVTightMVA6"		    ));
      

      nBranches_->tau_againstMuonLoose3                          .push_back(boostedTau.tauID("againstMuonLoose3"			    ));
      nBranches_->tau_againstMuonTight3                          .push_back(boostedTau.tauID("againstMuonTight3"			    ));
      
      // nBranches_->tau_byPileupWeightedIsolationRaw3Hits          .push_back(boostedTau.tauID("byPileupWeightedIsolationRaw3Hits"           ));
      // nBranches_->tau_byLoosePileupWeightedIsolation3Hits        .push_back(boostedTau.tauID("byLoosePileupWeightedIsolation3Hits"         ));
      // nBranches_->tau_byMediumPileupWeightedIsolation3Hits       .push_back(boostedTau.tauID("byMediumPileupWeightedIsolation3Hits"        ));
      // nBranches_->tau_byTightPileupWeightedIsolation3Hits        .push_back(boostedTau.tauID("byTightPileupWeightedIsolation3Hits"         ));
      nBranches_->tau_byPhotonPtSumOutsideSignalCone             .push_back(boostedTau.tauID("byPhotonPtSumOutsideSignalCone"              ));
      nBranches_->tau_footprintCorrection                        .push_back(boostedTau.tauID("footprintCorrection"                         ));

   /*======================================================*/   
        
  }        
  nBranches_->tau_N =  taus_->size() + boostedTaus_->size();
  
}  
