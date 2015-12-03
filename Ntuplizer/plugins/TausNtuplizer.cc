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
      // nBranches_->tau_byLooseCombinedIsolationDeltaBetaCorr3Hits .push_back(tau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"  ));
      // nBranches_->tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.push_back(tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits" ));
      // nBranches_->tau_byTightCombinedIsolationDeltaBetaCorr3Hits .push_back(tau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"  ));
      // nBranches_->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits   .push_back(tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"    ));
      // nBranches_->tau_chargedIsoPtSum                            .push_back(tau.tauID("chargedIsoPtSum" 			    ));
      // nBranches_->tau_neutralIsoPtSum                            .push_back(tau.tauID("neutralIsoPtSum" 			    ));
      // nBranches_->tau_puCorrPtSum                                .push_back(tau.tauID("puCorrPtSum"				    ));
      // nBranches_->tau_byIsolationMVA3oldDMwoLTraw                .push_back(tau.tauID("byIsolationMVA3oldDMwoLTraw"		    )); 
      // nBranches_->tau_byVLooseIsolationMVA3oldDMwoLT             .push_back(tau.tauID("byVLooseIsolationMVA3oldDMwoLT"  	    ));
      // nBranches_->tau_byLooseIsolationMVA3oldDMwoLT              .push_back(tau.tauID("byLooseIsolationMVA3oldDMwoLT"		    ));
      // nBranches_->tau_byMediumIsolationMVA3oldDMwoLT             .push_back(tau.tauID("byMediumIsolationMVA3oldDMwoLT"  	    ));
      // nBranches_->tau_byTightIsolationMVA3oldDMwoLT              .push_back(tau.tauID("byTightIsolationMVA3oldDMwoLT"		    ));
      // nBranches_->tau_byVTightIsolationMVA3oldDMwoLT             .push_back(tau.tauID("byVTightIsolationMVA3oldDMwoLT"  	    ));
      // nBranches_->tau_byVVTightIsolationMVA3oldDMwoLT            .push_back(tau.tauID("byVVTightIsolationMVA3oldDMwoLT" 	    ));
      // nBranches_->tau_byIsolationMVA3oldDMwLTraw                 .push_back(tau.tauID("byIsolationMVA3oldDMwLTraw"		    ));
      // nBranches_->tau_byVLooseIsolationMVA3oldDMwLT              .push_back(tau.tauID("byVLooseIsolationMVA3oldDMwLT"		    ));
      // nBranches_->tau_byLooseIsolationMVA3oldDMwLT               .push_back(tau.tauID("byLooseIsolationMVA3oldDMwLT"		    ));
      // nBranches_->tau_byMediumIsolationMVA3oldDMwLT              .push_back(tau.tauID("byMediumIsolationMVA3oldDMwLT"		    ));
      // nBranches_->tau_byTightIsolationMVA3oldDMwLT               .push_back(tau.tauID("byTightIsolationMVA3oldDMwLT"		    ));
      // nBranches_->tau_byVTightIsolationMVA3oldDMwLT              .push_back(tau.tauID("byVTightIsolationMVA3oldDMwLT"		    ));
      // nBranches_->tau_byVVTightIsolationMVA3oldDMwLT             .push_back(tau.tauID("byVVTightIsolationMVA3oldDMwLT"  	    ));
      // nBranches_->tau_byIsolationMVA3newDMwoLTraw                .push_back(tau.tauID("byIsolationMVA3newDMwoLTraw"		    ));
      // nBranches_->tau_byVLooseIsolationMVA3newDMwoLT             .push_back(tau.tauID("byVLooseIsolationMVA3newDMwoLT"  	    ));
      // nBranches_->tau_byLooseIsolationMVA3newDMwoLT              .push_back(tau.tauID("byLooseIsolationMVA3newDMwoLT"		    ));
      // nBranches_->tau_byMediumIsolationMVA3newDMwoLT             .push_back(tau.tauID("byMediumIsolationMVA3newDMwoLT"  	    ));
      // nBranches_->tau_byTightIsolationMVA3newDMwoLT              .push_back(tau.tauID("byTightIsolationMVA3newDMwoLT"		    ));
      // nBranches_->tau_byVTightIsolationMVA3newDMwoLT             .push_back(tau.tauID("byVTightIsolationMVA3newDMwoLT"  	    ));
      // nBranches_->tau_byVVTightIsolationMVA3newDMwoLT            .push_back(tau.tauID("byVVTightIsolationMVA3newDMwoLT" 	    ));
      // nBranches_->tau_byIsolationMVA3newDMwLTraw                 .push_back(tau.tauID("byIsolationMVA3newDMwLTraw"		    ));
      // nBranches_->tau_byVLooseIsolationMVA3newDMwLT              .push_back(tau.tauID("byVLooseIsolationMVA3newDMwLT"		    ));
      // nBranches_->tau_byLooseIsolationMVA3newDMwLT               .push_back(tau.tauID("byLooseIsolationMVA3newDMwLT"		    ));
      // nBranches_->tau_byMediumIsolationMVA3newDMwLT              .push_back(tau.tauID("byMediumIsolationMVA3newDMwLT"		    ));
      // nBranches_->tau_byTightIsolationMVA3newDMwLT               .push_back(tau.tauID("byTightIsolationMVA3newDMwLT"		    ));
      // nBranches_->tau_byVTightIsolationMVA3newDMwLT              .push_back(tau.tauID("byVTightIsolationMVA3newDMwLT"		    ));
      // nBranches_->tau_byVVTightIsolationMVA3newDMwLT             .push_back(tau.tauID("byVVTightIsolationMVA3newDMwLT"  	    ));
      // nBranches_->tau_againstElectronLoose                       .push_back(tau.tauID("againstElectronLoose"			    ));
      // nBranches_->tau_againstElectronMedium                      .push_back(tau.tauID("againstElectronMedium"			    ));
      // nBranches_->tau_againstElectronTight                       .push_back(tau.tauID("againstElectronTight"			    ));
      // nBranches_->tau_againstElectronMVA5raw                     .push_back(tau.tauID("againstElectronMVA5raw"  		    ));
      // nBranches_->tau_againstElectronMVA5category                .push_back(tau.tauID("againstElectronMVA5category"		    ));
      // nBranches_->tau_againstElectronVLooseMVA5                  .push_back(tau.tauID("againstElectronVLooseMVA5"		    ));
      // nBranches_->tau_againstElectronLooseMVA5                   .push_back(tau.tauID("againstElectronLooseMVA5"		    ));
      // nBranches_->tau_againstElectronMediumMVA5                  .push_back(tau.tauID("againstElectronMediumMVA5"		    ));
      // nBranches_->tau_againstElectronTightMVA5                   .push_back(tau.tauID("againstElectronTightMVA5"		    ));
      // nBranches_->tau_againstElectronVTightMVA5                  .push_back(tau.tauID("againstElectronVTightMVA5"		    ));
      // nBranches_->tau_againstMuonLoose                           .push_back(tau.tauID("againstMuonLoose"			    ));
      // nBranches_->tau_againstMuonMedium                          .push_back(tau.tauID("againstMuonMedium"			    ));
      // nBranches_->tau_againstMuonTight                           .push_back(tau.tauID("againstMuonTight"			    ));
      // nBranches_->tau_againstMuonLoose2                          .push_back(tau.tauID("againstMuonLoose2"			    ));
      // nBranches_->tau_againstMuonMedium2                         .push_back(tau.tauID("againstMuonMedium2"			    ));
      // nBranches_->tau_againstMuonTight2                          .push_back(tau.tauID("againstMuonTight2"			    ));
      // nBranches_->tau_againstMuonLoose3                          .push_back(tau.tauID("againstMuonLoose3"			    ));
      // nBranches_->tau_againstMuonTight3                          .push_back(tau.tauID("againstMuonTight3"			    ));
      // nBranches_->tau_againstMuonMVAraw                          .push_back(tau.tauID("againstMuonMVAraw"			    ));
      // nBranches_->tau_againstMuonLooseMVA                        .push_back(tau.tauID("againstMuonLooseMVA"			    ));
      // nBranches_->tau_againstMuonMediumMVA                       .push_back(tau.tauID("againstMuonMediumMVA"			    ));
      // nBranches_->tau_againstMuonTightMVA                        .push_back(tau.tauID("againstMuonTightMVA"			    ));  
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
      nBranches_->tau_decayModeFindingNewDMs                     .push_back(boostedTau.tauID("decayModeFindingNewDMs"		       ));
      nBranches_->tau_decayModeFinding	                     .push_back(boostedTau.tauID("decayModeFinding"			       ));
      // nBranches_->tau_byLooseCombinedIsolationDeltaBetaCorr3Hits .push_back(boostedTau.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"  ));
      // nBranches_->tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.push_back(boostedTau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits" ));
      // nBranches_->tau_byTightCombinedIsolationDeltaBetaCorr3Hits .push_back(boostedTau.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"  ));
      // nBranches_->tau_byCombinedIsolationDeltaBetaCorrRaw3Hits   .push_back(boostedTau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"    ));
      // nBranches_->tau_chargedIsoPtSum                            .push_back(boostedTau.tauID("chargedIsoPtSum"			       ));
      // nBranches_->tau_neutralIsoPtSum                            .push_back(boostedTau.tauID("neutralIsoPtSum"			       ));
      // nBranches_->tau_puCorrPtSum                                .push_back(boostedTau.tauID("puCorrPtSum"  			       ));
      // nBranches_->tau_byIsolationMVA3oldDMwoLTraw                .push_back(boostedTau.tauID("byIsolationMVA3oldDMwoLTraw"  	       )); 
      // nBranches_->tau_byVLooseIsolationMVA3oldDMwoLT             .push_back(boostedTau.tauID("byVLooseIsolationMVA3oldDMwoLT"	       ));
      // nBranches_->tau_byLooseIsolationMVA3oldDMwoLT              .push_back(boostedTau.tauID("byLooseIsolationMVA3oldDMwoLT"	       ));
      // nBranches_->tau_byMediumIsolationMVA3oldDMwoLT             .push_back(boostedTau.tauID("byMediumIsolationMVA3oldDMwoLT"	       ));
      // nBranches_->tau_byTightIsolationMVA3oldDMwoLT              .push_back(boostedTau.tauID("byTightIsolationMVA3oldDMwoLT"	       ));
      // nBranches_->tau_byVTightIsolationMVA3oldDMwoLT             .push_back(boostedTau.tauID("byVTightIsolationMVA3oldDMwoLT"	       ));
      // nBranches_->tau_byVVTightIsolationMVA3oldDMwoLT            .push_back(boostedTau.tauID("byVVTightIsolationMVA3oldDMwoLT"	       ));
      // nBranches_->tau_byIsolationMVA3oldDMwLTraw                 .push_back(boostedTau.tauID("byIsolationMVA3oldDMwLTraw"		       ));
      // nBranches_->tau_byVLooseIsolationMVA3oldDMwLT              .push_back(boostedTau.tauID("byVLooseIsolationMVA3oldDMwLT"	       ));
      // nBranches_->tau_byLooseIsolationMVA3oldDMwLT               .push_back(boostedTau.tauID("byLooseIsolationMVA3oldDMwLT" 	       ));
      // nBranches_->tau_byMediumIsolationMVA3oldDMwLT              .push_back(boostedTau.tauID("byMediumIsolationMVA3oldDMwLT"	       ));
      // nBranches_->tau_byTightIsolationMVA3oldDMwLT               .push_back(boostedTau.tauID("byTightIsolationMVA3oldDMwLT" 	       ));
      // nBranches_->tau_byVTightIsolationMVA3oldDMwLT              .push_back(boostedTau.tauID("byVTightIsolationMVA3oldDMwLT"	       ));
      // nBranches_->tau_byVVTightIsolationMVA3oldDMwLT             .push_back(boostedTau.tauID("byVVTightIsolationMVA3oldDMwLT"	       ));
      // nBranches_->tau_byIsolationMVA3newDMwoLTraw                .push_back(boostedTau.tauID("byIsolationMVA3newDMwoLTraw"  	       ));
      // nBranches_->tau_byVLooseIsolationMVA3newDMwoLT             .push_back(boostedTau.tauID("byVLooseIsolationMVA3newDMwoLT"	       ));
      // nBranches_->tau_byLooseIsolationMVA3newDMwoLT              .push_back(boostedTau.tauID("byLooseIsolationMVA3newDMwoLT"	       ));
      // nBranches_->tau_byMediumIsolationMVA3newDMwoLT             .push_back(boostedTau.tauID("byMediumIsolationMVA3newDMwoLT"	       ));
      // nBranches_->tau_byTightIsolationMVA3newDMwoLT              .push_back(boostedTau.tauID("byTightIsolationMVA3newDMwoLT"	       ));
      // nBranches_->tau_byVTightIsolationMVA3newDMwoLT             .push_back(boostedTau.tauID("byVTightIsolationMVA3newDMwoLT"	       ));
      // nBranches_->tau_byVVTightIsolationMVA3newDMwoLT            .push_back(boostedTau.tauID("byVVTightIsolationMVA3newDMwoLT"	       ));
      // nBranches_->tau_byIsolationMVA3newDMwLTraw                 .push_back(boostedTau.tauID("byIsolationMVA3newDMwLTraw"		       ));
      // nBranches_->tau_byVLooseIsolationMVA3newDMwLT              .push_back(boostedTau.tauID("byVLooseIsolationMVA3newDMwLT"	       ));
      // nBranches_->tau_byLooseIsolationMVA3newDMwLT               .push_back(boostedTau.tauID("byLooseIsolationMVA3newDMwLT" 	       ));
      // nBranches_->tau_byMediumIsolationMVA3newDMwLT              .push_back(boostedTau.tauID("byMediumIsolationMVA3newDMwLT"	       ));
      // nBranches_->tau_byTightIsolationMVA3newDMwLT               .push_back(boostedTau.tauID("byTightIsolationMVA3newDMwLT" 	       ));
      // nBranches_->tau_byVTightIsolationMVA3newDMwLT              .push_back(boostedTau.tauID("byVTightIsolationMVA3newDMwLT"	       ));
      // nBranches_->tau_byVVTightIsolationMVA3newDMwLT             .push_back(boostedTau.tauID("byVVTightIsolationMVA3newDMwLT"	       ));
      // nBranches_->tau_againstElectronLoose                       .push_back(boostedTau.tauID("againstElectronLoose" 		       ));
      // nBranches_->tau_againstElectronMedium                      .push_back(boostedTau.tauID("againstElectronMedium"		       ));
      // nBranches_->tau_againstElectronTight                       .push_back(boostedTau.tauID("againstElectronTight" 		       ));
      // nBranches_->tau_againstElectronMVA5raw                     .push_back(boostedTau.tauID("againstElectronMVA5raw"		       ));
      // nBranches_->tau_againstElectronMVA5category                .push_back(boostedTau.tauID("againstElectronMVA5category"  	       ));
      // nBranches_->tau_againstElectronVLooseMVA5                  .push_back(boostedTau.tauID("againstElectronVLooseMVA5"		       ));
      // nBranches_->tau_againstElectronLooseMVA5                   .push_back(boostedTau.tauID("againstElectronLooseMVA5"		       ));
      // nBranches_->tau_againstElectronMediumMVA5                  .push_back(boostedTau.tauID("againstElectronMediumMVA5"		       ));
      // nBranches_->tau_againstElectronTightMVA5                   .push_back(boostedTau.tauID("againstElectronTightMVA5"		       ));
      // nBranches_->tau_againstElectronVTightMVA5                  .push_back(boostedTau.tauID("againstElectronVTightMVA5"		     ));				
      // nBranches_->tau_againstMuonLoose                           .push_back(boostedTau.tauID("againstMuonLoose"			       ));
      // nBranches_->tau_againstMuonMedium                          .push_back(boostedTau.tauID("againstMuonMedium"			       ));
      // nBranches_->tau_againstMuonTight                           .push_back(boostedTau.tauID("againstMuonTight"			       ));
      // nBranches_->tau_againstMuonLoose2                          .push_back(boostedTau.tauID("againstMuonLoose2"			       ));
      // nBranches_->tau_againstMuonMedium2                         .push_back(boostedTau.tauID("againstMuonMedium2"			       ));
      // nBranches_->tau_againstMuonTight2                          .push_back(boostedTau.tauID("againstMuonTight2"			       ));
      // nBranches_->tau_againstMuonLoose3                          .push_back(boostedTau.tauID("againstMuonLoose3"			       ));
      // nBranches_->tau_againstMuonTight3                          .push_back(boostedTau.tauID("againstMuonTight3"			       ));
      // nBranches_->tau_againstMuonMVAraw                          .push_back(boostedTau.tauID("againstMuonMVAraw"			       ));
      // nBranches_->tau_againstMuonLooseMVA                        .push_back(boostedTau.tauID("againstMuonLooseMVA"  		       ));
      // nBranches_->tau_againstMuonMediumMVA                       .push_back(boostedTau.tauID("againstMuonMediumMVA" 		       ));
      // nBranches_->tau_againstMuonTightMVA                        .push_back(boostedTau.tauID("againstMuonTightMVA"  		       ));   
      /*======================================================*/                          
      
 
  
        
  }        
  nBranches_->tau_N =  taus_->size() + boostedTaus_->size();
  
}  
