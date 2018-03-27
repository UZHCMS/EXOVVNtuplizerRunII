#include "../interface/TausNtuplizer.h"

#include <cmath>

#include "TVector3.h"
#include "TLorentzVector.h"



//TauIdMVAAuxiliaries clusterVariables_;

//===================================================================================================================        
TausNtuplizer::TausNtuplizer( edm::EDGetTokenT<pat::TauCollection> tauToken,edm::EDGetTokenT<pat::TauCollection> tauBoostedTauToken,
			      edm::EDGetTokenT<double> rhoToken, 
			      edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken,
			      edm::EDGetTokenT<reco::VertexCollection> verticeToken,
			      NtupleBranches* nBranches,
			      std::map< std::string, bool >& runFlags )
  : CandidateNtuplizer( nBranches )
  , tauInputToken_ (tauToken)
  , tauBoostedTauInputToken_ ( tauBoostedTauToken)
  , rhoToken_	       	    ( rhoToken  )
  , packedpfcandidatesToken_(packedpfcandidatesToken)
  , verticeToken_     	    ( verticeToken  )	 
  , doBoostedTaus_     	    ( runFlags["doBoostedTaus"]  )	 
 
{
 
}

//===================================================================================================================
TausNtuplizer::~TausNtuplizer( void )
{

}

//===================================================================================================================
void TausNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup // , std::map< std::string, bool >& runFlags 
				  ){
 
  bool doMultipleTauMVA_versions =true;

  if ( !doBoostedTaus_ ) return;
  event.getByToken( tauBoostedTauInputToken_ , boostedTaus_  ); 
  event.getByToken( rhoToken_	 	   , rho_      );
  event.getByToken( verticeToken_ 	   , vertices_ );
  event.getByToken(  tauInputToken_ , taus_ ); 
  event.getByToken( packedpfcandidatesToken_	 	   , packedpfcandidates_      );

//  int n_total = 0;
//  for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){
//    pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
//    
//    float pf_pt = pf.pt();
//    if(pf_pt < 0.5) continue;
//
//    n_total ++;
//
//  }
//
//  nBranches_->tau_n_total.push_back(n_total);

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


      // YT added : 2 Aug 2017

      int n_ch = 0;
      int n_nh = 0;
     //      int n_h_f = 0;
      //      int n_em_f = 0;
      int n_gamma = 0;
      //      int n_e = 0;
      //      int n_mu = 0;

//      std::vector<int> associated_pdgId;
//      std::vector<float> associated_pt;
//      std::vector<float> associated_eta;
//      std::vector<float> associated_phi;
//      std::vector<float> associated_dr;

//      associated_pdgId.clear();
//      associated_pt.clear();
//      associated_eta.clear();
//      associated_phi.clear();
//      associated_dr.clear();
      
      for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){
	pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
    
	float pf_pt = pf.pt();
	if(pf_pt < 0.5) continue;

	float pf_phi = pf.phi();
	float pf_eta = pf.eta();

	float res = tau.phi() - pf_phi;
	while(res > TMath::Pi()){
	  res -= 2*TMath::Pi();
	}
	while(res < -TMath::Pi()){
	  res += 2*TMath::Pi();
	}


	float deta = tau.eta() - pf_eta;
	float deltaR = TMath::Sqrt(TMath::Power(deta,2) + TMath::Power(res,2));

	if(deltaR > 0.5) continue;

	
	int pf_pdgId = abs(pf.pdgId());
	
	if(pf_pdgId==211) n_ch++;
	if(pf_pdgId==130) n_nh++;
	if(pf_pdgId==22) n_gamma++;
	//	if(pf_pdgId==11) n_e++;
	//	if(pf_pdgId==13) n_mu++;
	//	if(pf_pdgId==1) n_h_f++;
	//	if(pf_pdgId==2) n_em_f++;

//	associated_pdgId.push_back(pf_pdgId);
//	associated_pt.push_back(pf_pt);
//	associated_eta.push_back(pf_eta);
//	associated_phi.push_back(pf_phi);
//	associated_dr.push_back(deltaR);
	
	//    std::cout << "YUTA: index=" << ii << ", pt=" << pf.pt() << ", pdg ID=" << pf.pdgId() << std::endl;
      }
      
      //      std::cout << "YUTA: " << n_ch << " " << n_nh << " " << n_gamma << " " << n_e << " " << n_mu << " " << n_h_f << " " << n_em_f << std::endl;

//      nBranches_->tau_associated_pdgId                .push_back(associated_pdgId);
//      nBranches_->tau_associated_pt                .push_back(associated_pt);
//      nBranches_->tau_associated_eta                .push_back(associated_eta);
//      nBranches_->tau_associated_phi                .push_back(associated_phi);
//      nBranches_->tau_associated_dr                .push_back(associated_dr);

      nBranches_->tau_n_ch.push_back(n_ch);
      nBranches_->tau_n_nh.push_back(n_nh);
      //      nBranches_->tau_n_h_f.push_back(n_h_f);
      //      nBranches_->tau_n_em_f.push_back(n_em_f);
      nBranches_->tau_n_gamma.push_back(n_gamma);
      //      nBranches_->tau_n_e.push_back(n_e);
      //      nBranches_->tau_n_mu.push_back(n_mu);


      // YT added : 17 Aug 2016
      pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
      nBranches_->tau_dz	      	      .push_back(fabs(packedLeadTauCand->dz()));


      int tauDecayMode = tau.decayMode();
      float decayDistX = tau.flightLength().x();
      float decayDistY = tau.flightLength().y();
      float decayDistZ = tau.flightLength().z();
      float decayDistMag = std::sqrt(decayDistX*decayDistX + decayDistY*decayDistY + decayDistZ*decayDistZ);
      float nPhoton = (float)clusterVariables_.tau_n_photons_total(tau);
      float ptWeightedDetaStrip = clusterVariables_.tau_pt_weighted_deta_strip(tau, tauDecayMode);
      float ptWeightedDphiStrip = clusterVariables_.tau_pt_weighted_dphi_strip(tau, tauDecayMode);
      float ptWeightedDrSignal = clusterVariables_.tau_pt_weighted_dr_signal(tau, tauDecayMode);
      float ptWeightedDrIsolation = clusterVariables_.tau_pt_weighted_dr_iso(tau, tauDecayMode);
      float leadingTrackChi2 = tau.leadingTrackNormChi2();
      float leadingTrackPt = packedLeadTauCand->pt();
      float eRatio = clusterVariables_.tau_Eratio(tau);
      float dxy_Sig = tau.dxy_Sig();
      float ip3d = tau.ip3d();
      float ip3d_Sig = tau.ip3d_Sig();
      bool hasSecondaryVertex = tau.hasSecondaryVertex();
      float flightLenthSig = tau.flightLengthSig();


      nBranches_->tau_nPhoton		     .push_back(nPhoton);
      nBranches_->tau_ptWeightedDetaStrip    .push_back(ptWeightedDetaStrip);
      nBranches_->tau_ptWeightedDphiStrip    .push_back(ptWeightedDphiStrip);
      nBranches_->tau_ptWeightedDrSignal     .push_back(ptWeightedDrSignal);
      nBranches_->tau_ptWeightedDrIsolation  .push_back(ptWeightedDrIsolation);
      nBranches_->tau_leadingTrackChi2       .push_back(leadingTrackChi2);
      nBranches_->tau_leadingTrackPt         .push_back(leadingTrackPt);
      nBranches_->tau_eRatio                 .push_back(eRatio);
      nBranches_->tau_dxy_Sig                .push_back(dxy_Sig);
      nBranches_->tau_ip3d                   .push_back(ip3d);
      nBranches_->tau_ip3d_Sig               .push_back(ip3d_Sig);
      nBranches_->tau_hasSecondaryVertex     .push_back(hasSecondaryVertex);
      //      nBranches_->tau_decayDistMag_x         .push_back(decayDistX);
      //      nBranches_->tau_decayDistMag_y         .push_back(decayDistY);
      //      nBranches_->tau_decayDistMag_z         .push_back(decayDistZ);
      nBranches_->tau_decayDistMag           .push_back(decayDistMag);
      nBranches_->tau_flightLenthSig         .push_back(flightLenthSig);
      

//      std::cout << "DBG1 : " << ptWeightedDetaStrip << " " << ptWeightedDphiStrip <<  " " << ptWeightedDrSignal << " " << ptWeightedDrIsolation << " "  << eRatio << std::endl;
//      std::cout << "DBG2: " << nPhoton << " " << dxy_Sig << " " << ip3d<< " " << ip3d_Sig << " " << hasSecondaryVertex << " " << flightLenthSig << std::endl;
//      std::cout << "DBG3 : " << decayDistX << " " << decayDistY << " " << decayDistZ << " " << decayDistMag << " " << leadingTrackChi2 << std::endl;

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
      //double rho = *(rho_.product());
      // float  deltaR = 0.3;
      // double energy = TMath::Pi()*deltaR*deltaR*rho;      
      // nBranches_->tau_pfRhoCorrRelIso03.push_back((tau.chargedHadronIso() + std::max(0.,tau.neutralHadronIso() + tau.photonIso() - energy))/tau.pt());      
      // nBranches_->tau_pfRhoCorrRelIso03Boost.push_back((tau.userIsolation(pat::PfChargedHadronIso) + std::max(0.,tau.userIsolation(pat::PfNeutralHadronIso) + tau.userIsolation(pat::PfGammaIso) - energy))/tau.pt()); 

      // deltaR = 0.4;
      // energy = TMath::Pi()*deltaR*deltaR*rho;      
      // nBranches_->tau_pfRhoCorrRelIso04.push_back((tau.chargedHadronIso() + std::max(0.,tau.neutralHadronIso() + tau.photonIso() - energy))/tau.pt());
      // nBranches_->tau_pfRhoCorrRelIso04Boost.push_back((tau.userIsolation(pat::PfChargedHadronIso) + std::max(0.,tau.userIsolation(pat::PfNeutralHadronIso) + tau.userIsolation(pat::PfGammaIso) - energy))/tau.pt()); 
            
      // nBranches_->tau_pfDeltaCorrRelIso     .push_back((tau.chargedHadronIso() + std::max(0., tau.neutralHadronIso() + tau.photonIso() - 0.5*tau.puChargedHadronIso()))/tau.pt());
      // nBranches_->tau_pfRelIso              .push_back((tau.chargedHadronIso() + tau.neutralHadronIso()+ tau.photonIso())/tau.pt()) ; 
      // nBranches_->tau_photonIso             .push_back(tau.photonIso());
      // nBranches_->tau_neutralHadIso         .push_back(tau.neutralHadronIso());
      // nBranches_->tau_chargedHadIso         .push_back(tau.chargedHadronIso());

      nBranches_->tau_photonPtSumOutsideSignalCone               .push_back(tau.tauID("photonPtSumOutsideSignalCone"));
      // nBranches_->tau_trackIso              .push_back(tau.trackIso());      
      // nBranches_->tau_pfDeltaCorrRelIsoBoost.push_back((tau.userIsolation(pat::PfChargedHadronIso) + std::max(0., tau.userIsolation(pat::PfNeutralHadronIso) + tau.userIsolation(pat::PfGammaIso) - 0.5*tau.userIsolation(pat::PfPUChargedHadronIso)))/tau.pt());
      // nBranches_->tau_pfRelIsoBoost         .push_back((tau.userIsolation(pat::PfChargedHadronIso) + tau.userIsolation(pat::PfNeutralHadronIso)+ tau.userIsolation(pat::PfGammaIso))/tau.pt()) ; 
      // nBranches_->tau_photonIsoBoost        .push_back(tau.userIsolation(pat::PfGammaIso));
      // nBranches_->tau_neutralHadIsoBoost    .push_back(tau.userIsolation(pat::PfNeutralHadronIso));
      // nBranches_->tau_chargedHadIsoBoost    .push_back(tau.userIsolation(pat::PfChargedHadronIso));
      // std::cout << " tau.userIsolation(pat::PfChargedHadronIso) "<< tau.userIsolation(pat::PfChargedHadronIso) <<" tau.userIsolation(pat::PfNeutralHadronIso) " << tau.userIsolation(pat::PfNeutralHadronIso) << " tau.userIsolation(pat::PfGammaIso) "<< std::endl;
 
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
      //      nBranches_->tau_footprintCorrectiondR03                    .push_back(tau.tauID("footprintCorrectiondR03"                         ));
      nBranches_->tau_neutralIsoPtSumdR03                        .push_back(tau.tauID("neutralIsoPtSumdR03" 			    ));
      //      nBranches_->tau_neutralIsoPtSumWeight                      .push_back(tau.tauID("neutralIsoPtSumWeight" 			    ));
      //      nBranches_->tau_neutralIsoPtSumWeightdR03                  .push_back(tau.tauID("neutralIsoPtSumWeightdR03" 			    ));
      nBranches_->tau_photonPtSumOutsideSignalConedR03           .push_back(tau.tauID("photonPtSumOutsideSignalConedR03"               ));

      nBranches_->tau_byIsolationMVArun2v1DBdR03oldDMwLTraw      .push_back(tau.tauID("byIsolationMVArun2v2DBoldDMdR0p3wLTrawNew"));
      nBranches_->tau_byIsolationMVArun2v1DBnewDMwLTraw          .push_back(tau.tauID("byIsolationMVArun2v2DBnewDMwLTrawNew"));
      nBranches_->tau_byIsolationMVArun2v1DBoldDMwLTraw          .push_back(tau.tauID("byIsolationMVArun2v1DBoldDMwLTrawNew"				    ));
      //      nBranches_->tau_byIsolationMVArun2v1DBoldDMwoLTraw          .push_back(tau.tauID("byIsolationMVArun2v1DBoldDMwoLTraw"				    ));
      //      nBranches_->tau_byIsolationMVArun2v1PWdR03oldDMwLTraw      .push_back(tau.tauID("byIsolationMVArun2v1PWdR03oldDMwLTraw"				    ));
      nBranches_->tau_byIsolationMVArun2v1PWnewDMwLTraw          .push_back(tau.tauID("byIsolationMVArun2v1PWnewDMwLTraw"				    ));
      //      nBranches_->tau_byIsolationMVArun2v1PWoldDMwLTraw          .push_back(tau.tauID("byIsolationMVArun2v1PWoldDMwLTraw"				    ));
      nBranches_->tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT    .push_back(tau.tauID("byLooseIsolationMVArun2v2DBoldDMdR0p3wLTNew"				    ));
      nBranches_->tau_byLooseIsolationMVArun2v1DBnewDMwLT        .push_back(tau.tauID("byLooseIsolationMVArun2v2DBnewDMwLTNew"				    ));
      nBranches_->tau_byLooseIsolationMVArun2v1DBoldDMwLT        .push_back(tau.tauID("byLooseIsolationMVArun2v1DBoldDMwLTNew"			    ));
      //      nBranches_->tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT    .push_back(tau.tauID("byLooseIsolationMVArun2v1PWdR03oldDMwLT"				    ));
      nBranches_->tau_byLooseIsolationMVArun2v1PWnewDMwLT        .push_back(	tau.tauID("byLooseIsolationMVArun2v1PWnewDMwLT"			    ));
      //      nBranches_->tau_byLooseIsolationMVArun2v1PWoldDMwLT        .push_back(	tau.tauID("byLooseIsolationMVArun2v1PWoldDMwLT"			    ));
      
      nBranches_->tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT    .push_back(tau.tauID("byMediumIsolationMVArun2v2DBoldDMdR0p3wLTNew"				    ));
      nBranches_->tau_byMediumIsolationMVArun2v1DBnewDMwLT        .push_back(	tau.tauID("byMediumIsolationMVArun2v2DBnewDMwLTNew"			    ));
      nBranches_->tau_byMediumIsolationMVArun2v1DBoldDMwLT        .push_back(tau.tauID("byMediumIsolationMVArun2v1DBoldDMwLTNew"				    ));
      //      nBranches_->tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT    .push_back(tau.tauID("byMediumIsolationMVArun2v1PWdR03oldDMwLT"				    ));
      nBranches_->tau_byMediumIsolationMVArun2v1PWnewDMwLT        .push_back(tau.tauID("byMediumIsolationMVArun2v1PWnewDMwLT"				    ));
      //      nBranches_->tau_byMediumIsolationMVArun2v1PWoldDMwLT        .push_back(tau.tauID("byMediumIsolationMVArun2v1PWoldDMwLT"				    ));

    
      nBranches_->tau_byTightIsolationMVArun2v1DBdR03oldDMwLT     .push_back(tau.tauID("byTightIsolationMVArun2v2DBoldDMdR0p3wLTNew"				    ));
      nBranches_->tau_byTightIsolationMVArun2v1DBnewDMwLT         .push_back(tau.tauID("byTightIsolationMVArun2v2DBnewDMwLTNew"			    ));
      nBranches_->tau_byTightIsolationMVArun2v1DBoldDMwLT         .push_back(tau.tauID("byTightIsolationMVArun2v1DBoldDMwLTNew"				    ));
      //      nBranches_->tau_byTightIsolationMVArun2v1PWdR03oldDMwLT     .push_back(tau.tauID("byTightIsolationMVArun2v1PWdR03oldDMwLT"				    ));
      nBranches_->tau_byTightIsolationMVArun2v1PWnewDMwLT         .push_back(tau.tauID("byTightIsolationMVArun2v1PWnewDMwLT"			    ));
      //      nBranches_->tau_byTightIsolationMVArun2v1PWoldDMwLT         .push_back(tau.tauID("byTightIsolationMVArun2v1PWoldDMwLT"			    ));
      nBranches_->tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT    .push_back(tau.tauID("byVLooseIsolationMVArun2v2DBoldDMdR0p3wLTNew"				    ));

      nBranches_->tau_byVLooseIsolationMVArun2v1DBnewDMwLT        .push_back(	tau.tauID("byVLooseIsolationMVArun2v2DBnewDMwLTNew"			    ));
      nBranches_->tau_byVLooseIsolationMVArun2v1DBoldDMwLT        .push_back(tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLTNew"			    ));
      nBranches_->tau_byVVLooseIsolationMVArun2v1DBoldDMwLT       .push_back(tau.tauID("byVVLooseIsolationMVArun2v1DBoldDMwLTNew"                       ));
           
      //     nBranches_->tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT    .push_back(tau.tauID("byVLooseIsolationMVArun2v1PWdR03oldDMwLT"			    ));
      nBranches_->tau_byVLooseIsolationMVArun2v1PWnewDMwLT        .push_back(tau.tauID("byVLooseIsolationMVArun2v1PWnewDMwLT"			    ));
      //      nBranches_->tau_byVLooseIsolationMVArun2v1PWoldDMwLT        .push_back(tau.tauID("byVLooseIsolationMVArun2v1PWoldDMwLT"			    ));
      nBranches_->tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT    .push_back(tau.tauID("byVTightIsolationMVArun2v2DBoldDMdR0p3wLTNew"			    ));
      nBranches_->tau_byVTightIsolationMVArun2v1DBnewDMwLT        .push_back(tau.tauID("byVTightIsolationMVArun2v2DBnewDMwLTNew"			    ));
      nBranches_->tau_byVTightIsolationMVArun2v1DBoldDMwLT        .push_back(tau.tauID("byVTightIsolationMVArun2v1DBoldDMwLTNew"			    ));

      //      nBranches_->tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT     .push_back(tau.tauID("byVTightIsolationMVArun2v1PWdR03oldDMwLT"			    ));
      nBranches_->tau_byVTightIsolationMVArun2v1PWnewDMwLT         .push_back(tau.tauID("byVTightIsolationMVArun2v1PWnewDMwLT"			    ));
      //      nBranches_->tau_byVTightIsolationMVArun2v1PWoldDMwLT         .push_back(tau.tauID("byVTightIsolationMVArun2v1PWoldDMwLT"			    ));
      nBranches_->tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT    .push_back(tau.tauID("byVVTightIsolationMVArun2v2DBoldDMdR0p3wLTNew"			    ));
      nBranches_->tau_byVVTightIsolationMVArun2v1DBnewDMwLT        .push_back(tau.tauID("byVVTightIsolationMVArun2v2DBnewDMwLTNew"			    ));
      nBranches_->tau_byVVTightIsolationMVArun2v1DBoldDMwLT        .push_back(tau.tauID("byVVTightIsolationMVArun2v1DBoldDMwLTNew"			    ));
      //      nBranches_->tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT    .push_back(tau.tauID("byVVTightIsolationMVArun2v1PWdR03oldDMwLT"			    ));
      nBranches_->tau_byVVTightIsolationMVArun2v1PWnewDMwLT        .push_back(tau.tauID("byVVTightIsolationMVArun2v1PWnewDMwLT"			    ));

      //      nBranches_->tau_byVVTightIsolationMVArun2v1PWoldDMwLT        .push_back(tau.tauID("byVVTightIsolationMVArun2v1PWoldDMwLT"			    ));
      // if  (runFlags["doMultipleTauMVA_versions"]){
      nBranches_->tau_byIsolationMVArun2v2DBoldDMwLTraw.push_back(tau.tauID("byIsolationMVArun2v2DBoldDMwLTrawNew"));
      nBranches_->tau_byVVLooseIsolationMVArun2v2DBoldDMwLT.push_back(tau.tauID("byVVLooseIsolationMVArun2v2DBoldDMwLTNew"));
      nBranches_->tau_byVLooseIsolationMVArun2v2DBoldDMwLT.push_back(tau.tauID("byVLooseIsolationMVArun2v2DBoldDMwLTNew"));
      nBranches_->tau_byLooseIsolationMVArun2v2DBoldDMwLT.push_back(tau.tauID("byLooseIsolationMVArun2v2DBoldDMwLTNew"));
      nBranches_->tau_byMediumIsolationMVArun2v2DBoldDMwLT.push_back(tau.tauID("byMediumIsolationMVArun2v2DBoldDMwLTNew"));
      nBranches_->tau_byTightIsolationMVArun2v2DBoldDMwLT.push_back(tau.tauID("byTightIsolationMVArun2v2DBoldDMwLTNew"));
      nBranches_->tau_byVTightIsolationMVArun2v2DBoldDMwLT.push_back(tau.tauID("byVTightIsolationMVArun2v2DBoldDMwLTNew")); 
      nBranches_->tau_byVVTightIsolationMVArun2v2DBoldDMwLT.push_back(tau.tauID("byVVTightIsolationMVArun2v2DBoldDMwLTNew"));

    

      // }

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
      //      nBranches_->tau_footprintCorrection                        .push_back(tau.tauID("footprintCorrection"                         ));

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


      // YT added : 2 Aug 2017

      int n_ch = 0;
      int n_nh = 0;
      //      int n_h_f = 0;
      //      int n_em_f = 0;
      int n_gamma = 0;
      //      int n_e = 0;
      //      int n_mu = 0;

//      std::vector<int> associated_pdgId;
//      std::vector<float> associated_pt;
//      std::vector<float> associated_eta;
//      std::vector<float> associated_phi;
//      std::vector<float> associated_dr;
//
//      associated_pdgId.clear();
//      associated_pt.clear();
//      associated_eta.clear();
//      associated_phi.clear();
//      associated_dr.clear();
      
      for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){
	pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
    
	float pf_pt = pf.pt();
	if(pf_pt < 0.5) continue;

	float pf_phi = pf.phi();
	float pf_eta = pf.eta();

	float res = boostedTau.phi() - pf_phi;
	while(res > TMath::Pi()){
	  res -= 2*TMath::Pi();
	}
	while(res < -TMath::Pi()){
	  res += 2*TMath::Pi();
	}


	float deta = boostedTau.eta() - pf_eta;
	float deltaR = TMath::Sqrt(TMath::Power(deta,2) + TMath::Power(res,2));

	if(deltaR > 0.5) continue;

	
	int pf_pdgId = abs(pf.pdgId());
	
	if(pf_pdgId==211) n_ch++;
	if(pf_pdgId==130) n_nh++;
	if(pf_pdgId==22) n_gamma++;
	//	if(pf_pdgId==11) n_e++;
	//	if(pf_pdgId==13) n_mu++;
	//	if(pf_pdgId==1) n_h_f++;
	//	if(pf_pdgId==2) n_em_f++;

//	associated_pdgId.push_back(pf_pdgId);
//	associated_pt.push_back(pf_pt);
//	associated_eta.push_back(pf_eta);
//	associated_phi.push_back(pf_phi);
//	associated_dr.push_back(deltaR);
	
	//    std::cout << "YUTA: index=" << ii << ", pt=" << pf.pt() << ", pdg ID=" << pf.pdgId() << std::endl;
      }
      
      //      std::cout << "YUTA: " << n_ch << " " << n_nh << " " << n_gamma << " " << n_e << " " << n_mu << " " << n_h_f << " " << n_em_f << std::endl;

//      nBranches_->tau_associated_pdgId                .push_back(associated_pdgId);
//      nBranches_->tau_associated_pt                .push_back(associated_pt);
//      nBranches_->tau_associated_eta                .push_back(associated_eta);
//      nBranches_->tau_associated_phi                .push_back(associated_phi);
//      nBranches_->tau_associated_dr                .push_back(associated_dr);

      nBranches_->tau_n_ch.push_back(n_ch);
      nBranches_->tau_n_nh.push_back(n_nh);
      //      nBranches_->tau_n_h_f.push_back(n_h_f);
      //      nBranches_->tau_n_em_f.push_back(n_em_f);
      nBranches_->tau_n_gamma.push_back(n_gamma);
      //      nBranches_->tau_n_e.push_back(n_e);
      //      nBranches_->tau_n_mu.push_back(n_mu);

      // YT added : 17 Aug 2016
      pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(boostedTau.leadChargedHadrCand().get());
      nBranches_->tau_dz	      	      .push_back(fabs(packedLeadTauCand->dz()));

      int tauDecayMode = boostedTau.decayMode();
      float decayDistX = boostedTau.flightLength().x();
      float decayDistY = boostedTau.flightLength().y();
      float decayDistZ = boostedTau.flightLength().z();
      float decayDistMag = std::sqrt(decayDistX*decayDistX + decayDistY*decayDistY + decayDistZ*decayDistZ);

      float nPhoton = (float)clusterVariables_.tau_n_photons_total(boostedTau);
      float ptWeightedDetaStrip = clusterVariables_.tau_pt_weighted_deta_strip(boostedTau, tauDecayMode);
      float ptWeightedDphiStrip = clusterVariables_.tau_pt_weighted_dphi_strip(boostedTau, tauDecayMode);
      float ptWeightedDrSignal = clusterVariables_.tau_pt_weighted_dr_signal(boostedTau, tauDecayMode);
      float ptWeightedDrIsolation = clusterVariables_.tau_pt_weighted_dr_iso(boostedTau, tauDecayMode);
      float leadingTrackChi2 = boostedTau.leadingTrackNormChi2();
      float leadingTrackPt = packedLeadTauCand->pt();
      float eRatio = clusterVariables_.tau_Eratio(boostedTau);
      float dxy_Sig = boostedTau.dxy_Sig();
      float ip3d = boostedTau.ip3d();
      float ip3d_Sig = boostedTau.ip3d_Sig();
      bool hasSecondaryVertex = boostedTau.hasSecondaryVertex();
      float flightLenthSig = boostedTau.flightLengthSig();


      nBranches_->tau_nPhoton		     .push_back(nPhoton);
      nBranches_->tau_ptWeightedDetaStrip    .push_back(ptWeightedDetaStrip);
      nBranches_->tau_ptWeightedDphiStrip    .push_back(ptWeightedDphiStrip);
      nBranches_->tau_ptWeightedDrSignal     .push_back(ptWeightedDrSignal);
      nBranches_->tau_ptWeightedDrIsolation  .push_back(ptWeightedDrIsolation);
      nBranches_->tau_leadingTrackChi2       .push_back(leadingTrackChi2);
      nBranches_->tau_leadingTrackPt         .push_back(leadingTrackPt);
      nBranches_->tau_eRatio                 .push_back(eRatio);
      nBranches_->tau_dxy_Sig                .push_back(dxy_Sig);
      nBranches_->tau_ip3d                   .push_back(ip3d);
      nBranches_->tau_ip3d_Sig               .push_back(ip3d_Sig);
      nBranches_->tau_hasSecondaryVertex     .push_back(hasSecondaryVertex);
      //      nBranches_->tau_decayDistMag_x         .push_back(decayDistX);
      //      nBranches_->tau_decayDistMag_y         .push_back(decayDistY);
      //      nBranches_->tau_decayDistMag_z         .push_back(decayDistZ);
      nBranches_->tau_decayDistMag           .push_back(decayDistMag);
      nBranches_->tau_flightLenthSig         .push_back(flightLenthSig);


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
      // double rho = *(rho_.product());
      // float  deltaR = 0.3;
      // double energy = TMath::Pi()*deltaR*deltaR*rho;      
      // nBranches_->tau_pfRhoCorrRelIso03.push_back((boostedTau.chargedHadronIso() + std::max(0.,boostedTau.neutralHadronIso() + boostedTau.photonIso() - energy))/boostedTau.pt());      
      // nBranches_->tau_pfRhoCorrRelIso03Boost.push_back((boostedTau.userIsolation(pat::PfChargedHadronIso) + std::max(0.,boostedTau.userIsolation(pat::PfNeutralHadronIso) + boostedTau.userIsolation(pat::PfGammaIso) - energy))/boostedTau.pt()); 

      // deltaR = 0.4;
      // energy = TMath::Pi()*deltaR*deltaR*rho;      
      // nBranches_->tau_pfRhoCorrRelIso04.push_back((boostedTau.chargedHadronIso() + std::max(0.,boostedTau.neutralHadronIso() + boostedTau.photonIso() - energy))/boostedTau.pt());
      // nBranches_->tau_pfRhoCorrRelIso04Boost.push_back((boostedTau.userIsolation(pat::PfChargedHadronIso) + std::max(0.,boostedTau.userIsolation(pat::PfNeutralHadronIso) + boostedTau.userIsolation(pat::PfGammaIso) - energy))/boostedTau.pt()); 
            
      // nBranches_->tau_pfDeltaCorrRelIso     .push_back((boostedTau.chargedHadronIso() + std::max(0., boostedTau.neutralHadronIso() + boostedTau.photonIso() - 0.5*boostedTau.puChargedHadronIso()))/boostedTau.pt());
      // nBranches_->tau_pfRelIso              .push_back((boostedTau.chargedHadronIso() + boostedTau.neutralHadronIso()+ boostedTau.photonIso())/boostedTau.pt()) ; 
      // nBranches_->tau_photonIso             .push_back(boostedTau.photonIso());
      // nBranches_->tau_neutralHadIso         .push_back(boostedTau.neutralHadronIso());
      // nBranches_->tau_chargedHadIso         .push_back(boostedTau.chargedHadronIso());
      nBranches_->tau_photonPtSumOutsideSignalCone               .push_back(boostedTau.tauID("photonPtSumOutsideSignalCone"));
      // nBranches_->tau_trackIso              .push_back(boostedTau.trackIso());            
      // nBranches_->tau_pfDeltaCorrRelIsoBoost.push_back((boostedTau.userIsolation(pat::PfChargedHadronIso) + std::max(0., boostedTau.userIsolation(pat::PfNeutralHadronIso) + boostedTau.userIsolation(pat::PfGammaIso) - 0.5*boostedTau.userIsolation(pat::PfPUChargedHadronIso)))/boostedTau.pt());
      // nBranches_->tau_pfRelIsoBoost         .push_back((boostedTau.userIsolation(pat::PfChargedHadronIso) + boostedTau.userIsolation(pat::PfNeutralHadronIso)+ boostedTau.userIsolation(pat::PfGammaIso))/boostedTau.pt()) ; 
      // nBranches_->tau_photonIsoBoost        .push_back(boostedTau.userIsolation(pat::PfGammaIso));
      // nBranches_->tau_neutralHadIsoBoost    .push_back(boostedTau.userIsolation(pat::PfNeutralHadronIso));
      // nBranches_->tau_chargedHadIsoBoost    .push_back(boostedTau.userIsolation(pat::PfChargedHadronIso));

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
      //      nBranches_->tau_footprintCorrectiondR03                    .push_back(boostedTau.tauID("footprintCorrectiondR03"                         ));
      nBranches_->tau_neutralIsoPtSumdR03                        .push_back(boostedTau.tauID("neutralIsoPtSumdR03" 			    ));
      //      nBranches_->tau_neutralIsoPtSumWeight                      .push_back(boostedTau.tauID("neutralIsoPtSumWeight" 			    ));
      //      nBranches_->tau_neutralIsoPtSumWeightdR03                  .push_back(boostedTau.tauID("neutralIsoPtSumWeightdR03" 			    ));
      nBranches_->tau_photonPtSumOutsideSignalConedR03           .push_back(boostedTau.tauID("photonPtSumOutsideSignalConedR03"               ));

      nBranches_->tau_byIsolationMVArun2v1DBdR03oldDMwLTraw      .push_back(boostedTau.tauID("byIsolationMVArun2v1DBdR03oldDMwLTraw"));
      nBranches_->tau_byIsolationMVArun2v1DBnewDMwLTraw          .push_back(boostedTau.tauID("byIsolationMVArun2v1DBnewDMwLTraw"));
      nBranches_->tau_byIsolationMVArun2v1DBoldDMwLTraw          .push_back(boostedTau.tauID("byIsolationMVArun2v1DBoldDMwLTraw"				    ));
      //      nBranches_->tau_byIsolationMVArun2v1DBoldDMwoLTraw          .push_back(boostedTau.tauID("byIsolationMVArun2v1DBoldDMwoLTraw"				    ));
      //      nBranches_->tau_byIsolationMVArun2v1PWdR03oldDMwLTraw      .push_back(boostedTau.tauID("byIsolationMVArun2v1PWdR03oldDMwLTraw"				    ));
      nBranches_->tau_byIsolationMVArun2v1PWnewDMwLTraw          .push_back(boostedTau.tauID("byIsolationMVArun2v1PWnewDMwLTraw"				    ));
      //      nBranches_->tau_byIsolationMVArun2v1PWoldDMwLTraw          .push_back(boostedTau.tauID("byIsolationMVArun2v1PWoldDMwLTraw"				    ));
      nBranches_->tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT    .push_back(boostedTau.tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT"				    ));
      nBranches_->tau_byLooseIsolationMVArun2v1DBnewDMwLT        .push_back(boostedTau.tauID("byLooseIsolationMVArun2v1DBnewDMwLT"				    ));
      nBranches_->tau_byLooseIsolationMVArun2v1DBoldDMwLT        .push_back(boostedTau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT"			    ));
      //      nBranches_->tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT    .push_back(boostedTau.tauID("byLooseIsolationMVArun2v1PWdR03oldDMwLT"				    ));
      nBranches_->tau_byLooseIsolationMVArun2v1PWnewDMwLT        .push_back(	boostedTau.tauID("byLooseIsolationMVArun2v1PWnewDMwLT"			    ));
      //      nBranches_->tau_byLooseIsolationMVArun2v1PWoldDMwLT        .push_back(	boostedTau.tauID("byLooseIsolationMVArun2v1PWoldDMwLT"			    ));
       nBranches_->tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT    .push_back(boostedTau.tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT"				    ));
      nBranches_->tau_byMediumIsolationMVArun2v1DBnewDMwLT        .push_back(	boostedTau.tauID("byMediumIsolationMVArun2v1DBnewDMwLT"			    ));
      nBranches_->tau_byMediumIsolationMVArun2v1DBoldDMwLT        .push_back(boostedTau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT"				    ));
      //      nBranches_->tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT    .push_back(boostedTau.tauID("byMediumIsolationMVArun2v1PWdR03oldDMwLT"				    ));
      nBranches_->tau_byMediumIsolationMVArun2v1PWnewDMwLT        .push_back(boostedTau.tauID("byMediumIsolationMVArun2v1PWnewDMwLT"				    ));
      //      nBranches_->tau_byMediumIsolationMVArun2v1PWoldDMwLT        .push_back(boostedTau.tauID("byMediumIsolationMVArun2v1PWoldDMwLT"				    ));

     
      nBranches_->tau_byTightIsolationMVArun2v1DBdR03oldDMwLT     .push_back(boostedTau.tauID("byTightIsolationMVArun2v1DBdR03oldDMwLT"				    ));
      nBranches_->tau_byTightIsolationMVArun2v1DBnewDMwLT         .push_back(boostedTau.tauID("byTightIsolationMVArun2v1DBnewDMwLT"			    ));
      nBranches_->tau_byTightIsolationMVArun2v1DBoldDMwLT         .push_back(boostedTau.tauID("byTightIsolationMVArun2v1DBoldDMwLT"				    ));
      //      nBranches_->tau_byTightIsolationMVArun2v1PWdR03oldDMwLT     .push_back(boostedTau.tauID("byTightIsolationMVArun2v1PWdR03oldDMwLT"				    ));
      nBranches_->tau_byTightIsolationMVArun2v1PWnewDMwLT         .push_back(boostedTau.tauID("byTightIsolationMVArun2v1PWnewDMwLT"			    ));
      //      nBranches_->tau_byTightIsolationMVArun2v1PWoldDMwLT         .push_back(boostedTau.tauID("byTightIsolationMVArun2v1PWoldDMwLT"			    ));
      nBranches_->tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT    .push_back(boostedTau.tauID("byVLooseIsolationMVArun2v1DBdR03oldDMwLT"				    ));

      nBranches_->tau_byVLooseIsolationMVArun2v1DBnewDMwLT        .push_back(	boostedTau.tauID("byVLooseIsolationMVArun2v1DBnewDMwLT"			    ));
      nBranches_->tau_byVLooseIsolationMVArun2v1DBoldDMwLT        .push_back(boostedTau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"));
      nBranches_->tau_byVVLooseIsolationMVArun2v1DBoldDMwLT       .push_back(false);//boostedTau.tauID("byVVLooseIsolationMVArun2v1DBoldDMwLT"));// latest training just after re-running the standard tau ID. WP not there yet in the CMSSW release.
     
      //      nBranches_->tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT    .push_back(boostedTau.tauID("byVLooseIsolationMVArun2v1PWdR03oldDMwLT"			    ));
      nBranches_->tau_byVLooseIsolationMVArun2v1PWnewDMwLT        .push_back(boostedTau.tauID("byVLooseIsolationMVArun2v1PWnewDMwLT"			    ));
      //      nBranches_->tau_byVLooseIsolationMVArun2v1PWoldDMwLT        .push_back(boostedTau.tauID("byVLooseIsolationMVArun2v1PWoldDMwLT"			    ));
      nBranches_->tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT    .push_back(boostedTau.tauID("byVTightIsolationMVArun2v1DBdR03oldDMwLT"			    ));
      nBranches_->tau_byVTightIsolationMVArun2v1DBnewDMwLT        .push_back(boostedTau.tauID("byVTightIsolationMVArun2v1DBnewDMwLT"			    ));
      nBranches_->tau_byVTightIsolationMVArun2v1DBoldDMwLT        .push_back(boostedTau.tauID("byVTightIsolationMVArun2v1DBoldDMwLT"			    ));

      //      nBranches_->tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT     .push_back(boostedTau.tauID("byVTightIsolationMVArun2v1PWdR03oldDMwLT"			    ));
      nBranches_->tau_byVTightIsolationMVArun2v1PWnewDMwLT         .push_back(boostedTau.tauID("byVTightIsolationMVArun2v1PWnewDMwLT"			    ));
      //      nBranches_->tau_byVTightIsolationMVArun2v1PWoldDMwLT         .push_back(boostedTau.tauID("byVTightIsolationMVArun2v1PWoldDMwLT"			    ));
      nBranches_->tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT    .push_back(boostedTau.tauID("byVVTightIsolationMVArun2v1DBdR03oldDMwLT"			    ));
      nBranches_->tau_byVVTightIsolationMVArun2v1DBnewDMwLT        .push_back(boostedTau.tauID("byVVTightIsolationMVArun2v1DBnewDMwLT"			    ));
      nBranches_->tau_byVVTightIsolationMVArun2v1DBoldDMwLT        .push_back(boostedTau.tauID("byVVTightIsolationMVArun2v1DBoldDMwLT"			    ));
      //      nBranches_->tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT    .push_back(boostedTau.tauID("byVVTightIsolationMVArun2v1PWdR03oldDMwLT"			    ));
      nBranches_->tau_byVVTightIsolationMVArun2v1PWnewDMwLT        .push_back(boostedTau.tauID("byVVTightIsolationMVArun2v1PWnewDMwLT"			    ));

      //      nBranches_->tau_byVVTightIsolationMVArun2v1PWoldDMwLT        .push_back(boostedTau.tauID("byVVTightIsolationMVArun2v1PWoldDMwLT"			    ));
      // if(runFlags["doMultipleTauMVA_versions"]){
      nBranches_->tau_byIsolationMVArun2v2DBoldDMwLTraw.push_back(-1);
      nBranches_->tau_byVVLooseIsolationMVArun2v2DBoldDMwLT.push_back(-1);
      nBranches_->tau_byVLooseIsolationMVArun2v2DBoldDMwLT.push_back(-1);
      nBranches_->tau_byLooseIsolationMVArun2v2DBoldDMwLT.push_back(-1);
      nBranches_->tau_byMediumIsolationMVArun2v2DBoldDMwLT.push_back(-1);
      nBranches_->tau_byTightIsolationMVArun2v2DBoldDMwLT.push_back(-1);
      nBranches_->tau_byVTightIsolationMVArun2v2DBoldDMwLT.push_back(-1);
      nBranches_->tau_byVVTightIsolationMVArun2v2DBoldDMwLT.push_back(-1);
      //}
      // nBranches_->tau_byIsolationMVArun2v1DBoldDMwLTrawNew.push_back(boostedTau.tauID("byIsolationMVArun2v1DBoldDMwLTrawNew"));
      // nBranches_->tau_byVLooseIsolationMVArun2v1DBoldDMwLTNew.push_back(boostedTau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLTNew"));
      // nBranches_->tau_byLooseIsolationMVArun2v1DBoldDMwLTNew.push_back(boostedTau.tauID("byLooseIsolationMVArun2v1DBoldDMwLTNew"));
      // nBranches_->tau_byMediumIsolationMVArun2v1DBoldDMwLTNew.push_back(boostedTau.tauID("byMediumIsolationMVArun2v1DBoldDMwLTNew"));
      // nBranches_->tau_byTightIsolationMVArun2v1DBoldDMwLTNew.push_back(boostedTau.tauID("byTightIsolationMVArun2v1DBoldDMwLTNew"));
      // nBranches_->tau_byVTightIsolationMVArun2v1DBoldDMwLTNew.push_back(boostedTau.tauID("byVTightIsolationMVArun2v1DBoldDMwLTNew"));
      // nBranches_->tau_byVVTightIsolationMVArun2v1DBoldDMwLTNew.push_back(boostedTau.tauID("byVVTightIsolationMVArun2v1DBoldDMwLTNew"));

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
      //      nBranches_->tau_footprintCorrection                        .push_back(boostedTau.tauID("footprintCorrection"                         ));

   /*======================================================*/   
        
  }        
  nBranches_->tau_N =  taus_->size() + boostedTaus_->size();
  
}  
