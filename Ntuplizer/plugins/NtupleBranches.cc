#include "../interface/NtupleBranches.h"

//===================================================================================================================        
NtupleBranches::NtupleBranches( std::map< std::string, bool >& runFlags, TTree* tree )
   : tree_( tree )
{
   branch( runFlags );
}

//===================================================================================================================
NtupleBranches::~NtupleBranches( void )
{
}

//===================================================================================================================      
void NtupleBranches::branch( std::map< std::string, bool >& runFlags ){


  if ( runFlags["runOnMC"] ){
    if ( runFlags["doGenParticles"] ){
      /** genParticles */
      tree_->Branch( "genParticle_N"	     , &genParticle_N	       );
      tree_->Branch( "genParticle_pt"	     , &genParticle_pt	       ); 
//      tree_->Branch( "genParticle_px"	     , &genParticle_px	       ); 
//      tree_->Branch( "genParticle_py"	     , &genParticle_py	       ); 
//      tree_->Branch( "genParticle_pz"	     , &genParticle_pz	       ); 
//      tree_->Branch( "genParticle_e" 	     , &genParticle_e	       ); 
      tree_->Branch( "genParticle_eta"	     , &genParticle_eta        ); 
      tree_->Branch( "genParticle_phi"	     , &genParticle_phi        ); 
      tree_->Branch( "genParticle_mass"	     , &genParticle_mass       ); 
      tree_->Branch( "genParticle_pdgId"     , &genParticle_pdgId      );
      tree_->Branch( "genParticle_status"    , &genParticle_status     );
      tree_->Branch( "genParticle_isPrompt"  , &genParticle_isPrompt   );
      tree_->Branch( "genParticle_isDirectPromptTauDecayProduct"  , &genParticle_isDirectPromptTauDecayProduct);
      tree_->Branch( "genParticle_isDirectHardProcessTauDecayProductFinalState"  , &genParticle_isDirectHardProcessTauDecayProductFinalState);
      tree_->Branch( "genParticle_fromHardProcessFinalState"  , &genParticle_fromHardProcessFinalState   );
      tree_->Branch( "genParticle_mother"    , &genParticle_mother     );
      tree_->Branch( "genParticle_nMoth"     , &genParticle_nMoth      );
      tree_->Branch( "genParticle_nDau"	     , &genParticle_nDau       ); 
      tree_->Branch( "genParticle_dau"	     , &genParticle_dau        );
      tree_->Branch( "genParticle_tauvispt"	     , &genParticle_tauvispt        );
      tree_->Branch( "genParticle_tauviseta"	     , &genParticle_tauviseta        );
      tree_->Branch( "genParticle_tauvisphi"	     , &genParticle_tauvisphi       );
      tree_->Branch( "genParticle_tauvismass"	     , &genParticle_tauvismass        );
      tree_->Branch( "genParticle_taudecay"	     , &genParticle_taudecay        );


    } //doGenParticles
    
    if ( runFlags["doGenEvent"] ){
      /** generator info */
      tree_->Branch( "lheV_pt"	             , &lheV_pt                ); 
      tree_->Branch( "lheHT"	             , &lheHT                  ); 
      tree_->Branch( "lheNj"	             , &lheNj                  );
      tree_->Branch( "lheNb"	             , &lheNb                  );
      tree_->Branch( "lheNl"	             , &lheNl                  );
      tree_->Branch( "lheV_mass"           , &lheV_mass              ); 
      tree_->Branch( "genWeight"	         , &genWeight              );
      tree_->Branch( "genFacWeightUp"	     , &genFacWeightUp         );
      tree_->Branch( "genFacWeightDown"	   , &genFacWeightDown       );
      tree_->Branch( "genRenWeightUp"	     , &genRenWeightUp         );
      tree_->Branch( "genRenWeightDown"	   , &genRenWeightDown       );
      tree_->Branch( "genFacRenWeightUp"	 , &genFacRenWeightUp      );
      tree_->Branch( "genFacRenWeightDown" , &genFacRenWeightDown    );
      tree_->Branch( "qScale"	             , &qScale                 );
      tree_->Branch( "PDF_rms"	           , &PDF_rms                );
      tree_->Branch( "PDF_x"	             , &PDF_x                  );
      tree_->Branch( "PDF_xPDF"	           , &PDF_xPDF               );
      tree_->Branch( "PDF_id"	             , &PDF_id                 );

    } //doGenEvent
  } //runOnMC
  
  if ( runFlags["doElectrons"] ){
    /** electrons */
    tree_->Branch( "el_N"  		     	, &el_N 	               );
    tree_->Branch( "el_pdgId"  		     	, &el_pdgId	               );
    tree_->Branch( "el_charge"		     	, &el_charge	               );
    tree_->Branch( "el_e" 		     	, &el_e 	               );
    tree_->Branch( "el_eta"		     	, &el_eta	               );
    tree_->Branch( "el_phi"		     	, &el_phi	               );
    tree_->Branch( "el_mass"		     	, &el_mass	               );
    tree_->Branch( "el_pt"		     	, &el_pt	               );
    tree_->Branch( "el_et"		     	, &el_et	               );
    tree_->Branch( "el_superCluster_eta"     	, &el_superCluster_eta         );
    tree_->Branch( "el_pfRhoCorrRelIso03"    	, &el_pfRhoCorrRelIso03        );
    tree_->Branch( "el_pfRhoCorrRelIso04"    	, &el_pfRhoCorrRelIso04        );
    tree_->Branch( "el_pfDeltaCorrRelIso"    	, &el_pfDeltaCorrRelIso        );
    // tree_->Branch( "el_pfRelIso"  	     	, &el_pfRelIso  	       );
    tree_->Branch( "el_photonIso" 	     	, &el_photonIso 	       );
    tree_->Branch( "el_neutralHadIso"	     	, &el_neutralHadIso	       );
    tree_->Branch( "el_chargedHadIso"	     	, &el_chargedHadIso	       );
    // tree_->Branch( "el_trackIso"	      	, &el_trackIso  	       );
    tree_->Branch( "el_passConversionVeto"      , &el_passConversionVeto       );
    tree_->Branch( "el_full5x5_sigmaIetaIeta"   , &el_full5x5_sigmaIetaIeta    );
    tree_->Branch( "el_dEtaIn"	                , &el_dEtaIn                   );
    tree_->Branch( "el_dPhiIn"                  , &el_dPhiIn                   );
    tree_->Branch( "el_hOverE"   	        , &el_hOverE                   );
    tree_->Branch( "el_relIsoWithDBeta"         , &el_relIsoWithDBeta          );
    tree_->Branch( "el_ooEmooP"	                , &el_ooEmooP                  );
    tree_->Branch( "el_expectedMissingInnerHits", &el_expectedMissingInnerHits );  
    tree_->Branch( "el_d0"                      , &el_d0                       );
    tree_->Branch( "el_dz"		        , &el_dz                       );  
    tree_->Branch( "el_d0_allvertices"          , &el_d0_allvertices           );
    tree_->Branch( "el_dz_allvertices"          , &el_dz_allvertices           );  
    tree_->Branch( "el_dr03EcalRecHitSumEt"     , &el_dr03EcalRecHitSumEt      );
    tree_->Branch( "el_dr03HcalDepth1TowerSumEt", &el_dr03HcalDepth1TowerSumEt );
    tree_->Branch( "el_rho"                     , &el_rho                      );
    tree_->Branch( "el_ecalDriven"              , &el_ecalDriven               );
    tree_->Branch( "el_dEtaInSeed"              , &el_dEtaInSeed               );
    tree_->Branch( "el_full5x5_e2x5Max"         , &el_full5x5_e2x5Max          );
    tree_->Branch( "el_full5x5_e5x5"            , &el_full5x5_e5x5             );
    tree_->Branch( "el_full5x5_e1x5"            , &el_full5x5_e1x5             );
    tree_->Branch( "el_full5x5_r9"              , &el_full5x5_r9               );
    tree_->Branch( "el_dr03TkSumPt"             , &el_dr03TkSumPt              );
    tree_->Branch( "el_superCluster_e"          , &el_superCluster_e           );
    tree_->Branch( "el_hadronicOverEm"          , &el_hadronicOverEm           );  
    tree_->Branch( "el_seedEnergy"              , &el_seedEnergy               );  
    tree_->Branch( "el_isVetoElectron"	        , &el_isVetoElectron           );
    tree_->Branch( "el_isLooseElectron"	        , &el_isLooseElectron          );  
    tree_->Branch( "el_isMediumElectron"	, &el_isMediumElectron         );
    tree_->Branch( "el_isTightElectron"         , &el_isTightElectron          );  
    tree_->Branch( "el_isHeepElectron"	        , &el_isHeepElectron           );
    tree_->Branch( "el_isHltElectron"	        , &el_isHltElectron            );
    tree_->Branch( "el_isMVAMediumElectron"	, &el_isMVAMediumElectron      );
    tree_->Branch( "el_isMVATightElectron"      , &el_isMVATightElectron       );  
    tree_->Branch( "el_MVAscore"                , &el_MVAscore                 );  
    tree_->Branch( "el_MVAcategory"             , &el_MVAcategory              );  
//    tree_->Branch( "el_isVetoElectronBoosted"	, &el_isVetoElectronBoosted    );
//    tree_->Branch( "el_isMediumElectronBoosted"	, &el_isMediumElectronBoosted  );
//    tree_->Branch( "el_isTightElectronBoosted"  , &el_isTightElectronBoosted   );  
//    tree_->Branch( "el_isHeepElectronBoosted"   , &el_isHeepElectronBoosted    );
//    tree_->Branch( "el_isLooseElectronBoosted"  , &el_isLooseElectronBoosted   );  
    tree_->Branch( "el_isVetoElectronWithoutIPandIsolation"	, &el_isVetoElectronWithoutIPandIsolation    );
    tree_->Branch( "el_isMediumElectronWithoutIPandIsolation"	, &el_isMediumElectronWithoutIPandIsolation  );
    tree_->Branch( "el_isTightElectronWithoutIPandIsolation"  , &el_isTightElectronWithoutIPandIsolation   );  
    tree_->Branch( "el_isHeepElectronWithoutIPandIsolation"   , &el_isHeepElectronWithoutIPandIsolation    );
    tree_->Branch( "el_isLooseElectronWithoutIPandIsolation"  , &el_isLooseElectronWithoutIPandIsolation   );  
    // tree_->Branch( "el_pfRhoCorrRelIso03Boost"  , &el_pfRhoCorrRelIso03Boost   );
    // tree_->Branch( "el_pfRhoCorrRelIso04Boost"  , &el_pfRhoCorrRelIso04Boost   );
    // tree_->Branch( "el_pfDeltaCorrRelIsoBoost"  , &el_pfDeltaCorrRelIsoBoost   );
    // tree_->Branch( "el_pfRelIsoBoost"  	        , &el_pfRelIsoBoost	       );
    // tree_->Branch( "el_photonIsoBoost" 	        , &el_photonIsoBoost	       );
    // tree_->Branch( "el_neutralHadIsoBoost"      , &el_neutralHadIsoBoost       );
    // tree_->Branch( "el_chargedHadIsoBoost"      , &el_chargedHadIsoBoost       );  
    tree_->Branch( "el_SemileptonicPFIso"       , &el_SemileptonicPFIso        );
    // tree_->Branch( "el_SemileptonicCorrPFIso"   , &el_SemileptonicCorrPFIso    );
  } //doElectrons
  
  if ( runFlags["doMuons"] ){
    /** muons */
    tree_->Branch( "mu_N"  		      , &mu_N			   );
    tree_->Branch( "mu_pdgId"  		      , &mu_pdgId		   );
    tree_->Branch( "mu_charge"		      , &mu_charge		   );
    tree_->Branch( "mu_e" 		      , &mu_e			   );
    tree_->Branch( "mu_eta"		      , &mu_eta 		   );
    tree_->Branch( "mu_phi"		      , &mu_phi 		   );
    tree_->Branch( "mu_mass"		      , &mu_mass		   );
    tree_->Branch( "mu_pt"		      , &mu_pt  	  	   );
    tree_->Branch( "mu_isHighPtMuon"	      , &mu_isHighPtMuon      	   );
    tree_->Branch( "mu_isTightMuon"	      , &mu_isTightMuon       	   );
    tree_->Branch( "mu_isMediumMuon"	      , &mu_isMediumMuon    	   );
    tree_->Branch( "mu_isMediumMuonGH"	      , &mu_isMediumMuonGH    	   );
    tree_->Branch( "mu_isLooseMuon"	      , &mu_isLooseMuon       	   );
    tree_->Branch( "mu_isPFMuon"	      , &mu_isPFMuon	      	   );
    tree_->Branch( "mu_isSoftMuon"	      , &mu_isSoftMuon             );
    tree_->Branch( "mu_isGlobalMuon"	      , &mu_isGlobalMuon    	   );
    tree_->Branch( "mu_isTrackerMuon"	      , &mu_isTrackerMuon    	   );
    tree_->Branch( "mu_isTrackerHighPtMuon"	      , &mu_isTrackerHighPtMuon    	   );
    // tree_->Branch( "mu_pfRhoCorrRelIso03"     , &mu_pfRhoCorrRelIso03 	   );
    // tree_->Branch( "mu_pfRhoCorrRelIso04"     , &mu_pfRhoCorrRelIso04 	   );
    tree_->Branch( "mu_pfDeltaCorrRelIso"     , &mu_pfDeltaCorrRelIso 	   );
    // tree_->Branch( "mu_pfRelIso"  	      , &mu_pfRelIso	      	   );
    tree_->Branch( "mu_photonIso" 	      , &mu_photonIso	      	   );
    tree_->Branch( "mu_neutralHadIso"	      , &mu_neutralHadIso     	   );
    tree_->Branch( "mu_chargedHadIso"	      , &mu_chargedHadIso     	   );
    tree_->Branch( "mu_trackIso"	      , &mu_trackIso	      	   );
    tree_->Branch( "mu_trackCorrIso"	      , &mu_trackCorrIso	   );
    tree_->Branch( "mu_d0"                    , &mu_d0  	      	   );
    tree_->Branch( "mu_dz"                    , &mu_dz  	      	   );
    tree_->Branch( "mu_d0_allvertices"        , &mu_d0_allvertices  	   );
    tree_->Branch( "mu_dz_allvertices"        , &mu_dz_allvertices    	   );
    tree_->Branch( "mu_innerTrack_pt"	        , &mu_innerTrack_pt      	   );
    tree_->Branch( "mu_bestTrack_pt"	        , &mu_bestTrack_pt      	   );
    tree_->Branch( "mu_bestTrack_ptErr"	      , &mu_bestTrack_ptErr        );
    tree_->Branch( "mu_tunePTrack_pt"	        , &mu_tunePTrack_pt      	   );
    tree_->Branch( "mu_tunePTrack_ptErr"	    , &mu_tunePTrack_ptErr        );
    // tree_->Branch( "mu_pfRhoCorrRelIso03Boost", &mu_pfRhoCorrRelIso03Boost );
    // tree_->Branch( "mu_pfRhoCorrRelIso04Boost", &mu_pfRhoCorrRelIso04Boost );
    // tree_->Branch( "mu_pfDeltaCorrRelIsoBoost", &mu_pfDeltaCorrRelIsoBoost );
    // tree_->Branch( "mu_pfRelIsoBoost"  	      , &mu_pfRelIsoBoost	   );
    // tree_->Branch( "mu_photonIsoBoost" 	      , &mu_photonIsoBoost	   );
    // tree_->Branch( "mu_neutralHadIsoBoost"    , &mu_neutralHadIsoBoost     );
    // tree_->Branch( "mu_chargedHadIsoBoost"    , &mu_chargedHadIsoBoost     );  
    tree_->Branch( "mu_normChi2"  	      , &mu_normChi2	    	   );
    tree_->Branch( "mu_trackerHits"	      , &mu_trackerHits     	   );
    tree_->Branch( "mu_matchedStations"	      , &mu_matchedStations 	   );
    tree_->Branch( "mu_pixelHits" 	      , &mu_pixelHits	    	   );
    tree_->Branch( "mu_globalHits"	      , &mu_globalHits        	   );
    tree_->Branch( "mu_SemileptonicPFIso"     , &mu_SemileptonicPFIso	   );
    // tree_->Branch( "mu_SemileptonicCorrPFIso" , &mu_SemileptonicCorrPFIso  );
  } //doMuons
  
  if ( runFlags["doTaus"] ){
    /** taus */
    tree_->Branch( "tau_N"  		     	 , &tau_N		       );
    tree_->Branch( "tau_pdgId"               	 , &tau_pdgId		       );
    tree_->Branch( "tau_charge"		     	 , &tau_charge  	       );
    tree_->Branch( "tau_e" 		     	 , &tau_e		       );
    tree_->Branch( "tau_eta"		     	 , &tau_eta		       );
    tree_->Branch( "tau_phi"		     	 , &tau_phi		       );
    tree_->Branch( "tau_mass"		     	 , &tau_mass		       );
    tree_->Branch( "tau_pt"		     	 , &tau_pt		       );
    tree_->Branch( "tau_d0"                  	 , &tau_d0		       );
    tree_->Branch( "tau_dz"                  	 , &tau_dz		       );
    
    // YT added
//    tree_->Branch( "tau_associated_pdgId"    	 , &tau_associated_pdgId       );
//    tree_->Branch( "tau_associated_pt"    	 , &tau_associated_pt       );
//    tree_->Branch( "tau_associated_eta"    	 , &tau_associated_eta       );
//    tree_->Branch( "tau_associated_phi"    	 , &tau_associated_phi       );
//    tree_->Branch( "tau_associated_dr"    	 , &tau_associated_dr       );

//    tree_->Branch( "tau_n_total"    	 , &tau_n_total       );
    tree_->Branch( "tau_n_ch"    	 , &tau_n_ch       );
    tree_->Branch( "tau_n_nh"    	 , &tau_n_nh       );
    //    tree_->Branch( "tau_n_h_f"    	 , &tau_n_h_f       );
    //    tree_->Branch( "tau_n_em_f"    	 , &tau_n_em_f       );
    tree_->Branch( "tau_n_gamma"    	 , &tau_n_gamma       );
    //    tree_->Branch( "tau_n_e"    	 , &tau_n_e       );
    //    tree_->Branch( "tau_n_mu"    	 , &tau_n_mu       );

    // tree_->Branch( "tau_pfRhoCorrRelIso03"   	 , &tau_pfRhoCorrRelIso03      );
    // tree_->Branch( "tau_pfRhoCorrRelIso04"   	 , &tau_pfRhoCorrRelIso04      );
    // tree_->Branch( "tau_pfDeltaCorrRelIso"   	 , &tau_pfDeltaCorrRelIso      );
    // tree_->Branch( "tau_pfRelIso"  	     	 , &tau_pfRelIso	       );
    // tree_->Branch( "tau_photonIso" 	     	 , &tau_photonIso	       );
    // tree_->Branch( "tau_neutralHadIso"	     	 , &tau_neutralHadIso	       );
    // tree_->Branch( "tau_chargedHadIso"	     	 , &tau_chargedHadIso	       );
    tree_->Branch( "tau_photonPtSumOutsideSignalCone"      , &tau_photonPtSumOutsideSignalCone     );  
    // tree_->Branch( "tau_trackIso"	     	 , &tau_trackIso	       );
    // tree_->Branch( "tau_pfRhoCorrRelIso03Boost"  , &tau_pfRhoCorrRelIso03Boost );
    // tree_->Branch( "tau_pfRhoCorrRelIso04Boost"  , &tau_pfRhoCorrRelIso04Boost );
    // tree_->Branch( "tau_pfDeltaCorrRelIsoBoost"  , &tau_pfDeltaCorrRelIsoBoost );
    // tree_->Branch( "tau_pfRelIsoBoost"  	 , &tau_pfRelIsoBoost	       );
    // tree_->Branch( "tau_photonIsoBoost" 	 , &tau_photonIsoBoost         );
    // tree_->Branch( "tau_neutralHadIsoBoost"      , &tau_neutralHadIsoBoost     );
    // tree_->Branch( "tau_chargedHadIsoBoost"      , &tau_chargedHadIsoBoost     );  
    tree_->Branch( "tau_TauType"		 , &tau_TauType		       );
    tree_->Branch( "tau_decayMode"		 , &tau_decayMode	       ); // YT added
    tree_->Branch( "tau_chargedPionPt"		 , &tau_chargedPionPt	       ); // YT added
    tree_->Branch( "tau_neutralPionPt"		 , &tau_neutralPionPt	       ); // YT added

    // YT added
    tree_->Branch( "tau_nPhoton"                 , &tau_nPhoton                );
    tree_->Branch( "tau_ptWeightedDetaStrip"     , &tau_ptWeightedDetaStrip    );
    tree_->Branch( "tau_ptWeightedDphiStrip"     , &tau_ptWeightedDphiStrip    );
    tree_->Branch( "tau_ptWeightedDrSignal"      , &tau_ptWeightedDrSignal     );
    tree_->Branch( "tau_ptWeightedDrIsolation"   , &tau_ptWeightedDrIsolation  );
    tree_->Branch( "tau_leadingTrackChi2"        , &tau_leadingTrackChi2       );
    tree_->Branch( "tau_leadingTrackPt"          , &tau_leadingTrackPt         );
    tree_->Branch( "tau_eRatio"                  , &tau_eRatio                 );
    tree_->Branch( "tau_dxy_Sig"                 , &tau_dxy_Sig                );
    tree_->Branch( "tau_ip3d"                    , &tau_ip3d                   );
    tree_->Branch( "tau_ip3d_Sig"                , &tau_ip3d_Sig               );
    tree_->Branch( "tau_hasSecondaryVertex"      , &tau_hasSecondaryVertex     );
    //    tree_->Branch( "tau_decayDistMag_x"          , &tau_decayDistMag_x         );
    //    tree_->Branch( "tau_decayDistMag_y"          , &tau_decayDistMag_y         );
    //    tree_->Branch( "tau_decayDistMag_z"          , &tau_decayDistMag_z         );
    tree_->Branch( "tau_decayDistMag"            , &tau_decayDistMag           );
    tree_->Branch( "tau_flightLenthSig"          , &tau_flightLenthSig         );
    // YT added end
  
    
    /** tau discriminants */
    tree_->Branch( "tau_decayModeFindingNewDMs"	              , &tau_decayModeFindingNewDMs			 );
    tree_->Branch( "tau_decayModeFinding"  			      , &tau_decayModeFinding				 );
    tree_->Branch( "tau_byLooseCombinedIsolationDeltaBetaCorr3Hits" , &tau_byLooseCombinedIsolationDeltaBetaCorr3Hits  );
    tree_->Branch( "tau_byMediumCombinedIsolationDeltaBetaCorr3Hits", &tau_byMediumCombinedIsolationDeltaBetaCorr3Hits );
    tree_->Branch( "tau_byTightCombinedIsolationDeltaBetaCorr3Hits" , &tau_byTightCombinedIsolationDeltaBetaCorr3Hits  );
    tree_->Branch( "tau_byCombinedIsolationDeltaBetaCorrRaw3Hits"   , &tau_byCombinedIsolationDeltaBetaCorrRaw3Hits	 );
    tree_->Branch( "tau_chargedIsoPtSum"			      , &tau_chargedIsoPtSum			         );
    tree_->Branch( "tau_neutralIsoPtSum"			      , &tau_neutralIsoPtSum			         );
    tree_->Branch( "tau_puCorrPtSum"				      , &tau_puCorrPtSum 				 );
     

    tree_->Branch( "tau_chargedIsoPtSumdR03"                       , &tau_chargedIsoPtSumdR03);
    //      tree_->Branch( "tau_footprintCorrectiondR03"                    , &tau_footprintCorrectiondR03);
    tree_->Branch( "tau_neutralIsoPtSumdR03"                        , &tau_neutralIsoPtSumdR03);
    //      tree_->Branch( "tau_neutralIsoPtSumWeight"                      , &tau_neutralIsoPtSumWeight);
    //      tree_->Branch( "tau_neutralIsoPtSumWeightdR03"                  , &tau_neutralIsoPtSumWeightdR03);
    tree_->Branch( "tau_photonPtSumOutsideSignalConedR03"           , &tau_photonPtSumOutsideSignalConedR03);

    tree_->Branch( "tau_byIsolationMVArun2v1DBdR03oldDMwLTraw"      , &tau_byIsolationMVArun2v1DBdR03oldDMwLTraw);
    tree_->Branch( "tau_byIsolationMVArun2v1DBnewDMwLTraw"          , &tau_byIsolationMVArun2v1DBnewDMwLTraw);
    tree_->Branch( "tau_byIsolationMVArun2v1DBoldDMwLTraw"          , &tau_byIsolationMVArun2v1DBoldDMwLTraw);
    tree_->Branch( "tau_byIsolationMVArun2v1PWnewDMwLTraw"          , &tau_byIsolationMVArun2v1PWnewDMwLTraw);
    tree_->Branch( "tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT"    , &tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT);
    tree_->Branch( "tau_byLooseIsolationMVArun2v1DBnewDMwLT"        , &tau_byLooseIsolationMVArun2v1DBnewDMwLT);
    tree_->Branch( "tau_byLooseIsolationMVArun2v1DBoldDMwLT"        , &tau_byLooseIsolationMVArun2v1DBoldDMwLT);
    tree_->Branch( "tau_byLooseIsolationMVArun2v1PWnewDMwLT"        , &tau_byLooseIsolationMVArun2v1PWnewDMwLT);
    tree_->Branch( "tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT"    , &tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT);
    tree_->Branch( "tau_byMediumIsolationMVArun2v1DBnewDMwLT"        , &tau_byMediumIsolationMVArun2v1DBnewDMwLT);
    tree_->Branch( "tau_byMediumIsolationMVArun2v1DBoldDMwLT"        , &tau_byMediumIsolationMVArun2v1DBoldDMwLT);
    
    tree_->Branch( "tau_byMediumIsolationMVArun2v1PWnewDMwLT"        , &tau_byMediumIsolationMVArun2v1PWnewDMwLT);
     
    tree_->Branch( "tau_byTightIsolationMVArun2v1DBdR03oldDMwLT"     , &tau_byTightIsolationMVArun2v1DBdR03oldDMwLT);
    tree_->Branch( "tau_byTightIsolationMVArun2v1DBnewDMwLT"         , &tau_byTightIsolationMVArun2v1DBnewDMwLT);
    tree_->Branch( "tau_byTightIsolationMVArun2v1DBoldDMwLT"         , &tau_byTightIsolationMVArun2v1DBoldDMwLT);
     
    tree_->Branch( "tau_byTightIsolationMVArun2v1PWnewDMwLT"         , &tau_byTightIsolationMVArun2v1PWnewDMwLT);
    tree_->Branch( "tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT"    , &tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT);

    tree_->Branch( "tau_byVLooseIsolationMVArun2v1DBnewDMwLT"        , &tau_byVLooseIsolationMVArun2v1DBnewDMwLT);
    tree_->Branch( "tau_byVLooseIsolationMVArun2v1DBoldDMwLT"        , &tau_byVLooseIsolationMVArun2v1DBoldDMwLT);
    tree_->Branch( "tau_byVVLooseIsolationMVArun2v1DBoldDMwLT"       , &tau_byVVLooseIsolationMVArun2v1DBoldDMwLT);
      
   
    tree_->Branch( "tau_byVLooseIsolationMVArun2v1PWnewDMwLT"        , &tau_byVLooseIsolationMVArun2v1PWnewDMwLT);
    tree_->Branch( "tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT"    , &tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT);
    tree_->Branch( "tau_byVTightIsolationMVArun2v1DBnewDMwLT"        , &tau_byVTightIsolationMVArun2v1DBnewDMwLT);
    tree_->Branch( "tau_byVTightIsolationMVArun2v1DBoldDMwLT"        , &tau_byVTightIsolationMVArun2v1DBoldDMwLT);

    
    tree_->Branch( "tau_byVTightIsolationMVArun2v1PWnewDMwLT"         , &tau_byVTightIsolationMVArun2v1PWnewDMwLT);
    
    tree_->Branch( "tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT"    , &tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT);
    tree_->Branch( "tau_byVVTightIsolationMVArun2v1DBnewDMwLT"        , &tau_byVVTightIsolationMVArun2v1DBnewDMwLT);
    tree_->Branch( "tau_byVVTightIsolationMVArun2v1DBoldDMwLT"        , &tau_byVVTightIsolationMVArun2v1DBoldDMwLT);
      
    tree_->Branch( "tau_byVVTightIsolationMVArun2v1PWnewDMwLT"        , &tau_byVVTightIsolationMVArun2v1PWnewDMwLT);

      
    // add multimple variables to save multiple MVA versions
    // if ( runFlags["doMultipleTauMVAversions"] )
    tree_->Branch("tau_byIsolationMVArun2v2DBoldDMwLTraw", &tau_byIsolationMVArun2v2DBoldDMwLTraw);
    tree_->Branch("tau_byVVLooseIsolationMVArun2v2DBoldDMwLT", &tau_byVVLooseIsolationMVArun2v2DBoldDMwLT);
    tree_->Branch("tau_byVLooseIsolationMVArun2v2DBoldDMwLT", &tau_byVLooseIsolationMVArun2v2DBoldDMwLT);
    tree_->Branch("tau_byLooseIsolationMVArun2v2DBoldDMwLT", &tau_byLooseIsolationMVArun2v2DBoldDMwLT);
    tree_->Branch("tau_byMediumIsolationMVArun2v2DBoldDMwLT", &tau_byMediumIsolationMVArun2v2DBoldDMwLT);
    tree_->Branch("tau_byTightIsolationMVArun2v2DBoldDMwLT", &tau_byTightIsolationMVArun2v2DBoldDMwLT);
    tree_->Branch("tau_byVTightIsolationMVArun2v2DBoldDMwLT", &tau_byVTightIsolationMVArun2v2DBoldDMwLT);
    tree_->Branch("tau_byVVTightIsolationMVArun2v2DBoldDMwLT", &tau_byVVTightIsolationMVArun2v2DBoldDMwLT);
      
    // }

    tree_->Branch( "tau_againstElectronMVA6raw"                     , &tau_againstElectronMVA6raw);
    tree_->Branch( "tau_againstElectronMVA6category"                , &tau_againstElectronMVA6category);
    tree_->Branch( "tau_againstElectronVLooseMVA6"                  , &tau_againstElectronVLooseMVA6);
    tree_->Branch( "tau_againstElectronLooseMVA6"                   , &tau_againstElectronLooseMVA6);
    tree_->Branch( "tau_againstElectronMediumMVA6"                  , &tau_againstElectronMediumMVA6);
    tree_->Branch( "tau_againstElectronTightMVA6"                   , &tau_againstElectronTightMVA6);
    tree_->Branch( "tau_againstElectronVTightMVA6"                  , &tau_againstElectronVTightMVA6);
      

    tree_->Branch( "tau_againstMuonLoose3"                          , &tau_againstMuonLoose3);
    tree_->Branch( "tau_againstMuonTight3"                          , &tau_againstMuonTight3); 
      
    tree_->Branch( "tau_byPhotonPtSumOutsideSignalCone"             , &tau_byPhotonPtSumOutsideSignalCone);
    //      tree_->Branch( "tau_footprintCorrection"                        , &tau_footprintCorrection);


    //doTausBoosted
  } //doTaus
      

 

  
 
  
  if (runFlags["doTriggerDecisions"]) {
    /** HLT trigger decisions */
    tree_->Branch("HLT_isFired", &HLT_isFired );
  }

  
  if (runFlags["doTriggerObjects"]) {
    /** HLT trigger objects */
    tree_->Branch("triggerObject_pt"		, &triggerObject_pt		);
    tree_->Branch("triggerObject_eta"		, &triggerObject_eta		);
    tree_->Branch("triggerObject_phi"		, &triggerObject_phi	        );
    tree_->Branch("triggerObject_mass"		, &triggerObject_mass	        );
    tree_->Branch("triggerObject_lastname"	, &triggerObject_lastname	);
    tree_->Branch("triggerObject_filterLabels"	, &triggerObject_filterLabels	);
    tree_->Branch("triggerObject_firedTrigger"	, &triggerObject_firedTrigger	);
    tree_->Branch("triggerObject_filterIDs"	, &triggerObject_filterIDs	);

  } //doTriggerObjects
  
  if (runFlags["doHltFilters"]) {
    /** HLT filter decisions */
    tree_->Branch("passFilter_HBHE"                 ,&passFilter_HBHE_                ,"passFilter_HBHE_/O");
    tree_->Branch("passFilter_HBHELoose"            ,&passFilter_HBHELoose_	          ,"passFilter_HBHELoose_/O");
    tree_->Branch("passFilter_HBHETight"            ,&passFilter_HBHETight_	          ,"passFilter_HBHETight_/O");
    tree_->Branch("passFilter_HBHEIso"              ,&passFilter_HBHEIso_	            ,"passFilter_HBHEIso_/O");
    tree_->Branch("passFilter_CSCHalo"              ,&passFilter_CSCHalo_             ,"passFilter_CSCHalo_/O");
    tree_->Branch("passFilter_CSCTightHalo2015"     ,&passFilter_CSCTightHalo2015_    ,"passFilter_CSCTightHalo2015_/O");
    tree_->Branch("passFilter_HCALlaser"            ,&passFilter_HCALlaser_           ,"passFilter_HCALlaser_/O");
    tree_->Branch("passFilter_ECALDeadCell"         ,&passFilter_ECALDeadCell_        ,"passFilter_ECALDeadCell_/O");
    tree_->Branch("passFilter_GoodVtx"              ,&passFilter_GoodVtx_             ,"passFilter_GoodVtx_/O");
    tree_->Branch("passFilter_TrkFailure"           ,&passFilter_TrkFailure_          ,"passFilter_TrkFailure_/O");
    tree_->Branch("passFilter_EEBadSc"              ,&passFilter_EEBadSc_             ,"passFilter_EEBadSc_/O");
    tree_->Branch("passFilter_ECALlaser"            ,&passFilter_ECALlaser_           ,"passFilter_ECALlaser_/O");
    tree_->Branch("passFilter_TrkPOG"               ,&passFilter_TrkPOG_              ,"passFilter_TrkPOG_/O");
    tree_->Branch("passFilter_TrkPOG_manystrip"     ,&passFilter_TrkPOG_manystrip_    ,"passFilter_TrkPOG_manystrip_/O");
    tree_->Branch("passFilter_TrkPOG_toomanystrip"  ,&passFilter_TrkPOG_toomanystrip_ ,"passFilter_TrkPOG_toomanystrip_/O");
    tree_->Branch("passFilter_TrkPOG_logError"      ,&passFilter_TrkPOG_logError_     ,"passFilter_TrkPOG_logError_/O");
    tree_->Branch("passFilter_METFilters"           ,&passFilter_METFilters_          ,"passFilter_METFilters_/O");
    
     //NEW FOR ICHEP
    tree_->Branch("passFilter_CSCTightHaloTrkMuUnvetoFilter", &passFilter_CSCTightHaloTrkMuUnvetoFilter_   ,"passFilter_CSCTightHaloTrkMuUnvetoFilter_/O");
    tree_->Branch("passFilter_globalTightHalo2016"          , &passFilter_globalTightHalo2016_             ,"passFilter_globalTightHalo2016_/O");
    tree_->Branch("passFilter_globalSuperTightHalo2016"          , &passFilter_globalSuperTightHalo2016_             ,"passFilter_globalSuperTightHalo2016_/O");
    tree_->Branch("passFilter_HcalStripHalo"                , &passFilter_HcalStripHalo_                   ,"passFilter_HcalStripHalo_/O");
    tree_->Branch("passFilter_chargedHadronTrackResolution" , &passFilter_chargedHadronTrackResolution_    ,"passFilter_chargedHadronTrackResolution_/O");
    tree_->Branch("passFilter_muonBadTrack"                 , &passFilter_muonBadTrack_                    ,"passFilter_muonBadTrack_/O");
    tree_->Branch("flag_badMuons"                 , &flag_badMuons_                    ,"flag_badMuons_/O");
    tree_->Branch("flag_duplicateMuons"                 , &flag_duplicateMuons_                    ,"flag_duplicateMuons_/O");
    tree_->Branch("flag_nobadMuons"                 , &flag_nobadMuons_                    ,"flag_nobadMuons_/O");
    tree_->Branch("passFilter_ecalBadCalib_"    ,&passFilter_ecalBadCalib_, "passFilter_ecalBadCalib_/O");
  } //do HltFilters

  if (runFlags["doMissingEt"]) {
    /** MET */
    tree_->Branch( "rho", &rho );
    tree_->Branch("METraw_et"		        , &METraw_et	     );
    tree_->Branch("METraw_phi"		        , &METraw_phi	     ); 
    tree_->Branch("METraw_sumEt"		, &METraw_sumEt	     );   
    tree_->Branch("MET_corrPx"		        , &MET_corrPx	     ); 
    tree_->Branch("MET_corrPy"		        , &MET_corrPy	     );   
    tree_->Branch("MET_et"	                , &MET_et  	     ); 
    tree_->Branch("MET_phi"	                , &MET_phi           );
    tree_->Branch("MET_puppi_et"	        , &MET_puppi_et      ); 
    tree_->Branch("MET_puppi_phi"               , &MET_puppi_phi     );
    tree_->Branch("MET_sumEt"	                , &MET_sumEt 	     ); 
    tree_->Branch("MET_JetEnUp"	                , &MET_JetEnUp 	     ); 
    tree_->Branch("MET_JetEnDown"	                , &MET_JetEnDown 	     ); 
    tree_->Branch("MET_JetResUp"	                , &MET_JetResUp 	     ); 
    tree_->Branch("MET_JetResDown"	                , &MET_JetResDown 	     ); 
    tree_->Branch("MET_UnclusteredEnUp"	                , &MET_UnclusteredEnUp 	     ); 
    tree_->Branch("MET_UnclusteredEnDown"	                , &MET_UnclusteredEnDown 	     ); 
    
  } //doMissingEt

  if ( runFlags["doMETSVFIT"] ){
    /** MET SVift*/
    tree_->Branch( "MET_significance"                                 , &MET_significance );
    tree_->Branch( "MET_cov00"                                        , &MET_cov00 );
    tree_->Branch( "MET_cov10"                                        , &MET_cov10 );
    tree_->Branch( "MET_cov11"                                        , &MET_cov11 );
  }

  if ( runFlags["doMVAMET"] ){
    /** MET SVift*/
    tree_->Branch("MET_Nmva"	                , &MET_Nmva 	     ); 
    tree_->Branch("MET_mva_et"	                , &MET_mva_et        ); 
    tree_->Branch("MET_mva_phi"                 , &MET_mva_phi       );
    tree_->Branch( "MET_mva_cov00"                                        , &MET_mva_cov00 );
    tree_->Branch( "MET_mva_cov10"                                        , &MET_mva_cov10 );
    tree_->Branch( "MET_mva_cov11"                                        , &MET_mva_cov11 );
    tree_->Branch( "MET_mva_recoil_pt"                                        , &MET_mva_recoil_pt );
    tree_->Branch( "MET_mva_recoil_eta"                                        , &MET_mva_recoil_eta );
    tree_->Branch( "MET_mva_recoil_phi"                                        , &MET_mva_recoil_phi );
    tree_->Branch( "MET_mva_recoil_pdgId"                                        , &MET_mva_recoil_pdgId );

  }

  
  /*------------- ------EVENT infos-----------------------------*/
  tree_->Branch("EVENT_event"	 , &EVENT_event     );
  tree_->Branch("EVENT_run"	 , &EVENT_run	    );
  tree_->Branch("EVENT_lumiBlock", &EVENT_lumiBlock );
  
  if (runFlags["runOnMC"]) {
    if (runFlags["doPileUp"]) {
      /*--------------------------PU infos--------------------------*/
      tree_->Branch("nPuVtxTrue", &nPuVtxTrue );	
      tree_->Branch("nPuVtx"    , &nPuVtx     );
      tree_->Branch("bX"	, &bX	      );
    } //doPileUp
  } //runOnMC
  
  if (runFlags["doVertices"]) {  
    /*--------------------------PV infos--------------------------*/
    tree_->Branch("PV_N"     , &PV_N      );
    tree_->Branch("PV_filter", &PV_filter );
    tree_->Branch("PV_chi2"  , &PV_chi2   );
    tree_->Branch("PV_ndof"  , &PV_ndof   );
    tree_->Branch("PV_rho"   , &PV_rho    );
    tree_->Branch("PV_z"     , &PV_z      );
    tree_->Branch("BeamSpot_x0" , &BeamSpot_x0     );
    tree_->Branch("BeamSpot_y0" , &BeamSpot_y0     );
    tree_->Branch("BeamSpot_z0" , &BeamSpot_z0     );


  }
  if (runFlags["doJpsiMu"]){
    tree_->Branch("IsJpsiMu" , &IsJpsiMu  );
    tree_->Branch("IsJpsiEle", &IsJpsiEle );

    tree_->Branch("Jpsi_mu1_pt", &Jpsi_mu1_pt );
    tree_->Branch("Jpsi_mu1_eta", &Jpsi_mu1_eta );
    tree_->Branch("Jpsi_mu1_phi", &Jpsi_mu1_phi );
    tree_->Branch("Jpsi_mu1_ch", &Jpsi_mu1_ch );
    tree_->Branch("Jpsi_mu1_isLoose"  , &Jpsi_mu1_isLoose   );
    tree_->Branch("Jpsi_mu1_isTight"  , &Jpsi_mu1_isTight   );
    tree_->Branch("Jpsi_mu1_isPF"     , &Jpsi_mu1_isPF      );
    tree_->Branch("Jpsi_mu1_isGlobal" , &Jpsi_mu1_isGlobal  );
    tree_->Branch("Jpsi_mu1_isTracker", &Jpsi_mu1_isTracker );
    tree_->Branch("Jpsi_mu1_isSoft"   , &Jpsi_mu1_isSoft    );




    tree_->Branch("Jpsi_mu2_pt", &Jpsi_mu2_pt );
    tree_->Branch("Jpsi_mu2_eta", &Jpsi_mu2_eta );
    tree_->Branch("Jpsi_mu2_phi", &Jpsi_mu2_phi );
    tree_->Branch("Jpsi_mu2_ch", &Jpsi_mu2_ch );
    tree_->Branch("Jpsi_mu2_isLoose"  , &Jpsi_mu2_isLoose   );
    tree_->Branch("Jpsi_mu2_isTight"  , &Jpsi_mu2_isTight   );
    tree_->Branch("Jpsi_mu2_isPF"     , &Jpsi_mu2_isPF      );
    tree_->Branch("Jpsi_mu2_isGlobal" , &Jpsi_mu2_isGlobal  );
    tree_->Branch("Jpsi_mu2_isTracker", &Jpsi_mu2_isTracker );
    tree_->Branch("Jpsi_mu2_isSoft"   , &Jpsi_mu2_isSoft    );

    tree_->Branch("Jpsi_mu3_pt", &Jpsi_mu3_pt );
    tree_->Branch("Jpsi_mu3_eta", &Jpsi_mu3_eta );
    tree_->Branch("Jpsi_mu3_phi", &Jpsi_mu3_phi );
    tree_->Branch("Jpsi_mu3_ch", &Jpsi_mu3_ch );
    tree_->Branch("Jpsi_mu3_isLoose"  , &Jpsi_mu3_isLoose   );
    tree_->Branch("Jpsi_mu3_isTight"  , &Jpsi_mu3_isTight   );
    tree_->Branch("Jpsi_mu3_isPF"     , &Jpsi_mu3_isPF      );
    tree_->Branch("Jpsi_mu3_isGlobal" , &Jpsi_mu3_isGlobal  );
    tree_->Branch("Jpsi_mu3_isTracker", &Jpsi_mu3_isTracker );
    tree_->Branch("Jpsi_mu3_isSoft"   , &Jpsi_mu3_isSoft    );
    tree_->Branch("Jpsi_mu3_x"   , &Jpsi_mu3_x    );
    tree_->Branch("Jpsi_mu3_y"   , &Jpsi_mu3_y    );
    tree_->Branch("Jpsi_mu3_z"   , &Jpsi_mu3_z    );

    tree_->Branch("Jpsi_mu3_isopt03", &Jpsi_mu3_isopt03 );
    tree_->Branch("Jpsi_mu3_isopt04", &Jpsi_mu3_isopt04 );
    tree_->Branch("Jpsi_mu3_isopt05", &Jpsi_mu3_isopt05 );
    tree_->Branch("Jpsi_mu3_isopt06", &Jpsi_mu3_isopt06 );
    tree_->Branch("Jpsi_mu3_isopt07", &Jpsi_mu3_isopt07 );
    tree_->Branch("Jpsi_dr_mu3pf"    , &Jpsi_dr_mu3pf     );



    tree_->Branch("Jpsi_dx", &Jpsi_dx );
    tree_->Branch("Jpsi_dy", &Jpsi_dy );
    tree_->Branch("Jpsi_dz", &Jpsi_dz );
    tree_->Branch("Jpsi_pt", &Jpsi_pt );
    tree_->Branch("Jpsi_eta", &Jpsi_eta );
    tree_->Branch("Jpsi_phi", &Jpsi_phi );
    tree_->Branch("Jpsi_mass", &Jpsi_mass );
    tree_->Branch("Jpsi_vtxprob", &Jpsi_vtxprob );
    tree_->Branch("Jpsi_vtxz", &Jpsi_vtxz );

    tree_->Branch("Jpsi_trimu_dx", &Jpsi_trimu_dx );
    tree_->Branch("Jpsi_trimu_dy", &Jpsi_trimu_dy );
    tree_->Branch("Jpsi_trimu_dz", &Jpsi_trimu_dz );
    tree_->Branch("Jpsi_trimu_pt", &Jpsi_trimu_pt );
    tree_->Branch("Jpsi_trimu_eta", &Jpsi_trimu_eta );
    tree_->Branch("Jpsi_trimu_phi", &Jpsi_trimu_phi );
    tree_->Branch("Jpsi_trimu_mass", &Jpsi_trimu_mass );
    tree_->Branch("Jpsi_trimu_vtxprob", &Jpsi_trimu_vtxprob );
    tree_->Branch("Jpsi_trimu_vtxz", &Jpsi_trimu_vtxz );

    tree_->Branch("Jpsi_PV_x", &Jpsi_PV_x );
    tree_->Branch("Jpsi_PV_y", &Jpsi_PV_y );
    tree_->Branch("Jpsi_PV_z", &Jpsi_PV_z );

    //flight part
    tree_->Branch("Jpsi_flightSig3D", &Jpsi_flightSig3D);
    tree_->Branch("Jpsi_flightLength3D", &Jpsi_flightLength3D);
    tree_->Branch("Jpsi_flightLengthErr3D", &Jpsi_flightLengthErr3D);
    tree_->Branch("Jpsi_flightSig2D", &Jpsi_flightSig2D);
    tree_->Branch("Jpsi_flightLength2D", &Jpsi_flightLength2D);
    tree_->Branch("Jpsi_flightLengthErr2D", &Jpsi_flightLengthErr2D);

    tree_->Branch("Jpsi_trimu_flightSig3D", &Jpsi_trimu_flightSig3D);
    tree_->Branch("Jpsi_trimu_flightLength3D", &Jpsi_trimu_flightLength3D);
    tree_->Branch("Jpsi_trimu_flightLengthErr3D", &Jpsi_trimu_flightLengthErr3D);
    tree_->Branch("Jpsi_trimu_flightSig2D", &Jpsi_trimu_flightSig2D);
    tree_->Branch("Jpsi_trimu_flightLength2D", &Jpsi_trimu_flightLength2D);
    tree_->Branch("Jpsi_trimu_flightLengthErr2D", &Jpsi_trimu_flightLengthErr2D);

  }
}

//=================================================================================================================== 
void NtupleBranches::reset( void ){

  /** genParticle */
  genParticle_N = 0;
  genParticle_pt.clear();
//  genParticle_px.clear();
//  genParticle_py.clear();
//  genParticle_pz.clear();
//  genParticle_e.clear();
  genParticle_eta.clear();
  genParticle_phi.clear();
  genParticle_mass.clear();
  genParticle_pdgId.clear();
  genParticle_isPrompt.clear();
  genParticle_isDirectPromptTauDecayProduct.clear();
  genParticle_fromHardProcessFinalState.clear();
  genParticle_isDirectHardProcessTauDecayProductFinalState.clear();
  genParticle_status.clear();
  genParticle_mother.clear();
  genParticle_nMoth.clear();
  genParticle_nDau.clear();
  genParticle_dau.clear();
  genParticle_tauvispt.clear();
  genParticle_tauviseta.clear();
  genParticle_tauvisphi.clear();
  genParticle_tauvismass.clear();
  genParticle_taudecay.clear();
  
  /** generator info */
  genWeight   = 0;
  qScale      = 0;
  genFacWeightUp       = 0;
  genFacWeightDown     = 0;
  genRenWeightUp       = 0;
  genRenWeightDown     = 0;
  genFacRenWeightUp    = 0;
  genFacRenWeightDown  = 0;
  PDF_rms = 0;
  PDF_id.clear();  
  PDF_x.clear();	
  PDF_xPDF.clear();
  lheV_pt = 0;
  lheHT = 0;
  lheNj = 0;
  lheNb = 0;
  lheV_mass = 0;
    
  /** electrons */
  el_N        = 0;
  el_pdgId.clear();
  el_charge.clear();
  el_e.clear();
  el_eta.clear();
  el_phi.clear();
  el_mass.clear();
  el_pt.clear();
  el_et.clear();
  el_superCluster_eta.clear();
  el_pfRhoCorrRelIso03.clear();
  el_pfRhoCorrRelIso04.clear();
  el_pfDeltaCorrRelIso.clear();
  // el_pfRelIso.clear();
  el_photonIso.clear();
  el_neutralHadIso.clear();
  el_chargedHadIso.clear();
  // el_trackIso.clear();
  // el_pfRhoCorrRelIso03Boost.clear();
  // el_pfRhoCorrRelIso04Boost.clear();
  // el_pfDeltaCorrRelIsoBoost.clear();
  // el_pfRelIsoBoost.clear();
  // el_photonIsoBoost.clear();
  // el_neutralHadIsoBoost.clear();
  // el_chargedHadIsoBoost.clear();
  el_passConversionVeto.clear();
  el_full5x5_sigmaIetaIeta.clear();
  el_dEtaIn.clear();
  el_dPhiIn.clear();
  el_hOverE.clear();
  el_relIsoWithDBeta.clear();
  el_ooEmooP.clear();
  el_expectedMissingInnerHits.clear();
  el_d0.clear();
  el_dz.clear();
  el_d0_allvertices.clear();
  el_dz_allvertices.clear();
  el_dr03EcalRecHitSumEt.clear();
  el_dr03HcalDepth1TowerSumEt.clear();
  el_rho.clear();
  el_ecalDriven.clear();
  el_dEtaInSeed.clear();
  el_full5x5_e2x5Max.clear();
  el_full5x5_e5x5.clear();
  el_full5x5_e1x5.clear();
  el_full5x5_r9.clear();
  el_dr03TkSumPt.clear();
  el_superCluster_e.clear();
  el_hadronicOverEm.clear();
  el_seedEnergy.clear();
  el_isVetoElectron.clear();
  el_isLooseElectron.clear();
  el_isMediumElectron.clear();
  el_isTightElectron.clear();
  el_isHeepElectron.clear();
  el_isHltElectron.clear();
  el_isMVAMediumElectron.clear();
  el_isMVATightElectron.clear();
  el_MVAscore.clear();
  el_MVAcategory.clear();
//  el_isVetoElectronBoosted.clear();
//  el_isMediumElectronBoosted.clear();
//  el_isTightElectronBoosted.clear();
//  el_isHeepElectronBoosted.clear();
//  el_isLooseElectronBoosted.clear();
  el_isVetoElectronWithoutIPandIsolation.clear();
  el_isMediumElectronWithoutIPandIsolation.clear();
  el_isTightElectronWithoutIPandIsolation.clear();
  el_isHeepElectronWithoutIPandIsolation.clear();
  el_isLooseElectronWithoutIPandIsolation.clear();
  //  el_pfRhoCorrRelIso03Boost.clear();
  //  el_pfRhoCorrRelIso04Boost.clear();
  //  el_pfDeltaCorrRelIsoBoost.clear();
  //  el_pfRelIsoBoost.clear();
  //  el_photonIsoBoost.clear();
  //  el_neutralHadIsoBoost.clear();
  //  el_chargedHadIsoBoost.clear();
  el_SemileptonicPFIso.clear();
  //el_SemileptonicCorrPFIso.clear();
  
  /** muons */
  mu_N        = 0;
  mu_pdgId.clear();
  mu_charge.clear();
  mu_e.clear();
  mu_eta.clear();
  mu_phi.clear();
  mu_mass.clear();
  mu_pt.clear();
  mu_isHighPtMuon.clear();
  mu_isTightMuon.clear();
  mu_isMediumMuon.clear();
  mu_isMediumMuonGH.clear();
  mu_isLooseMuon.clear();
  mu_isPFMuon.clear();
  mu_isSoftMuon.clear();
  mu_isGlobalMuon.clear();
  mu_isTrackerMuon.clear();
  mu_isTrackerHighPtMuon.clear();
  // mu_pfRhoCorrRelIso03.clear();
  // mu_pfRhoCorrRelIso04.clear();
  mu_pfDeltaCorrRelIso.clear();
  // mu_pfRelIso.clear();
  mu_photonIso.clear();
  mu_neutralHadIso.clear();
  mu_chargedHadIso.clear();
  mu_trackIso.clear();
  mu_trackCorrIso.clear();
  mu_d0.clear();
  mu_dz.clear();
  mu_d0_allvertices.clear();
  mu_dz_allvertices.clear();
  mu_innerTrack_pt.clear();
  mu_bestTrack_pt.clear();
  mu_bestTrack_ptErr.clear();
  mu_tunePTrack_pt.clear();
  mu_tunePTrack_ptErr.clear();
  // mu_pfRhoCorrRelIso03Boost.clear();
  // mu_pfRhoCorrRelIso04Boost.clear();
  // mu_pfDeltaCorrRelIsoBoost.clear();
  // mu_pfRelIsoBoost.clear();
  // mu_photonIsoBoost.clear();
  // mu_neutralHadIsoBoost.clear();
  // mu_chargedHadIsoBoost.clear();
  mu_normChi2.clear();
  mu_trackerHits.clear();
  mu_matchedStations.clear();
  mu_pixelHits.clear();
  mu_globalHits.clear();
  mu_SemileptonicPFIso.clear();
  //mu_SemileptonicCorrPFIso.clear();

  /** taus */
  tau_N       = 0;
  tau_pdgId.clear();
  tau_charge.clear();
  tau_e.clear();
  tau_eta.clear();
  tau_phi.clear();
  tau_mass.clear();
  tau_pt.clear();
  tau_d0.clear();
  tau_dz.clear();
 
  // YT added
//  tau_associated_pdgId.clear();
//  tau_associated_pt.clear();
//  tau_associated_eta.clear();
//  tau_associated_phi.clear();
//  tau_associated_dr.clear();
//  tau_n_total.clear();
  tau_n_ch.clear();
  tau_n_nh.clear();
  //  tau_n_h_f.clear();
  //  tau_n_em_f.clear();
  tau_n_gamma.clear();
  //  tau_n_e.clear();
  //  tau_n_mu.clear();
  //  tau_n_total.clear();

  // tau_pfRhoCorrRelIso03.clear();
  // tau_pfRhoCorrRelIso04.clear();
  // tau_pfDeltaCorrRelIso.clear();
  // tau_pfRelIso.clear();
  // tau_photonIso.clear();
  // tau_neutralHadIso.clear();
  // tau_chargedHadIso.clear();
  tau_photonPtSumOutsideSignalCone.clear();
  // tau_trackIso.clear();
  // tau_pfRhoCorrRelIso03Boost.clear();
  // tau_pfRhoCorrRelIso04Boost.clear();
  // tau_pfDeltaCorrRelIsoBoost.clear();
  // tau_pfRelIsoBoost.clear();
  // tau_photonIsoBoost.clear();
  // tau_neutralHadIsoBoost.clear();
  // tau_chargedHadIsoBoost.clear();
  tau_TauType.clear();
  tau_decayMode.clear();
  tau_chargedPionPt.clear();
  tau_neutralPionPt.clear();

  tau_nPhoton.clear();  
  tau_ptWeightedDetaStrip.clear();
  tau_ptWeightedDphiStrip.clear();
  tau_ptWeightedDrSignal.clear();
  tau_ptWeightedDrIsolation.clear();
  tau_leadingTrackChi2.clear();
  tau_leadingTrackPt.clear();
  tau_eRatio.clear();
  tau_dxy_Sig.clear();
  tau_ip3d.clear();
  tau_ip3d_Sig.clear();
  tau_hasSecondaryVertex.clear();
  //  tau_decayDistMag_x.clear();
  //  tau_decayDistMag_y.clear();
  //  tau_decayDistMag_z.clear();
  tau_decayDistMag.clear();
  tau_flightLenthSig.clear();


  /** tau discriminants */
  tau_decayModeFindingNewDMs.clear();
  tau_decayModeFinding.clear();
  tau_byLooseCombinedIsolationDeltaBetaCorr3Hits.clear();
  tau_byMediumCombinedIsolationDeltaBetaCorr3Hits.clear();
  tau_byTightCombinedIsolationDeltaBetaCorr3Hits.clear();
  tau_byCombinedIsolationDeltaBetaCorrRaw3Hits.clear();
  tau_chargedIsoPtSum.clear();
  tau_neutralIsoPtSum.clear();
  tau_puCorrPtSum.clear();

  
  tau_chargedIsoPtSumdR03.clear();
  
  tau_neutralIsoPtSumdR03.clear();
 
  tau_photonPtSumOutsideSignalConedR03.clear();

  tau_byIsolationMVArun2v1DBdR03oldDMwLTraw.clear();
  tau_byIsolationMVArun2v1DBnewDMwLTraw.clear();
  tau_byIsolationMVArun2v1DBoldDMwLTraw.clear();
   tau_byIsolationMVArun2v1PWnewDMwLTraw.clear();
  
  tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT.clear();
  tau_byLooseIsolationMVArun2v1DBnewDMwLT.clear();
  tau_byLooseIsolationMVArun2v1DBoldDMwLT.clear();

  tau_byLooseIsolationMVArun2v1PWnewDMwLT.clear();
  
  
  tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT.clear();
  tau_byMediumIsolationMVArun2v1DBnewDMwLT.clear();
  tau_byMediumIsolationMVArun2v1DBoldDMwLT.clear();
 
  tau_byMediumIsolationMVArun2v1PWnewDMwLT.clear();
 

  tau_byTightIsolationMVArun2v1DBdR03oldDMwLT.clear();
  tau_byTightIsolationMVArun2v1DBnewDMwLT.clear();
  tau_byTightIsolationMVArun2v1DBoldDMwLT.clear();
 
  tau_byTightIsolationMVArun2v1PWnewDMwLT.clear();
 
  tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT.clear();

  tau_byVLooseIsolationMVArun2v1DBnewDMwLT.clear();
  tau_byVLooseIsolationMVArun2v1DBoldDMwLT.clear();
  tau_byVVLooseIsolationMVArun2v1DBoldDMwLT.clear();
   
 
  tau_byVLooseIsolationMVArun2v1PWnewDMwLT.clear();
 
  tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT.clear();
  tau_byVTightIsolationMVArun2v1DBnewDMwLT.clear();
  tau_byVTightIsolationMVArun2v1DBoldDMwLT.clear();

  
  tau_byVTightIsolationMVArun2v1PWnewDMwLT.clear();
 
  tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT.clear();
  tau_byVVTightIsolationMVArun2v1DBnewDMwLT.clear();
  tau_byVVTightIsolationMVArun2v1DBoldDMwLT.clear();
 
  tau_byVVTightIsolationMVArun2v1PWnewDMwLT.clear();

  

 
  tau_byIsolationMVArun2v2DBoldDMwLTraw.clear();
  tau_byVVLooseIsolationMVArun2v2DBoldDMwLT.clear();
  tau_byVLooseIsolationMVArun2v2DBoldDMwLT.clear();
  tau_byLooseIsolationMVArun2v2DBoldDMwLT.clear();
  tau_byMediumIsolationMVArun2v2DBoldDMwLT.clear();
  tau_byTightIsolationMVArun2v2DBoldDMwLT.clear();
  tau_byVTightIsolationMVArun2v2DBoldDMwLT.clear();
  tau_byVVTightIsolationMVArun2v2DBoldDMwLT.clear();
  

  tau_againstElectronMVA6raw.clear();
  tau_againstElectronMVA6category.clear();
  tau_againstElectronVLooseMVA6.clear();
  tau_againstElectronLooseMVA6.clear();
  tau_againstElectronMediumMVA6.clear();
  tau_againstElectronTightMVA6.clear();
  tau_againstElectronVTightMVA6.clear();
      

  tau_againstMuonLoose3.clear();
  tau_againstMuonTight3.clear(); 
      
  tau_byPhotonPtSumOutsideSignalCone.clear();
  //  tau_footprintCorrection.clear();



  
  /** energy density */
  rho = 0;
  
 

     
  

  /** HLT trigger decisions */
  HLT_isFired.clear();
  
  /** HLT trigger objects */
  triggerObject_pt.clear();
  triggerObject_eta.clear();
  triggerObject_phi.clear();
  triggerObject_mass.clear();
  triggerObject_lastname.clear();
  triggerObject_filterIDs.clear();
  triggerObject_filterLabels.clear();
  triggerObject_firedTrigger.clear();

  /** HLT filter decisions */
  passFilter_HBHE_                  = false;
  passFilter_HBHELoose_             = false;
  passFilter_HBHETight_             = false;
  passFilter_HBHEIso_               = false;
  passFilter_CSCHalo_               = false;
  passFilter_CSCTightHalo2015_      = false;
  passFilter_HCALlaser_             = false;
  passFilter_ECALDeadCell_          = false;
  passFilter_GoodVtx_               = false;
  passFilter_TrkFailure_            = false;
  passFilter_EEBadSc_               = false;
  passFilter_ECALlaser_             = false;
  passFilter_TrkPOG_                = false;
  passFilter_TrkPOG_manystrip_      = false;
  passFilter_TrkPOG_toomanystrip_   = false;
  passFilter_TrkPOG_logError_       = false;
  passFilter_METFilters_            = false;
   //NEW FOR ICHEP
  passFilter_CSCTightHaloTrkMuUnvetoFilter_   = false;
  passFilter_globalTightHalo2016_             = false;
  passFilter_globalSuperTightHalo2016_             = false;
  passFilter_HcalStripHalo_                   = false;
  passFilter_chargedHadronTrackResolution_    = false;
  passFilter_muonBadTrack_                    = false;
  flag_badMuons_                    = false;
  flag_duplicateMuons_              = false;
  flag_nobadMuons_                  = false;

  /** MET */
  METraw_et.clear();	 
  METraw_phi.clear();
  METraw_sumEt.clear();
  MET_corrPx.clear();
  MET_corrPy.clear();
  MET_et.clear();
  MET_phi.clear();
  MET_puppi_et.clear();
  MET_puppi_phi.clear();

  MET_sumEt.clear();
  MET_T1Uncertainty.clear();
  
  MET_JetEnUp.clear();
  MET_JetEnDown.clear();
  MET_JetResUp.clear();
  MET_JetResDown.clear();
  MET_UnclusteredEnUp.clear();
  MET_UnclusteredEnDown.clear();

  /** MET SVift*/
  MET_significance.clear();
  MET_cov00.clear();
  MET_cov10.clear();
  MET_cov11.clear();
  MET_mva_et.clear();
  MET_mva_phi.clear();
  MET_mva_cov00.clear();
  MET_mva_cov10.clear();
  MET_mva_cov11.clear();
  MET_mva_recoil_pt.clear();
  MET_mva_recoil_eta.clear();
  MET_mva_recoil_phi.clear();
  MET_mva_recoil_pdgId.clear();
  MET_Nmva.clear();

  /*------------------------EVENT infos-------------------------*/    
  EVENT_event = 0;
  EVENT_run = 0;
  EVENT_lumiBlock = 0;

  /*--------------------------PV infos--------------------------*/
  PV_N = 0;
  PV_filter = true;
  PV_chi2.clear();
  PV_ndof.clear();
  PV_rho.clear();
  PV_z.clear();
  BeamSpot_x0.clear();
  BeamSpot_y0.clear();
  BeamSpot_z0.clear();
  /*--------------------------PU infos--------------------------*/  			       
  nPuVtxTrue.clear();
  nPuVtx.clear();
  bX.clear();

  /*-------------------------JPSI infos--------------------------*/ 
  IsJpsiMu.clear();
  IsJpsiEle.clear();

  Jpsi_mu1_pt.clear();
  Jpsi_mu1_eta.clear();
  Jpsi_mu1_phi.clear();
  Jpsi_mu1_ch.clear();
  Jpsi_mu1_isLoose.clear();
  Jpsi_mu1_isTight.clear();
  Jpsi_mu1_isPF.clear();
  Jpsi_mu1_isGlobal.clear();
  Jpsi_mu1_isTracker.clear();
  Jpsi_mu1_isSoft.clear();

  Jpsi_mu2_pt.clear();
  Jpsi_mu2_eta.clear();
  Jpsi_mu2_phi.clear();
  Jpsi_mu2_ch.clear();
  Jpsi_mu2_isLoose.clear();
  Jpsi_mu2_isTight.clear();
  Jpsi_mu2_isPF.clear();
  Jpsi_mu2_isGlobal.clear();
  Jpsi_mu2_isTracker.clear();
  Jpsi_mu2_isSoft.clear();

  Jpsi_mu3_pt.clear();
  Jpsi_mu3_eta.clear();
  Jpsi_mu3_phi.clear();
  Jpsi_mu3_ch.clear();
  Jpsi_mu3_isLoose.clear();
  Jpsi_mu3_isTight.clear();
  Jpsi_mu3_isPF.clear();
  Jpsi_mu3_isGlobal.clear();
  Jpsi_mu3_isTracker.clear();
  Jpsi_mu3_isSoft.clear();
  Jpsi_mu3_x.clear();
  Jpsi_mu3_y.clear();
  Jpsi_mu3_z.clear();


  Jpsi_mu3_isopt03.clear();
  Jpsi_mu3_isopt04.clear();
  Jpsi_mu3_isopt05.clear();
  Jpsi_mu3_isopt06.clear();
  Jpsi_mu3_isopt07.clear();


  Jpsi_dx.clear();
  Jpsi_dy.clear();
  Jpsi_dz.clear();
  Jpsi_pt.clear();
  Jpsi_eta.clear();
  Jpsi_phi.clear();
  Jpsi_mass.clear();
  Jpsi_vtxprob.clear();
  Jpsi_vtxz.clear();

  Jpsi_trimu_dx.clear();
  Jpsi_trimu_dy.clear();
  Jpsi_trimu_dz.clear();
  Jpsi_trimu_pt.clear();
  Jpsi_trimu_eta.clear();
  Jpsi_trimu_phi.clear();
  Jpsi_trimu_mass.clear();
  Jpsi_trimu_vtxprob.clear();
  Jpsi_trimu_vtxz.clear();

  Jpsi_PV_x.clear();
  Jpsi_PV_y.clear();
  Jpsi_PV_z.clear();

  Jpsi_flightSig3D.clear();
  Jpsi_flightLength3D.clear();
  Jpsi_flightLengthErr3D.clear();
  Jpsi_flightSig2D.clear();
  Jpsi_flightLength2D.clear();
  Jpsi_flightLengthErr2D.clear();
  Jpsi_trimu_flightSig3D.clear(); 
  Jpsi_trimu_flightLength3D.clear();
  Jpsi_trimu_flightLengthErr3D.clear();
  Jpsi_trimu_flightSig2D.clear(); 
  Jpsi_trimu_flightLength2D.clear();
  Jpsi_trimu_flightLengthErr2D.clear();
  Jpsi_dr_mu3pf.clear();
 
} 
