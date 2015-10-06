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
      tree_->Branch( "genParticle_px"	     , &genParticle_px	       ); 
      tree_->Branch( "genParticle_py"	     , &genParticle_py	       ); 
      tree_->Branch( "genParticle_pz"	     , &genParticle_pz	       ); 
      tree_->Branch( "genParticle_e" 	     , &genParticle_e	       ); 
      tree_->Branch( "genParticle_eta"	     , &genParticle_eta        ); 
      tree_->Branch( "genParticle_phi"	     , &genParticle_phi        ); 
      tree_->Branch( "genParticle_mass"	     , &genParticle_mass       ); 
      tree_->Branch( "genParticle_pdgId"     , &genParticle_pdgId      );
      tree_->Branch( "genParticle_status"    , &genParticle_status     );
      tree_->Branch( "genParticle_mother"    , &genParticle_mother     );
      tree_->Branch( "genParticle_nMoth"     , &genParticle_nMoth      );
      tree_->Branch( "genParticle_nDau"	     , &genParticle_nDau       ); 
      tree_->Branch( "genParticle_dau"	     , &genParticle_dau        );
    } //doGenParticles
    
    if ( runFlags["doGenEvent"] ){
      /** generator info */
      tree_->Branch( "lheV_pt"	             , &lheV_pt                ); 
      tree_->Branch( "lheHT"	             , &lheHT                  ); 
      tree_->Branch( "lheNj"	             , &lheNj                  ); 
      tree_->Branch( "genWeight"	     , &genWeight              ); 
      tree_->Branch( "qScale"	             , &qScale                 );
      tree_->Branch( "PDF_x"	             , &PDF_x                  );
      tree_->Branch( "PDF_xPDF"	             , &PDF_xPDF               );
      tree_->Branch( "PDF_id"	             , &PDF_id                 );
    } //doGenEvent
  } //runOnMC
  
  if ( runFlags["doElectrons"] ){
    /** electrons */
    tree_->Branch( "el_N"  		     , &el_N	            );
    tree_->Branch( "el_pdgId"  		     , &el_pdgId            );
    tree_->Branch( "el_charge"		     , &el_charge           );
    tree_->Branch( "el_e" 		     , &el_e	            );
    tree_->Branch( "el_eta"		     , &el_eta 	            );
    tree_->Branch( "el_phi"		     , &el_phi 	            );
    tree_->Branch( "el_mass"		     , &el_mass	            );
    tree_->Branch( "el_pt"		     , &el_pt  	            );
    tree_->Branch( "el_et"		     , &el_et  	            );
    tree_->Branch( "el_superCluster_eta"     , &el_superCluster_eta );
  
    tree_->Branch( "el_pfRhoCorrRelIso03"    , &el_pfRhoCorrRelIso03 );
    tree_->Branch( "el_pfRhoCorrRelIso04"    , &el_pfRhoCorrRelIso04 );
    tree_->Branch( "el_pfDeltaCorrRelIso"    , &el_pfDeltaCorrRelIso );
    tree_->Branch( "el_pfRelIso"  	     , &el_pfRelIso	     );
    tree_->Branch( "el_photonIso" 	     , &el_photonIso	     );
    tree_->Branch( "el_neutralHadIso"	     , &el_neutralHadIso     );
    tree_->Branch( "el_chargedHadIso"	     , &el_chargedHadIso     );
    tree_->Branch( "el_trackIso"	     , &el_trackIso	     );
  
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
  
    tree_->Branch( "el_dr03EcalRecHitSumEt"     , &el_dr03EcalRecHitSumEt      );
    tree_->Branch( "el_dr03HcalDepth1TowerSumEt", &el_dr03HcalDepth1TowerSumEt );
    tree_->Branch( "el_rho"                     , &el_rho                      );
    tree_->Branch( "el_ecalDriven"              , &el_ecalDriven               );
    tree_->Branch( "el_dEtaInSeed"              , &el_dEtaInSeed               );
    tree_->Branch( "el_full5x5_e2x5Max"         , &el_full5x5_e2x5Max          );
    tree_->Branch( "el_full5x5_e5x5"            , &el_full5x5_e5x5             );
    tree_->Branch( "el_full5x5_e1x5"            , &el_full5x5_e1x5             );
    tree_->Branch( "el_dr03TkSumPt"             , &el_dr03TkSumPt              );
    tree_->Branch( "el_superCluster_e"          , &el_superCluster_e           );
    tree_->Branch( "el_hadronicOverEm"          , &el_hadronicOverEm           );
  
    tree_->Branch( "el_isVetoElectron"	        , &el_isVetoElectron           );
    tree_->Branch( "el_isMediumElectron"	, &el_isMediumElectron         );
    tree_->Branch( "el_isTightElectron"         , &el_isTightElectron          );  
    tree_->Branch( "el_isHeepElectron"	        , &el_isHeepElectron           );
    tree_->Branch( "el_isHeep51Electron"	, &el_isHeep51Electron         );
    tree_->Branch( "el_isLooseElectron"	        , &el_isLooseElectron          );
  
    tree_->Branch( "el_isVetoElectronBoosted"	, &el_isVetoElectronBoosted           );
    tree_->Branch( "el_isMediumElectronBoosted"	, &el_isMediumElectronBoosted         );
    tree_->Branch( "el_isTightElectronBoosted"  , &el_isTightElectronBoosted          );  
    tree_->Branch( "el_isHeepElectronBoosted"   , &el_isHeepElectronBoosted           );
    tree_->Branch( "el_isHeep51ElectronBoosted"	, &el_isHeep51ElectronBoosted         );
    tree_->Branch( "el_isLooseElectronBoosted"  , &el_isLooseElectronBoosted          );
  
    tree_->Branch( "el_pfRhoCorrRelIso03Boost"  , &el_pfRhoCorrRelIso03Boost   );
    tree_->Branch( "el_pfRhoCorrRelIso04Boost"  , &el_pfRhoCorrRelIso04Boost   );
    tree_->Branch( "el_pfDeltaCorrRelIsoBoost"  , &el_pfDeltaCorrRelIsoBoost   );
    tree_->Branch( "el_pfRelIsoBoost"  	        , &el_pfRelIsoBoost	       );
    tree_->Branch( "el_photonIsoBoost" 	        , &el_photonIsoBoost	       );
    tree_->Branch( "el_neutralHadIsoBoost"      , &el_neutralHadIsoBoost       );
    tree_->Branch( "el_chargedHadIsoBoost"      , &el_chargedHadIsoBoost       );
  
    tree_->Branch( "el_SemileptonicPFIso"       , &el_SemileptonicPFIso      );
    tree_->Branch( "el_SemileptonicCorrPFIso"   , &el_SemileptonicCorrPFIso  );
  } //doElectrons
  
  if ( runFlags["doMuons"] ){
    /** muons */
    tree_->Branch( "mu_N"  		     , &mu_N	       );
    tree_->Branch( "mu_pdgId"  		     , &mu_pdgId       );
    tree_->Branch( "mu_charge"		     , &mu_charge      );
    tree_->Branch( "mu_e" 		     , &mu_e	       );
    tree_->Branch( "mu_eta"		     , &mu_eta 	       );
    tree_->Branch( "mu_phi"		     , &mu_phi 	       );
    tree_->Branch( "mu_mass"		     , &mu_mass	       );
    tree_->Branch( "mu_pt"		     , &mu_pt  	       );
  
    tree_->Branch( "mu_isHighPtMuon"	     , &mu_isHighPtMuon       );
    tree_->Branch( "mu_isTightMuon"	     , &mu_isTightMuon        );
    tree_->Branch( "mu_isLooseMuon"	     , &mu_isLooseMuon        );
    tree_->Branch( "mu_isPFMuon"	     , &mu_isPFMuon           );
    tree_->Branch( "mu_isSoftMuon"	     , &mu_isSoftMuon	      );
  
    tree_->Branch( "mu_pfRhoCorrRelIso03"    , &mu_pfRhoCorrRelIso03  );
    tree_->Branch( "mu_pfRhoCorrRelIso04"    , &mu_pfRhoCorrRelIso04  );
    tree_->Branch( "mu_pfDeltaCorrRelIso"    , &mu_pfDeltaCorrRelIso  );
    tree_->Branch( "mu_pfRelIso"  	     , &mu_pfRelIso	      );
    tree_->Branch( "mu_photonIso" 	     , &mu_photonIso	      );
    tree_->Branch( "mu_neutralHadIso"	     , &mu_neutralHadIso      );
    tree_->Branch( "mu_chargedHadIso"	     , &mu_chargedHadIso      );
    tree_->Branch( "mu_trackIso"	     , &mu_trackIso	      );
    tree_->Branch( "mu_d0"                   , &mu_d0                 );
    tree_->Branch( "mu_bestTrack_pt"	     , &mu_bestTrack_pt       );
    tree_->Branch( "mu_bestTrack_ptErr"	     , &mu_bestTrack_ptErr    );
  
    tree_->Branch( "mu_pfRhoCorrRelIso03Boost", &mu_pfRhoCorrRelIso03Boost );
    tree_->Branch( "mu_pfRhoCorrRelIso04Boost", &mu_pfRhoCorrRelIso04Boost );
    tree_->Branch( "mu_pfDeltaCorrRelIsoBoost", &mu_pfDeltaCorrRelIsoBoost );
    tree_->Branch( "mu_pfRelIsoBoost"  	      , &mu_pfRelIsoBoost	   );
    tree_->Branch( "mu_photonIsoBoost" 	      , &mu_photonIsoBoost	   );
    tree_->Branch( "mu_neutralHadIsoBoost"    , &mu_neutralHadIsoBoost     );
    tree_->Branch( "mu_chargedHadIsoBoost"    , &mu_chargedHadIsoBoost     );
  
    tree_->Branch( "mu_normChi2"  	      , &mu_normChi2	    );
    tree_->Branch( "mu_isGlobalMuon"	      , &mu_isGlobalMuon    );
    tree_->Branch( "mu_trackerHits"	      , &mu_trackerHits     );
    tree_->Branch( "mu_matchedStations"	      , &mu_matchedStations );
    tree_->Branch( "mu_pixelHits" 	      , &mu_pixelHits	    );
    tree_->Branch( "mu_globalHits"	      , &mu_globalHits      );
  
    tree_->Branch( "mu_SemileptonicPFIso"       , &mu_SemileptonicPFIso      );
    tree_->Branch( "mu_SemileptonicCorrPFIso"   , &mu_SemileptonicCorrPFIso  );
  } //doMuons
  
  if ( runFlags["doTaus"] ){
    /** taus */
    tree_->Branch( "tau_N"  		     , &tau_N		       );
    tree_->Branch( "tau_pdgId"               , &tau_pdgId              );
    tree_->Branch( "tau_charge"		     , &tau_charge	       );
    tree_->Branch( "tau_e" 		     , &tau_e		       );
    tree_->Branch( "tau_eta"		     , &tau_eta 	       );
    tree_->Branch( "tau_phi"		     , &tau_phi 	       );
    tree_->Branch( "tau_mass"		     , &tau_mass	       );
    tree_->Branch( "tau_pt"		     , &tau_pt  	       );
    
    tree_->Branch( "tau_pfRhoCorrRelIso03"   , &tau_pfRhoCorrRelIso03  );
    tree_->Branch( "tau_pfRhoCorrRelIso04"   , &tau_pfRhoCorrRelIso04  );
    tree_->Branch( "tau_pfDeltaCorrRelIso"   , &tau_pfDeltaCorrRelIso  );
    tree_->Branch( "tau_pfRelIso"  	     , &tau_pfRelIso	       );
    tree_->Branch( "tau_photonIso" 	     , &tau_photonIso	       );
    tree_->Branch( "tau_neutralHadIso"	     , &tau_neutralHadIso      );
    tree_->Branch( "tau_chargedHadIso"	     , &tau_chargedHadIso      );
    tree_->Branch( "tau_trackIso"	     , &tau_trackIso	       );

    tree_->Branch( "tau_d0"                  , &tau_d0                 );

    tree_->Branch( "tau_pfRhoCorrRelIso03Boost"  , &tau_pfRhoCorrRelIso03Boost );
    tree_->Branch( "tau_pfRhoCorrRelIso04Boost"  , &tau_pfRhoCorrRelIso04Boost );
    tree_->Branch( "tau_pfDeltaCorrRelIsoBoost"  , &tau_pfDeltaCorrRelIsoBoost );
    tree_->Branch( "tau_pfRelIsoBoost"  	 , &tau_pfRelIsoBoost	       );
    tree_->Branch( "tau_photonIsoBoost" 	 , &tau_photonIsoBoost         );
    tree_->Branch( "tau_neutralHadIsoBoost"      , &tau_neutralHadIsoBoost     );
    tree_->Branch( "tau_chargedHadIsoBoost"      , &tau_chargedHadIsoBoost     );
  
    tree_->Branch( "tau_TauType"		 , &tau_TauType		       );
  
    if ( runFlags["doTausBoosted"] ){
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
      tree_->Branch( "tau_byIsolationMVA3oldDMwoLTraw"		      , &tau_byIsolationMVA3oldDMwoLTraw 		 );
      tree_->Branch( "tau_byVLooseIsolationMVA3oldDMwoLT"	      , &tau_byVLooseIsolationMVA3oldDMwoLT	         );
      tree_->Branch( "tau_byLooseIsolationMVA3oldDMwoLT"	      , &tau_byLooseIsolationMVA3oldDMwoLT	         );
      tree_->Branch( "tau_byMediumIsolationMVA3oldDMwoLT"	      , &tau_byMediumIsolationMVA3oldDMwoLT	         );
      tree_->Branch( "tau_byTightIsolationMVA3oldDMwoLT"	      , &tau_byTightIsolationMVA3oldDMwoLT	         );
      tree_->Branch( "tau_byVTightIsolationMVA3oldDMwoLT"	      , &tau_byVTightIsolationMVA3oldDMwoLT	         );
      tree_->Branch( "tau_byVVTightIsolationMVA3oldDMwoLT"	      , &tau_byVVTightIsolationMVA3oldDMwoLT	         );
      tree_->Branch( "tau_byIsolationMVA3oldDMwLTraw"		      , &tau_byIsolationMVA3oldDMwLTraw  		 );
      tree_->Branch( "tau_byVLooseIsolationMVA3oldDMwLT"	      , &tau_byVLooseIsolationMVA3oldDMwLT		 );
      tree_->Branch( "tau_byLooseIsolationMVA3oldDMwLT"		      , &tau_byLooseIsolationMVA3oldDMwLT		 );
      tree_->Branch( "tau_byMediumIsolationMVA3oldDMwLT"	      , &tau_byMediumIsolationMVA3oldDMwLT		 );
      tree_->Branch( "tau_byTightIsolationMVA3oldDMwLT"		      , &tau_byTightIsolationMVA3oldDMwLT		 );
      tree_->Branch( "tau_byVTightIsolationMVA3oldDMwLT"	      , &tau_byVTightIsolationMVA3oldDMwLT	         );
      tree_->Branch( "tau_byVVTightIsolationMVA3oldDMwLT"	      , &tau_byVVTightIsolationMVA3oldDMwLT	         );
      tree_->Branch( "tau_byIsolationMVA3newDMwoLTraw"		      , &tau_byIsolationMVA3newDMwoLTraw 		 );
      tree_->Branch( "tau_byVLooseIsolationMVA3newDMwoLT"	      , &tau_byVLooseIsolationMVA3newDMwoLT	         );
      tree_->Branch( "tau_byLooseIsolationMVA3newDMwoLT"	      , &tau_byLooseIsolationMVA3newDMwoLT	         );
      tree_->Branch( "tau_byMediumIsolationMVA3newDMwoLT"	      , &tau_byMediumIsolationMVA3newDMwoLT	         );
      tree_->Branch( "tau_byTightIsolationMVA3newDMwoLT"	      , &tau_byTightIsolationMVA3newDMwoLT	         );
      tree_->Branch( "tau_byVTightIsolationMVA3newDMwoLT"	      , &tau_byVTightIsolationMVA3newDMwoLT	         );
      tree_->Branch( "tau_byVVTightIsolationMVA3newDMwoLT"	      , &tau_byVVTightIsolationMVA3newDMwoLT	         );
      tree_->Branch( "tau_byIsolationMVA3newDMwLTraw"		      , &tau_byIsolationMVA3newDMwLTraw  		 );
      tree_->Branch( "tau_byVLooseIsolationMVA3newDMwLT"	      , &tau_byVLooseIsolationMVA3newDMwLT		 );
      tree_->Branch( "tau_byLooseIsolationMVA3newDMwLT"		      , &tau_byLooseIsolationMVA3newDMwLT		 );
      tree_->Branch( "tau_byMediumIsolationMVA3newDMwLT"	      , &tau_byMediumIsolationMVA3newDMwLT		 );
      tree_->Branch( "tau_byTightIsolationMVA3newDMwLT"		      , &tau_byTightIsolationMVA3newDMwLT		 );
      tree_->Branch( "tau_byVTightIsolationMVA3newDMwLT"	      , &tau_byVTightIsolationMVA3newDMwLT	         );
      tree_->Branch( "tau_byVVTightIsolationMVA3newDMwLT"	      , &tau_byVVTightIsolationMVA3newDMwLT	         );
      tree_->Branch( "tau_againstElectronLoose"			      , &tau_againstElectronLoose			 );
      tree_->Branch( "tau_againstElectronMedium"		      , &tau_againstElectronMedium			 );
      tree_->Branch( "tau_againstElectronTight"			      , &tau_againstElectronTight			 );
      tree_->Branch( "tau_againstElectronMVA5raw"		      , &tau_againstElectronMVA5raw			 );
      tree_->Branch( "tau_againstElectronMVA5category"		      , &tau_againstElectronMVA5category 		 );
      tree_->Branch( "tau_againstElectronVLooseMVA5" 		      , &tau_againstElectronVLooseMVA5			 );
      tree_->Branch( "tau_againstElectronLooseMVA5"  		      , &tau_againstElectronLooseMVA5			 );
      tree_->Branch( "tau_againstElectronMediumMVA5" 		      , &tau_againstElectronMediumMVA5			 );
      tree_->Branch( "tau_againstElectronTightMVA5"  		      , &tau_againstElectronTightMVA5			 );
      tree_->Branch( "tau_againstElectronVTightMVA5" 		      , &tau_againstElectronVTightMVA5			 );
      tree_->Branch( "tau_againstMuonLoose"  			      , &tau_againstMuonLoose				 );
      tree_->Branch( "tau_againstMuonMedium" 			      , &tau_againstMuonTight				 );
      tree_->Branch( "tau_againstMuonTight"  			      , &tau_againstMuonTight				 );
      tree_->Branch( "tau_againstMuonLoose2" 			      , &tau_againstMuonLoose2				 );
      tree_->Branch( "tau_againstMuonMedium2"			      , &tau_againstMuonMedium2  			 );
      tree_->Branch( "tau_againstMuonTight2" 			      , &tau_againstMuonLoose3				 );
      tree_->Branch( "tau_againstMuonLoose3" 			      , &tau_againstMuonLoose3				 );
      tree_->Branch( "tau_againstMuonTight3" 			      , &tau_againstMuonTight3				 );
      tree_->Branch( "tau_againstMuonMVAraw" 			      , &tau_againstMuonMVAraw				 );
      tree_->Branch( "tau_againstMuonLooseMVA"			      , &tau_againstMuonLooseMVA 			 );
      tree_->Branch( "tau_againstMuonMediumMVA"			      , &tau_againstMuonMediumMVA			 );
      tree_->Branch( "tau_againstMuonTightMVA"			      , &tau_againstMuonTightMVA 			 );
    } //doTausBoosted
  } //doTaus
      

  if (runFlags["doAK4Jets"] || runFlags["doAK8Jets"]) {
    /** energy density */
    tree_->Branch( "rho", &rho );
  }

  if (runFlags["doAK4Jets"]) {
    /** AK4 jets */

    tree_->Branch( "jetAK4_N"		    , &jetAK4_N 	 );
    tree_->Branch( "jetAK4_pt"		    , &jetAK4_pt	 );
    tree_->Branch( "jetAK4_eta"		    , &jetAK4_eta	 );
    tree_->Branch( "jetAK4_mass"	    , &jetAK4_mass	 );
    tree_->Branch( "jetAK4_phi"		    , &jetAK4_phi	 );
    tree_->Branch( "jetAK4_e"		    , &jetAK4_e 	 );
    tree_->Branch( "jetAK4_jec"		    , &jetAK4_jec 	 );
    tree_->Branch( "jetAK4_IDLoose"	    , &jetAK4_IDLoose	 );
    tree_->Branch( "jetAK4_IDTight"	    , &jetAK4_IDTight	 );
    tree_->Branch( "jetAK4_muf" 	    , &jetAK4_muf	 );
    tree_->Branch( "jetAK4_phf" 	    , &jetAK4_phf	 );
    tree_->Branch( "jetAK4_emf" 	    , &jetAK4_emf	 );
    tree_->Branch( "jetAK4_nhf" 	    , &jetAK4_nhf	 );
    tree_->Branch( "jetAK4_chf" 	    , &jetAK4_chf	 );
    tree_->Branch( "jetAK4_area" 	    , &jetAK4_area	 );
    tree_->Branch( "jetAK4_cm"     	    , &jetAK4_cm	 );
    tree_->Branch( "jetAK4_nm"     	    , &jetAK4_nm	 );
    tree_->Branch( "jetAK4_che"     	    , &jetAK4_che	 );
    tree_->Branch( "jetAK4_ne"     	    , &jetAK4_ne	 );
    tree_->Branch( "jetAK4_hf_hf"  	    , &jetAK4_hf_hf      );
    tree_->Branch( "jetAK4_hf_emf"  	    , &jetAK4_hf_emf     );
    tree_->Branch( "jetAK4_hof"  	    , &jetAK4_hof        );
    tree_->Branch( "jetAK4_chm"  	    , &jetAK4_chm        );
    tree_->Branch( "jetAK4_neHadMult"  	    , &jetAK4_neHadMult  );
    tree_->Branch( "jetAK4_phoMult"  	    , &jetAK4_phoMult    );
    tree_->Branch( "jetAK4_nemf"  	    , &jetAK4_nemf       );
    tree_->Branch( "jetAK4_cemf"  	    , &jetAK4_cemf       );
    tree_->Branch( "jetAK4_charge" 	    , &jetAK4_charge     );
  //tree_->Branch( "jetAK4_ssv"             , &jetAK4_ssv        );
    tree_->Branch( "jetAK4_cisv"	    , &jetAK4_cisv	 );
  //tree_->Branch( "jetAK4_tchp"            , &jetAK4_tchp       );
  //tree_->Branch( "jetAK4_tche"            , &jetAK4_tche       );
  //tree_->Branch( "jetAK4_jp"              , &jetAK4_jp         );
  //tree_->Branch( "jetAK4_jbp"             , &jetAK4_jbp        );
    tree_->Branch( "jetAK4_vtxMass"	    , &jetAK4_vtxMass    );
    tree_->Branch( "jetAK4_vtxNtracks"      , &jetAK4_vtxNtracks );
    tree_->Branch( "jetAK4_vtx3DVal"	    , &jetAK4_vtx3DVal   );
    tree_->Branch( "jetAK4_vtx3DSig"	    , &jetAK4_vtx3DSig   );
  //tree_->Branch( "jetAK4_nSVs"	    , &jetAK4_nSVs	 );
    if ( runFlags["runOnMC"] ){
      tree_->Branch( "jetAK4_partonFlavour" , &jetAK4_partonFlavour);
      tree_->Branch( "jetAK4_hadronFlavour" , &jetAK4_hadronFlavour);
      tree_->Branch( "jetAK4_genParton_pdgID" , &jetAK4_genParton_pdgID);
      tree_->Branch( "jetAK4_nbHadrons"       , &jetAK4_nbHadrons);
      tree_->Branch( "jetAK4_ncHadrons"       , &jetAK4_ncHadrons);
    }

  } //doAK4Jets
        
  if (runFlags["doAK8Jets"]) {
    /** AK8 jets */

    tree_->Branch( "jetAK8_N"		     , &jetAK8_N 		 );
    tree_->Branch( "jetAK8_pt"		     , &jetAK8_pt		 );
    tree_->Branch( "jetAK8_eta"		     , &jetAK8_eta		 );
    tree_->Branch( "jetAK8_mass"	     , &jetAK8_mass		 );
    tree_->Branch( "jetAK8_phi"		     , &jetAK8_phi		 );
    tree_->Branch( "jetAK8_e"		     , &jetAK8_e 		 );
    tree_->Branch( "jetAK8_jec"		     , &jetAK8_jec 		 );
    tree_->Branch( "jetAK8_IDLoose"	     , &jetAK8_IDLoose		 );
    tree_->Branch( "jetAK8_IDTight"	     , &jetAK8_IDTight           );
    tree_->Branch( "jetAK8_muf" 	     , &jetAK8_muf		 );
    tree_->Branch( "jetAK8_phf" 	     , &jetAK8_phf		 );
    tree_->Branch( "jetAK8_emf" 	     , &jetAK8_emf		 );
    tree_->Branch( "jetAK8_nhf" 	     , &jetAK8_nhf		 );
    tree_->Branch( "jetAK8_chf" 	     , &jetAK8_chf		 );
    tree_->Branch( "jetAK8_area" 	     , &jetAK8_area		 );
    tree_->Branch( "jetAK8_cm"     	     , &jetAK8_cm		 );
    tree_->Branch( "jetAK8_nm"     	     , &jetAK8_nm		 );
    tree_->Branch( "jetAK8_che"     	     , &jetAK8_che		 );
    tree_->Branch( "jetAK8_ne"     	     , &jetAK8_ne		 );
    tree_->Branch( "jetAK8_hf_hf"  	     , &jetAK8_hf_hf             );
    tree_->Branch( "jetAK8_hf_emf"  	     , &jetAK8_hf_emf            );
    tree_->Branch( "jetAK8_hof"  	     , &jetAK8_hof               );
    tree_->Branch( "jetAK8_chm"  	     , &jetAK8_chm               );
    tree_->Branch( "jetAK8_neHadMult"  	     , &jetAK8_neHadMult         );
    tree_->Branch( "jetAK8_phoMult"  	     , &jetAK8_phoMult           );
    tree_->Branch( "jetAK8_nemf"  	     , &jetAK8_nemf              );
    tree_->Branch( "jetAK8_cemf"  	     , &jetAK8_cemf              );
    tree_->Branch( "jetAK8_charge" 	     , &jetAK8_charge    	 );
    if ( runFlags["runOnMC"] ){
      tree_->Branch( "jetAK8_partonFlavour" , &jetAK8_partonFlavour);
      tree_->Branch( "jetAK8_hadronFlavour" , &jetAK8_hadronFlavour);
      tree_->Branch( "jetAK8_genParton_pdgID" , &jetAK8_genParton_pdgID);
      tree_->Branch( "jetAK8_nbHadrons"       , &jetAK8_nbHadrons);
      tree_->Branch( "jetAK8_ncHadrons"       , &jetAK8_ncHadrons);
    }
    
    if (runFlags["doHbbTag"]) {
      tree_->Branch( "jetAK8_Hbbtag"	     , &jetAK8_Hbbtag		 );
    }
    //tree_->Branch( "jetAK8_ssv"            , &jetAK8_ssv               );
    tree_->Branch( "jetAK8_csv"		     , &jetAK8_csv		 );
  //tree_->Branch( "jetAK8_tchp"             , &jetAK8_tchp              );
  //tree_->Branch( "jetAK8_tche"             , &jetAK8_tche              );
  //tree_->Branch( "jetAK8_jp"               , &jetAK8_jp                );
  //tree_->Branch( "jetAK8_jbp"              , &jetAK8_jbp               );
    tree_->Branch( "jetAK8_tau1"	     , &jetAK8_tau1		 );
    tree_->Branch( "jetAK8_tau2"	     , &jetAK8_tau2      	 );
    tree_->Branch( "jetAK8_tau3"	     , &jetAK8_tau3	    	 );
    tree_->Branch( "jetAK8_pruned_mass"      , &jetAK8_pruned_mass       );
    tree_->Branch( "jetAK8_pruned_massCorr"  , &jetAK8_pruned_massCorr	 );
    tree_->Branch( "jetAK8_pruned_jec"	     , &jetAK8_pruned_jec        );
    tree_->Branch( "jetAK8_softdrop_mass"    , &jetAK8_softdrop_mass     ); 
    tree_->Branch( "jetAK8_softdrop_massCorr", &jetAK8_softdrop_massCorr );
    tree_->Branch( "jetAK8_softdrop_jec"     , &jetAK8_softdrop_jec      );

    if (runFlags["doTrimming"]) {
      /*----------------------AK10 trimming ---------------------------*/   
     tree_->Branch( "jetAK10_trimmed_mass"    , &jetAK10_trimmed_mass     );
     tree_->Branch( "jetAK10_trimmed_massCorr", &jetAK10_trimmed_massCorr );
     tree_->Branch( "jetAK10_trimmed_jec"     , &jetAK10_trimmed_jec      );
     tree_->Branch( "jetAK10_ecf1"	     , &jetAK10_ecf1		 );
     tree_->Branch( "jetAK10_ecf2"	     , &jetAK10_ecf2      	 );
     tree_->Branch( "jetAK10_ecf3"	     , &jetAK10_ecf3	    	 );
    }

    if (runFlags["doPuppi"]) {
      /*----------------------PUPPI ---------------------------*/   
     tree_->Branch( "jetAK8_puppi_pruned_mass"      , &jetAK8_puppi_pruned_mass       );
     tree_->Branch( "jetAK8_puppi_pruned_massCorr"  , &jetAK8_puppi_pruned_massCorr	 );
     tree_->Branch( "jetAK8_puppi_pruned_jec"	     , &jetAK8_puppi_pruned_jec        );
     tree_->Branch( "jetAK8_puppi_softdrop_mass"    , &jetAK8_puppi_softdrop_mass     ); 
     tree_->Branch( "jetAK8_puppi_softdrop_massCorr", &jetAK8_puppi_softdrop_massCorr );
     tree_->Branch( "jetAK8_puppi_softdrop_jec"     , &jetAK8_puppi_softdrop_jec      );
     tree_->Branch( "jetAK8_puppi_tau1"	     , &jetAK8_puppi_tau1		 );
     tree_->Branch( "jetAK8_puppi_tau2"	     , &jetAK8_puppi_tau2      	 );
     tree_->Branch( "jetAK8_puppi_tau3"	     , &jetAK8_puppi_tau3	    	 );
    }
      // /*----------------------Softdrop AK8 subjets---------------------------*/
      tree_->Branch( "subjetAK8_softdrop_N"      , &subjetAK8_softdrop_N           );
      tree_->Branch( "subjetAK8_softdrop_pt"     , &subjetAK8_softdrop_pt      );
      tree_->Branch( "subjetAK8_softdrop_eta"    , &subjetAK8_softdrop_eta     );
      tree_->Branch( "subjetAK8_softdrop_mass"   , &subjetAK8_softdrop_mass    );
      tree_->Branch( "subjetAK8_softdrop_phi"    , &subjetAK8_softdrop_phi     );
      tree_->Branch( "subjetAK8_softdrop_e"      , &subjetAK8_softdrop_e       );
      tree_->Branch( "subjetAK8_softdrop_charge" , &subjetAK8_softdrop_charge  );
      tree_->Branch( "subjetAK8_softdrop_partonFlavour", &subjetAK8_softdrop_partonFlavour );
      tree_->Branch( "subjetAK8_softdrop_hadronFlavour", &subjetAK8_softdrop_hadronFlavour );
    //tree_->Branch( "subjetAK8_softdrop_ssv"    , &subjetAK8_softdrop_ssv     );
      tree_->Branch( "subjetAK8_softdrop_csv"    , &subjetAK8_softdrop_csv     );
    //tree_->Branch( "subjetAK8_softdrop_tchp"   , &subjetAK8_softdrop_tchp    );
    //tree_->Branch( "subjetAK8_softdrop_tche"   , &subjetAK8_softdrop_tche    );
    //tree_->Branch( "subjetAK8_softdrop_jp"     , &subjetAK8_softdrop_jp      );
    //tree_->Branch( "subjetAK8_softdrop_jbp"    , &subjetAK8_softdrop_jbp     );
    
    // /*----------------------AK8 jets pruned-----------------------*/
    // tree_->Branch( "njetsAK8pruned"   , &njetsAK8pruned   );
    // tree_->Branch( "jetAK8pruned_pt"   , &jetAK8pruned_pt   );
    // tree_->Branch( "jetAK8pruned_eta"   , &jetAK8pruned_eta   );
    // tree_->Branch( "jetAK8pruned_mass"   , &jetAK8pruned_mass    );
    // tree_->Branch( "jetAK8pruned_phi"   , &jetAK8pruned_phi   );
    // tree_->Branch( "jetAK8pruned_e"   , &jetAK8pruned_e   );
    // tree_->Branch( "jetAK8pruned_flavour", &jetAK8pruned_flavour );
    // tree_->Branch( "jetAK8pruned_charge" , &jetAK8pruned_charge  );
    // tree_->Branch( "jetAK8pruned_ssv"   , &jetAK8pruned_ssv   );
    // tree_->Branch( "jetAK8pruned_csv"   , &jetAK8pruned_csv   );
    // tree_->Branch( "jetAK8pruned_tchp"   , &jetAK8pruned_tchp    );
    // tree_->Branch( "jetAK8pruned_tche"   , &jetAK8pruned_tche    );
    // tree_->Branch( "jetAK8pruned_jp"   , &jetAK8pruned_jp   );
    // tree_->Branch( "jetAK8pruned_jbp"   , &jetAK8pruned_jbp   );
    // tree_->Branch( "jetAK8pruned_nSVs"   , &jetAK8pruned_nSVs    );
  
    // /*----------------------AK8 jets softdrop-----------------------*/
    //tree_->Branch( "njetsAK8softdrop"	     , &njetsAK8softdrop	);
    //tree_->Branch( "jetAK8softdrop_pt"     , &jetAK8softdrop_pt       );
    //tree_->Branch( "jetAK8softdrop_eta"    , &jetAK8softdrop_eta      );
    //tree_->Branch( "jetAK8softdrop_mass"   , &jetAK8softdrop_mass     );
    //tree_->Branch( "jetAK8softdrop_phi"    , &jetAK8softdrop_phi      );
    //tree_->Branch( "jetAK8softdrop_e"	     , &jetAK8softdrop_e        );
    //tree_->Branch( "jetAK8softdrop_flavour", &jetAK8softdrop_flavour  );
    //tree_->Branch( "jetAK8softdrop_charge" , &jetAK8softdrop_charge   );
    //tree_->Branch( "jetAK8softdrop_ssv"    , &jetAK8softdrop_ssv      );
    //tree_->Branch( "jetAK8softdrop_csv"    , &jetAK8softdrop_csv      );
    //tree_->Branch( "jetAK8softdrop_tchp"   , &jetAK8softdrop_tchp     );
    //tree_->Branch( "jetAK8softdrop_tche"   , &jetAK8softdrop_tche     );
    //tree_->Branch( "jetAK8softdrop_jp"     , &jetAK8softdrop_jp       );
    //tree_->Branch( "jetAK8softdrop_jbp"    , &jetAK8softdrop_jbp      );
    //tree_->Branch( "jetAK8softdrop_nSVs"   , &jetAK8softdrop_nSVs     );
	  
    if (runFlags["doPrunedSubjets"]) {
      /*----------------------Pruned AK8 subjets---------------------------*/   
      tree_->Branch( "subjetAK8_pruned_N"	 , &subjetAK8_pruned_N	     );
      tree_->Branch( "subjetAK8_pruned_pt"       , &subjetAK8_pruned_pt	     );
      tree_->Branch( "subjetAK8_pruned_eta"      , &subjetAK8_pruned_eta     );
      tree_->Branch( "subjetAK8_pruned_mass"     , &subjetAK8_pruned_mass    );
      tree_->Branch( "subjetAK8_pruned_phi"      , &subjetAK8_pruned_phi     );
      tree_->Branch( "subjetAK8_pruned_e"        , &subjetAK8_pruned_e	     );
      tree_->Branch( "subjetAK8_pruned_charge"   , &subjetAK8_pruned_charge  );
      tree_->Branch( "subjetAK8_pruned_partonFlavour"  , &subjetAK8_pruned_partonFlavour );
      tree_->Branch( "subjetAK8_pruned_hadronFlavour"  , &subjetAK8_pruned_hadronFlavour );
    //tree_->Branch( "subjetAK8_pruned_ssv"      , &subjetAK8_pruned_ssv     );
      tree_->Branch( "subjetAK8_pruned_csv"      , &subjetAK8_pruned_csv     );
    //tree_->Branch( "subjetAK8_pruned_tchp"     , &subjetAK8_pruned_tchp    );
    //tree_->Branch( "subjetAK8_pruned_tche"     , &subjetAK8_pruned_tche    );
    //tree_->Branch( "subjetAK8_pruned_jp"       , &subjetAK8_pruned_jp      );
    //tree_->Branch( "subjetAK8_pruned_jbp"      , &subjetAK8_pruned_jbp     );
  
    } //doPruning
  } //doAK8Jets
  
  if (runFlags["runOnMC"]) {
    if (runFlags["doGenJets"]) {
      /*-------------------------AK4 GenJets---------------------------*/	 
      tree_->Branch( "genJetAK4_N"	    , &genJetAK4_N 	  );
      tree_->Branch( "genJetAK4_pt"         , &genJetAK4_pt	  );
      tree_->Branch( "genJetAK4_eta"	    , &genJetAK4_eta	  );
      tree_->Branch( "genJetAK4_mass"	    , &genJetAK4_mass	  );
      tree_->Branch( "genJetAK4_phi"	    , &genJetAK4_phi	  );
      tree_->Branch( "genJetAK4_e"	    , &genJetAK4_e 	  );
      tree_->Branch( "genJetNoNuAK4_pt"	    , &genJetNoNuAK4_pt	  );
      tree_->Branch( "genJetNoNuAK4_mass"   , &genJetNoNuAK4_mass );
      tree_->Branch( "genJetNoNuAK4_e"	    , &genJetNoNuAK4_e    );
      
      tree_->Branch( "genJetAK8_N"	     , &genJetAK8_N 	       );
      tree_->Branch( "genJetAK8_pt"	     , &genJetAK8_pt	       );
      tree_->Branch( "genJetAK8_eta"	     , &genJetAK8_eta	       );
      tree_->Branch( "genJetAK8_mass"	     , &genJetAK8_mass	       );
      tree_->Branch( "genJetAK8_phi"	     , &genJetAK8_phi	       );
      tree_->Branch( "genJetAK8_e"	     , &genJetAK8_e 	       );
      tree_->Branch( "genJetAK8_prunedmass"  , &genJetAK8_prunedmass   );
      tree_->Branch( "genJetAK8_softdropmass", &genJetAK8_softdropmass );
      
    } //runOnMC
  } //doGenJets

  if (runFlags["doTriggerDecisions"]) {
    /** HLT trigger decisions */
    tree_->Branch("HLT_isFired", &HLT_isFired );
  }

  
  if (runFlags["doTriggerObjects"]) {
    /** HLT trigger objects */
    tree_->Branch("triggerObject_pt"		, &triggerObject_pt		);
    tree_->Branch("triggerObject_eta"		, &triggerObject_eta		);
    tree_->Branch("triggerObject_phi"		, &triggerObject_phi	        );
    tree_->Branch("triggerObject_mass"		, &triggerObject_mass		);
    tree_->Branch("triggerObject_filterIDs"	, &triggerObject_filterIDs	);
    tree_->Branch("triggerObject_firedTrigger"	, &triggerObject_firedTrigger	);
  } //doTriggerObjects
  
  if (runFlags["doHltFilters"]) {
    /** HLT filter decisions */
    tree_->Branch("passFilter_HBHE"                 ,&passFilter_HBHE_                ,"passFilter_HBHE_/O");
    tree_->Branch("passFilter_HBHELoose"            ,&passFilter_HBHELoose_	      ,"passFilter_HBHELoose_/O");
    tree_->Branch("passFilter_HBHETight"            ,&passFilter_HBHETight_	      ,"passFilter_HBHETight_/O");
    tree_->Branch("passFilter_CSCHalo"              ,&passFilter_CSCHalo_             ,"passFilter_CSCHalo_/O");
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
  } //do HltFilters

  if (runFlags["doMissingEt"]) {
    /** MET */
    tree_->Branch("METraw_et"		        , &METraw_et	     );
    tree_->Branch("METraw_phi"		        , &METraw_phi	     ); 
    tree_->Branch("METraw_sumEt"		, &METraw_sumEt	     );   
    tree_->Branch("MET_corrPx"		        , &MET_corrPx	     ); 
    tree_->Branch("MET_corrPy"		        , &MET_corrPy	     );   
    tree_->Branch("MET_et"	                , &MET_et  	     ); 
    tree_->Branch("MET_phi"	                , &MET_phi           );
    tree_->Branch("MET_sumEt"	                , &MET_sumEt 	     ); 
    //tree_->Branch("METdefault_et"	        , &METdefault_et     ); 
    //tree_->Branch("METdefault_sumEt"	        , &METdefault_sumEt  ); 
    //tree_->Branch("METdefault_phi"	        , &METdefault_phi    );
    //tree_->Branch("MET_T1Uncertainty"	        , &MET_T1Uncertainty );
  } //doMissingEt


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
  }
  
}

//=================================================================================================================== 
void NtupleBranches::reset( void ){

  /** genParticle */
  genParticle_N = 0;
  genParticle_pt.clear();
  genParticle_px.clear();
  genParticle_py.clear();
  genParticle_pz.clear();
  genParticle_e.clear();
  genParticle_eta.clear();
  genParticle_phi.clear();
  genParticle_mass.clear();
  genParticle_pdgId.clear();
  genParticle_status.clear();
  genParticle_mother.clear();
  genParticle_nMoth.clear();
  genParticle_nDau.clear();
  genParticle_dau.clear();
  
  /** generator info */
  lheV_pt     = 0;
  lheHT       = 0;
  lheNj       = 0;
  genWeight   = 0;
  qScale      = 0;
  PDF_id.clear();  
  PDF_x.clear();	
  PDF_xPDF.clear();
    
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
  el_pfRelIso.clear();
  el_photonIso.clear();
  el_neutralHadIso.clear();
  el_chargedHadIso.clear();
  el_trackIso.clear();


  el_pfRhoCorrRelIso03Boost.clear();
  el_pfRhoCorrRelIso04Boost.clear();
  el_pfDeltaCorrRelIsoBoost.clear();
  el_pfRelIsoBoost.clear();
  el_photonIsoBoost.clear();
  el_neutralHadIsoBoost.clear();
  el_chargedHadIsoBoost.clear();

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

  el_dr03EcalRecHitSumEt.clear();
  el_dr03HcalDepth1TowerSumEt.clear();
  el_rho.clear();
  el_ecalDriven.clear();
  el_dEtaInSeed.clear();
  el_full5x5_e2x5Max.clear();
  el_full5x5_e5x5.clear();
  el_full5x5_e1x5.clear();
  el_dr03TkSumPt.clear();
  el_superCluster_e.clear();
  el_hadronicOverEm.clear();

  el_isVetoElectron.clear();
  el_isMediumElectron.clear();
  el_isTightElectron.clear();
  el_isHeepElectron.clear();
  el_isHeep51Electron.clear();
  el_isLooseElectron.clear();
  el_isVetoElectronBoosted.clear();
  el_isMediumElectronBoosted.clear();
  el_isTightElectronBoosted.clear();
  el_isHeepElectronBoosted.clear();
  el_isHeep51ElectronBoosted.clear();
  el_isLooseElectronBoosted.clear();

  el_pfRhoCorrRelIso03Boost.clear();
  el_pfRhoCorrRelIso04Boost.clear();
  el_pfDeltaCorrRelIsoBoost.clear();
  el_pfRelIsoBoost.clear();
  el_photonIsoBoost.clear();
  el_neutralHadIsoBoost.clear();
  el_chargedHadIsoBoost.clear();

  el_SemileptonicPFIso.clear();
  el_SemileptonicCorrPFIso.clear();
  
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
  mu_isLooseMuon.clear();
  mu_isPFMuon.clear();
  mu_isSoftMuon.clear();

  mu_pfRhoCorrRelIso03.clear();
  mu_pfRhoCorrRelIso04.clear();
  mu_pfDeltaCorrRelIso.clear();
  mu_pfRelIso.clear();
  mu_photonIso.clear();
  mu_neutralHadIso.clear();
  mu_chargedHadIso.clear();
  mu_trackIso.clear();
  mu_d0.clear();
  mu_dz.clear();
  mu_bestTrack_pt.clear();
  mu_bestTrack_ptErr.clear();

  mu_pfRhoCorrRelIso03Boost.clear();
  mu_pfRhoCorrRelIso04Boost.clear();
  mu_pfDeltaCorrRelIsoBoost.clear();
  mu_pfRelIsoBoost.clear();
  mu_photonIsoBoost.clear();
  mu_neutralHadIsoBoost.clear();
  mu_chargedHadIsoBoost.clear();

  mu_normChi2.clear();
  mu_isGlobalMuon.clear();
  mu_trackerHits.clear();
  mu_matchedStations.clear();
  mu_pixelHits.clear();
  mu_globalHits.clear();

  mu_SemileptonicPFIso.clear();
  mu_SemileptonicCorrPFIso.clear();


  /** taus */
  tau_N       = 0;
  tau_pdgId.clear();
  tau_charge.clear();
  tau_e.clear();
  tau_eta.clear();
  tau_phi.clear();
  tau_mass.clear();
  tau_pt.clear();

  tau_pfRhoCorrRelIso03.clear();
  tau_pfRhoCorrRelIso04.clear();
  tau_pfDeltaCorrRelIso.clear();
  tau_pfRelIso.clear();
  tau_photonIso.clear();
  tau_neutralHadIso.clear();
  tau_chargedHadIso.clear();
  tau_trackIso.clear();
  tau_d0.clear();

  tau_pfRhoCorrRelIso03Boost.clear();
  tau_pfRhoCorrRelIso04Boost.clear();
  tau_pfDeltaCorrRelIsoBoost.clear();
  tau_pfRelIsoBoost.clear();
  tau_photonIsoBoost.clear();
  tau_neutralHadIsoBoost.clear();
  tau_chargedHadIsoBoost.clear();

  tau_TauType.clear();

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
  tau_byIsolationMVA3oldDMwoLTraw.clear();
  tau_byVLooseIsolationMVA3oldDMwoLT.clear();
  tau_byLooseIsolationMVA3oldDMwoLT.clear();
  tau_byMediumIsolationMVA3oldDMwoLT.clear();
  tau_byTightIsolationMVA3oldDMwoLT.clear();
  tau_byVTightIsolationMVA3oldDMwoLT.clear();
  tau_byVVTightIsolationMVA3oldDMwoLT.clear();
  tau_byIsolationMVA3oldDMwLTraw.clear();
  tau_byVLooseIsolationMVA3oldDMwLT.clear();
  tau_byLooseIsolationMVA3oldDMwLT.clear();
  tau_byMediumIsolationMVA3oldDMwLT.clear();
  tau_byTightIsolationMVA3oldDMwLT.clear();
  tau_byVTightIsolationMVA3oldDMwLT.clear();
  tau_byVVTightIsolationMVA3oldDMwLT.clear();
  tau_byIsolationMVA3newDMwoLTraw.clear();
  tau_byVLooseIsolationMVA3newDMwoLT.clear();
  tau_byLooseIsolationMVA3newDMwoLT.clear();
  tau_byMediumIsolationMVA3newDMwoLT.clear();
  tau_byTightIsolationMVA3newDMwoLT.clear();
  tau_byVTightIsolationMVA3newDMwoLT.clear();
  tau_byVVTightIsolationMVA3newDMwoLT.clear();
  tau_byIsolationMVA3newDMwLTraw.clear();
  tau_byVLooseIsolationMVA3newDMwLT.clear();
  tau_byLooseIsolationMVA3newDMwLT.clear();
  tau_byMediumIsolationMVA3newDMwLT.clear();
  tau_byTightIsolationMVA3newDMwLT.clear();
  tau_byVTightIsolationMVA3newDMwLT.clear();
  tau_byVVTightIsolationMVA3newDMwLT.clear();
  tau_againstElectronLoose.clear();
  tau_againstElectronMedium.clear();
  tau_againstElectronTight.clear();
  tau_againstElectronMVA5raw.clear();
  tau_againstElectronMVA5category.clear();
  tau_againstElectronVLooseMVA5.clear();
  tau_againstElectronLooseMVA5.clear();
  tau_againstElectronMediumMVA5.clear();
  tau_againstElectronTightMVA5.clear();
  tau_againstElectronVTightMVA5.clear();
  tau_againstMuonLoose.clear();
  tau_againstMuonMedium.clear();
  tau_againstMuonTight.clear();
  tau_againstMuonLoose2.clear();
  tau_againstMuonMedium2.clear();
  tau_againstMuonTight2.clear();
  tau_againstMuonLoose3.clear();
  tau_againstMuonTight3.clear();
  tau_againstMuonMVAraw.clear();
  tau_againstMuonLooseMVA.clear();
  tau_againstMuonMediumMVA.clear();
  tau_againstMuonTightMVA.clear();
  
  /** energy density */
  rho = 0;
  
  /** AK4 jets */
  jetAK4_N    = 0;
  jetAK4_pt.clear();
  jetAK4_eta.clear();
  jetAK4_mass.clear();
  jetAK4_phi.clear();
  jetAK4_e.clear();
  jetAK4_jec.clear();
  //jetAK4_jecUp.clear();
  //jetAK4_jecDown.clear(); 
  jetAK4_IDTight.clear();
  jetAK4_IDLoose.clear();
  jetAK4_muf.clear();
  jetAK4_phf.clear();
  jetAK4_emf.clear();
  jetAK4_nhf.clear();
  jetAK4_chf.clear();
  jetAK4_area.clear();
  jetAK4_cm.clear();
  jetAK4_nm.clear();
  jetAK4_che.clear();
  jetAK4_ne.clear();
  jetAK4_hf_hf.clear();
  jetAK4_hf_emf.clear();
  jetAK4_hof.clear();
  jetAK4_chm.clear();
  jetAK4_neHadMult.clear();
  jetAK4_phoMult.clear();
  jetAK4_nemf.clear();
  jetAK4_cemf.clear();
  jetAK4_charge.clear();
  jetAK4_partonFlavour.clear();
  jetAK4_hadronFlavour.clear();
  jetAK4_genParton_pdgID.clear();
  jetAK4_nbHadrons.clear();
  jetAK4_ncHadrons.clear();
  //jetAK4_ssv.clear();
  jetAK4_cisv.clear();
  //jetAK4_tchp.clear();
  //jetAK4_tche.clear();
  //jetAK4_jp.clear();
  //jetAK4_jbp.clear();
  jetAK4_vtxMass.clear();
  jetAK4_vtxNtracks.clear();
  jetAK4_vtx3DVal.clear();
  jetAK4_vtx3DSig.clear();
  //jetAK4_nSVs.clear();
  
  jetAK8_N = 0;
  jetAK8_pt.clear();
  jetAK8_eta.clear();
  jetAK8_mass.clear();
  jetAK8_phi.clear();
  jetAK8_e.clear();
  jetAK8_jec.clear();
  //jetAK8_jecUp.clear();
  //jetAK8_jecDown.clear();
  jetAK8_IDLoose.clear();
  jetAK8_IDTight.clear();
  jetAK8_muf.clear();
  jetAK8_phf.clear();
  jetAK8_emf.clear();
  jetAK8_nhf.clear();
  jetAK8_chf.clear();    
  jetAK8_area.clear();
  jetAK8_cm.clear();
  jetAK8_nm.clear();
  jetAK8_che.clear();
  jetAK8_ne.clear();
  jetAK8_hf_hf.clear();
  jetAK8_hf_emf.clear();
  jetAK8_hof.clear();
  jetAK8_chm.clear();
  jetAK8_neHadMult.clear();
  jetAK8_phoMult.clear();
  jetAK8_nemf.clear();
  jetAK8_cemf.clear();
  jetAK8_charge.clear();
  jetAK8_partonFlavour.clear();
  jetAK8_hadronFlavour.clear();
  jetAK8_genParton_pdgID.clear();
  jetAK8_nbHadrons.clear();
  jetAK8_ncHadrons.clear();
  jetAK8_Hbbtag.clear();
  //jetAK8_ssv.clear();
  jetAK8_csv.clear();    
  //jetAK8_tchp.clear();
  //jetAK8_tche.clear();
  //jetAK8_jp.clear();
  //jetAK8_jbp.clear();
  jetAK8_tau1.clear();
  jetAK8_tau2.clear();
  jetAK8_tau3.clear();    
  jetAK8_pruned_mass.clear();
  jetAK8_softdrop_mass.clear();
  jetAK10_trimmed_mass.clear();
  jetAK8_pruned_massCorr.clear();
  jetAK8_softdrop_massCorr.clear();
  jetAK10_trimmed_massCorr.clear();
  jetAK8_pruned_jec.clear();
  jetAK8_softdrop_jec.clear();
  jetAK10_trimmed_jec.clear();
  jetAK10_ecf1.clear();
  jetAK10_ecf2.clear();
  jetAK10_ecf3.clear();    
  jetAK8_puppi_tau1.clear();
  jetAK8_puppi_tau2.clear();
  jetAK8_puppi_tau3.clear();    
  jetAK8_puppi_pruned_mass.clear();
  jetAK8_puppi_softdrop_mass.clear();
  jetAK8_puppi_pruned_massCorr.clear();
  jetAK8_puppi_softdrop_massCorr.clear();
  jetAK8_puppi_pruned_jec.clear();
  jetAK8_puppi_softdrop_jec.clear();
  //jetAK8_trimmed_mass.clear();
  //jetAK8_filtered_mass.clear();
  //jetAK8_nSubJets.clear();

  /** AK8 jets pruned */
  //njetsAK8_pruned = 0;
  //jetAK8_pruned_pt.clear();
  //jetAK8_pruned_eta.clear();
  //jetAK8_pruned_mass.clear();
  //jetAK8_pruned_phi.clear();
  //jetAK8_pruned_e.clear();
  //jetAK8_pruned_charge.clear();
  //jetAK8_pruned_flavour.clear();
  //jetAK8_pruned_ssv.clear();
  //jetAK8_pruned_csv.clear();
  //jetAK8_pruned_tchp.clear();
  //jetAK8_pruned_tche.clear();
  //jetAK8_pruned_jp.clear();
  //jetAK8_pruned_jbp.clear();
  //jetAK8_pruned_nSVs.clear();

  /** AK8 jets softdrop */
  //njetsAK8_softdrop = 0;
  //jetAK8_softdrop_pt.clear();
  //jetAK8_softdrop_eta.clear();
  //jetAK8_softdrop_mass.clear();
  //jetAK8_softdrop_phi.clear();
  //jetAK8_softdrop_e.clear();
  //jetAK8_softdrop_charge.clear();
  //jetAK8_softdrop_flavour.clear();
  //jetAK8_softdrop_ssv.clear();
  //jetAK8_softdrop_csv.clear();
  //jetAK8_softdrop_tchp.clear();
  //jetAK8_softdrop_tche.clear();
  //jetAK8_softdrop_jp.clear();
  //jetAK8_softdrop_jbp.clear();
  //jetAK8_softdrop_nSVs.clear();

  /** pruned AK8 subjets  */
  subjetAK8_pruned_N.clear();
  subjetAK8_pruned_pt.clear();
  subjetAK8_pruned_eta.clear();
  subjetAK8_pruned_mass.clear();
  subjetAK8_pruned_phi.clear();
  subjetAK8_pruned_e.clear();
  subjetAK8_pruned_charge.clear();
  subjetAK8_pruned_partonFlavour.clear();
  subjetAK8_pruned_hadronFlavour.clear();
  //subjetAK8_pruned_ssv.clear();
  subjetAK8_pruned_csv.clear();    
  //subjetAK8_pruned_tchp.clear();
  //subjetAK8_pruned_tche.clear();
  //subjetAK8_pruned_jp.clear();
  //subjetAK8_pruned_jbp.clear();

  /** softdrop AK8 subjets */
  subjetAK8_softdrop_N.clear()         ;
  subjetAK8_softdrop_pt.clear();
  subjetAK8_softdrop_eta.clear();
  subjetAK8_softdrop_mass.clear();
  subjetAK8_softdrop_phi.clear();
  subjetAK8_softdrop_e.clear();
  subjetAK8_softdrop_charge.clear();
  subjetAK8_softdrop_partonFlavour.clear();
  subjetAK8_softdrop_hadronFlavour.clear();
  //subjetAK8_softdrop_ssv.clear();
  subjetAK8_softdrop_csv.clear();
  //subjetAK8_softdrop_tchp.clear();
  //subjetAK8_softdrop_tche.clear();
  //subjetAK8_softdrop_jp.clear();
  //subjetAK8_softdrop_jbp.clear();

  /** AK4 genJets*/
  genJetAK4_N = 0;
  genJetAK4_pt.clear();
  genJetAK4_eta.clear();
  genJetAK4_mass.clear();
  genJetAK4_phi.clear();
  genJetAK4_e.clear();
  genJetNoNuAK4_pt.clear();
  genJetNoNuAK4_mass.clear();
  genJetNoNuAK4_e.clear();
  
  genJetAK8_N = 0;
  genJetAK8_pt              .clear();
  genJetAK8_eta             .clear();
  genJetAK8_mass            .clear();
  genJetAK8_phi             .clear();
  genJetAK8_e	              .clear();
  genJetAK8_prunedmass      .clear();
  genJetAK8_softdropmass    .clear();

  /** HLT trigger decisions */
  HLT_isFired.clear();
  
  /** HLT trigger objects */
  triggerObject_pt.clear();
  triggerObject_eta.clear();
  triggerObject_phi.clear();
  triggerObject_mass.clear();
  triggerObject_filterIDs.clear();
  triggerObject_firedTrigger.clear();

  /** HLT filter decisions */
  passFilter_HBHE_                  = false;
  passFilter_HBHELoose_             = false;
  passFilter_HBHETight_             = false;
  passFilter_CSCHalo_               = false;
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

  /** MET */
  METraw_et.clear();	 
  METraw_phi.clear();
  METraw_sumEt.clear();
  MET_corrPx.clear();
  MET_corrPy.clear();
  MET_et.clear();
  MET_phi.clear();
  MET_sumEt.clear();
  MET_T1Uncertainty.clear();

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
  
  /*--------------------------PU infos--------------------------*/  			       
  nPuVtxTrue.clear();
  nPuVtx.clear();
  bX.clear();



} 
