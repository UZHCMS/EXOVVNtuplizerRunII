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

  if (runFlags["doJpsiEle"]){
    tree_->Branch("IsJpsiEle", &IsJpsiEle );
  }

  if (runFlags["doJpsiMu"]){
    tree_->Branch("IsJpsiMu" , &IsJpsiMu  );

    tree_->Branch("JpsiMu_nCandidates", &JpsiMu_nCandidates );

    tree_->Branch("JpsiMu_mu1_pt", &JpsiMu_mu1_pt );
    tree_->Branch("JpsiMu_mu1_eta", &JpsiMu_mu1_eta );
    tree_->Branch("JpsiMu_mu1_phi", &JpsiMu_mu1_phi );
    tree_->Branch("JpsiMu_mu1_mass", &JpsiMu_mu1_mass );
    tree_->Branch("JpsiMu_mu1_unfit_pt", &JpsiMu_mu1_unfit_pt );
    tree_->Branch("JpsiMu_mu1_unfit_eta", &JpsiMu_mu1_unfit_eta );
    tree_->Branch("JpsiMu_mu1_unfit_phi", &JpsiMu_mu1_unfit_phi );
    tree_->Branch("JpsiMu_mu1_unfit_mass", &JpsiMu_mu1_unfit_mass );
    tree_->Branch("JpsiMu_mu1_q", &JpsiMu_mu1_q );
    tree_->Branch("JpsiMu_mu1_isLoose"  , &JpsiMu_mu1_isLoose   );
    tree_->Branch("JpsiMu_mu1_isTight"  , &JpsiMu_mu1_isTight   );
    tree_->Branch("JpsiMu_mu1_isPF"     , &JpsiMu_mu1_isPF      );
    tree_->Branch("JpsiMu_mu1_isGlobal" , &JpsiMu_mu1_isGlobal  );
    tree_->Branch("JpsiMu_mu1_isTracker", &JpsiMu_mu1_isTracker );
    tree_->Branch("JpsiMu_mu1_isSoft"   , &JpsiMu_mu1_isSoft    );
    tree_->Branch("JpsiMu_mu1_vx"   , &JpsiMu_mu1_vx    );
    tree_->Branch("JpsiMu_mu1_vy"   , &JpsiMu_mu1_vy    );
    tree_->Branch("JpsiMu_mu1_vz"   , &JpsiMu_mu1_vz    );
    tree_->Branch("JpsiMu_mu1_iso"   , &JpsiMu_mu1_iso    );
    tree_->Branch("JpsiMu_mu1_dbiso"   , &JpsiMu_mu1_dbiso    );

    tree_->Branch("JpsiMu_mu2_pt", &JpsiMu_mu2_pt );
    tree_->Branch("JpsiMu_mu2_eta", &JpsiMu_mu2_eta );
    tree_->Branch("JpsiMu_mu2_phi", &JpsiMu_mu2_phi );
    tree_->Branch("JpsiMu_mu2_mass", &JpsiMu_mu2_mass );
    tree_->Branch("JpsiMu_mu2_unfit_pt", &JpsiMu_mu2_unfit_pt );
    tree_->Branch("JpsiMu_mu2_unfit_eta", &JpsiMu_mu2_unfit_eta );
    tree_->Branch("JpsiMu_mu2_unfit_phi", &JpsiMu_mu2_unfit_phi );
    tree_->Branch("JpsiMu_mu2_unfit_mass", &JpsiMu_mu2_unfit_mass );
    tree_->Branch("JpsiMu_mu2_q", &JpsiMu_mu2_q );
    tree_->Branch("JpsiMu_mu2_isLoose"  , &JpsiMu_mu2_isLoose   );
    tree_->Branch("JpsiMu_mu2_isTight"  , &JpsiMu_mu2_isTight   );
    tree_->Branch("JpsiMu_mu2_isPF"     , &JpsiMu_mu2_isPF      );
    tree_->Branch("JpsiMu_mu2_isGlobal" , &JpsiMu_mu2_isGlobal  );
    tree_->Branch("JpsiMu_mu2_isTracker", &JpsiMu_mu2_isTracker );
    tree_->Branch("JpsiMu_mu2_isSoft"   , &JpsiMu_mu2_isSoft    );
    tree_->Branch("JpsiMu_mu2_vx"   , &JpsiMu_mu2_vx    );
    tree_->Branch("JpsiMu_mu2_vy"   , &JpsiMu_mu2_vy    );
    tree_->Branch("JpsiMu_mu2_vz"   , &JpsiMu_mu2_vz    );
    tree_->Branch("JpsiMu_mu2_iso"   , &JpsiMu_mu2_iso    );
    tree_->Branch("JpsiMu_mu2_dbiso"   , &JpsiMu_mu2_dbiso    );

    tree_->Branch("JpsiMu_mu3_pt", &JpsiMu_mu3_pt );
    tree_->Branch("JpsiMu_mu3_eta", &JpsiMu_mu3_eta );
    tree_->Branch("JpsiMu_mu3_phi", &JpsiMu_mu3_phi );
    tree_->Branch("JpsiMu_mu3_mass", &JpsiMu_mu3_mass );
    tree_->Branch("JpsiMu_mu3_unfit_pt", &JpsiMu_mu3_unfit_pt );
    tree_->Branch("JpsiMu_mu3_unfit_eta", &JpsiMu_mu3_unfit_eta );
    tree_->Branch("JpsiMu_mu3_unfit_phi", &JpsiMu_mu3_unfit_phi );
    tree_->Branch("JpsiMu_mu3_unfit_mass", &JpsiMu_mu3_unfit_mass );
    tree_->Branch("JpsiMu_mu3_doca2mu1", &JpsiMu_mu3_doca2mu1 );
    tree_->Branch("JpsiMu_mu3_doca2mu2", &JpsiMu_mu3_doca2mu2 );
    tree_->Branch("JpsiMu_mu3_q", &JpsiMu_mu3_q );
    tree_->Branch("JpsiMu_mu3_isLoose"  , &JpsiMu_mu3_isLoose   );
    tree_->Branch("JpsiMu_mu3_isTight"  , &JpsiMu_mu3_isTight   );
    tree_->Branch("JpsiMu_mu3_isPF"     , &JpsiMu_mu3_isPF      );
    tree_->Branch("JpsiMu_mu3_isGlobal" , &JpsiMu_mu3_isGlobal  );
    tree_->Branch("JpsiMu_mu3_isTracker", &JpsiMu_mu3_isTracker );
    tree_->Branch("JpsiMu_mu3_isSoft"   , &JpsiMu_mu3_isSoft    );
    tree_->Branch("JpsiMu_mu3_vx"   , &JpsiMu_mu3_vx    );
    tree_->Branch("JpsiMu_mu3_vy"   , &JpsiMu_mu3_vy    );
    tree_->Branch("JpsiMu_mu3_vz"   , &JpsiMu_mu3_vz    );
    tree_->Branch("JpsiMu_mu3_iso"   , &JpsiMu_mu3_iso    );
    tree_->Branch("JpsiMu_mu3_dbiso"   , &JpsiMu_mu3_dbiso    );

    tree_->Branch("JpsiMu_PV_vx", &JpsiMu_PV_vx );
    tree_->Branch("JpsiMu_PV_vy", &JpsiMu_PV_vy );
    tree_->Branch("JpsiMu_PV_vz", &JpsiMu_PV_vz );

    tree_->Branch("JpsiMu_bbPV_vx", &JpsiMu_bbPV_vx );
    tree_->Branch("JpsiMu_bbPV_vy", &JpsiMu_bbPV_vy );
    tree_->Branch("JpsiMu_bbPV_vz", &JpsiMu_bbPV_vz );

    tree_->Branch("JpsiMu_bbPV_refit_vx", &JpsiMu_bbPV_vx );
    tree_->Branch("JpsiMu_bbPV_refit_vy", &JpsiMu_bbPV_vy );
    tree_->Branch("JpsiMu_bbPV_refit_vz", &JpsiMu_bbPV_vz );

    tree_->Branch("JpsiMu_genPV_vx", &JpsiMu_genPV_vx );
    tree_->Branch("JpsiMu_genPV_vy", &JpsiMu_genPV_vy );
    tree_->Branch("JpsiMu_genPV_vz", &JpsiMu_genPV_vz );

    tree_->Branch("JpsiMu_Jpsi_pt", &JpsiMu_Jpsi_pt );
    tree_->Branch("JpsiMu_Jpsi_eta", &JpsiMu_Jpsi_eta );
    tree_->Branch("JpsiMu_Jpsi_phi", &JpsiMu_Jpsi_phi );
    tree_->Branch("JpsiMu_Jpsi_mass", &JpsiMu_Jpsi_mass );
    tree_->Branch("JpsiMu_Jpsi_vprob", &JpsiMu_Jpsi_vprob );
    tree_->Branch("JpsiMu_Jpsi_lip", &JpsiMu_Jpsi_lip);
    tree_->Branch("JpsiMu_Jpsi_lips", &JpsiMu_Jpsi_lips);
    tree_->Branch("JpsiMu_Jpsi_pvip", &JpsiMu_Jpsi_pvip);
    tree_->Branch("JpsiMu_Jpsi_pvips", &JpsiMu_Jpsi_pvips);
    tree_->Branch("JpsiMu_Jpsi_fl3d", &JpsiMu_Jpsi_fl3d);
    tree_->Branch("JpsiMu_Jpsi_fls3d", &JpsiMu_Jpsi_fls3d);
    tree_->Branch("JpsiMu_Jpsi_alpha", &JpsiMu_Jpsi_alpha);
    tree_->Branch("JpsiMu_Jpsi_maxdoca", &JpsiMu_Jpsi_maxdoca);
    tree_->Branch("JpsiMu_Jpsi_mindoca", &JpsiMu_Jpsi_mindoca);
    tree_->Branch("JpsiMu_Jpsi_vx", &JpsiMu_Jpsi_vx );
    tree_->Branch("JpsiMu_Jpsi_vy", &JpsiMu_Jpsi_vy );
    tree_->Branch("JpsiMu_Jpsi_vz", &JpsiMu_Jpsi_vz );
    tree_->Branch("JpsiMu_Jpsi_unfit_pt", &JpsiMu_Jpsi_unfit_pt );
    tree_->Branch("JpsiMu_Jpsi_unfit_mass", &JpsiMu_Jpsi_unfit_mass );
    tree_->Branch("JpsiMu_Jpsi_unfit_vprob", &JpsiMu_Jpsi_unfit_vprob );
    tree_->Branch("JpsiMu_Jpsi_unfit_vx", &JpsiMu_Jpsi_unfit_vx );
    tree_->Branch("JpsiMu_Jpsi_unfit_vy", &JpsiMu_Jpsi_unfit_vy );
    tree_->Branch("JpsiMu_Jpsi_unfit_vz", &JpsiMu_Jpsi_unfit_vz );


    tree_->Branch("JpsiMu_B_pt", &JpsiMu_B_pt );
    tree_->Branch("JpsiMu_B_eta", &JpsiMu_B_eta );
    tree_->Branch("JpsiMu_B_phi", &JpsiMu_B_phi );
    tree_->Branch("JpsiMu_B_mass", &JpsiMu_B_mass );
    tree_->Branch("JpsiMu_B_vprob", &JpsiMu_B_vprob );
    tree_->Branch("JpsiMu_B_lip", &JpsiMu_B_lip);
    tree_->Branch("JpsiMu_B_lips", &JpsiMu_B_lips);
    tree_->Branch("JpsiMu_B_pvip", &JpsiMu_B_pvip);
    tree_->Branch("JpsiMu_B_pvips", &JpsiMu_B_pvips);
    tree_->Branch("JpsiMu_B_fl3d", &JpsiMu_B_fl3d);
    tree_->Branch("JpsiMu_B_fls3d", &JpsiMu_B_fls3d);
    tree_->Branch("JpsiMu_B_alpha", &JpsiMu_B_alpha);
    tree_->Branch("JpsiMu_B_maxdoca", &JpsiMu_B_maxdoca);
    tree_->Branch("JpsiMu_B_mindoca", &JpsiMu_B_mindoca);
    tree_->Branch("JpsiMu_B_vx", &JpsiMu_B_vx );
    tree_->Branch("JpsiMu_B_vy", &JpsiMu_B_vy );
    tree_->Branch("JpsiMu_B_vz", &JpsiMu_B_vz );
    tree_->Branch("JpsiMu_B_iso", &JpsiMu_B_iso);
    tree_->Branch("JpsiMu_B_iso_ntracks", &JpsiMu_B_iso_ntracks );
    tree_->Branch("JpsiMu_B_iso_mindoca", &JpsiMu_B_iso_mindoca );
    tree_->Branch("JpsiMu_B_unfit_pt", &JpsiMu_B_unfit_pt );
    tree_->Branch("JpsiMu_B_unfit_mass", &JpsiMu_B_unfit_mass );
    tree_->Branch("JpsiMu_B_unfit_vprob", &JpsiMu_B_unfit_vprob );
    tree_->Branch("JpsiMu_B_unfit_vx", &JpsiMu_B_unfit_vx );
    tree_->Branch("JpsiMu_B_unfit_vy", &JpsiMu_B_unfit_vy );
    tree_->Branch("JpsiMu_B_unfit_vz", &JpsiMu_B_unfit_vz );

    tree_->Branch("JpsiMu_ngenmuons", &JpsiMu_ngenmuons);
    tree_->Branch("JpsiMu_isgenmatched", &JpsiMu_isgenmatched);
    tree_->Branch("JpsiMu_mu3_isgenmatched", &JpsiMu_mu3_isgenmatched);

  }



  if (runFlags["doJpsiTau"]){
    tree_->Branch("IsJpsiTau", &IsJpsiTau );

    tree_->Branch("JpsiTau_nCandidates", &JpsiTau_nCandidates );

    tree_->Branch("JpsiTau_mu1_pt", &JpsiTau_mu1_pt );
    tree_->Branch("JpsiTau_mu1_eta", &JpsiTau_mu1_eta );
    tree_->Branch("JpsiTau_mu1_phi", &JpsiTau_mu1_phi );
    tree_->Branch("JpsiTau_mu1_mass", &JpsiTau_mu1_mass );
    tree_->Branch("JpsiTau_mu1_unfit_pt", &JpsiTau_mu1_unfit_pt );
    tree_->Branch("JpsiTau_mu1_unfit_eta", &JpsiTau_mu1_unfit_eta );
    tree_->Branch("JpsiTau_mu1_unfit_phi", &JpsiTau_mu1_unfit_phi );
    tree_->Branch("JpsiTau_mu1_unfit_mass", &JpsiTau_mu1_unfit_mass );
    tree_->Branch("JpsiTau_mu1_q", &JpsiTau_mu1_q );
    tree_->Branch("JpsiTau_mu1_isLoose"  , &JpsiTau_mu1_isLoose   );
    tree_->Branch("JpsiTau_mu1_isTight"  , &JpsiTau_mu1_isTight   );
    tree_->Branch("JpsiTau_mu1_isPF"     , &JpsiTau_mu1_isPF      );
    tree_->Branch("JpsiTau_mu1_isGlobal" , &JpsiTau_mu1_isGlobal  );
    tree_->Branch("JpsiTau_mu1_isTracker", &JpsiTau_mu1_isTracker );
    tree_->Branch("JpsiTau_mu1_isSoft"   , &JpsiTau_mu1_isSoft    );
    tree_->Branch("JpsiTau_mu1_vx"   , &JpsiTau_mu1_vx    );
    tree_->Branch("JpsiTau_mu1_vy"   , &JpsiTau_mu1_vy    );
    tree_->Branch("JpsiTau_mu1_vz"   , &JpsiTau_mu1_vz    );
    tree_->Branch("JpsiTau_mu1_iso"   , &JpsiTau_mu1_iso    );
    tree_->Branch("JpsiTau_mu1_dbiso"   , &JpsiTau_mu1_dbiso    );

    tree_->Branch("JpsiTau_mu2_pt", &JpsiTau_mu2_pt );
    tree_->Branch("JpsiTau_mu2_eta", &JpsiTau_mu2_eta );
    tree_->Branch("JpsiTau_mu2_phi", &JpsiTau_mu2_phi );
    tree_->Branch("JpsiTau_mu2_mass", &JpsiTau_mu2_mass );
    tree_->Branch("JpsiTau_mu2_unfit_pt", &JpsiTau_mu2_unfit_pt );
    tree_->Branch("JpsiTau_mu2_unfit_eta", &JpsiTau_mu2_unfit_eta );
    tree_->Branch("JpsiTau_mu2_unfit_phi", &JpsiTau_mu2_unfit_phi );
    tree_->Branch("JpsiTau_mu2_unfit_mass", &JpsiTau_mu2_unfit_mass );
    tree_->Branch("JpsiTau_mu2_q", &JpsiTau_mu2_q );
    tree_->Branch("JpsiTau_mu2_isLoose"  , &JpsiTau_mu2_isLoose   );
    tree_->Branch("JpsiTau_mu2_isTight"  , &JpsiTau_mu2_isTight   );
    tree_->Branch("JpsiTau_mu2_isPF"     , &JpsiTau_mu2_isPF      );
    tree_->Branch("JpsiTau_mu2_isGlobal" , &JpsiTau_mu2_isGlobal  );
    tree_->Branch("JpsiTau_mu2_isTracker", &JpsiTau_mu2_isTracker );
    tree_->Branch("JpsiTau_mu2_isSoft"   , &JpsiTau_mu2_isSoft    );
    tree_->Branch("JpsiTau_mu2_vx"   , &JpsiTau_mu2_vx    );
    tree_->Branch("JpsiTau_mu2_vy"   , &JpsiTau_mu2_vy    );
    tree_->Branch("JpsiTau_mu2_vz"   , &JpsiTau_mu2_vz    );
    tree_->Branch("JpsiTau_mu2_iso"   , &JpsiTau_mu2_iso    );
    tree_->Branch("JpsiTau_mu2_dbiso"   , &JpsiTau_mu2_dbiso    );

    tree_->Branch("JpsiTau_tau_pt", &JpsiTau_tau_pt );
    tree_->Branch("JpsiTau_tau_eta", &JpsiTau_tau_eta );
    tree_->Branch("JpsiTau_tau_phi", &JpsiTau_tau_phi );
    tree_->Branch("JpsiTau_tau_mass", &JpsiTau_tau_mass );
    tree_->Branch("JpsiTau_tau_q", &JpsiTau_tau_q );
    tree_->Branch("JpsiTau_tau_vx"   , &JpsiTau_tau_vx    );
    tree_->Branch("JpsiTau_tau_vy"   , &JpsiTau_tau_vy    );
    tree_->Branch("JpsiTau_tau_vz"   , &JpsiTau_tau_vz    );
    tree_->Branch("JpsiTau_tau_iso"   , &JpsiTau_tau_iso    );

    tree_->Branch("JpsiTau_PV_vx", &JpsiTau_PV_vx );
    tree_->Branch("JpsiTau_PV_vy", &JpsiTau_PV_vy );
    tree_->Branch("JpsiTau_PV_vz", &JpsiTau_PV_vz );

    tree_->Branch("JpsiTau_bbPV_vx", &JpsiTau_bbPV_vx );
    tree_->Branch("JpsiTau_bbPV_vy", &JpsiTau_bbPV_vy );
    tree_->Branch("JpsiTau_bbPV_vz", &JpsiTau_bbPV_vz );

    tree_->Branch("JpsiTau_bbPV_refit_vx", &JpsiTau_bbPV_vx );
    tree_->Branch("JpsiTau_bbPV_refit_vy", &JpsiTau_bbPV_vy );
    tree_->Branch("JpsiTau_bbPV_refit_vz", &JpsiTau_bbPV_vz );

    tree_->Branch("JpsiTau_genPV_vx", &JpsiTau_genPV_vx );
    tree_->Branch("JpsiTau_genPV_vy", &JpsiTau_genPV_vy );
    tree_->Branch("JpsiTau_genPV_vz", &JpsiTau_genPV_vz );

    tree_->Branch("JpsiTau_Jpsi_pt", &JpsiTau_Jpsi_pt );
    tree_->Branch("JpsiTau_Jpsi_eta", &JpsiTau_Jpsi_eta );
    tree_->Branch("JpsiTau_Jpsi_phi", &JpsiTau_Jpsi_phi );
    tree_->Branch("JpsiTau_Jpsi_mass", &JpsiTau_Jpsi_mass );
    tree_->Branch("JpsiTau_Jpsi_vprob", &JpsiTau_Jpsi_vprob );
    tree_->Branch("JpsiTau_Jpsi_lip", &JpsiTau_Jpsi_lip);
    tree_->Branch("JpsiTau_Jpsi_lips", &JpsiTau_Jpsi_lips);
    tree_->Branch("JpsiTau_Jpsi_pvip", &JpsiTau_Jpsi_pvip);
    tree_->Branch("JpsiTau_Jpsi_pvips", &JpsiTau_Jpsi_pvips);
    tree_->Branch("JpsiTau_Jpsi_fl3d", &JpsiTau_Jpsi_fl3d);
    tree_->Branch("JpsiTau_Jpsi_fls3d", &JpsiTau_Jpsi_fls3d);
    tree_->Branch("JpsiTau_Jpsi_alpha", &JpsiTau_Jpsi_alpha);
    tree_->Branch("JpsiTau_Jpsi_maxdoca", &JpsiTau_Jpsi_maxdoca);
    tree_->Branch("JpsiTau_Jpsi_mindoca", &JpsiTau_Jpsi_mindoca);
    tree_->Branch("JpsiTau_Jpsi_vx", &JpsiTau_Jpsi_vx );
    tree_->Branch("JpsiTau_Jpsi_vy", &JpsiTau_Jpsi_vy );
    tree_->Branch("JpsiTau_Jpsi_vz", &JpsiTau_Jpsi_vz );
    tree_->Branch("JpsiTau_Jpsi_unfit_pt", &JpsiTau_Jpsi_unfit_pt );
    tree_->Branch("JpsiTau_Jpsi_unfit_mass", &JpsiTau_Jpsi_unfit_mass );
    tree_->Branch("JpsiTau_Jpsi_unfit_vprob", &JpsiTau_Jpsi_unfit_vprob );
    tree_->Branch("JpsiTau_Jpsi_unfit_vx", &JpsiTau_Jpsi_unfit_vx );
    tree_->Branch("JpsiTau_Jpsi_unfit_vy", &JpsiTau_Jpsi_unfit_vy );
    tree_->Branch("JpsiTau_Jpsi_unfit_vz", &JpsiTau_Jpsi_unfit_vz );


    tree_->Branch("JpsiTau_B_pt", &JpsiTau_B_pt );
    tree_->Branch("JpsiTau_B_eta", &JpsiTau_B_eta );
    tree_->Branch("JpsiTau_B_phi", &JpsiTau_B_phi );
    tree_->Branch("JpsiTau_B_mass", &JpsiTau_B_mass );
    tree_->Branch("JpsiTau_B_vprob", &JpsiTau_B_vprob );
    tree_->Branch("JpsiTau_B_lip", &JpsiTau_B_lip);
    tree_->Branch("JpsiTau_B_lips", &JpsiTau_B_lips);
    tree_->Branch("JpsiTau_B_pvip", &JpsiTau_B_pvip);
    tree_->Branch("JpsiTau_B_pvips", &JpsiTau_B_pvips);
    tree_->Branch("JpsiTau_B_fl3d", &JpsiTau_B_fl3d);
    tree_->Branch("JpsiTau_B_fls3d", &JpsiTau_B_fls3d);
    tree_->Branch("JpsiTau_B_alpha", &JpsiTau_B_alpha);
    tree_->Branch("JpsiTau_B_maxdoca", &JpsiTau_B_maxdoca);
    tree_->Branch("JpsiTau_B_mindoca", &JpsiTau_B_mindoca);
    tree_->Branch("JpsiTau_B_vx", &JpsiTau_B_vx );
    tree_->Branch("JpsiTau_B_vy", &JpsiTau_B_vy );
    tree_->Branch("JpsiTau_B_vz", &JpsiTau_B_vz );
    tree_->Branch("JpsiTau_B_iso", &JpsiTau_B_iso);
    tree_->Branch("JpsiTau_B_iso_ntracks", &JpsiTau_B_iso_ntracks );
    tree_->Branch("JpsiTau_B_iso_mindoca", &JpsiTau_B_iso_mindoca );

    tree_->Branch("JpsiTau_ngenmuons", &JpsiTau_ngenmuons);
    tree_->Branch("JpsiTau_isgenmatched", &JpsiTau_isgenmatched);
    tree_->Branch("JpsiTau_isgen3", &JpsiTau_isgen3);
    tree_->Branch("JpsiTau_isgen3matched", &JpsiTau_isgen3matched);

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

  /** energy density */
  rho = 0;
  
 

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
  IsJpsiTau.clear();

  JpsiMu_nCandidates.clear();

  JpsiMu_mu1_pt.clear();
  JpsiMu_mu1_eta.clear();
  JpsiMu_mu1_phi.clear();
  JpsiMu_mu1_mass.clear();
  JpsiMu_mu1_unfit_pt.clear();
  JpsiMu_mu1_unfit_eta.clear();
  JpsiMu_mu1_unfit_phi.clear();
  JpsiMu_mu1_unfit_mass.clear();
  JpsiMu_mu1_q.clear();
  JpsiMu_mu1_isLoose.clear();
  JpsiMu_mu1_isTight.clear();
  JpsiMu_mu1_isPF.clear();
  JpsiMu_mu1_isGlobal.clear();
  JpsiMu_mu1_isTracker.clear();
  JpsiMu_mu1_isSoft.clear();
  JpsiMu_mu1_vx.clear();
  JpsiMu_mu1_vy.clear();
  JpsiMu_mu1_vz.clear();
  JpsiMu_mu1_iso.clear();
  JpsiMu_mu1_dbiso.clear();

  JpsiMu_mu2_pt.clear();
  JpsiMu_mu2_eta.clear();
  JpsiMu_mu2_phi.clear();
  JpsiMu_mu2_mass.clear();
  JpsiMu_mu2_unfit_pt.clear();
  JpsiMu_mu2_unfit_eta.clear();
  JpsiMu_mu2_unfit_phi.clear();
  JpsiMu_mu2_unfit_mass.clear();
  JpsiMu_mu2_q.clear();
  JpsiMu_mu2_isLoose.clear();
  JpsiMu_mu2_isTight.clear();
  JpsiMu_mu2_isPF.clear();
  JpsiMu_mu2_isGlobal.clear();
  JpsiMu_mu2_isTracker.clear();
  JpsiMu_mu2_isSoft.clear();
  JpsiMu_mu2_vx.clear();
  JpsiMu_mu2_vy.clear();
  JpsiMu_mu2_vz.clear();
  JpsiMu_mu2_iso.clear();
  JpsiMu_mu2_dbiso.clear();

  JpsiMu_mu3_pt.clear();
  JpsiMu_mu3_eta.clear();
  JpsiMu_mu3_phi.clear();
  JpsiMu_mu3_mass.clear();
  JpsiMu_mu3_unfit_pt.clear();
  JpsiMu_mu3_unfit_eta.clear();
  JpsiMu_mu3_unfit_phi.clear();
  JpsiMu_mu3_unfit_mass.clear();
  JpsiMu_mu3_doca2mu1.clear();
  JpsiMu_mu3_doca2mu2.clear();
  JpsiMu_mu3_q.clear();
  JpsiMu_mu3_isLoose.clear();
  JpsiMu_mu3_isTight.clear();
  JpsiMu_mu3_isPF.clear();
  JpsiMu_mu3_isGlobal.clear();
  JpsiMu_mu3_isTracker.clear();
  JpsiMu_mu3_isSoft.clear();
  JpsiMu_mu3_vx.clear();
  JpsiMu_mu3_vy.clear();
  JpsiMu_mu3_vz.clear();
  JpsiMu_mu3_iso.clear();
  JpsiMu_mu3_dbiso.clear();

  JpsiMu_Jpsi_pt.clear();
  JpsiMu_Jpsi_eta.clear();
  JpsiMu_Jpsi_phi.clear();
  JpsiMu_Jpsi_mass.clear();
  JpsiMu_Jpsi_vprob.clear();
  JpsiMu_Jpsi_lip.clear();
  JpsiMu_Jpsi_lips.clear();
  JpsiMu_Jpsi_pvip.clear();
  JpsiMu_Jpsi_pvips.clear();
  JpsiMu_Jpsi_fl3d.clear();
  JpsiMu_Jpsi_fls3d.clear();
  JpsiMu_Jpsi_alpha.clear();
  JpsiMu_Jpsi_maxdoca.clear();
  JpsiMu_Jpsi_mindoca.clear();
  JpsiMu_Jpsi_vx.clear();
  JpsiMu_Jpsi_vy.clear();
  JpsiMu_Jpsi_vz.clear();
  JpsiMu_Jpsi_unfit_pt.clear();
  JpsiMu_Jpsi_unfit_mass.clear();
  JpsiMu_Jpsi_unfit_vprob.clear();
  JpsiMu_Jpsi_unfit_vx.clear();
  JpsiMu_Jpsi_unfit_vy.clear();
  JpsiMu_Jpsi_unfit_vz.clear();

  JpsiMu_B_pt.clear();
  JpsiMu_B_eta.clear();
  JpsiMu_B_phi.clear();
  JpsiMu_B_mass.clear();
  JpsiMu_B_vprob.clear();
  JpsiMu_B_lip.clear();
  JpsiMu_B_lips.clear();
  JpsiMu_B_pvip.clear();
  JpsiMu_B_pvips.clear();
  JpsiMu_B_fl3d.clear();
  JpsiMu_B_fls3d.clear();
  JpsiMu_B_alpha.clear();
  JpsiMu_B_maxdoca.clear();
  JpsiMu_B_mindoca.clear();
  JpsiMu_B_vx.clear();
  JpsiMu_B_vy.clear();
  JpsiMu_B_vz.clear();
  JpsiMu_B_iso.clear();
  JpsiMu_B_iso_ntracks.clear();
  JpsiMu_B_iso_mindoca.clear();
  JpsiMu_B_unfit_pt.clear();
  JpsiMu_B_unfit_mass.clear();
  JpsiMu_B_unfit_vprob.clear();
  JpsiMu_B_unfit_vx.clear();
  JpsiMu_B_unfit_vy.clear();
  JpsiMu_B_unfit_vz.clear();

  JpsiMu_PV_vx.clear();
  JpsiMu_PV_vy.clear();
  JpsiMu_PV_vz.clear();

  JpsiMu_bbPV_vx.clear();
  JpsiMu_bbPV_vy.clear();
  JpsiMu_bbPV_vz.clear();

  JpsiMu_bbPV_refit_vx.clear();
  JpsiMu_bbPV_refit_vy.clear();
  JpsiMu_bbPV_refit_vz.clear();

  JpsiMu_genPV_vx.clear();
  JpsiMu_genPV_vy.clear();
  JpsiMu_genPV_vz.clear();

  JpsiMu_ngenmuons.clear();
  JpsiMu_isgenmatched.clear();
  JpsiMu_mu3_isgenmatched.clear();






  JpsiTau_nCandidates.clear();

  JpsiTau_mu1_pt.clear();
  JpsiTau_mu1_eta.clear();
  JpsiTau_mu1_phi.clear();
  JpsiTau_mu1_mass.clear();
  JpsiTau_mu1_unfit_pt.clear();
  JpsiTau_mu1_unfit_eta.clear();
  JpsiTau_mu1_unfit_phi.clear();
  JpsiTau_mu1_unfit_mass.clear();
  JpsiTau_mu1_q.clear();
  JpsiTau_mu1_isLoose.clear();
  JpsiTau_mu1_isTight.clear();
  JpsiTau_mu1_isPF.clear();
  JpsiTau_mu1_isGlobal.clear();
  JpsiTau_mu1_isTracker.clear();
  JpsiTau_mu1_isSoft.clear();
  JpsiTau_mu1_vx.clear();
  JpsiTau_mu1_vy.clear();
  JpsiTau_mu1_vz.clear();
  JpsiTau_mu1_iso.clear();
  JpsiTau_mu1_dbiso.clear();

  JpsiTau_mu2_pt.clear();
  JpsiTau_mu2_eta.clear();
  JpsiTau_mu2_phi.clear();
  JpsiTau_mu2_mass.clear();
  JpsiTau_mu2_unfit_pt.clear();
  JpsiTau_mu2_unfit_eta.clear();
  JpsiTau_mu2_unfit_phi.clear();
  JpsiTau_mu2_unfit_mass.clear();
  JpsiTau_mu2_q.clear();
  JpsiTau_mu2_isLoose.clear();
  JpsiTau_mu2_isTight.clear();
  JpsiTau_mu2_isPF.clear();
  JpsiTau_mu2_isGlobal.clear();
  JpsiTau_mu2_isTracker.clear();
  JpsiTau_mu2_isSoft.clear();
  JpsiTau_mu2_vx.clear();
  JpsiTau_mu2_vy.clear();
  JpsiTau_mu2_vz.clear();
  JpsiTau_mu2_iso.clear();
  JpsiTau_mu2_dbiso.clear();

  JpsiTau_tau_pt.clear();
  JpsiTau_tau_eta.clear();
  JpsiTau_tau_phi.clear();
  JpsiTau_tau_mass.clear();
  JpsiTau_tau_q.clear();
  JpsiTau_tau_vx.clear();
  JpsiTau_tau_vy.clear();
  JpsiTau_tau_vz.clear();
  JpsiTau_tau_iso.clear();

  JpsiTau_Jpsi_pt.clear();
  JpsiTau_Jpsi_eta.clear();
  JpsiTau_Jpsi_phi.clear();
  JpsiTau_Jpsi_mass.clear();
  JpsiTau_Jpsi_vprob.clear();
  JpsiTau_Jpsi_lip.clear();
  JpsiTau_Jpsi_lips.clear();
  JpsiTau_Jpsi_pvip.clear();
  JpsiTau_Jpsi_pvips.clear();
  JpsiTau_Jpsi_fl3d.clear();
  JpsiTau_Jpsi_fls3d.clear();
  JpsiTau_Jpsi_alpha.clear();
  JpsiTau_Jpsi_maxdoca.clear();
  JpsiTau_Jpsi_mindoca.clear();
  JpsiTau_Jpsi_vx.clear();
  JpsiTau_Jpsi_vy.clear();
  JpsiTau_Jpsi_vz.clear();
  JpsiTau_Jpsi_unfit_pt.clear();
  JpsiTau_Jpsi_unfit_mass.clear();
  JpsiTau_Jpsi_unfit_vprob.clear();
  JpsiTau_Jpsi_unfit_vx.clear();
  JpsiTau_Jpsi_unfit_vy.clear();
  JpsiTau_Jpsi_unfit_vz.clear();

  JpsiTau_B_pt.clear();
  JpsiTau_B_eta.clear();
  JpsiTau_B_phi.clear();
  JpsiTau_B_mass.clear();
  JpsiTau_B_vprob.clear();
  JpsiTau_B_lip.clear();
  JpsiTau_B_lips.clear();
  JpsiTau_B_pvip.clear();
  JpsiTau_B_pvips.clear();
  JpsiTau_B_fl3d.clear();
  JpsiTau_B_fls3d.clear();
  JpsiTau_B_alpha.clear();
  JpsiTau_B_maxdoca.clear();
  JpsiTau_B_mindoca.clear();
  JpsiTau_B_vx.clear();
  JpsiTau_B_vy.clear();
  JpsiTau_B_vz.clear();
  JpsiTau_B_iso.clear();
  JpsiTau_B_iso_ntracks.clear();
  JpsiTau_B_iso_mindoca.clear();

  JpsiTau_PV_vx.clear();
  JpsiTau_PV_vy.clear();
  JpsiTau_PV_vz.clear();

  JpsiTau_bbPV_vx.clear();
  JpsiTau_bbPV_vy.clear();
  JpsiTau_bbPV_vz.clear();

  JpsiTau_bbPV_refit_vx.clear();
  JpsiTau_bbPV_refit_vy.clear();
  JpsiTau_bbPV_refit_vz.clear();

  JpsiTau_genPV_vx.clear();
  JpsiTau_genPV_vy.clear();
  JpsiTau_genPV_vz.clear();

  JpsiTau_ngenmuons.clear();
  JpsiTau_isgenmatched.clear();

  JpsiTau_isgen3.clear();
  JpsiTau_isgen3matched.clear();






 
} 


void NtupleBranches::LabelHistograms( std::map< std::string, bool >& runFlags ){
    if ( runFlags["doCutFlow"] ){
        std::vector bins_string = {"Precut", "p_{T} & #eta (#mu_{1})","p_{T} & #eta (#mu_{2})","Trigger","Trig matched #mu_{1}","Trig matched #mu_{2}","J/#psi mass", "J/#psi kinematic","p_{T}(#mu_{3})","B kinematic"};
        for (int i=0; i< cutflow_perevt->GetNbinsX(); i++){
            cutflow_perevt->GetXaxis()->SetBinLabel(i+1, bins_string[i]);
        }
    }

    if ( runFlags["doGenHist"] ){
        std::vector bins_string = { "#mu", "#pi^{0}", "#pi^{#pm}","#rho^{0}","#rho^{#pm}","#eta","#eta^{`}","#omega","#phi","K^{0}","K^{#pm}","K^{*0}","K^{*#pm}","D^{#pm}","D^{0}","#eta_{c}","#eta_{b}","#Upsilon"};
        for (int i=0; i< genParticle_Bdau_X_id->GetNbinsX(); i++){
            genParticle_Bdau_X_id->GetXaxis()->SetBinLabel(i+1, bins_string[i]);
        }       
    }
}
