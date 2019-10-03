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
  if (runFlags["doJpsiMu"]){
    tree_->Branch("IsJpsiMu" , &IsJpsiMu  );
    tree_->Branch("IsJpsiEle", &IsJpsiEle );

    tree_->Branch("Jpsi_mu1_pt", &Jpsi_mu1_pt );
    tree_->Branch("Jpsi_mu1_eta", &Jpsi_mu1_eta );
    tree_->Branch("Jpsi_mu1_phi", &Jpsi_mu1_phi );
    tree_->Branch("Jpsi_mu1_q", &Jpsi_mu1_q );
    tree_->Branch("Jpsi_mu1_isLoose"  , &Jpsi_mu1_isLoose   );
    tree_->Branch("Jpsi_mu1_isTight"  , &Jpsi_mu1_isTight   );
    tree_->Branch("Jpsi_mu1_isPF"     , &Jpsi_mu1_isPF      );
    tree_->Branch("Jpsi_mu1_isGlobal" , &Jpsi_mu1_isGlobal  );
    tree_->Branch("Jpsi_mu1_isTracker", &Jpsi_mu1_isTracker );
    tree_->Branch("Jpsi_mu1_isSoft"   , &Jpsi_mu1_isSoft    );
    tree_->Branch("Jpsi_mu1_vx"   , &Jpsi_mu1_vx    );
    tree_->Branch("Jpsi_mu1_vy"   , &Jpsi_mu1_vy    );
    tree_->Branch("Jpsi_mu1_vz"   , &Jpsi_mu1_vz    );




    tree_->Branch("Jpsi_mu2_pt", &Jpsi_mu2_pt );
    tree_->Branch("Jpsi_mu2_eta", &Jpsi_mu2_eta );
    tree_->Branch("Jpsi_mu2_phi", &Jpsi_mu2_phi );
    tree_->Branch("Jpsi_mu2_q", &Jpsi_mu2_q );
    tree_->Branch("Jpsi_mu2_isLoose"  , &Jpsi_mu2_isLoose   );
    tree_->Branch("Jpsi_mu2_isTight"  , &Jpsi_mu2_isTight   );
    tree_->Branch("Jpsi_mu2_isPF"     , &Jpsi_mu2_isPF      );
    tree_->Branch("Jpsi_mu2_isGlobal" , &Jpsi_mu2_isGlobal  );
    tree_->Branch("Jpsi_mu2_isTracker", &Jpsi_mu2_isTracker );
    tree_->Branch("Jpsi_mu2_isSoft"   , &Jpsi_mu2_isSoft    );
    tree_->Branch("Jpsi_mu2_vx"   , &Jpsi_mu2_vx    );
    tree_->Branch("Jpsi_mu2_vy"   , &Jpsi_mu2_vy    );
    tree_->Branch("Jpsi_mu2_vz"   , &Jpsi_mu2_vz    );


    tree_->Branch("Jpsi_mu3_pt", &Jpsi_mu3_pt );
    tree_->Branch("Jpsi_mu3_eta", &Jpsi_mu3_eta );
    tree_->Branch("Jpsi_mu3_phi", &Jpsi_mu3_phi );
    tree_->Branch("Jpsi_mu3_q", &Jpsi_mu3_q );
    tree_->Branch("Jpsi_mu3_isLoose"  , &Jpsi_mu3_isLoose   );
    tree_->Branch("Jpsi_mu3_isTight"  , &Jpsi_mu3_isTight   );
    tree_->Branch("Jpsi_mu3_isPF"     , &Jpsi_mu3_isPF      );
    tree_->Branch("Jpsi_mu3_isGlobal" , &Jpsi_mu3_isGlobal  );
    tree_->Branch("Jpsi_mu3_isTracker", &Jpsi_mu3_isTracker );
    tree_->Branch("Jpsi_mu3_isSoft"   , &Jpsi_mu3_isSoft    );
    tree_->Branch("Jpsi_mu3_vx"   , &Jpsi_mu3_vx    );
    tree_->Branch("Jpsi_mu3_vy"   , &Jpsi_mu3_vy    );
    tree_->Branch("Jpsi_mu3_vz"   , &Jpsi_mu3_vz    );

    tree_->Branch("Jpsi_mu3_isopt03", &Jpsi_mu3_isopt03 );
    tree_->Branch("Jpsi_mu3_isopt04", &Jpsi_mu3_isopt04 );
    tree_->Branch("Jpsi_mu3_isopt05", &Jpsi_mu3_isopt05 );
    tree_->Branch("Jpsi_mu3_isopt06", &Jpsi_mu3_isopt06 );
    tree_->Branch("Jpsi_mu3_isopt07", &Jpsi_mu3_isopt07 );
    tree_->Branch("Jpsi_dr_mu3pf"    , &Jpsi_dr_mu3pf     );

    tree_->Branch("Jpsi_PV_vx", &Jpsi_PV_vx );
    tree_->Branch("Jpsi_PV_vy", &Jpsi_PV_vy );
    tree_->Branch("Jpsi_PV_vz", &Jpsi_PV_vz );

    tree_->Branch("Jpsi_bbPV_vx", &Jpsi_bbPV_vx );
    tree_->Branch("Jpsi_bbPV_vy", &Jpsi_bbPV_vy );
    tree_->Branch("Jpsi_bbPV_vz", &Jpsi_bbPV_vz );

    tree_->Branch("Jpsi_genPV_vx", &Jpsi_genPV_vx );
    tree_->Branch("Jpsi_genPV_vy", &Jpsi_genPV_vy );
    tree_->Branch("Jpsi_genPV_vz", &Jpsi_genPV_vz );

    tree_->Branch("Jpsi_pt", &Jpsi_pt );
    tree_->Branch("Jpsi_eta", &Jpsi_eta );
    tree_->Branch("Jpsi_phi", &Jpsi_phi );
    tree_->Branch("Jpsi_mass", &Jpsi_mass );
    tree_->Branch("Jpsi_vprob", &Jpsi_vprob );
    tree_->Branch("Jpsi_lip", &Jpsi_lip);
    tree_->Branch("Jpsi_lips", &Jpsi_lips);
    tree_->Branch("Jpsi_pvip", &Jpsi_pvip);
    tree_->Branch("Jpsi_pvips", &Jpsi_pvips);
    tree_->Branch("Jpsi_fl3d", &Jpsi_fl3d);
    tree_->Branch("Jpsi_fls3d", &Jpsi_fls3d);
    tree_->Branch("Jpsi_alpha", &Jpsi_alpha);
    tree_->Branch("Jpsi_maxdoca", &Jpsi_maxdoca);
    tree_->Branch("Jpsi_mindoca", &Jpsi_mindoca);
    tree_->Branch("Jpsi_vx", &Jpsi_vx );
    tree_->Branch("Jpsi_vy", &Jpsi_vy );
    tree_->Branch("Jpsi_vz", &Jpsi_vz );
    tree_->Branch("Jpsi_unfitpt", &Jpsi_unfitpt );
    tree_->Branch("Jpsi_unfitmass", &Jpsi_unfitmass );
    tree_->Branch("Jpsi_unfitvprob", &Jpsi_unfitvprob );
    tree_->Branch("Jpsi_unfit_vx", &Jpsi_unfit_vx );
    tree_->Branch("Jpsi_unfit_vy", &Jpsi_unfit_vy );
    tree_->Branch("Jpsi_unfit_vz", &Jpsi_unfit_vz );


    tree_->Branch("Jpsi_trimu_pt", &Jpsi_trimu_pt );
    tree_->Branch("Jpsi_trimu_eta", &Jpsi_trimu_eta );
    tree_->Branch("Jpsi_trimu_phi", &Jpsi_trimu_phi );
    tree_->Branch("Jpsi_trimu_mass", &Jpsi_trimu_mass );
    tree_->Branch("Jpsi_trimu_vprob", &Jpsi_trimu_vprob );
    tree_->Branch("Jpsi_trimu_lip", &Jpsi_trimu_lip);
    tree_->Branch("Jpsi_trimu_lips", &Jpsi_trimu_lips);
    tree_->Branch("Jpsi_trimu_pvip", &Jpsi_trimu_pvip);
    tree_->Branch("Jpsi_trimu_pvips", &Jpsi_trimu_pvips);
    tree_->Branch("Jpsi_trimu_fl3d", &Jpsi_trimu_fl3d);
    tree_->Branch("Jpsi_trimu_fls3d", &Jpsi_trimu_fls3d);
    tree_->Branch("Jpsi_trimu_alpha", &Jpsi_trimu_alpha);
    tree_->Branch("Jpsi_trimu_maxdoca", &Jpsi_trimu_maxdoca);
    tree_->Branch("Jpsi_trimu_mindoca", &Jpsi_trimu_mindoca);
    tree_->Branch("Jpsi_trimu_vx", &Jpsi_trimu_vx );
    tree_->Branch("Jpsi_trimu_vy", &Jpsi_trimu_vy );
    tree_->Branch("Jpsi_trimu_vz", &Jpsi_trimu_vz );
    tree_->Branch("Jpsi_trimu_unfitpt", &Jpsi_trimu_unfitpt );
    tree_->Branch("Jpsi_trimu_unfitmass", &Jpsi_trimu_unfitmass );
    tree_->Branch("Jpsi_trimu_unfitvprob", &Jpsi_trimu_unfitvprob );
    tree_->Branch("Jpsi_trimu_unfit_vx", &Jpsi_trimu_unfit_vx );
    tree_->Branch("Jpsi_trimu_unfit_vy", &Jpsi_trimu_unfit_vy );
    tree_->Branch("Jpsi_trimu_unfit_vz", &Jpsi_trimu_unfit_vz );

    tree_->Branch("Jpsi_ngenmuons", &Jpsi_ngenmuons);
    tree_->Branch("Jpsi_isgenmatched", &Jpsi_isgenmatched);
    tree_->Branch("Jpsi_mu3_isgenmatched", &Jpsi_mu3_isgenmatched);

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

  Jpsi_mu1_pt.clear();
  Jpsi_mu1_eta.clear();
  Jpsi_mu1_phi.clear();
  Jpsi_mu1_q.clear();
  Jpsi_mu1_isLoose.clear();
  Jpsi_mu1_isTight.clear();
  Jpsi_mu1_isPF.clear();
  Jpsi_mu1_isGlobal.clear();
  Jpsi_mu1_isTracker.clear();
  Jpsi_mu1_isSoft.clear();
  Jpsi_mu1_vx.clear();
  Jpsi_mu1_vy.clear();
  Jpsi_mu1_vz.clear();


  Jpsi_mu2_pt.clear();
  Jpsi_mu2_eta.clear();
  Jpsi_mu2_phi.clear();
  Jpsi_mu2_q.clear();
  Jpsi_mu2_isLoose.clear();
  Jpsi_mu2_isTight.clear();
  Jpsi_mu2_isPF.clear();
  Jpsi_mu2_isGlobal.clear();
  Jpsi_mu2_isTracker.clear();
  Jpsi_mu2_isSoft.clear();
  Jpsi_mu2_vx.clear();
  Jpsi_mu2_vy.clear();
  Jpsi_mu2_vz.clear();

  Jpsi_mu3_pt.clear();
  Jpsi_mu3_eta.clear();
  Jpsi_mu3_phi.clear();
  Jpsi_mu3_q.clear();
  Jpsi_mu3_isLoose.clear();
  Jpsi_mu3_isTight.clear();
  Jpsi_mu3_isPF.clear();
  Jpsi_mu3_isGlobal.clear();
  Jpsi_mu3_isTracker.clear();
  Jpsi_mu3_isSoft.clear();
  Jpsi_mu3_vx.clear();
  Jpsi_mu3_vy.clear();
  Jpsi_mu3_vz.clear();


  Jpsi_mu3_isopt03.clear();
  Jpsi_mu3_isopt04.clear();
  Jpsi_mu3_isopt05.clear();
  Jpsi_mu3_isopt06.clear();
  Jpsi_mu3_isopt07.clear();


  Jpsi_pt.clear();
  Jpsi_eta.clear();
  Jpsi_phi.clear();
  Jpsi_mass.clear();
  Jpsi_vprob.clear();
  Jpsi_lip.clear();
  Jpsi_lips.clear();
  Jpsi_pvip.clear();
  Jpsi_pvips.clear();
  Jpsi_fl3d.clear();
  Jpsi_fls3d.clear();
  Jpsi_alpha.clear();
  Jpsi_maxdoca.clear();
  Jpsi_mindoca.clear();
  Jpsi_vx.clear();
  Jpsi_vy.clear();
  Jpsi_vz.clear();
  Jpsi_unfitpt.clear();
  Jpsi_unfitmass.clear();
  Jpsi_unfitvprob.clear();
  Jpsi_unfit_vx.clear();
  Jpsi_unfit_vy.clear();
  Jpsi_unfit_vz.clear();

  Jpsi_trimu_pt.clear();
  Jpsi_trimu_eta.clear();
  Jpsi_trimu_phi.clear();
  Jpsi_trimu_mass.clear();
  Jpsi_trimu_vprob.clear();
  Jpsi_trimu_lip.clear();
  Jpsi_trimu_lips.clear();
  Jpsi_trimu_pvip.clear();
  Jpsi_trimu_pvips.clear();
  Jpsi_trimu_fl3d.clear();
  Jpsi_trimu_fls3d.clear();
  Jpsi_trimu_alpha.clear();
  Jpsi_trimu_maxdoca.clear();
  Jpsi_trimu_mindoca.clear();
  Jpsi_trimu_vx.clear();
  Jpsi_trimu_vy.clear();
  Jpsi_trimu_vz.clear();
  Jpsi_trimu_unfitpt.clear();
  Jpsi_trimu_unfitmass.clear();
  Jpsi_trimu_unfitvprob.clear();
  Jpsi_trimu_unfit_vx.clear();
  Jpsi_trimu_unfit_vy.clear();
  Jpsi_trimu_unfit_vz.clear();

  Jpsi_PV_vx.clear();
  Jpsi_PV_vy.clear();
  Jpsi_PV_vz.clear();

  Jpsi_bbPV_vx.clear();
  Jpsi_bbPV_vy.clear();
  Jpsi_bbPV_vz.clear();

  Jpsi_genPV_vx.clear();
  Jpsi_genPV_vy.clear();
  Jpsi_genPV_vz.clear();

  Jpsi_dr_mu3pf.clear();

  Jpsi_ngenmuons.clear();
  Jpsi_isgenmatched.clear();
  Jpsi_mu3_isgenmatched.clear();

 
} 
