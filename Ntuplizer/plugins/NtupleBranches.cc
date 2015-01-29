#include "../interface/NtupleBranches.h"

//===================================================================================================================        
NtupleBranches::NtupleBranches( TTree* tree )
   : tree_( tree )
{
   branch();

}

//===================================================================================================================
NtupleBranches::~NtupleBranches( void )
{
}

//===================================================================================================================      
void NtupleBranches::branch( void ){

  /*----------------------gen particles-------------------------*/
  tree_->Branch( "genParticle_pt"	     , &genParticle_pt	       ); 
  tree_->Branch( "genParticle_px"	     , &genParticle_px	       ); 
  tree_->Branch( "genParticle_py"	     , &genParticle_py	       ); 
  tree_->Branch( "genParticle_pz"	     , &genParticle_pz	       ); 
  tree_->Branch( "genParticle_e" 	     , &genParticle_e	       ); 
  tree_->Branch( "genParticle_eta"	     , &genParticle_eta        ); 
  tree_->Branch( "genParticle_phi"	     , &genParticle_phi        ); 
  tree_->Branch( "genParticle_mass"	     , &genParticle_mass       ); 
  tree_->Branch( "genParticle_pdgId"	     , &genParticle_pdgId      );
  tree_->Branch( "genParticle_status"	     , &genParticle_status     );
  tree_->Branch( "genParticle_mother"	     , &genParticle_mother     );
  tree_->Branch( "genParticle_nMoth"	     , &genParticle_nMoth      );
  tree_->Branch( "genParticle_nDau"	     , &genParticle_nDau       ); 
  tree_->Branch( "genParticle_dau"	     , &genParticle_dau        ); 
  tree_->Branch( "lheV_pt"	             , &lheV_pt                ); 
  tree_->Branch( "lheHT"	             , &lheHT                  ); 
  tree_->Branch( "lheNj"	             , &lheNj                  ); 
  
  /*-------------------------leptons----------------------------*/
  tree_->Branch( "nlep"  		     , &nlep		       );
  tree_->Branch( "lep_type"		     , &lep_type	       );
  tree_->Branch( "lep_charge"		     , &lep_charge	       );
  tree_->Branch( "lep_e" 		     , &lep_e		       );
  tree_->Branch( "lep_eta"		     , &lep_eta 	       );
  tree_->Branch( "lep_mass"		     , &lep_mass	       );
  tree_->Branch( "lep_pt"		     , &lep_pt  	       );
  tree_->Branch( "lep_phi"		     , &lep_phi 	       );
  tree_->Branch( "lep_isHEEP"	             , &lep_isHEEP             );
  tree_->Branch( "lep_isHighPtMuon"	     , &lep_isHighPtMuon       );
  tree_->Branch( "lep_pfRhoCorrRelIso03"     , &lep_pfRhoCorrRelIso03  );
  tree_->Branch( "lep_pfRhoCorrRelIso04"     , &lep_pfRhoCorrRelIso04  );
  tree_->Branch( "lep_pfDeltaCorrRelIso"     , &lep_pfDeltaCorrRelIso  );
  tree_->Branch( "lep_pfRelIso"  	     , &lep_pfRelIso	       );
  tree_->Branch( "lep_photonIso" 	     , &lep_photonIso	       );
  tree_->Branch( "lep_neutralHadIso"	     , &lep_neutralHadIso      );
  tree_->Branch( "lep_chargedHadIso"	     , &lep_chargedHadIso      );
  tree_->Branch( "lep_trackIso"	             , &lep_trackIso	       );
      
  /*-------------------------AK5 jets---------------------------*/	 
  tree_->Branch( "njetsAK5"		    , &njetsAK5 	    );
  tree_->Branch( "jetAK5_pt"		    , &jetAK5_pt	    );
  tree_->Branch( "jetAK5_eta"		    , &jetAK5_eta	    );
  tree_->Branch( "jetAK5_mass"		    , &jetAK5_mass	    );
  tree_->Branch( "jetAK5_phi"		    , &jetAK5_phi	    );
  tree_->Branch( "jetAK5_e"		    , &jetAK5_e 	    );
  tree_->Branch( "jetAK5_jec"		    , &jetAK5_jec	    );
  tree_->Branch( "jetAK5_jecUp"		    , &jetAK5_jecUp	    );
  tree_->Branch( "jetAK5_jecDown"	    , &jetAK5_jecDown       );
  tree_->Branch( "jetAK5_IDLoose"	    , &jetAK5_IDLoose	    );
  tree_->Branch( "jetAK5_cm"     	    , &jetAK5_cm	    );
  tree_->Branch( "jetAK5_nm"     	    , &jetAK5_nm	    );
  tree_->Branch( "jetAK5_che"     	    , &jetAK5_che	    );
  tree_->Branch( "jetAK5_ne"     	    , &jetAK5_ne	    );
  tree_->Branch( "jetAK5_charge" 	    , &jetAK5_charge        );
  tree_->Branch( "jetAK5_flavour"	    , &jetAK5_flavour       );
  tree_->Branch( "jetAK5_ssv"		    , &jetAK5_ssv	    );
  tree_->Branch( "jetAK5_csv"		    , &jetAK5_csv	    );
  tree_->Branch( "jetAK5_tchp"		    , &jetAK5_tchp	    );
  tree_->Branch( "jetAK5_tche"		    , &jetAK5_tche	    );
  tree_->Branch( "jetAK5_jp"		    , &jetAK5_jp	    );
  tree_->Branch( "jetAK5_jbp"		    , &jetAK5_jbp	    );
  tree_->Branch( "jetAK5_nSVs"		    , &jetAK5_nSVs	    );
        
  /*-------------------------CA8 jets---------------------------*/    
  tree_->Branch( "njetsCA8"		    , &njetsCA8 	    );
  tree_->Branch( "jetCA8_pt"		    , &jetCA8_pt	    );
  tree_->Branch( "jetCA8_eta"		    , &jetCA8_eta	    );
  tree_->Branch( "jetCA8_mass"		    , &jetCA8_mass	    );
  tree_->Branch( "jetCA8_phi"		    , &jetCA8_phi	    );
  tree_->Branch( "jetCA8_e"		    , &jetCA8_e 	    );
  tree_->Branch( "jetCA8_jec"		    , &jetCA8_jec 	    );
  tree_->Branch( "jetCA8_jecUp"		    , &jetCA8_jecUp	    );
  tree_->Branch( "jetCA8_jecDown"	    , &jetCA8_jecDown       );  
  tree_->Branch( "jetCA8_IDLoose"	    , &jetCA8_IDLoose	    );
  tree_->Branch( "jetCA8_cm"     	    , &jetCA8_cm	    );
  tree_->Branch( "jetCA8_nm"     	    , &jetCA8_nm	    );
  tree_->Branch( "jetCA8_che"     	    , &jetCA8_che	    );
  tree_->Branch( "jetCA8_ne"     	    , &jetCA8_ne	    );
  tree_->Branch( "jetCA8_charge" 	    , &jetCA8_charge        );
  tree_->Branch( "jetCA8_flavour"	    , &jetCA8_flavour       );
  tree_->Branch( "jetCA8_ssv"		    , &jetCA8_ssv	    );
  tree_->Branch( "jetCA8_csv"		    , &jetCA8_csv	    );
  tree_->Branch( "jetCA8_tchp"		    , &jetCA8_tchp	    );
  tree_->Branch( "jetCA8_tche"		    , &jetCA8_tche	    );
  tree_->Branch( "jetCA8_jp"		    , &jetCA8_jp	    );
  tree_->Branch( "jetCA8_jbp"		    , &jetCA8_jbp	    );
  tree_->Branch( "jetCA8_nSVs"		    , &jetCA8_nSVs	    );
  tree_->Branch( "jetCA8_tau1"		    , &jetCA8_tau1	    );
  tree_->Branch( "jetCA8_tau2"		    , &jetCA8_tau2	    );
  tree_->Branch( "jetCA8_tau3"		    , &jetCA8_tau3	    );
  
  /*----------------------CA8 jets pruned-----------------------*/    
  tree_->Branch( "njetsCA8pruned"	    , &njetsCA8pruned	     );
  tree_->Branch( "jetCA8pruned_pt"	    , &jetCA8pruned_pt       );
  tree_->Branch( "jetCA8pruned_eta"	    , &jetCA8pruned_eta      );
  tree_->Branch( "jetCA8pruned_mass"	    , &jetCA8pruned_mass     );
  tree_->Branch( "jetCA8pruned_phi"	    , &jetCA8pruned_phi      );
  tree_->Branch( "jetCA8pruned_e"	    , &jetCA8pruned_e	     );
  tree_->Branch( "jetCA8pruned_jec"	    , &jetCA8pruned_jec      );
  tree_->Branch( "jetCA8pruned_jecUp"	    , &jetCA8pruned_jecUp    );
  tree_->Branch( "jetCA8pruned_jecDown"	    , &jetCA8pruned_jecDown  );
  tree_->Branch( "jetCA8pruned_flavour"     , &jetCA8pruned_flavour  );
  tree_->Branch( "jetCA8pruned_charge"	    , &jetCA8pruned_charge   );
  tree_->Branch( "jetCA8pruned_ssv"	    , &jetCA8pruned_ssv      );
  tree_->Branch( "jetCA8pruned_csv"	    , &jetCA8pruned_csv      );
  tree_->Branch( "jetCA8pruned_tchp"	    , &jetCA8pruned_tchp     );
  tree_->Branch( "jetCA8pruned_tche"	    , &jetCA8pruned_tche     );
  tree_->Branch( "jetCA8pruned_jp"	    , &jetCA8pruned_jp       );
  tree_->Branch( "jetCA8pruned_jbp"	    , &jetCA8pruned_jbp      );
  tree_->Branch( "jetCA8pruned_nSVs"	    , &jetCA8pruned_nSVs     ); 
	  
  /*----------------------CA8 subjets---------------------------*/    
  tree_->Branch( "nsubjets"		    , &nsubjets 	       );
  tree_->Branch( "subjetCA8pruned_pt"	    , &subjetCA8pruned_pt      );
  tree_->Branch( "subjetCA8pruned_eta"	    , &subjetCA8pruned_eta     );
  tree_->Branch( "subjetCA8pruned_mass"     , &subjetCA8pruned_mass    );
  tree_->Branch( "subjetCA8pruned_phi"	    , &subjetCA8pruned_phi     );
  tree_->Branch( "subjetCA8pruned_e"	    , &subjetCA8pruned_e       );
  tree_->Branch( "subjetCA8pruned_charge"   , &subjetCA8pruned_charge  );
  tree_->Branch( "subjetCA8pruned_flavour"  , &subjetCA8pruned_flavour );
  tree_->Branch( "subjetCA8pruned_ssv"	    , &subjetCA8pruned_ssv     );
  tree_->Branch( "subjetCA8pruned_csv"	    , &subjetCA8pruned_csv     );
  tree_->Branch( "subjetCA8pruned_tchp"     , &subjetCA8pruned_tchp    );
  tree_->Branch( "subjetCA8pruned_tche"     , &subjetCA8pruned_tche    );
  tree_->Branch( "subjetCA8pruned_jp"	    , &subjetCA8pruned_jp      );
  tree_->Branch( "subjetCA8pruned_jbp"	    , &subjetCA8pruned_jbp     );
  
  /*----------------------HLT trigger---------------------------*/	  
  tree_->Branch("isFired_HLT_Mu40_eta2p1_v9"             , &isFired_HLT_Mu40_eta2p1_v9              );
  tree_->Branch("isFired_HLT_Mu40_eta2p1_v10"            , &isFired_HLT_Mu40_eta2p1_v10             );
  tree_->Branch("isFired_HLT_Mu40_eta2p1_v11"            , &isFired_HLT_Mu40_eta2p1_v11             );  
  tree_->Branch("isFired_HLT_Ele80_CaloIdVT_TrkIdT_v8"   , &isFired_HLT_Ele80_CaloIdVT_TrkIdT_v8    );
  tree_->Branch("isFired_HLT_Ele80_CaloIdVT_TrkIdT_v9"   , &isFired_HLT_Ele80_CaloIdVT_TrkIdT_v9    );
  tree_->Branch("isFired_HLT_Ele80_CaloIdVT_TrkIdT_v10"  , &isFired_HLT_Ele80_CaloIdVT_TrkIdT_v10   );
  tree_->Branch("isFired_HLT_Ele80_CaloIdVT_GsfTrkIdT_v1", &isFired_HLT_Ele80_CaloIdVT_GsfTrkIdT_v1 );
  tree_->Branch("isFired_HLT_Ele80_CaloIdVT_GsfTrkIdT_v2", &isFired_HLT_Ele80_CaloIdVT_GsfTrkIdT_v2 );
  tree_->Branch("isFired_HLT_PFHT650_v5"		 , &isFired_HLT_PFHT650_v5		    );
  tree_->Branch("isFired_HLT_PFHT650_v6"		 , &isFired_HLT_PFHT650_v6		    );
  tree_->Branch("isFired_HLT_PFHT650_v7"		 , &isFired_HLT_PFHT650_v7		    );
  tree_->Branch("isFired_HLT_PFHT650_v8"		 , &isFired_HLT_PFHT650_v8		    );
  tree_->Branch("isFired_HLT_PFHT650_v9"		 , &isFired_HLT_PFHT650_v9		    );
  tree_->Branch("isFired_HLT_PFNoPUHT650_v1"             , &isFired_HLT_PFNoPUHT650_v1		    );
  tree_->Branch("isFired_HLT_PFNoPUHT650_v3"             , &isFired_HLT_PFNoPUHT650_v3		    );
  tree_->Branch("isFired_HLT_PFNoPUHT650_v4"             , &isFired_HLT_PFNoPUHT650_v4		    );
  tree_->Branch("isFired_HLT_PFJet320_v3"		 , &isFired_HLT_PFJet320_v3		    );
  tree_->Branch("isFired_HLT_PFJet320_v4"		 , &isFired_HLT_PFJet320_v4		    );
  tree_->Branch("isFired_HLT_PFJet320_v5"		 , &isFired_HLT_PFJet320_v5		    );
  tree_->Branch("isFired_HLT_PFJet320_v6"		 , &isFired_HLT_PFJet320_v6		    );
  tree_->Branch("isFired_HLT_PFJet320_v8"		 , &isFired_HLT_PFJet320_v8		    );
  tree_->Branch("isFired_HLT_PFJet320_v9"  		 , &isFired_HLT_PFJet320_v9  		    );

  /*-------------------------MET--------------------------------*/
  tree_->Branch("METraw_et"		        , &METraw_et	     );
  tree_->Branch("METraw_phi"		        , &METraw_phi	     ); 	
  tree_->Branch("MET_et"	                , &MET_et  	     ); 
  tree_->Branch("MET_phi"	                , &MET_phi           );

  /*------------- ------EVENT infos-----------------------------*/
  tree_->Branch("EVENT_event"		        , &EVENT_event	     );
  tree_->Branch("EVENT_run"		        , &EVENT_run	     );
  tree_->Branch("EVENT_lumiBlock"	        , &EVENT_lumiBlock   );
  
  /*--------------------------PU infos--------------------------*/  			         
  tree_->Branch("nPuVtxTrue"                    , &nPuVtxTrue        );    
  tree_->Branch("nPuVtx"	                , &nPuVtx	     );
  tree_->Branch("bX"	    		        , &bX	  	     );
  
  /*--------------------------PV infos--------------------------*/
  tree_->Branch("nPVs"                          , &nPVs              );
  
}

//=================================================================================================================== 
void NtupleBranches::reset( void ){

  lheV_pt     = 0;
  lheHT       = 0;
  lheNj       = 0;
  nlep        = 0;
  njetsAK5    = 0;
  njetsCA8    = 0;  
  njetsCA8pruned = 0;
  isFired_HLT_Mu40_eta2p1_v9  = false;
  isFired_HLT_Mu40_eta2p1_v10 = false;
  isFired_HLT_Mu40_eta2p1_v11 = false;
  isFired_HLT_Ele80_CaloIdVT_TrkIdT_v8  = false;
  isFired_HLT_Ele80_CaloIdVT_TrkIdT_v9  = false;
  isFired_HLT_Ele80_CaloIdVT_TrkIdT_v10 = false;
  isFired_HLT_Ele80_CaloIdVT_GsfTrkIdT_v1 = false;
  isFired_HLT_Ele80_CaloIdVT_GsfTrkIdT_v2 = false;
  isFired_HLT_PFHT650_v5 = false;
  isFired_HLT_PFHT650_v6 = false;
  isFired_HLT_PFHT650_v7 = false;
  isFired_HLT_PFHT650_v8 = false;
  isFired_HLT_PFHT650_v9 = false;
  isFired_HLT_PFNoPUHT650_v1 = false;
  isFired_HLT_PFNoPUHT650_v3 = false;
  isFired_HLT_PFNoPUHT650_v4 = false;
  isFired_HLT_PFJet320_v3 = false;
  isFired_HLT_PFJet320_v4 = false;
  isFired_HLT_PFJet320_v5 = false;
  isFired_HLT_PFJet320_v6 = false;
  isFired_HLT_PFJet320_v8 = false;
  isFired_HLT_PFJet320_v9 = false;
    
  /************************************/    
  genParticle_pt	    .clear();
  genParticle_px	    .clear();
  genParticle_py	    .clear();
  genParticle_pz	    .clear();
  genParticle_e 	    .clear();
  genParticle_eta	    .clear();
  genParticle_phi	    .clear();
  genParticle_mass	    .clear();
  genParticle_pdgId	    .clear();
  genParticle_status	    .clear();
  genParticle_mother   	    .clear();
  genParticle_nMoth	    .clear();
  genParticle_nDau	    .clear();
  genParticle_dau	    .clear();
  /************************************/
  lep_type		    .clear();
  lep_charge		    .clear();
  lep_e 		    .clear();
  lep_eta		    .clear();
  lep_mass		    .clear();
  lep_pt		    .clear();
  lep_phi		    .clear();
  lep_isHEEP	            .clear();
  lep_isHighPtMuon	    .clear();
  lep_pfRhoCorrRelIso03     .clear();
  lep_pfRhoCorrRelIso04     .clear();
  lep_pfDeltaCorrRelIso     .clear();
  lep_pfRelIso  	    .clear();
  lep_photonIso 	    .clear();
  lep_neutralHadIso	    .clear();
  lep_chargedHadIso	    .clear();
  lep_trackIso	            .clear();  	  
  /************************************/
  jetAK5_pt		    .clear();
  jetAK5_eta		    .clear();
  jetAK5_mass		    .clear();
  jetAK5_phi		    .clear();
  jetAK5_e		    .clear();
  jetAK5_jec                .clear();
  jetAK5_jecUp              .clear();
  jetAK5_jecDown            .clear();
  jetAK5_IDLoose            .clear();
  jetAK5_cm                 .clear();
  jetAK5_nm                 .clear();
  jetAK5_che                .clear();
  jetAK5_ne                 .clear();
  jetAK5_charge 	    .clear();
  jetAK5_flavour	    .clear();
  jetAK5_ssv		    .clear();
  jetAK5_csv		    .clear();
  jetAK5_tchp		    .clear();
  jetAK5_tche		    .clear();
  jetAK5_jp		    .clear();
  jetAK5_jbp		    .clear();
  jetAK5_nSVs		    .clear();
  /************************************/
  jetCA8_pt		    .clear();
  jetCA8_eta		    .clear();
  jetCA8_mass		    .clear();
  jetCA8_phi		    .clear();
  jetCA8_e		    .clear();
  jetCA8_jec                .clear();
  jetCA8_jecUp              .clear();
  jetCA8_jecDown            .clear();  
  jetCA8_IDLoose            .clear();
  jetCA8_cm                 .clear();
  jetCA8_nm                 .clear();
  jetCA8_che                .clear();
  jetCA8_ne                 .clear();
  jetCA8_charge 	    .clear();
  jetCA8_flavour	    .clear();
  jetCA8_ssv		    .clear();
  jetCA8_csv		    .clear();
  jetCA8_tchp		    .clear();
  jetCA8_tche		    .clear();
  jetCA8_jp		    .clear();
  jetCA8_jbp		    .clear();
  jetCA8_nSVs		    .clear(); 
  jetCA8_tau1		    .clear();
  jetCA8_tau2		    .clear();
  jetCA8_tau3		    .clear();
  /************************************/
  jetCA8pruned_pt	    .clear();
  jetCA8pruned_eta	    .clear();
  jetCA8pruned_mass	    .clear();
  jetCA8pruned_phi	    .clear();
  jetCA8pruned_e	    .clear();
  jetCA8pruned_jec          .clear();
  jetCA8pruned_jecUp        .clear();
  jetCA8pruned_jecDown      .clear();
  jetCA8pruned_flavour      .clear();
  jetCA8pruned_charge	    .clear();
  jetCA8pruned_ssv	    .clear();
  jetCA8pruned_csv	    .clear();
  jetCA8pruned_tchp	    .clear();
  jetCA8pruned_tche	    .clear();
  jetCA8pruned_jp	    .clear();
  jetCA8pruned_jbp	    .clear();
  jetCA8pruned_nSVs	    .clear(); 
  /************************************/
  nsubjets		    .clear();
  subjetCA8pruned_pt	    .clear();
  subjetCA8pruned_eta	    .clear();
  subjetCA8pruned_mass      .clear();
  subjetCA8pruned_phi	    .clear();
  subjetCA8pruned_e	    .clear();
  subjetCA8pruned_charge    .clear();
  subjetCA8pruned_flavour   .clear();
  subjetCA8pruned_ssv	    .clear();
  subjetCA8pruned_csv	    .clear();
  subjetCA8pruned_tchp      .clear();
  subjetCA8pruned_tche      .clear();
  subjetCA8pruned_jp	    .clear();
  subjetCA8pruned_jbp	    .clear();
  /************************************/
  METraw_et		    .clear();
  METraw_phi		    .clear();
  MET_et	            .clear();
  MET_phi	            .clear();
  /************************************/  
  nPuVtxTrue                .clear();
  nPuVtx       		    .clear();
  bX	     		    .clear();

} 
