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
  tree_->Branch( "genWeight"	             , &genWeight              ); 
  tree_->Branch( "qScale"	             , &qScale                 ); 
  
  /*-------------------------leptons----------------------------*/
  tree_->Branch( "nlep"  		     , &nlep		       );
  tree_->Branch( "lep_type"		     , &lep_type	       );
  tree_->Branch( "lep_charge"		     , &lep_charge	       );
  tree_->Branch( "lep_e" 		     , &lep_e		       );
  tree_->Branch( "lep_eta"		     , &lep_eta 	       );
  tree_->Branch( "lep_etaSC"		     , &lep_etaSC 	       );
  tree_->Branch( "lep_mass"		     , &lep_mass	       );
  tree_->Branch( "lep_pt"		     , &lep_pt  	       );
  tree_->Branch( "lep_phi"		     , &lep_phi 	       );
  tree_->Branch( "lep_isHighPtMuon"	     , &lep_isHighPtMuon       );
  tree_->Branch( "lep_isTightMuon"	     , &lep_isTightMuon        );
  tree_->Branch( "lep_isLooseMuon"	     , &lep_isLooseMuon        );
  tree_->Branch( "lep_pfRhoCorrRelIso03"     , &lep_pfRhoCorrRelIso03  );
  tree_->Branch( "lep_pfRhoCorrRelIso04"     , &lep_pfRhoCorrRelIso04  );
  tree_->Branch( "lep_pfDeltaCorrRelIso"     , &lep_pfDeltaCorrRelIso  );
  tree_->Branch( "lep_pfRelIso"  	     , &lep_pfRelIso	       );
  tree_->Branch( "lep_photonIso" 	     , &lep_photonIso	       );
  tree_->Branch( "lep_neutralHadIso"	     , &lep_neutralHadIso      );
  tree_->Branch( "lep_chargedHadIso"	     , &lep_chargedHadIso      );
  tree_->Branch( "lep_trackIso"	             , &lep_trackIso	       );  
  tree_->Branch( "lep_passConversionVeto"      , &lep_passConversionVeto       );
  tree_->Branch( "lep_full5x5_sigmaIetaIeta"   , &lep_full5x5_sigmaIetaIeta    );
  tree_->Branch( "lep_dEtaIn"	               , &lep_dEtaIn                   );
  tree_->Branch( "lep_dPhiIn"                  , &lep_dPhiIn                   );
  tree_->Branch( "lep_hOverE"   	       , &lep_hOverE                   );
  tree_->Branch( "lep_relIsoWithDBeta"         , &lep_relIsoWithDBeta          );
  tree_->Branch( "lep_ooEmooP"	               , &lep_ooEmooP                  );
  tree_->Branch( "lep_d0"                      , &lep_d0                       );
  tree_->Branch( "lep_dz"		       , &lep_dz                       );
  tree_->Branch( "lep_expectedMissingInnerHits", &lep_expectedMissingInnerHits );
  tree_->Branch( "lep_isVetoElectron"	       , &lep_isVetoElectron           );
  tree_->Branch( "lep_isMediumElectron"	       , &lep_isMediumElectron         );
  tree_->Branch( "lep_isTightElectron"         , &lep_isTightElectron          );  
  tree_->Branch( "lep_isHeepElectron"	       , &lep_isHeepElectron           );
  tree_->Branch( "lep_isLooseElectron"	       , &lep_isLooseElectron          );
  tree_->Branch( "lep_isSoftMuon"	       , &lep_isSoftMuon	     );
  tree_->Branch( "lep_pfRhoCorrRelIso03Boost"  , &lep_pfRhoCorrRelIso03Boost );
  tree_->Branch( "lep_pfRhoCorrRelIso04Boost"  , &lep_pfRhoCorrRelIso04Boost );
  tree_->Branch( "lep_pfDeltaCorrRelIsoBoost"  , &lep_pfDeltaCorrRelIsoBoost );
  tree_->Branch( "lep_pfRelIsoBoost"  	       , &lep_pfRelIsoBoost	     );
  tree_->Branch( "lep_photonIsoBoost" 	       , &lep_photonIsoBoost	     );
  tree_->Branch( "lep_neutralHadIsoBoost"      , &lep_neutralHadIsoBoost     );
  tree_->Branch( "lep_chargedHadIsoBoost"      , &lep_chargedHadIsoBoost     );
  tree_->Branch( "lep_normChi2"  	       , &lep_normChi2  	     );
  tree_->Branch( "lep_isGlobalMuon"	       , &lep_isGlobalMuon	     );
  tree_->Branch( "lep_trackerHits"	       , &lep_trackerHits	     );
  tree_->Branch( "lep_matchedStations"	       , &lep_matchedStations	     );
  tree_->Branch( "lep_pixelHits" 	       , &lep_pixelHits 	     );
  tree_->Branch( "lep_globalHits"	       , &lep_globalHits	     );
  tree_->Branch( "lep_TauType"		       , &lep_TauType		     );
  tree_->Branch( "lep_SemileptonicPFIso"       , &lep_SemileptonicPFIso      );
  tree_->Branch( "lep_SemileptonicCorrPFIso"   , &lep_SemileptonicCorrPFIso  );
  
 /*-------------------------Tau Discriminant-------------------*/	
  tree_->Branch( "decayModeFindingNewDMs"		      , &decayModeFindingNewDMs 		     );
  tree_->Branch( "decayModeFinding"			      , &decayModeFinding			     ); 
  tree_->Branch( "byLooseCombinedIsolationDeltaBetaCorr3Hits" , &byLooseCombinedIsolationDeltaBetaCorr3Hits  );
  tree_->Branch( "byMediumCombinedIsolationDeltaBetaCorr3Hits", &byMediumCombinedIsolationDeltaBetaCorr3Hits );
  tree_->Branch( "byTightCombinedIsolationDeltaBetaCorr3Hits" , &byTightCombinedIsolationDeltaBetaCorr3Hits  );
  tree_->Branch( "byCombinedIsolationDeltaBetaCorrRaw3Hits"   , &byCombinedIsolationDeltaBetaCorrRaw3Hits    );
  tree_->Branch( "chargedIsoPtSum"			      , &chargedIsoPtSum			     );
  tree_->Branch( "neutralIsoPtSum"			      , &neutralIsoPtSum			     );
  tree_->Branch( "puCorrPtSum"				      , &puCorrPtSum				     );
  tree_->Branch( "byIsolationMVA3oldDMwoLTraw"		      , &byIsolationMVA3oldDMwoLTraw		     );
  tree_->Branch( "byVLooseIsolationMVA3oldDMwoLT"	      , &byVLooseIsolationMVA3oldDMwoLT 	     );
  tree_->Branch( "byLooseIsolationMVA3oldDMwoLT" 	      , &byLooseIsolationMVA3oldDMwoLT  	     );
  tree_->Branch( "byMediumIsolationMVA3oldDMwoLT"	      , &byMediumIsolationMVA3oldDMwoLT 	     );
  tree_->Branch( "byTightIsolationMVA3oldDMwoLT" 	      , &byTightIsolationMVA3oldDMwoLT  	     );
  tree_->Branch( "byVTightIsolationMVA3oldDMwoLT"	      , &byVTightIsolationMVA3oldDMwoLT 	     );
  tree_->Branch( "byVVTightIsolationMVA3oldDMwoLT"	      , &byVVTightIsolationMVA3oldDMwoLT	     );
  tree_->Branch( "byIsolationMVA3oldDMwLTraw"		      , &byIsolationMVA3oldDMwLTraw		     );
  tree_->Branch( "byVLooseIsolationMVA3oldDMwLT" 	      , &byVLooseIsolationMVA3oldDMwLT  	     );
  tree_->Branch( "byLooseIsolationMVA3oldDMwLT"  	      , &byLooseIsolationMVA3oldDMwLT		     );
  tree_->Branch( "byMediumIsolationMVA3oldDMwLT" 	      , &byMediumIsolationMVA3oldDMwLT  	     );
  tree_->Branch( "byTightIsolationMVA3oldDMwLT"  	      , &byTightIsolationMVA3oldDMwLT		     );
  tree_->Branch( "byVTightIsolationMVA3oldDMwLT" 	      , &byVTightIsolationMVA3oldDMwLT  	     );
  tree_->Branch( "byVVTightIsolationMVA3oldDMwLT"	      , &byVVTightIsolationMVA3oldDMwLT 	     );
  tree_->Branch( "byIsolationMVA3newDMwoLTraw"		      , &byIsolationMVA3newDMwoLTraw		     );
  tree_->Branch( "byVLooseIsolationMVA3newDMwoLT"	      , &byVLooseIsolationMVA3newDMwoLT 	     );
  tree_->Branch( "byLooseIsolationMVA3newDMwoLT" 	      , &byLooseIsolationMVA3newDMwoLT  	     );
  tree_->Branch( "byMediumIsolationMVA3newDMwoLT"	      , &byMediumIsolationMVA3newDMwoLT 	     );
  tree_->Branch( "byTightIsolationMVA3newDMwoLT" 	      , &byTightIsolationMVA3newDMwoLT  	     );
  tree_->Branch( "byVTightIsolationMVA3newDMwoLT"	      , &byVTightIsolationMVA3newDMwoLT 	     );
  tree_->Branch( "byVVTightIsolationMVA3newDMwoLT"	      , &byVVTightIsolationMVA3newDMwoLT	     );
  tree_->Branch( "byIsolationMVA3newDMwLTraw"		      , &byIsolationMVA3newDMwLTraw		     );
  tree_->Branch( "byVLooseIsolationMVA3newDMwLT" 	      , &byVLooseIsolationMVA3newDMwLT  	     );
  tree_->Branch( "byLooseIsolationMVA3newDMwLT"  	      , &byLooseIsolationMVA3newDMwLT		     );
  tree_->Branch( "byMediumIsolationMVA3newDMwLT" 	      , &byMediumIsolationMVA3newDMwLT  	     );
  tree_->Branch( "byTightIsolationMVA3newDMwLT"  	      , &byTightIsolationMVA3newDMwLT		     );
  tree_->Branch( "byVTightIsolationMVA3newDMwLT" 	      , &byVTightIsolationMVA3newDMwLT  	     );
  tree_->Branch( "byVVTightIsolationMVA3newDMwLT"	      , &byVVTightIsolationMVA3newDMwLT 	     );
  tree_->Branch( "againstElectronLoose"  		      , &againstElectronLoose  		             );
  tree_->Branch( "againstElectronMedium" 		      , &againstElectronMedium  		     );
  tree_->Branch( "againstElectronTight"  		      , &againstElectronTight			     );
  tree_->Branch( "againstElectronMVA5raw"		      , &againstElectronMVA5raw 		     );
  tree_->Branch( "againstElectronMVA5category"		      , &againstElectronMVA5category		     );
  tree_->Branch( "againstElectronVLooseMVA5"		      , &againstElectronVLooseMVA5		     );
  tree_->Branch( "againstElectronLooseMVA5"		      , &againstElectronLooseMVA5		     );
  tree_->Branch( "againstElectronMediumMVA5"		      , &againstElectronMediumMVA5		     );
  tree_->Branch( "againstElectronTightMVA5"		      , &againstElectronTightMVA5		     );
  tree_->Branch( "againstElectronVTightMVA5"		      , &againstElectronVTightMVA5		     );
  tree_->Branch( "againstMuonLoose"			      , &againstMuonLoose			     );
  tree_->Branch( "againstMuonMedium"			      , &againstMuonTight			     );
  tree_->Branch( "againstMuonTight"			      , &againstMuonTight			     );
  tree_->Branch( "againstMuonLoose2"			      , &againstMuonLoose2			     );
  tree_->Branch( "againstMuonMedium2"			      , &againstMuonMedium2			     );
  tree_->Branch( "againstMuonTight2"			      , &againstMuonLoose3			     );
  tree_->Branch( "againstMuonLoose3"			      , &againstMuonLoose3			     );
  tree_->Branch( "againstMuonTight3"			      , &againstMuonTight3			     );
  tree_->Branch( "againstMuonMVAraw"			      , &againstMuonMVAraw			     );
  tree_->Branch( "againstMuonLooseMVA"			      , &againstMuonLooseMVA			     );
  tree_->Branch( "againstMuonMediumMVA"  		      , &againstMuonMediumMVA			     ); 
  tree_->Branch( "againstMuonTightMVA"			      , &againstMuonTightMVA			     );
      
  /*-------------------------energy density---------------------------*/	 
  tree_->Branch( "rho", &rho );

  /*-------------------------AK4 jets---------------------------*/	 
  tree_->Branch( "njetsAK4"		    , &njetsAK4 	);
  tree_->Branch( "jetAK4_pt"		    , &jetAK4_pt	);
  tree_->Branch( "jetAK4_eta"		    , &jetAK4_eta	);
  tree_->Branch( "jetAK4_mass"		    , &jetAK4_mass	);
  tree_->Branch( "jetAK4_phi"		    , &jetAK4_phi	);
  tree_->Branch( "jetAK4_e"		    , &jetAK4_e 	);
  tree_->Branch( "jetAK4_jec"		    , &jetAK4_jec 	);
  tree_->Branch( "jetAK4_IDLoose"	    , &jetAK4_IDLoose	);
  tree_->Branch( "jetAK4_muf" 		    , &jetAK4_muf	);
  tree_->Branch( "jetAK4_phf" 		    , &jetAK4_phf	);
  tree_->Branch( "jetAK4_emf" 		    , &jetAK4_emf	);
  tree_->Branch( "jetAK4_nhf" 		    , &jetAK4_nhf	);
  tree_->Branch( "jetAK4_chf" 		    , &jetAK4_chf	);
  tree_->Branch( "jetAK4_area" 		    , &jetAK4_area	);
  tree_->Branch( "jetAK4_cm"     	    , &jetAK4_cm	);
  tree_->Branch( "jetAK4_nm"     	    , &jetAK4_nm	);
  tree_->Branch( "jetAK4_che"     	    , &jetAK4_che	);
  tree_->Branch( "jetAK4_ne"     	    , &jetAK4_ne	);
  tree_->Branch( "jetAK4_charge" 	    , &jetAK4_charge    );
  tree_->Branch( "jetAK4_flavour"	    , &jetAK4_flavour   );
  tree_->Branch( "jetAK4_ssv"		    , &jetAK4_ssv	);
  tree_->Branch( "jetAK4_cisv"		    , &jetAK4_cisv	);
  tree_->Branch( "jetAK4_tchp"		    , &jetAK4_tchp	);
  tree_->Branch( "jetAK4_tche"		    , &jetAK4_tche	);
  tree_->Branch( "jetAK4_jp"		    , &jetAK4_jp	);
  tree_->Branch( "jetAK4_jbp"		    , &jetAK4_jbp	);
  tree_->Branch( "jetAK4_vtxMass"	    , &jetAK4_vtxMass   );
  tree_->Branch( "jetAK4_vtxNtracks"        , &jetAK4_vtxNtracks);
  tree_->Branch( "jetAK4_vtx3DVal"	    , &jetAK4_vtx3DVal  );
  tree_->Branch( "jetAK4_vtx3DSig"	    , &jetAK4_vtx3DSig  );
  //tree_->Branch( "jetAK4_nSVs"	    , &jetAK4_nSVs	);
        
  /*-------------------------AK8 jets---------------------------*/    
  tree_->Branch( "njetsAK8"		    , &njetsAK8 		 );
  tree_->Branch( "jetAK8_pt"		    , &jetAK8_pt		 );
  tree_->Branch( "jetAK8_eta"		    , &jetAK8_eta		 );
  tree_->Branch( "jetAK8_mass"		    , &jetAK8_mass		 );
  tree_->Branch( "jetAK8_phi"		    , &jetAK8_phi		 );
  tree_->Branch( "jetAK8_e"		    , &jetAK8_e 		 );
  tree_->Branch( "jetAK8_jec"		    , &jetAK8_jec 		 );
  tree_->Branch( "jetAK8_IDLoose"	    , &jetAK8_IDLoose		 );
  tree_->Branch( "jetAK8_muf" 		    , &jetAK8_muf		 );
  tree_->Branch( "jetAK8_phf" 		    , &jetAK8_phf		 );
  tree_->Branch( "jetAK8_emf" 		    , &jetAK8_emf		 );
  tree_->Branch( "jetAK8_nhf" 		    , &jetAK8_nhf		 );
  tree_->Branch( "jetAK8_chf" 		    , &jetAK8_chf		 );
  tree_->Branch( "jetAK8_area" 		    , &jetAK8_area		 );
  tree_->Branch( "jetAK8_cm"     	    , &jetAK8_cm		 );
  tree_->Branch( "jetAK8_nm"     	    , &jetAK8_nm		 );
  tree_->Branch( "jetAK8_che"     	    , &jetAK8_che		 );
  tree_->Branch( "jetAK8_ne"     	    , &jetAK8_ne		 );
  tree_->Branch( "jetAK8_charge" 	    , &jetAK8_charge    	 );
  tree_->Branch( "jetAK8_flavour"	    , &jetAK8_flavour   	 );
  tree_->Branch( "jetAK8_Hbbtag"	    , &jetAK8_Hbbtag		 );
  tree_->Branch( "jetAK8_ssv"		    , &jetAK8_ssv		 );
  tree_->Branch( "jetAK8_csv"		    , &jetAK8_csv		 );
  tree_->Branch( "jetAK8_tchp"		    , &jetAK8_tchp		 );
  tree_->Branch( "jetAK8_tche"		    , &jetAK8_tche		 );
  tree_->Branch( "jetAK8_jp"		    , &jetAK8_jp		 );
  tree_->Branch( "jetAK8_jbp"		    , &jetAK8_jbp		 );
  tree_->Branch( "jetAK8_tau1"		    , &jetAK8_tau1		 );
  tree_->Branch( "jetAK8_tau2"		    , &jetAK8_tau2      	 );
  tree_->Branch( "jetAK8_tau3"		    , &jetAK8_tau3	    	 );
  tree_->Branch( "jetAK8_prunedmass"        , &jetAK8_prunedmass   );
  tree_->Branch( "jetAK8_prunedmassCorr"    , &jetAK8_prunedmassCorr	 );
  tree_->Branch( "jetAK8pruned_jec"	    , &jetAK8pruned_jec          );
  tree_->Branch( "jetAK8_softdropmass"      , &jetAK8_softdropmass ); 
  tree_->Branch( "jetAK8_softdropmassCorr"  , &jetAK8_softdropmassCorr   );
  tree_->Branch( "jetAK8softdrop_jec"	    , &jetAK8softdrop_jec        );
  
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
  //tree_->Branch( "njetsAK8softdrop"	   , &njetsAK8softdrop	      );
  //tree_->Branch( "jetAK8softdrop_pt"     , &jetAK8softdrop_pt       );
  //tree_->Branch( "jetAK8softdrop_eta"    , &jetAK8softdrop_eta      );
  //tree_->Branch( "jetAK8softdrop_mass"   , &jetAK8softdrop_mass     );
  //tree_->Branch( "jetAK8softdrop_phi"    , &jetAK8softdrop_phi      );
  //tree_->Branch( "jetAK8softdrop_e"	   , &jetAK8softdrop_e        );
  //tree_->Branch( "jetAK8softdrop_flavour", &jetAK8softdrop_flavour  );
  //tree_->Branch( "jetAK8softdrop_charge" , &jetAK8softdrop_charge   );
  //tree_->Branch( "jetAK8softdrop_ssv"    , &jetAK8softdrop_ssv      );
  //tree_->Branch( "jetAK8softdrop_csv"    , &jetAK8softdrop_csv      );
  //tree_->Branch( "jetAK8softdrop_tchp"   , &jetAK8softdrop_tchp     );
  //tree_->Branch( "jetAK8softdrop_tche"   , &jetAK8softdrop_tche     );
  //tree_->Branch( "jetAK8softdrop_jp"     , &jetAK8softdrop_jp       );
  //tree_->Branch( "jetAK8softdrop_jbp"    , &jetAK8softdrop_jbp      );
  //tree_->Branch( "jetAK8softdrop_nSVs"   , &jetAK8softdrop_nSVs     );
	  
  /*----------------------Pruned AK8 subjets---------------------------*/    
  tree_->Branch( "nprunedsubjets"	    , &nprunedsubjets	       );
  tree_->Branch( "subjetAK8pruned_pt"	    , &subjetAK8pruned_pt      );
  tree_->Branch( "subjetAK8pruned_eta"	    , &subjetAK8pruned_eta     );
  tree_->Branch( "subjetAK8pruned_mass"     , &subjetAK8pruned_mass    );
  tree_->Branch( "subjetAK8pruned_phi"	    , &subjetAK8pruned_phi     );
  tree_->Branch( "subjetAK8pruned_e"	    , &subjetAK8pruned_e       );
  tree_->Branch( "subjetAK8pruned_charge"   , &subjetAK8pruned_charge  );
  tree_->Branch( "subjetAK8pruned_flavour"  , &subjetAK8pruned_flavour );
  tree_->Branch( "subjetAK8pruned_ssv"	    , &subjetAK8pruned_ssv     );
  tree_->Branch( "subjetAK8pruned_csv"	    , &subjetAK8pruned_csv     );
  tree_->Branch( "subjetAK8pruned_tchp"     , &subjetAK8pruned_tchp    );
  tree_->Branch( "subjetAK8pruned_tche"     , &subjetAK8pruned_tche    );
  tree_->Branch( "subjetAK8pruned_jp"	    , &subjetAK8pruned_jp      );
  tree_->Branch( "subjetAK8pruned_jbp"	    , &subjetAK8pruned_jbp     );
  
  // /*----------------------Softdrop AK8 subjets---------------------------*/
  tree_->Branch( "nsoftdropsubjets"         , &nsoftdropsubjets           );
  tree_->Branch( "subjetAK8softdrop_pt"     , &subjetAK8softdrop_pt      );
  tree_->Branch( "subjetAK8softdrop_eta"    , &subjetAK8softdrop_eta     );
  tree_->Branch( "subjetAK8softdrop_mass"   , &subjetAK8softdrop_mass    );
  tree_->Branch( "subjetAK8softdrop_phi"    , &subjetAK8softdrop_phi     );
  tree_->Branch( "subjetAK8softdrop_e"      , &subjetAK8softdrop_e       );
  tree_->Branch( "subjetAK8softdrop_charge" , &subjetAK8softdrop_charge  );
  tree_->Branch( "subjetAK8softdrop_flavour", &subjetAK8softdrop_flavour );
  tree_->Branch( "subjetAK8softdrop_ssv"    , &subjetAK8softdrop_ssv     );
  tree_->Branch( "subjetAK8softdrop_csv"    , &subjetAK8softdrop_csv     );
  tree_->Branch( "subjetAK8softdrop_tchp"   , &subjetAK8softdrop_tchp    );
  tree_->Branch( "subjetAK8softdrop_tche"   , &subjetAK8softdrop_tche    );
  tree_->Branch( "subjetAK8softdrop_jp"     , &subjetAK8softdrop_jp      );
  tree_->Branch( "subjetAK8softdrop_jbp"    , &subjetAK8softdrop_jbp     );
    
  /*-------------------------AK4 GenJets---------------------------*/	 
  tree_->Branch( "ngenJetsAK4"		    , &ngenJetsAK4 	  );
  tree_->Branch( "genJetAK4_pt"		    , &genJetAK4_pt	  );
  tree_->Branch( "genJetAK4_eta"	    , &genJetAK4_eta	  );
  tree_->Branch( "genJetAK4_mass"	    , &genJetAK4_mass	  );
  tree_->Branch( "genJetAK4_phi"	    , &genJetAK4_phi	  );
  tree_->Branch( "genJetAK4_e"		    , &genJetAK4_e 	  );
  tree_->Branch( "genJetNoNuAK4_pt"	    , &genJetNoNuAK4_pt	  );
  tree_->Branch( "genJetNoNuAK4_mass"	    , &genJetNoNuAK4_mass );
  tree_->Branch( "genJetNoNuAK4_e"	    , &genJetNoNuAK4_e    );
        
  /*----------------------HLT trigger---------------------------*/	  
  tree_->Branch("isFired_HLT_AK8PFJet360_TrimMass30_v1"			 , &isFired_HLT_AK8PFJet360_TrimMass30_v1		   );
  tree_->Branch("isFired_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v1"         , &isFired_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v1	   );
  tree_->Branch("isFired_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p41_v1", &isFired_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p41_v1);
  tree_->Branch("isFired_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v1"		 , &isFired_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v1	   );
  tree_->Branch("isFired_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v1"		 , &isFired_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v1	   );
  tree_->Branch("isFired_HLT_PFHT900_v1"				 , &isFired_HLT_PFHT900_v1				   );
  tree_->Branch("isFired_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v1"            , &isFired_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v1 	   );
  tree_->Branch("isFired_HLT_Ele32_eta2p1_WP75_Gsf_v1"                   , &isFired_HLT_Ele32_eta2p1_WP75_Gsf_v1		   );
  tree_->Branch("isFired_HLT_Ele105_CaloIdVT_GsfTrkIdT_v1"		 , &isFired_HLT_Ele105_CaloIdVT_GsfTrkIdT_v1		   );
  tree_->Branch("isFired_HLT_IsoMu24_eta2p1_v1"                          , &isFired_HLT_IsoMu24_eta2p1_v1			   );
  tree_->Branch("isFired_HLT_Mu45_eta2p1_v1"				 , &isFired_HLT_Mu45_eta2p1_v1  			   );
  
  tree_->Branch("triggerObject_pt"		, &triggerObject_pt		);
  tree_->Branch("triggerObject_eta"		, &triggerObject_eta		);
  tree_->Branch("triggerObject_phi"		, &triggerObject_phi	        );
  tree_->Branch("triggerObject_mass"		, &triggerObject_mass		);
  tree_->Branch("triggerObject_filterIDs"	, &triggerObject_filterIDs	);
  tree_->Branch("triggerObject_firedTrigger"	, &triggerObject_firedTrigger	);

  /*-------------------------MET--------------------------------*/
  tree_->Branch("METraw_et"		        , &METraw_et	     );
  tree_->Branch("METraw_phi"		        , &METraw_phi	     ); 
  tree_->Branch("METraw_sumEt"		        , &METraw_sumEt	     );   
  tree_->Branch("MET_corrPx"		        , &MET_corrPx	     ); 
  tree_->Branch("MET_corrPy"		        , &MET_corrPy	     );   
  tree_->Branch("MET_et"	                , &MET_et  	     ); 
  tree_->Branch("MET_phi"	                , &MET_phi           );
  tree_->Branch("MET_sumEt"	                , &MET_sumEt 	     ); 
  //tree_->Branch("METdefault_et"	        , &METdefault_et     ); 
  //tree_->Branch("METdefault_sumEt"	        , &METdefault_sumEt  ); 
  //tree_->Branch("METdefault_phi"	        , &METdefault_phi    );
  //tree_->Branch("MET_T1Uncertainty"	        , &MET_T1Uncertainty );

  /*------------- ------EVENT infos-----------------------------*/
  tree_->Branch("EVENT_event"	 , &EVENT_event     );
  tree_->Branch("EVENT_run"	 , &EVENT_run	    );
  tree_->Branch("EVENT_lumiBlock", &EVENT_lumiBlock );
  
  /*--------------------------PU infos--------------------------*/  			         
  tree_->Branch("nPuVtxTrue", &nPuVtxTrue );	
  tree_->Branch("nPuVtx"    , &nPuVtx	  );
  tree_->Branch("bX"	    , &bX	  );
  
  /*--------------------------PV infos--------------------------*/
  tree_->Branch("nPVs", &nPVs);
  
}

//=================================================================================================================== 
void NtupleBranches::reset( void ){

  lheV_pt     = 0;
  lheHT       = 0;
  lheNj       = 0;
  genWeight   = 0;
  qScale      = 0;
  nlep        = 0;
  njetsAK4    = 0;
  njetsAK8    = 0;  
  ngenJetsAK4 = 0;

  //njetsAK8pruned = 0;
  //njetsAK8softdrop = 0;
  isFired_HLT_AK8PFJet360_TrimMass30_v1			  = false;
  isFired_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v1	  = false;
  isFired_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p41_v1 = false;
  isFired_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v1	 	  = false;
  isFired_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v1	 	  = false;
  isFired_HLT_PFHT900_v1			 	  = false;
  isFired_HLT_IsoMu24_eta2p1_v1			 	  = false;
  isFired_HLT_Mu45_eta2p1_v1			 	  = false;
  isFired_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v1    	  = false;
  isFired_HLT_Ele32_eta2p1_WP75_Gsf_v1           	  = false;
  isFired_HLT_Ele105_CaloIdVT_GsfTrkIdT_v1	   	  = false;
    
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
  genParticle_mother        .clear();
  genParticle_nMoth	    .clear();
  genParticle_nDau	    .clear();
  genParticle_dau	    .clear();
  /************************************/
  genJetAK4_pt              .clear();
  genJetAK4_eta             .clear();
  genJetAK4_mass            .clear();
  genJetAK4_phi             .clear();
  genJetAK4_e	            .clear();
  genJetNoNuAK4_pt          .clear();
  genJetNoNuAK4_mass        .clear();
  genJetNoNuAK4_e           .clear();
  /************************************/
  lep_type		      .clear();
  lep_charge		      .clear();
  lep_e 		      .clear();
  lep_eta		      .clear();
  lep_etaSC		      .clear();
  lep_mass		      .clear();
  lep_pt		      .clear();
  lep_phi		      .clear();
  lep_isHighPtMuon	      .clear();
  lep_isTightMuon	      .clear();
  lep_isLooseMuon	      .clear();
  lep_pfRhoCorrRelIso03	      .clear();
  lep_pfRhoCorrRelIso04       .clear();
  lep_pfDeltaCorrRelIso       .clear();
  lep_pfRelIso  	      .clear();
  lep_photonIso 	      .clear();
  lep_neutralHadIso	      .clear();
  lep_chargedHadIso	      .clear();
  lep_trackIso	              .clear(); 
  lep_passConversionVeto      .clear();
  lep_full5x5_sigmaIetaIeta   .clear();
  lep_dEtaIn		      .clear();
  lep_dPhiIn		      .clear();
  lep_hOverE		      .clear();
  lep_relIsoWithDBeta	      .clear();
  lep_ooEmooP		      .clear();
  lep_d0		      .clear();
  lep_dz		      .clear();
  lep_expectedMissingInnerHits.clear();
  lep_isVetoElectron	      .clear();
  lep_isMediumElectron	      .clear();
  lep_isTightElectron	      .clear();
  lep_isLooseElectron	      .clear();
  lep_isHeepElectron	      .clear();
  lep_TauType               .clear();
  lep_isSoftMuon	    .clear();
  lep_pfRhoCorrRelIso03Boost.clear();
  lep_pfRhoCorrRelIso04Boost.clear();
  lep_pfDeltaCorrRelIsoBoost.clear();
  lep_pfRelIsoBoost         .clear();    
  lep_photonIsoBoost        .clear();
  lep_neutralHadIsoBoost    .clear();
  lep_chargedHadIsoBoost    .clear();   
  lep_SemileptonicPFIso     .clear();
  lep_SemileptonicCorrPFIso .clear();
  lep_normChi2  	    .clear();
  lep_isGlobalMuon	    .clear();
  lep_trackerHits	    .clear();
  lep_matchedStations	    .clear();
  lep_pixelHits 	    .clear();
  lep_globalHits	    .clear();	
  /************************************/ 
  /*Tau discriminats*/
  decayModeFindingNewDMs		     .clear();      
  decayModeFinding			     .clear();
  byLooseCombinedIsolationDeltaBetaCorr3Hits .clear();
  byMediumCombinedIsolationDeltaBetaCorr3Hits.clear();
  byTightCombinedIsolationDeltaBetaCorr3Hits .clear();
  byCombinedIsolationDeltaBetaCorrRaw3Hits   .clear();
  chargedIsoPtSum			     .clear();
  neutralIsoPtSum			     .clear();
  puCorrPtSum				     .clear();
  byIsolationMVA3oldDMwoLTraw		     .clear();
  byVLooseIsolationMVA3oldDMwoLT	     .clear();
  byLooseIsolationMVA3oldDMwoLT 	     .clear();
  byMediumIsolationMVA3oldDMwoLT	     .clear();
  byTightIsolationMVA3oldDMwoLT 	     .clear();
  byVTightIsolationMVA3oldDMwoLT	     .clear();
  byVVTightIsolationMVA3oldDMwoLT	     .clear();
  byIsolationMVA3oldDMwLTraw		     .clear();
  byVLooseIsolationMVA3oldDMwLT 	     .clear();
  byLooseIsolationMVA3oldDMwLT  	     .clear();
  byMediumIsolationMVA3oldDMwLT 	     .clear();
  byTightIsolationMVA3oldDMwLT  	     .clear();
  byVTightIsolationMVA3oldDMwLT 	     .clear();
  byVVTightIsolationMVA3oldDMwLT	     .clear();
  byIsolationMVA3newDMwoLTraw		     .clear();
  byVLooseIsolationMVA3newDMwoLT	     .clear();
  byLooseIsolationMVA3newDMwoLT 	     .clear();
  byMediumIsolationMVA3newDMwoLT	     .clear();
  byTightIsolationMVA3newDMwoLT 	     .clear();
  byVTightIsolationMVA3newDMwoLT	     .clear();
  byVVTightIsolationMVA3newDMwoLT	     .clear();
  byIsolationMVA3newDMwLTraw		     .clear();
  byVLooseIsolationMVA3newDMwLT 	     .clear();
  byLooseIsolationMVA3newDMwLT  	     .clear();
  byMediumIsolationMVA3newDMwLT 	     .clear();
  byTightIsolationMVA3newDMwLT  	     .clear();
  byVTightIsolationMVA3newDMwLT 	     .clear();
  byVVTightIsolationMVA3newDMwLT	     .clear();
  againstElectronLoose  		     .clear();
  againstElectronMedium 		     .clear();
  againstElectronTight  		     .clear();
  againstElectronMVA5raw		     .clear();
  againstElectronMVA5category		     .clear();
  againstElectronVLooseMVA5		     .clear();
  againstElectronLooseMVA5		     .clear();
  againstElectronMediumMVA5		     .clear();
  againstElectronTightMVA5		     .clear();
  againstElectronVTightMVA5		     .clear();
  againstMuonLoose			     .clear();
  againstMuonTight			     .clear();
  againstMuonTight			     .clear();
  againstMuonLoose2			     .clear();
  againstMuonMedium2			     .clear();
  againstMuonLoose3			     .clear();
  againstMuonLoose3			     .clear();
  againstMuonTight3			     .clear();
  againstMuonMVAraw			     .clear();
  againstMuonLooseMVA			     .clear();
  againstMuonMediumMVA  		     .clear();
  againstMuonTightMVA			     .clear();

 /************************************/
  jetAK4_pt		    .clear();
  jetAK4_eta		    .clear();
  jetAK4_mass		    .clear();
  jetAK4_phi		    .clear();
  jetAK4_e		    .clear();
  jetAK4_jec                .clear();
  //jetAK4_jecUp              .clear();
  //jetAK4_jecDown            .clear();
  jetAK4_IDLoose            .clear();
  jetAK4_muf     	    .clear();
  jetAK4_phf     	    .clear();
  jetAK4_emf     	    .clear();
  jetAK4_nhf     	    .clear();
  jetAK4_chf     	    .clear();
  jetAK4_area        	    .clear();
  jetAK4_cm                 .clear();
  jetAK4_nm                 .clear();
  jetAK4_che                .clear();
  jetAK4_ne                 .clear();
  jetAK4_charge 	    .clear();
  jetAK4_flavour	    .clear();
  jetAK4_ssv		    .clear();
  jetAK4_cisv		    .clear();
  jetAK4_tchp		    .clear();
  jetAK4_tche		    .clear();
  jetAK4_jp		    .clear();
  jetAK4_jbp		    .clear();
  jetAK4_vtxMass	    .clear(); 
  jetAK4_vtxNtracks	    .clear(); 
  jetAK4_vtx3DVal	    .clear(); 
  jetAK4_vtx3DSig	    .clear(); 
  //jetAK4_nSVs		    .clear();
  /************************************/
  jetAK8_pt		    .clear();
  jetAK8_eta		    .clear();
  jetAK8_mass		    .clear();
  jetAK8_phi		    .clear();
  jetAK8_e		    .clear();
  jetAK8_jec                .clear();
  //jetAK8_jecUp 	     .clear();
  //jetAK8_jecDown	     .clear();
  jetAK8_IDLoose            .clear();
  jetAK8_muf     	    .clear();
  jetAK8_phf     	    .clear();
  jetAK8_emf     	    .clear();
  jetAK8_nhf     	    .clear();
  jetAK8_chf     	    .clear();
  jetAK8_area        	    .clear();
  jetAK8_cm                 .clear();
  jetAK8_nm                 .clear();
  jetAK8_che                .clear();
  jetAK8_ne                 .clear();
  jetAK8_charge 	    .clear();
  jetAK8_flavour	    .clear();
  jetAK8_Hbbtag		    .clear();
  jetAK8_ssv		    .clear();
  jetAK8_csv		    .clear();
  jetAK8_tchp		    .clear();
  jetAK8_tche		    .clear();
  jetAK8_jp		    .clear();
  jetAK8_jbp		    .clear();
  jetAK8_tau1		    .clear();
  jetAK8_tau2		    .clear();
  jetAK8_tau3		    .clear();
  jetAK8_prunedmass   .clear();
  jetAK8_softdropmass .clear();
  jetAK8_prunedmassCorr	    .clear();
  jetAK8_softdropmassCorr   .clear();
  jetAK8softdrop_jec	    .clear();
  jetAK8pruned_jec	    .clear();

  /************************************/
  // jetAK8pruned_pt     .clear();
  // jetAK8pruned_eta     .clear();
  // jetAK8pruned_mass     .clear();
  // jetAK8pruned_phi     .clear();
  // jetAK8pruned_e     .clear();
  // jetAK8pruned_flavour   .clear();
  // jetAK8pruned_charge    .clear();
  // jetAK8pruned_ssv     .clear();
  // jetAK8pruned_csv     .clear();
  // jetAK8pruned_tchp     .clear();
  // jetAK8pruned_tche     .clear();
  // jetAK8pruned_jp     .clear();
  // jetAK8pruned_jbp     .clear();
  /************************************/
  //jetAK8softdrop_pt	    .clear();
  //jetAK8softdrop_eta      .clear();
  //jetAK8softdrop_mass     .clear();
  //jetAK8softdrop_phi      .clear();
  //jetAK8softdrop_e	    .clear();
  //jetAK8softdrop_flavour  .clear();
  //jetAK8softdrop_charge   .clear();
  //jetAK8softdrop_ssv      .clear();
  //jetAK8softdrop_csv      .clear();
  //jetAK8softdrop_tchp     .clear();
  //jetAK8softdrop_tche     .clear();
  //jetAK8softdrop_jp	    .clear();
  //jetAK8softdrop_jbp      .clear();
  //jetAK8softdrop_nSVs	    .clear();
  /************************************/
  nprunedsubjets	    .clear();
  subjetAK8pruned_pt	    .clear();
  subjetAK8pruned_eta	    .clear();
  subjetAK8pruned_mass      .clear();
  subjetAK8pruned_phi	    .clear();
  subjetAK8pruned_e	    .clear();
  subjetAK8pruned_charge    .clear();
  subjetAK8pruned_flavour   .clear();
  subjetAK8pruned_ssv	    .clear();
  subjetAK8pruned_csv	    .clear();
  subjetAK8pruned_tchp      .clear();
  subjetAK8pruned_tche      .clear();
  subjetAK8pruned_jp	    .clear();
  subjetAK8pruned_jbp	    .clear();
  /************************************/
  nsoftdropsubjets         .clear();
  subjetAK8softdrop_pt     .clear();
  subjetAK8softdrop_eta    .clear();
  subjetAK8softdrop_mass   .clear();
  subjetAK8softdrop_phi    .clear();
  subjetAK8softdrop_e      .clear();
  subjetAK8softdrop_charge .clear();
  subjetAK8softdrop_flavour.clear();
  subjetAK8softdrop_ssv    .clear();
  subjetAK8softdrop_csv    .clear();
  subjetAK8softdrop_tchp   .clear();
  subjetAK8softdrop_tche   .clear();
  subjetAK8softdrop_jp     .clear();
  subjetAK8softdrop_jbp    .clear();
  /************************************/

  METraw_et		 .clear();
  METraw_phi		 .clear();
  METraw_sumEt		 .clear();
  MET_et	         .clear();
  MET_sumEt              .clear();
  MET_phi	         .clear();
  MET_corrPx	         .clear();
  MET_corrPy	         .clear();
  MET_T1Uncertainty	 .clear();

  /************************************/  
  nPuVtxTrue             .clear();
  nPuVtx       		 .clear();
  bX	     		 .clear();
  
  triggerObject_pt.clear();
  triggerObject_eta.clear();
  triggerObject_phi.clear();
  triggerObject_mass.clear();
  triggerObject_filterIDs.clear();
  triggerObject_firedTrigger.clear();

} 
