#ifndef NtupleBranches_H
#define NtupleBranches_H

#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
#include <DataFormats/MuonReco/interface/Muon.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/JetReco/interface/Jet.h>
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Framework/interface/Event.h"
#include "TMatrixD.h"
#include "TMatrix.h"

#include "TTree.h"
#include "TLorentzVector.h"

/*here we declare the input and output variables*/

class NtupleBranches {

public:
  NtupleBranches( std::map< std::string, bool >& runFlags, TTree* tree = 0);
  ~NtupleBranches( void );
   
  void branch( std::map< std::string, bool >& runFlags );
  void getEventByLabels( edm::EventBase const & event );
  void reset( void );
  void fillTree( void ){ tree_->Fill(); };
  void LabelHistograms (std::map< std::string, bool >& runFlags );
  //=================================================================================================================== 

  /* output histogram */
  /* Cutflow */
  TH1F* cutflow ;
  TH1F* bweight ;
  TH1F* nmuon ;
  TH1F* q2_nocut;
  TH1F* hammer_width;
  TH1F* hammer_width_lattice;
  
  /* Histogram for genParticles */
   TH1F* genParticle_Bdau_X_id;
   TH1F* genParticle_Bdau_X_pt; 
   TH1F* genParticle_Bdau_X_eta; 
   TH1F* genParticle_Bdau_X_phi; 
   TH1F* genParticle_Bdau_X_mass; 
   TH1F* genParticle_Bdau_mu1_pt;
   TH1F* genParticle_Bdau_mu1_eta;
   TH1F* genParticle_Bdau_mu1_phi; 
   TH1F* genParticle_Bdau_mu2_pt;
   TH1F* genParticle_Bdau_mu2_eta; 
   TH1F* genParticle_Bdau_mu2_phi; 
   TH1F* genParticle_Bdau_Jpsi_pt; 
   TH1F* genParticle_Bdau_Jpsi_eta; 
   TH1F* genParticle_Bdau_Jpsi_phi; 
   TH1F* genParticle_Bdau_Jpsi_mass;
   TH1F* genParticle_Bvis_pt; 
   TH1F* genParticle_Bvis_eta;
   TH1F* genParticle_Bvis_phi; 
   TH1F* genParticle_Bvis_mass;

  /* output tree variables*/
    
  /** genParticles */
  int                             genParticle_N;
 
  std::vector<float>              genParticle_pt       ;
//  std::vector<float>              genParticle_px       ;
//  std::vector<float>              genParticle_py       ;
//  std::vector<float>              genParticle_pz       ;
//  std::vector<float>              genParticle_e        ;
  std::vector<float>              genParticle_eta      ;
  std::vector<float>              genParticle_phi      ;
  std::vector<float>              genParticle_mass     ;
  std::vector<int  >              genParticle_pdgId    ;
  std::vector<int  >              genParticle_isPrompt ;
  std::vector<int  >              genParticle_isDirectPromptTauDecayProduct;
  std::vector<int  >              genParticle_fromHardProcessFinalState;
  std::vector<int  >              genParticle_isDirectHardProcessTauDecayProductFinalState;
  std::vector<int  >              genParticle_status   ;
  std::vector<int  >              genParticle_nDau     ;
  std::vector<int  >              genParticle_nMoth    ;
  std::vector<std::vector<int> >  genParticle_mother   ; 
  std::vector<std::vector<float>> genParticle_mother_pt;
  std::vector<std::vector<int> >  genParticle_dau      ;

  //  std::vector<int>  genParticle_pmother   ; 
  std::vector<std::vector<int>>  genParticle_pdgs   ; 
  std::vector<std::vector<int>>  genParticle_layers   ; 
  std::vector<std::vector<float>>  genParticle_ppt   ; 
  std::vector<std::vector<float>>  genParticle_peta   ; 
  std::vector<std::vector<float>>  genParticle_pphi   ; 
  std::vector<std::vector<int>>  genParticle_isfinal   ; 

  /** generator info */
  float                           lheV_pt              ;
  float                           lheHT                ;
  int                             lheNj                ;
  int                             lheNb                ;
  int                             lheNl                ;
  float                           lheV_mass            ;
  float                           genWeight            ;
  float                           genFacWeightUp       ;
  float                           genFacWeightDown     ;
  float                           genRenWeightUp       ;
  float                           genRenWeightDown     ;
  float                           genFacRenWeightUp    ;
  float                           genFacRenWeightDown  ;
  float                           qScale               ;
  float                           PDF_rms              ;
  std::vector<int  >              PDF_id               ;
  std::vector<float>              PDF_x                ;
  std::vector<float>              PDF_xPDF             ;

          
 


  
   
  
  /** HLT trigger decisions */
  std::map<std::string,bool> HLT_isFired;
	 
  /** HLT trigger objects */
  std::vector<float>  		    triggerObject_pt	      ;
  std::vector<float>  		    triggerObject_eta	      ;
  std::vector<float>  		    triggerObject_phi	      ;
  std::vector<float>  		    triggerObject_mass	      ;
  std::vector<std::string>  		    triggerObject_lastname    ;
  std::vector< std::vector<float> > triggerObject_filterIDs   ; // as defined in http://cmslxr.fnal.gov/lxr/source/DataFormats/HLTReco/interface/TriggerTypeDefs.h
  //  std::vector< std::vector<std::string> > triggerObject_filterLabels;
  std::map<std::string, std::vector<std::string> > triggerObject_filterLabels;
  std::vector< std::vector<int> >   triggerObject_firedTrigger; // as defined in plugins/TriggersNtuplizer.cc


  /** HLT filter decisions */
  bool passFilter_HBHE_;
  bool passFilter_HBHELoose_;
  bool passFilter_HBHETight_;
  bool passFilter_HBHEIso_;
  bool passFilter_CSCHalo_;
  bool passFilter_CSCTightHalo2015_;
  bool passFilter_HCALlaser_;
  bool passFilter_ECALDeadCell_;
  bool passFilter_GoodVtx_;
  bool passFilter_TrkFailure_;
  bool passFilter_EEBadSc_;
  bool passFilter_ECALlaser_;
  bool passFilter_TrkPOG_;
  bool passFilter_TrkPOG_manystrip_;
  bool passFilter_TrkPOG_toomanystrip_;
  bool passFilter_TrkPOG_logError_;
  bool passFilter_METFilters_;
  bool passFilter_CSCTightHaloTrkMuUnvetoFilter_   ;
  bool passFilter_globalTightHalo2016_             ;
  bool passFilter_globalSuperTightHalo2016_             ;
  bool passFilter_HcalStripHalo_                   ;
  bool passFilter_chargedHadronTrackResolution_    ;
  bool passFilter_muonBadTrack_                    ;
  bool flag_badMuons_;
  bool flag_duplicateMuons_;
  bool flag_nobadMuons_;
  bool passFilter_ecalBadCalib_;

  /** MET */
  /** energy density */
  float                             rho;
 
  std::vector<float>                METraw_et		      ;	 
  std::vector<float>                METraw_phi		      ;
  std::vector<float>  	            METraw_sumEt	      ;
  std::vector<float>  	      	    MET_corrPx  	      ;
  std::vector<float>  	            MET_corrPy  	      ;
  std::vector<float>  	      	    MET_et		      ;
  std::vector<float>  	      	    MET_phi		      ;
  std::vector<float>  	      	    MET_puppi_et		      ;
  std::vector<float>  	      	    MET_puppi_phi	      ;
  std::vector<float>  	      	    MET_mva_et		      ;
  std::vector<float>  	      	    MET_mva_phi	              ;
  std::vector<float>  	      	    MET_sumEt		      ;
  std::vector<float>  	      	    MET_T1Uncertainty	      ;
  
  std::vector<float>  	      	    MET_JetEnUp	      ;
  std::vector<float>  	      	    MET_JetEnDown	      ;
  std::vector<float>  	      	    MET_JetResUp	      ;
  std::vector<float>  	      	    MET_JetResDown	      ;
  std::vector<float>  	      	    MET_UnclusteredEnUp	      ;
  std::vector<float>  	      	    MET_UnclusteredEnDown	      ;
  

  /** MET SVift*/
  std::vector<float>                MET_significance	      ;	 
  std::vector<float>                MET_cov00	      ;	 
  std::vector<float>                MET_cov10	      ;	 
  std::vector<float>                MET_cov11	      ;	 
  std::vector<float>                MET_mva_cov00	      ;	 
  std::vector<float>                MET_mva_cov10	      ;	 
  std::vector<float>                MET_mva_cov11	      ;	 
  std::vector< std::vector<float> > MET_mva_recoil_pt;
  std::vector< std::vector<float> > MET_mva_recoil_eta;
  std::vector< std::vector<float> > MET_mva_recoil_phi;
  std::vector< std::vector<int> >   MET_mva_recoil_pdgId;
  std::vector<int>  	            MET_Nmva;

  /*------------------------EVENT infos-------------------------*/    
  int                               EVENT_event            ;
  int                               EVENT_run              ;
  int                               EVENT_lumiBlock        ;
  /*-------------------------JPSI VARIABLES---------------------*/
  std::vector<int  >                IsJpsiMu    ;   
  //  std::vector<int  >                IsBsTauTau   ;
  //  std::vector<int  >                IsBsTauTauFH   ;
  //  std::vector<int  >                IsBsDstarTauNu   ;

  std::vector<int  >                JpsiMu_nCandidates ;
  std::vector<float>                JpsiMu_mu1_pt      ;
  std::vector<float>                JpsiMu_mu1_eta     ;
  std::vector<float>                JpsiMu_mu1_phi     ;
  std::vector<float>                JpsiMu_mu1_mass     ;
  std::vector<float>                JpsiMu_mu1_unfit_pt      ;
  std::vector<float>                JpsiMu_mu1_unfit_eta     ;
  std::vector<float>                JpsiMu_mu1_unfit_phi     ;
  std::vector<float>                JpsiMu_mu1_unfit_mass     ;
  std::vector<int  >                JpsiMu_mu1_q      ;   
  std::vector<int  >                JpsiMu_mu1_isLoose   ;
  std::vector<int  >                JpsiMu_mu1_isTight   ;
  std::vector<int  >                JpsiMu_mu1_isPF      ;
  std::vector<int  >                JpsiMu_mu1_isGlobal  ;
  std::vector<int  >                JpsiMu_mu1_isTracker ;
  std::vector<int  >                JpsiMu_mu1_isSoft    ;
  std::vector<float>                JpsiMu_mu1_vx       ;
  std::vector<float>                JpsiMu_mu1_vy       ;
  std::vector<float>                JpsiMu_mu1_vz       ;
  std::vector<float>                JpsiMu_mu1_iso       ;
  std::vector<float>                JpsiMu_mu1_dbiso       ;

  std::vector<float>                JpsiMu_mu2_pt      ;
  std::vector<float>                JpsiMu_mu2_eta     ;
  std::vector<float>                JpsiMu_mu2_phi     ;
  std::vector<float>                JpsiMu_mu2_mass     ;
  std::vector<float>                JpsiMu_mu2_unfit_pt      ;
  std::vector<float>                JpsiMu_mu2_unfit_eta     ;
  std::vector<float>                JpsiMu_mu2_unfit_phi     ;
  std::vector<float>                JpsiMu_mu2_unfit_mass     ;
  std::vector<int  >                JpsiMu_mu2_q   ;   
  std::vector<int  >                JpsiMu_mu2_isLoose   ;
  std::vector<int  >                JpsiMu_mu2_isTight   ;
  std::vector<int  >                JpsiMu_mu2_isPF      ;
  std::vector<int  >                JpsiMu_mu2_isGlobal  ;
  std::vector<int  >                JpsiMu_mu2_isTracker ;
  std::vector<int  >                JpsiMu_mu2_isSoft    ;
  std::vector<float>                JpsiMu_mu2_vx       ;
  std::vector<float>                JpsiMu_mu2_vy       ;
  std::vector<float>                JpsiMu_mu2_vz       ;
  std::vector<float>                JpsiMu_mu2_iso       ;
  std::vector<float>                JpsiMu_mu2_dbiso       ;

  std::vector<float>                JpsiMu_mu3_pt      ;
  std::vector<float>                JpsiMu_mu3_eta     ;
  std::vector<float>                JpsiMu_mu3_phi     ;
  std::vector<float>                JpsiMu_mu3_mass     ;
  std::vector<float>                JpsiMu_mu3_unfit_pt      ;
  std::vector<float>                JpsiMu_mu3_unfit_eta     ;
  std::vector<float>                JpsiMu_mu3_unfit_phi     ;
  std::vector<float>                JpsiMu_mu3_unfit_mass     ;
  std::vector<float>                JpsiMu_mu3_doca2mu1;
  std::vector<float>                JpsiMu_mu3_doca2mu2;
  std::vector<int  >                JpsiMu_mu3_q   ;   
  std::vector<int  >                JpsiMu_mu3_isLoose   ;
  std::vector<int  >                JpsiMu_mu3_isTight   ;
  std::vector<int  >                JpsiMu_mu3_isPF      ;
  std::vector<int  >                JpsiMu_mu3_isGlobal  ;
  std::vector<int  >                JpsiMu_mu3_isTracker ;
  std::vector<int  >                JpsiMu_mu3_isSoft    ;
  std::vector<float>                JpsiMu_mu3_vx       ;
  std::vector<float>                JpsiMu_mu3_vy       ;
  std::vector<float>                JpsiMu_mu3_vz       ;
  std::vector<float>                JpsiMu_mu3_iso       ;
  std::vector<float>                JpsiMu_mu3_dbiso       ;

  std::vector<float>                JpsiMu_Jpsi_pt      ;
  std::vector<float>                JpsiMu_Jpsi_eta     ;
  std::vector<float>                JpsiMu_Jpsi_phi     ;
  std::vector<float>                JpsiMu_Jpsi_mass       ;
  std::vector<float>                JpsiMu_Jpsi_vprob    ;
  std::vector<float>                JpsiMu_Jpsi_lip;
  std::vector<float>                JpsiMu_Jpsi_lips;
  std::vector<float>                JpsiMu_Jpsi_pvip;
  std::vector<float>                JpsiMu_Jpsi_pvips;
  std::vector<float>                JpsiMu_Jpsi_fl3d;
  std::vector<float>                JpsiMu_Jpsi_fls3d;
  std::vector<float>                JpsiMu_Jpsi_alpha;
  std::vector<float>                JpsiMu_Jpsi_maxdoca;
  std::vector<float>                JpsiMu_Jpsi_mindoca;
  std::vector<float>                JpsiMu_Jpsi_vx      ;
  std::vector<float>                JpsiMu_Jpsi_vy      ;
  std::vector<float>                JpsiMu_Jpsi_vz      ;
  std::vector<float>                JpsiMu_Jpsi_unfit_pt      ;
  std::vector<float>                JpsiMu_Jpsi_unfit_mass       ;
  std::vector<float>                JpsiMu_Jpsi_unfit_vprob    ;
  std::vector<float>                JpsiMu_Jpsi_unfit_vx;
  std::vector<float>                JpsiMu_Jpsi_unfit_vy;
  std::vector<float>                JpsiMu_Jpsi_unfit_vz;


  std::vector<float>                JpsiMu_B_pt      ;
  std::vector<float>                JpsiMu_B_eta     ;
  std::vector<float>                JpsiMu_B_phi     ;
  std::vector<float>                JpsiMu_B_mass    ;
  std::vector<float>                JpsiMu_B_mcorr    ;
  std::vector<float>                JpsiMu_B_vprob ;
  std::vector<float>                JpsiMu_B_lip;
  std::vector<float>                JpsiMu_B_lips;
  std::vector<float>                JpsiMu_B_pvip;
  std::vector<float>                JpsiMu_B_pvips;
  std::vector<float>                JpsiMu_B_fl3d;
  std::vector<float>                JpsiMu_B_fls3d;
  std::vector<float>                JpsiMu_B_alpha;
  std::vector<float>                JpsiMu_B_maxdoca;
  std::vector<float>                JpsiMu_B_mindoca;
  std::vector<float>                JpsiMu_B_vx      ;
  std::vector<float>                JpsiMu_B_vy      ;
  std::vector<float>                JpsiMu_B_vz      ;
  std::vector<float>                JpsiMu_B_iso;
  std::vector<int  >                JpsiMu_B_iso_ntracks;
  std::vector<float>                JpsiMu_B_iso_mindoca;
  std::vector<float>                JpsiMu_B_unfit_pt      ;
  std::vector<float>                JpsiMu_B_unfit_mass       ;
  std::vector<float>                JpsiMu_B_unfit_vprob    ;
  std::vector<float>                JpsiMu_B_unfit_vx;
  std::vector<float>                JpsiMu_B_unfit_vy;
  std::vector<float>                JpsiMu_B_unfit_vz;
  std::vector<float>                JpsiMu_B_q2;
  std::vector<float>                JpsiMu_B_mm2;
  std::vector<float>                JpsiMu_B_ptmiss;
  std::vector<float>                JpsiMu_B_Es;
  std::vector<float>                JpsiMu_B_ptback;

  std::vector<float>                JpsiMu_PV_vx       ;
  std::vector<float>                JpsiMu_PV_vy       ;
  std::vector<float>                JpsiMu_PV_vz       ;

  std::vector<float>                JpsiMu_bbPV_vx       ;
  std::vector<float>                JpsiMu_bbPV_vy       ;
  std::vector<float>                JpsiMu_bbPV_vz       ;

  //  std::vector<float>                JpsiMu_bbPV_refit_vx       ;
  //  std::vector<float>                JpsiMu_bbPV_refit_vy       ;
  //  std::vector<float>                JpsiMu_bbPV_refit_vz       ;

  std::vector<float>                JpsiMu_genPV_vx       ;
  std::vector<float>                JpsiMu_genPV_vy       ;
  std::vector<float>                JpsiMu_genPV_vz       ;

  std::vector<int  >                JpsiMu_ngenmuons      ;
  std::vector<int  >                JpsiMu_isgenmatched;
  std::vector<int  >                JpsiMu_mu3_isgenmatched;
  std::vector<float> JpsiMu_q2_gen;
  std::vector<float> JpsiMu_B_pt_gen;
  std::vector<float> JpsiMu_B_eta_gen;
  std::vector<float> JpsiMu_B_phi_gen;
  std::vector<float> JpsiMu_B_mass_gen;

  std::vector<std::vector<float>> JpsiMu_hammer_ff;
  std::vector<float> JpsiMu_hammer_ebe;
  std::vector<std::vector<float>> JpsiMu_hammer_ebe_toy;
  //  std::vector<float> JpsiMu_hammer_ebe_up;
  //  std::vector<float> JpsiMu_hammer_ebe_down;
//  std::vector<float> JpsiMu_hammer_ebe_rate_up;
  //  std::vector<float> JpsiMu_hammer_ebe_rate_down;

//  std::vector<float> JpsiMu_hammer_ebe_a0_up;
//  std::vector<float> JpsiMu_hammer_ebe_a0_down;
//  std::vector<float> JpsiMu_hammer_ebe_a1_up;
//  std::vector<float> JpsiMu_hammer_ebe_a1_down;
//  std::vector<float> JpsiMu_hammer_ebe_a2_up;
//  std::vector<float> JpsiMu_hammer_ebe_a2_down;
//
//  std::vector<float> JpsiMu_hammer_ebe_b0_up;
//  std::vector<float> JpsiMu_hammer_ebe_b0_down;
//  std::vector<float> JpsiMu_hammer_ebe_b1_up;
//  std::vector<float> JpsiMu_hammer_ebe_b1_down;
//  std::vector<float> JpsiMu_hammer_ebe_b2_up;
//  std::vector<float> JpsiMu_hammer_ebe_b2_down;
//
//  std::vector<float> JpsiMu_hammer_ebe_c1_up;
//  std::vector<float> JpsiMu_hammer_ebe_c1_down;
//  std::vector<float> JpsiMu_hammer_ebe_c2_up;
//  std::vector<float> JpsiMu_hammer_ebe_c2_down;
//
//  std::vector<float> JpsiMu_hammer_ebe_d0_up;
//  std::vector<float> JpsiMu_hammer_ebe_d0_down;
//  std::vector<float> JpsiMu_hammer_ebe_d1_up;
//  std::vector<float> JpsiMu_hammer_ebe_d1_down;
//  std::vector<float> JpsiMu_hammer_ebe_d2_up;
//  std::vector<float> JpsiMu_hammer_ebe_d2_down;



  /** HLT trigger decisions for Jpsi */
  std::map<std::string,bool> HLT_BPH_isFired;

  ////////////////////////////////////////

  // event by event quantity ... 

  int                 JpsiTau_nCandidates ;

  bool JpsiTau_isJpsiMu;
  bool JpsiTau_isJpsiTau2Mu;

  float                JpsiTau_mu1_pt      ;
  float                JpsiTau_mu1_eta     ;
  float                JpsiTau_mu1_phi     ;
  float                JpsiTau_mu1_mass     ;
  int                  JpsiTau_mu1_q      ;   
  int                  JpsiTau_mu1_isLoose   ;
  int                  JpsiTau_mu1_isTight   ;
  int                  JpsiTau_mu1_isPF      ;
  int                  JpsiTau_mu1_isGlobal  ;
  int                  JpsiTau_mu1_isTracker ;
  int                  JpsiTau_mu1_isSoft    ;
  float                JpsiTau_mu1_vx       ;
  float                JpsiTau_mu1_vy       ;
  float                JpsiTau_mu1_vz       ;
  //  float                JpsiTau_mu1_iso       ;
  float                JpsiTau_mu1_dbiso       ;

  float                JpsiTau_mu2_pt      ;
  float                JpsiTau_mu2_eta     ;
  float                JpsiTau_mu2_phi     ;
  float                JpsiTau_mu2_mass     ;
  int                  JpsiTau_mu2_q   ;   
  int                  JpsiTau_mu2_isLoose   ;
  int                  JpsiTau_mu2_isTight   ;
  int                  JpsiTau_mu2_isPF      ;
  int                  JpsiTau_mu2_isGlobal  ;
  int                  JpsiTau_mu2_isTracker ;
  int                  JpsiTau_mu2_isSoft    ;
  float                JpsiTau_mu2_vx       ;
  float                JpsiTau_mu2_vy       ;
  float                JpsiTau_mu2_vz       ;
  //  float                JpsiTau_mu2_iso       ;
  float                JpsiTau_mu2_dbiso       ;


  float                JpsiTau_PV_vx       ;
  float                JpsiTau_PV_vy       ;
  float                JpsiTau_PV_vz       ;
  
  float                JpsiTau_bbPV_vx       ;
  float                JpsiTau_bbPV_vy       ;
  float                JpsiTau_bbPV_vz       ;
  float                JpsiTau_bbPV_chi2       ;
  float                JpsiTau_bbPV_ndof       ;
  float                JpsiTau_bbPV_rho       ;

  float                JpsiTau_Jpsi_pt      ;
  float                JpsiTau_Jpsi_eta     ;
  float                JpsiTau_Jpsi_phi     ;
  float                JpsiTau_Jpsi_mass       ;
  float                JpsiTau_Jpsi_vprob    ;
  float                JpsiTau_Jpsi_lip;
  float                JpsiTau_Jpsi_lips;
  float                JpsiTau_Jpsi_pvip;
  float                JpsiTau_Jpsi_pvips;
  float                JpsiTau_Jpsi_fl3d;
  float                JpsiTau_Jpsi_fls3d;
  float                JpsiTau_Jpsi_alpha;
  float                JpsiTau_Jpsi_maxdoca;
  float                JpsiTau_Jpsi_mindoca;
  float                JpsiTau_Jpsi_vx      ;
  float                JpsiTau_Jpsi_vy      ;
  float                JpsiTau_Jpsi_vz      ;
//  float                JpsiTau_Jpsi_unfit_pt      ;
//  float                JpsiTau_Jpsi_unfit_mass       ;
//  float                JpsiTau_Jpsi_unfit_vprob    ;
//  float                JpsiTau_Jpsi_unfit_vx;
//  float                JpsiTau_Jpsi_unfit_vy;
//  float                JpsiTau_Jpsi_unfit_vz;






  //  std::vector<float>                JpsiTau_tau_fullfit_pt      ;
  //  std::vector<float>                JpsiTau_tau_fullfit_eta     ;
  //  std::vector<float>                JpsiTau_tau_fullfit_phi     ;
  //  std::vector<float>                JpsiTau_tau_fullfit_mass     ;
  std::vector<float>                JpsiTau_tau_pt      ;
  std::vector<float>                JpsiTau_tau_eta     ;
  std::vector<float>                JpsiTau_tau_phi     ;
  std::vector<float>                JpsiTau_tau_mass     ;
  std::vector<int  >                JpsiTau_tau_q   ;   
  std::vector<float>                JpsiTau_tau_vx       ;
  std::vector<float>                JpsiTau_tau_vy       ;
  std::vector<float>                JpsiTau_tau_vz       ;

  std::vector<float>       JpsiTau_tau_max_dr_3prong;
  std::vector<float>       JpsiTau_tau_lip;
  std::vector<float>       JpsiTau_tau_lips;
  std::vector<float>       JpsiTau_tau_pvip;
  std::vector<float>       JpsiTau_tau_pvips;
  std::vector<float>       JpsiTau_tau_fl3d;
  std::vector<float>       JpsiTau_tau_fls3d;
  std::vector<float>       JpsiTau_tau_alpha;
  std::vector<float>       JpsiTau_tau_vprob;

  std::vector<float>       JpsiTau_tau_fl3d_wjpsi;
  std::vector<float>       JpsiTau_tau_fls3d_wjpsi;

//  std::vector<float>        JpsiTau_tau_dr1;
//  std::vector<float>        JpsiTau_tau_dr2;
//  std::vector<float>        JpsiTau_tau_dr3;
//  std::vector<float>        JpsiTau_tau_ptres1;
//  std::vector<float>        JpsiTau_tau_ptres2;
//  std::vector<float>        JpsiTau_tau_ptres3;
//  std::vector<int>         JpsiTau_tau_matched_ppdgId;
//  std::vector<float>       JpsiTau_tau_matched_gentaupt;
  std::vector<float>       JpsiTau_tau_sumofdnn; 
  std::vector<float>       JpsiTau_tau_sumofdnn_1prong; 
  std::vector<float>       JpsiTau_tau_sumofdnn_otherB; 
  std::vector<float>       JpsiTau_tau_sumofdnn_pu; 
  //  std::vector<float>       JpsiTau_tau_sumofdnn_old; 
  //  std::vector<float>       JpsiTau_tau_sumofdnn_others; 

//  std::vector<float>       JpsiTau_tau_pi1_dnn;
  //  std::vector<float>       JpsiTau_tau_pi2_dnn;
  //  std::vector<float>       JpsiTau_tau_pi3_dnn;

  //  std::vector<float>       JpsiTau_tau_pi1_doca;
  //  std::vector<float>       JpsiTau_tau_pi2_doca;
  //  std::vector<float>       JpsiTau_tau_pi3_doca;

  //  std::vector<int>       JpsiTau_tau_pi1_pv;
  //  std::vector<int>       JpsiTau_tau_pi2_pv;
  //  std::vector<int>       JpsiTau_tau_pi3_pv;


  std::vector<float>                JpsiTau_tau_rhomass1     ;
  std::vector<float>                JpsiTau_tau_rhomass2     ;
  std::vector<float>                JpsiTau_tau_rhomass_ss     ;
  std::vector<float>                JpsiTau_tau_rhopt1     ;
  std::vector<float>                JpsiTau_tau_rhopt2     ;

  std::vector<float>       JpsiTau_tau_pi1_pt;
  std::vector<float>       JpsiTau_tau_pi1_eta;
  std::vector<float>       JpsiTau_tau_pi1_phi;
  std::vector<float>       JpsiTau_tau_pi1_mass;
  std::vector<int>       JpsiTau_tau_pi1_q;

  std::vector<float> JpsiTau_tau_pi1_doca3d;
  std::vector<float> JpsiTau_tau_pi1_doca3de;
  std::vector<float> JpsiTau_tau_pi1_doca2d;
  std::vector<float> JpsiTau_tau_pi1_doca2de;
  //  std::vector<float> JpsiTau_tau_pi1_doca1d;
  //  std::vector<float> JpsiTau_tau_pi1_doca1de;
  //  std::vector<bool> JpsiTau_tau_pi1_isRight;
  std::vector<float> JpsiTau_tau_pi1_dz;
  std::vector<float> JpsiTau_tau_pi1_near_dz;
  std::vector<bool> JpsiTau_tau_pi1_isAssociate;
  std::vector<int> JpsiTau_tau_pi1_pvAssociationQuality;
  std::vector<bool> JpsiTau_tau_pi1_isBdecay;
  std::vector<int> JpsiTau_tau_pi1_isBdecaypdg;
  std::vector<int> JpsiTau_tau_pi1_isBdecayppdg;
  std::vector<bool> JpsiTau_tau_pi1_isSignal;
  std::vector<int> JpsiTau_tau_pi1_nprong;
  std::vector<int> JpsiTau_tau_pi1_nprong_pi0;

  std::vector<float> JpsiTau_tau_pi1_dnn;
  std::vector<float> JpsiTau_tau_pi1_dnn_1prong;
  std::vector<float> JpsiTau_tau_pi1_dnn_otherB;
  std::vector<float> JpsiTau_tau_pi1_dnn_pu;
  //  std::vector<float> JpsiTau_tau_pi1_dnn_old;

  std::vector<bool> JpsiTau_tau_pi1_trigMatch;
  std::vector<float> JpsiTau_tau_pi1_trigMatch_dr;

  std::vector<float>       JpsiTau_tau_pi2_pt;
  std::vector<float>       JpsiTau_tau_pi2_eta;
  std::vector<float>       JpsiTau_tau_pi2_phi;
  std::vector<float>       JpsiTau_tau_pi2_mass;
  std::vector<int>       JpsiTau_tau_pi2_q;

  std::vector<float> JpsiTau_tau_pi2_doca3d;
  std::vector<float> JpsiTau_tau_pi2_doca3de;
  std::vector<float> JpsiTau_tau_pi2_doca2d;
  std::vector<float> JpsiTau_tau_pi2_doca2de;
  //  std::vector<float> JpsiTau_tau_pi2_doca1d;
  //  std::vector<float> JpsiTau_tau_pi2_doca1de;
  //  std::vector<bool> JpsiTau_tau_pi2_isRight;
  std::vector<float> JpsiTau_tau_pi2_dz;
  std::vector<float> JpsiTau_tau_pi2_near_dz;
  std::vector<bool> JpsiTau_tau_pi2_isAssociate;
  std::vector<int> JpsiTau_tau_pi2_pvAssociationQuality;
  std::vector<bool> JpsiTau_tau_pi2_isBdecay;
  std::vector<int> JpsiTau_tau_pi2_isBdecaypdg;
  std::vector<int> JpsiTau_tau_pi2_isBdecayppdg;
  std::vector<bool> JpsiTau_tau_pi2_isSignal;
  std::vector<int> JpsiTau_tau_pi2_nprong;
  std::vector<int> JpsiTau_tau_pi2_nprong_pi0;
  std::vector<float> JpsiTau_tau_pi2_dnn;
  std::vector<float> JpsiTau_tau_pi2_dnn_1prong;
  std::vector<float> JpsiTau_tau_pi2_dnn_otherB;
  std::vector<float> JpsiTau_tau_pi2_dnn_pu;
  //  std::vector<float> JpsiTau_tau_pi2_dnn_old;

  std::vector<bool> JpsiTau_tau_pi2_trigMatch;
  std::vector<float> JpsiTau_tau_pi2_trigMatch_dr;


  std::vector<float>       JpsiTau_tau_pi3_pt;
  std::vector<float>       JpsiTau_tau_pi3_eta;
  std::vector<float>       JpsiTau_tau_pi3_phi;
  std::vector<float>       JpsiTau_tau_pi3_mass;
  std::vector<int>       JpsiTau_tau_pi3_q;

  std::vector<float> JpsiTau_tau_pi3_doca3d;
  std::vector<float> JpsiTau_tau_pi3_doca3de;
  std::vector<float> JpsiTau_tau_pi3_doca2d;
  std::vector<float> JpsiTau_tau_pi3_doca2de;
  //  std::vector<float> JpsiTau_tau_pi3_doca1d;
  //  std::vector<float> JpsiTau_tau_pi3_doca1de;
  //  std::vector<bool> JpsiTau_tau_pi3_isRight;
  std::vector<float> JpsiTau_tau_pi3_dz;
  std::vector<float> JpsiTau_tau_pi3_near_dz;
  std::vector<bool> JpsiTau_tau_pi3_isAssociate;
  std::vector<int> JpsiTau_tau_pi3_pvAssociationQuality;
  std::vector<bool> JpsiTau_tau_pi3_isBdecay;
  std::vector<int> JpsiTau_tau_pi3_isBdecaypdg;
  std::vector<int> JpsiTau_tau_pi3_isBdecayppdg;
  std::vector<bool> JpsiTau_tau_pi3_isSignal;
  std::vector<int> JpsiTau_tau_pi3_nprong;
  std::vector<int> JpsiTau_tau_pi3_nprong_pi0;

  std::vector<float> JpsiTau_tau_pi3_dnn;
  std::vector<float> JpsiTau_tau_pi3_dnn_1prong;
  std::vector<float> JpsiTau_tau_pi3_dnn_otherB;
  std::vector<float> JpsiTau_tau_pi3_dnn_pu;
  //  std::vector<float> JpsiTau_tau_pi3_dnn_old;

  std::vector<bool> JpsiTau_tau_pi3_trigMatch;
  std::vector<float> JpsiTau_tau_pi3_trigMatch_dr;


  std::vector<float> JpsiTau_tau_delta_chi2;
  std::vector<int> JpsiTau_tau_delta_n_ch;
  std::vector<int> JpsiTau_tau_delta_n_mu;
  std::vector<float> JpsiTau_tau_vweight;
  std::vector<float> JpsiTau_tau_refit_vx;
  std::vector<float> JpsiTau_tau_refit_vy;
  std::vector<float> JpsiTau_tau_refit_vz;
  std::vector<float> JpsiTau_tau_refit_chi2;
  std::vector<float> JpsiTau_tau_refit_ndof;
  std::vector<float> JpsiTau_tau_refit_rho;

  std::vector<std::vector<float>>                JpsiTau_tau_iso;
  std::vector<std::vector<int > >                JpsiTau_tau_iso_ntracks;
  std::vector<std::vector<float>>                JpsiTau_tau_iso_mindoca;


  std::vector<float>       JpsiTau_ptbal;
  std::vector<float>       JpsiTau_jpsi_tau_alpha;


  std::vector<float>                JpsiTau_B_pt      ;
  std::vector<float>                JpsiTau_B_eta     ;
  std::vector<float>                JpsiTau_B_phi     ;
  std::vector<float>                JpsiTau_B_mass    ;
  std::vector<float>                JpsiTau_B_mcorr    ;
  std::vector<float>                JpsiTau_B_vprob ;
  std::vector<float>                JpsiTau_B_lip;
  std::vector<float>                JpsiTau_B_lips;
  std::vector<float>                JpsiTau_B_pvip;
  std::vector<float>                JpsiTau_B_pvips;
  std::vector<float>                JpsiTau_B_fl3d;
  std::vector<float>                JpsiTau_B_fls3d;
  std::vector<float>                JpsiTau_B_alpha;
  std::vector<float>                JpsiTau_B_maxdoca;
  std::vector<float>                JpsiTau_B_mindoca;
  std::vector<float>                JpsiTau_B_vx      ;
  std::vector<float>                JpsiTau_B_vy      ;
  std::vector<float>                JpsiTau_B_vz      ;


  //  std::vector<float>                JpsiTau_B_iso_nocut;
  //  std::vector<int  >                JpsiTau_B_iso_ntracks_nocut;
  //  std::vector<float>                JpsiTau_B_iso_mindoca_nocut;

//  std::vector<float>                JpsiTau_B_unfit_pt      ;
//  std::vector<float>                JpsiTau_B_unfit_mass       ;
//  std::vector<float>                JpsiTau_B_unfit_vprob    ;
//  std::vector<float>                JpsiTau_B_unfit_vx;
//  std::vector<float>                JpsiTau_B_unfit_vy;
//  std::vector<float>                JpsiTau_B_unfit_vz;
  std::vector<float>                JpsiTau_B_q2;
  std::vector<float>                JpsiTau_B_mm2;
  std::vector<float>                JpsiTau_B_ptmiss;
  std::vector<float>                JpsiTau_B_Es;
  std::vector<float>                JpsiTau_B_ptback;

  std::vector<float>                JpsiTau_B_pt_simple;
  std::vector<float>                JpsiTau_B_eta_simple;
  std::vector<float>                JpsiTau_B_phi_simple;
  std::vector<float>                JpsiTau_B_mass_simple;
  std::vector<float>                JpsiTau_B_q2_simple;
  std::vector<float>                JpsiTau_B_mm2_simple;
  std::vector<float>                JpsiTau_B_ptmiss_simple;
  std::vector<float>                JpsiTau_B_Es_simple;
  std::vector<float>                JpsiTau_B_ptback_simple;


  //  std::vector<float>                JpsiTau_bbPV_refit_vx       ;
  //  std::vector<float>                JpsiTau_bbPV_refit_vy       ;
  //  std::vector<float>                JpsiTau_bbPV_refit_vz       ;

  
  // gen related quantity, event-by-event
  float                JpsiTau_genPV_vx       ;
  float                JpsiTau_genPV_vy       ;
  float                JpsiTau_genPV_vz       ;
  float                JpsiTau_genSV_vx       ;
  float                JpsiTau_genSV_vy       ;
  float                JpsiTau_genSV_vz       ;

  int                  JpsiTau_ngenmuons      ;
  int                  JpsiTau_isgenmatched;
  //  int                  JpsiTau_isgen3;
  //  int                  JpsiTau_isgen3matched;

  float genWeightBkgB;

  int JpsiTau_nch;
  //  int JpsiTau_nch_after_dnn;
  int JpsiTau_nch_before;
  //  int JpsiTau_nch_qr;

//  int JpsiTau_ngentau3;
//  int JpsiTau_ngentau;
//  float JpsiTau_gentaupt;
//  float JpsiTau_gentaueta;
//  float JpsiTau_gentauphi;
//  float JpsiTau_gentaumass;
//  float JpsiTau_gentaupt_bd;
//  float JpsiTau_gentaueta_bd;
//  float JpsiTau_gentauphi_bd;
//  float JpsiTau_gentaumass_bd;
//  int JpsiTau_gentaudm;
  float JpsiTau_q2_gen;
  int JpsiTau_nBc;
  float JpsiTau_B_pt_gen;
  float JpsiTau_B_eta_gen;
  float JpsiTau_B_phi_gen;
  float JpsiTau_B_mass_gen;

  std::vector<float> JpsiTau_gen_pion_pt;
  std::vector<float> JpsiTau_gen_pion_eta;
  std::vector<float> JpsiTau_gen_pion_phi;
  std::vector<bool> JpsiTau_gen_pion_matched;

  std::vector<float> JpsiTau_gen_tau_pt;
  std::vector<float> JpsiTau_gen_tau_eta;
  std::vector<float> JpsiTau_gen_tau_phi;
  std::vector<int> JpsiTau_gen_tau_nprong;
  std::vector<int> JpsiTau_gen_tau_nmatched;


  Int_t JpsiTau_st_nch;
  Int_t JpsiTau_st_nch_matched;
  //  Int_t JpsiTau_st_npi0;
  Int_t JpsiTau_st_n_charged_pions = 0;
  Int_t JpsiTau_st_n_neutral_pions = 0;
  Int_t JpsiTau_st_n_mu_decay = 0;
  Int_t JpsiTau_st_n_e_decay = 0;
  Int_t JpsiTau_st_n_occurance = 0;
  Int_t JpsiTau_st_decayid;
  Float_t JpsiTau_st_gentau_pt;
  Float_t JpsiTau_st_gentau_eta;
  Float_t JpsiTau_st_gentau_phi;
  Float_t JpsiTau_st_genjpsi_pt;
  Float_t JpsiTau_st_genjpsi_eta;
  Float_t JpsiTau_st_genjpsi_phi;

  Float_t JpsiTau_perEVT_mc;
  Float_t JpsiTau_perEVT_data;
  //  Float_t JpsiTau_perEVT_otherB;
  //  Float_t JpsiTau_perEVT_sig;
  //  Float_t JpsiTau_perEVT_leptonic;
  //  Float_t JpsiTau_perEVT_1prong;
 
  std::vector<int> JpsiTau_st_idx;

  std::vector<float> JpsiTau_st_doca3d;
  std::vector<float> JpsiTau_st_doca2d;
  std::vector<float> JpsiTau_st_doca3ds;
  std::vector<float> JpsiTau_st_doca2ds;
  std::vector<float> JpsiTau_st_doca3de;
  std::vector<float> JpsiTau_st_doca2de;

  std::vector<float> JpsiTau_st_doca3d_max;
  std::vector<float> JpsiTau_st_doca2d_max;
  std::vector<float> JpsiTau_st_doca3ds_max;
  std::vector<float> JpsiTau_st_doca2ds_max;
  std::vector<float> JpsiTau_st_doca3de_max;
  std::vector<float> JpsiTau_st_doca2de_max;
  std::vector<float> JpsiTau_st_dz_max;

  //  std::vector<bool> JpsiTau_st_isRight;
  std::vector<bool> JpsiTau_st_isBdecay;
  std::vector<bool> JpsiTau_st_isSignal;
  std::vector<int> JpsiTau_st_isBdecaypdg;
  std::vector<int> JpsiTau_st_isBdecayppdg;
  std::vector<int> JpsiTau_st_nprong;
  std::vector<int> JpsiTau_st_nprong_pi0;
  std::vector<float> JpsiTau_st_dz;
  std::vector<bool> JpsiTau_st_isAssociate;
  std::vector<float> JpsiTau_st_near_dz;
  std::vector<float> JpsiTau_st_dr_jpsi;
  std::vector<bool> JpsiTau_st_trigMatch;
  std::vector<float> JpsiTau_st_trigMatch_dr;
  std::vector<int> JpsiTau_st_pvAssociationQuality;
  std::vector<float> JpsiTau_st_pt;
  std::vector<float> JpsiTau_st_eta;
  std::vector<float> JpsiTau_st_phi;
  std::vector<int> JpsiTau_st_charge;
  std::vector<float> JpsiTau_st_mass;
  std::vector<float> JpsiTau_st_dnn;
  std::vector<float> JpsiTau_st_dnn_1prong;
  std::vector<float> JpsiTau_st_dnn_otherB;
  std::vector<float> JpsiTau_st_dnn_pu;
  //  std::vector<float> JpsiTau_st_dnn_old;
  std::vector<int> JpsiTau_st_matchidx;

//  std::vector<float> JpsiTau_ed_pfeta;
//  std::vector<float> JpsiTau_ed_pfphi;
//  std::vector<int> JpsiTau_ed_isRight;
//  std::vector<int> JpsiTau_ed_id;
//  std::vector<float> JpsiTau_ed_pfdnn;
//  std::vector<float> JpsiTau_ed_genpt;
    
//  std::vector<std::vector<float>> JpsiTau_hammer_ff;

  //  std::vector<std::vector<float>> JpsiTau_hammer_ebe_toy;
  //  std::vector<float> JpsiTau_hammer_ebe_up;
  //  std::vector<float> JpsiTau_hammer_ebe_down;
  //  std::vector<float> JpsiTau_hammer_ebe_rate_up;
//  std::vector<float> JpsiTau_hammer_ebe_rate_down;

  std::vector<float> JpsiTau_hammer_ebe;
  std::vector<float> JpsiTau_hammer_ebe_e0_up;
  std::vector<float> JpsiTau_hammer_ebe_e0_down;
  std::vector<float> JpsiTau_hammer_ebe_e1_up;
  std::vector<float> JpsiTau_hammer_ebe_e1_down;
  std::vector<float> JpsiTau_hammer_ebe_e2_up;
  std::vector<float> JpsiTau_hammer_ebe_e2_down;
  std::vector<float> JpsiTau_hammer_ebe_e3_up;
  std::vector<float> JpsiTau_hammer_ebe_e3_down;
  std::vector<float> JpsiTau_hammer_ebe_e4_up;
  std::vector<float> JpsiTau_hammer_ebe_e4_down;
  std::vector<float> JpsiTau_hammer_ebe_e5_up;
  std::vector<float> JpsiTau_hammer_ebe_e5_down;
  std::vector<float> JpsiTau_hammer_ebe_e6_up;
  std::vector<float> JpsiTau_hammer_ebe_e6_down;
  std::vector<float> JpsiTau_hammer_ebe_e7_up;
  std::vector<float> JpsiTau_hammer_ebe_e7_down;
  std::vector<float> JpsiTau_hammer_ebe_e8_up;
  std::vector<float> JpsiTau_hammer_ebe_e8_down;
  std::vector<float> JpsiTau_hammer_ebe_e9_up;
  std::vector<float> JpsiTau_hammer_ebe_e9_down;
  std::vector<float> JpsiTau_hammer_ebe_e10_up;
  std::vector<float> JpsiTau_hammer_ebe_e10_down;
  std::vector<float> JpsiTau_hammer_ebe_e11_up;
  std::vector<float> JpsiTau_hammer_ebe_e11_down;
  std::vector<float> JpsiTau_hammer_ebe_e12_up;
  std::vector<float> JpsiTau_hammer_ebe_e12_down;
  std::vector<float> JpsiTau_hammer_ebe_e13_up;
  std::vector<float> JpsiTau_hammer_ebe_e13_down;
  std::vector<float> JpsiTau_hammer_ebe_e14_up;
  std::vector<float> JpsiTau_hammer_ebe_e14_down;


  std::vector<float> JpsiTau_hammer_ebe_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e0_up_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e0_down_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e1_up_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e1_down_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e2_up_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e2_down_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e3_up_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e3_down_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e4_up_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e4_down_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e5_up_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e5_down_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e6_up_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e6_down_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e7_up_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e7_down_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e8_up_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e8_down_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e9_up_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e9_down_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e10_up_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e10_down_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e11_up_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e11_down_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e12_up_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e12_down_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e13_up_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e13_down_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e14_up_lattice;
  std::vector<float> JpsiTau_hammer_ebe_e14_down_lattice;







  ////////////////////////


















  ////////////////////////////////////////

  // event by event quantity ... 

  int                 JpsiK_nCandidates ;
  float                JpsiK_mu1_pt      ;
  float                JpsiK_mu1_eta     ;
  float                JpsiK_mu1_phi     ;
  float                JpsiK_mu1_mass     ;
  int                  JpsiK_mu1_q      ;   
  int                  JpsiK_mu1_isLoose   ;
  int                  JpsiK_mu1_isTight   ;
  int                  JpsiK_mu1_isPF      ;
  int                  JpsiK_mu1_isGlobal  ;
  int                  JpsiK_mu1_isTracker ;
  int                  JpsiK_mu1_isSoft    ;
  float                JpsiK_mu1_vx       ;
  float                JpsiK_mu1_vy       ;
  float                JpsiK_mu1_vz       ;
  //  float                JpsiK_mu1_iso       ;
  float                JpsiK_mu1_dbiso       ;

  float                JpsiK_mu2_pt      ;
  float                JpsiK_mu2_eta     ;
  float                JpsiK_mu2_phi     ;
  float                JpsiK_mu2_mass     ;
  int                  JpsiK_mu2_q   ;   
  int                  JpsiK_mu2_isLoose   ;
  int                  JpsiK_mu2_isTight   ;
  int                  JpsiK_mu2_isPF      ;
  int                  JpsiK_mu2_isGlobal  ;
  int                  JpsiK_mu2_isTracker ;
  int                  JpsiK_mu2_isSoft    ;
  float                JpsiK_mu2_vx       ;
  float                JpsiK_mu2_vy       ;
  float                JpsiK_mu2_vz       ;
  //  float                JpsiK_mu2_iso       ;
  float                JpsiK_mu2_dbiso       ;


  float                JpsiK_PV_vx       ;
  float                JpsiK_PV_vy       ;
  float                JpsiK_PV_vz       ;
  
  float                JpsiK_bbPV_vx       ;
  float                JpsiK_bbPV_vy       ;
  float                JpsiK_bbPV_vz       ;
  float                JpsiK_bbPV_chi2       ;
  float                JpsiK_bbPV_ndof       ;
  float                JpsiK_bbPV_rho       ;

  float                JpsiK_Jpsi_pt      ;
  float                JpsiK_Jpsi_eta     ;
  float                JpsiK_Jpsi_phi     ;
  float                JpsiK_Jpsi_mass       ;
  float                JpsiK_Jpsi_vprob    ;
  float                JpsiK_Jpsi_lip;
  float                JpsiK_Jpsi_lips;
  float                JpsiK_Jpsi_pvip;
  float                JpsiK_Jpsi_pvips;
  float                JpsiK_Jpsi_fl3d;
  float                JpsiK_Jpsi_fls3d;
  float                JpsiK_Jpsi_alpha;
  float                JpsiK_Jpsi_maxdoca;
  float                JpsiK_Jpsi_mindoca;
  float                JpsiK_Jpsi_vx      ;
  float                JpsiK_Jpsi_vy      ;
  float                JpsiK_Jpsi_vz      ;
//  float                JpsiK_Jpsi_unfit_pt      ;
//  float                JpsiK_Jpsi_unfit_mass       ;
//  float                JpsiK_Jpsi_unfit_vprob    ;
//  float                JpsiK_Jpsi_unfit_vx;
//  float                JpsiK_Jpsi_unfit_vy;
//  float                JpsiK_Jpsi_unfit_vz;





  std::vector<float>                JpsiK_B_pt      ;
  std::vector<float>                JpsiK_B_eta     ;
  std::vector<float>                JpsiK_B_phi     ;
  std::vector<float>                JpsiK_B_mass    ;
  std::vector<float>                JpsiK_B_mcorr    ;
  std::vector<float>                JpsiK_B_vprob ;
  std::vector<float>                JpsiK_B_lip;
  std::vector<float>                JpsiK_B_lips;
  std::vector<float>                JpsiK_B_pvip;
  std::vector<float>                JpsiK_B_pvips;
  std::vector<float>                JpsiK_B_fl3d;
  std::vector<float>                JpsiK_B_fls3d;
  std::vector<float>                JpsiK_B_alpha;
  std::vector<float>                JpsiK_B_maxdoca;
  std::vector<float>                JpsiK_B_mindoca;
  std::vector<float>                JpsiK_B_vx      ;
  std::vector<float>                JpsiK_B_vy      ;
  std::vector<float>                JpsiK_B_vz      ;


  std::vector<float>       JpsiK_pi_pt;
  std::vector<float>       JpsiK_pi_eta;
  std::vector<float>       JpsiK_pi_phi;
  std::vector<float>       JpsiK_pi_mass;
  std::vector<int>       JpsiK_pi_q;

  std::vector<float> JpsiK_pi_doca3d;
  std::vector<float> JpsiK_pi_doca3de;
  std::vector<float> JpsiK_pi_doca2d;
  std::vector<float> JpsiK_pi_doca2de;
  std::vector<float> JpsiK_pi_dz;
  std::vector<float> JpsiK_pi_near_dz;
  std::vector<bool> JpsiK_pi_isAssociate;
  std::vector<int> JpsiK_pi_pvAssociationQuality;
  std::vector<bool> JpsiK_pi_trigMatch;
  std::vector<float> JpsiK_pi_trigMatch_dr;


  //  std::vector<float>                JpsiK_B_iso_nocut;
  //  std::vector<int  >                JpsiK_B_iso_ntracks_nocut;
  //  std::vector<float>                JpsiK_B_iso_mindoca_nocut;

//  std::vector<float>                JpsiK_B_unfit_pt      ;
//  std::vector<float>                JpsiK_B_unfit_mass       ;
//  std::vector<float>                JpsiK_B_unfit_vprob    ;
//  std::vector<float>                JpsiK_B_unfit_vx;
//  std::vector<float>                JpsiK_B_unfit_vy;
//  std::vector<float>                JpsiK_B_unfit_vz;
  std::vector<float>                JpsiK_B_q2;
  std::vector<float>                JpsiK_B_mm2;
  std::vector<float>                JpsiK_B_ptmiss;
  std::vector<float>                JpsiK_B_Es;
  std::vector<float>                JpsiK_B_ptback;



  //  std::vector<float>                JpsiK_bbPV_refit_vx       ;
  //  std::vector<float>                JpsiK_bbPV_refit_vy       ;
  //  std::vector<float>                JpsiK_bbPV_refit_vz       ;

  
  // gen related quantity, event-by-event
  float                JpsiK_genPV_vx       ;
  float                JpsiK_genPV_vy       ;
  float                JpsiK_genPV_vz       ;
  float                JpsiK_genSV_vx       ;
  float                JpsiK_genSV_vy       ;
  float                JpsiK_genSV_vz       ;

  int                  JpsiK_ngenmuons      ;
  int                  JpsiK_isgenmatched;
  //  int                  JpsiK_isgen3;
  //  int                  JpsiK_isgen3matched;


  int JpsiK_nch;
  //  int JpsiK_nch_after_dnn;
  int JpsiK_nch_before;
  //  int JpsiK_nch_qr;

//  int JpsiK_ngentau3;
//  int JpsiK_ngentau;
//  float JpsiK_gentaupt;
//  float JpsiK_gentaueta;
//  float JpsiK_gentauphi;
//  float JpsiK_gentaumass;
//  float JpsiK_gentaupt_bd;
//  float JpsiK_gentaueta_bd;
//  float JpsiK_gentauphi_bd;
//  float JpsiK_gentaumass_bd;
//  int JpsiK_gentaudm;
  float JpsiK_q2_gen;
  int JpsiK_nBc;
  float JpsiK_B_pt_gen;
  float JpsiK_B_eta_gen;
  float JpsiK_B_phi_gen;
  float JpsiK_B_mass_gen;



  ///////////////////////////////////////////////////////






  //////////////////////////////////////// 
  // event by event quantity ... 

  int                 JpsiK_e_nCandidates ;
  float                JpsiK_e_mu1_pt      ;
  float                JpsiK_e_mu1_eta     ;
  float                JpsiK_e_mu1_phi     ;
  float                JpsiK_e_mu1_mass     ;
  int                  JpsiK_e_mu1_q      ;   
  int                  JpsiK_e_mu1_isLoose   ;
  int                  JpsiK_e_mu1_isTight   ;
  int                  JpsiK_e_mu1_isPF      ;
  int                  JpsiK_e_mu1_isGlobal  ;
  int                  JpsiK_e_mu1_isTracker ;
  int                  JpsiK_e_mu1_isSoft    ;
  float                JpsiK_e_mu1_vx       ;
  float                JpsiK_e_mu1_vy       ;
  float                JpsiK_e_mu1_vz       ;
  //  float                JpsiK_e_mu1_iso       ;
  float                JpsiK_e_mu1_dbiso       ;

  float                JpsiK_e_mu2_pt      ;
  float                JpsiK_e_mu2_eta     ;
  float                JpsiK_e_mu2_phi     ;
  float                JpsiK_e_mu2_mass     ;
  int                  JpsiK_e_mu2_q   ;   
  int                  JpsiK_e_mu2_isLoose   ;
  int                  JpsiK_e_mu2_isTight   ;
  int                  JpsiK_e_mu2_isPF      ;
  int                  JpsiK_e_mu2_isGlobal  ;
  int                  JpsiK_e_mu2_isTracker ;
  int                  JpsiK_e_mu2_isSoft    ;
  float                JpsiK_e_mu2_vx       ;
  float                JpsiK_e_mu2_vy       ;
  float                JpsiK_e_mu2_vz       ;
  //  float                JpsiK_e_mu2_iso       ;
  float                JpsiK_e_mu2_dbiso       ;


  float                JpsiK_e_PV_vx       ;
  float                JpsiK_e_PV_vy       ;
  float                JpsiK_e_PV_vz       ;
  
  float                JpsiK_e_bbPV_vx       ;
  float                JpsiK_e_bbPV_vy       ;
  float                JpsiK_e_bbPV_vz       ;
  float                JpsiK_e_bbPV_chi2       ;
  float                JpsiK_e_bbPV_ndof       ;
  float                JpsiK_e_bbPV_rho       ;

  float                JpsiK_e_Jpsi_pt      ;
  float                JpsiK_e_Jpsi_eta     ;
  float                JpsiK_e_Jpsi_phi     ;
  float                JpsiK_e_Jpsi_mass       ;
  float                JpsiK_e_Jpsi_vprob    ;
  float                JpsiK_e_Jpsi_lip;
  float                JpsiK_e_Jpsi_lips;
  float                JpsiK_e_Jpsi_pvip;
  float                JpsiK_e_Jpsi_pvips;
  float                JpsiK_e_Jpsi_fl3d;
  float                JpsiK_e_Jpsi_fls3d;
  float                JpsiK_e_Jpsi_alpha;
  float                JpsiK_e_Jpsi_maxdoca;
  float                JpsiK_e_Jpsi_mindoca;
  float                JpsiK_e_Jpsi_vx      ;
  float                JpsiK_e_Jpsi_vy      ;
  float                JpsiK_e_Jpsi_vz      ;
//  float                JpsiK_e_Jpsi_unfit_pt      ;
//  float                JpsiK_e_Jpsi_unfit_mass       ;
//  float                JpsiK_e_Jpsi_unfit_vprob    ;
//  float                JpsiK_e_Jpsi_unfit_vx;
//  float                JpsiK_e_Jpsi_unfit_vy;
//  float                JpsiK_e_Jpsi_unfit_vz;





  std::vector<float>                JpsiK_e_B_pt      ;
  std::vector<float>                JpsiK_e_B_eta     ;
  std::vector<float>                JpsiK_e_B_phi     ;
  std::vector<float>                JpsiK_e_B_mass    ;
  std::vector<float>                JpsiK_e_B_mcorr    ;
  std::vector<float>                JpsiK_e_B_vprob ;
  std::vector<float>                JpsiK_e_B_lip;
  std::vector<float>                JpsiK_e_B_lips;
  std::vector<float>                JpsiK_e_B_pvip;
  std::vector<float>                JpsiK_e_B_pvips;
  std::vector<float>                JpsiK_e_B_fl3d;
  std::vector<float>                JpsiK_e_B_fls3d;
  std::vector<float>                JpsiK_e_B_alpha;
  std::vector<float>                JpsiK_e_B_maxdoca;
  std::vector<float>                JpsiK_e_B_mindoca;
  std::vector<float>                JpsiK_e_B_vx      ;
  std::vector<float>                JpsiK_e_B_vy      ;
  std::vector<float>                JpsiK_e_B_vz      ;


  std::vector<float>       JpsiK_e_pi_pt;
  std::vector<float>       JpsiK_e_pi_eta;
  std::vector<float>       JpsiK_e_pi_phi;
  std::vector<float>       JpsiK_e_pi_mass;
  std::vector<int>       JpsiK_e_pi_q;

  std::vector<float> JpsiK_e_pi_doca3d;
  std::vector<float> JpsiK_e_pi_doca3de;
  std::vector<float> JpsiK_e_pi_doca2d;
  std::vector<float> JpsiK_e_pi_doca2de;
  std::vector<float> JpsiK_e_pi_dz;
  std::vector<float> JpsiK_e_pi_near_dz;
  std::vector<bool> JpsiK_e_pi_isAssociate;
  std::vector<int> JpsiK_e_pi_pvAssociationQuality;
  std::vector<bool> JpsiK_e_pi_trigMatch;
  std::vector<float> JpsiK_e_pi_trigMatch_dr;


  //  std::vector<float>                JpsiK_e_B_iso_nocut;
  //  std::vector<int  >                JpsiK_e_B_iso_ntracks_nocut;
  //  std::vector<float>                JpsiK_e_B_iso_mindoca_nocut;

//  std::vector<float>                JpsiK_e_B_unfit_pt      ;
//  std::vector<float>                JpsiK_e_B_unfit_mass       ;
//  std::vector<float>                JpsiK_e_B_unfit_vprob    ;
//  std::vector<float>                JpsiK_e_B_unfit_vx;
//  std::vector<float>                JpsiK_e_B_unfit_vy;
//  std::vector<float>                JpsiK_e_B_unfit_vz;
  std::vector<float>                JpsiK_e_B_q2;
  std::vector<float>                JpsiK_e_B_mm2;
  std::vector<float>                JpsiK_e_B_ptmiss;
  std::vector<float>                JpsiK_e_B_Es;
  std::vector<float>                JpsiK_e_B_ptback;



  //  std::vector<float>                JpsiK_e_bbPV_refit_vx       ;
  //  std::vector<float>                JpsiK_e_bbPV_refit_vy       ;
  //  std::vector<float>                JpsiK_e_bbPV_refit_vz       ;

  
  // gen related quantity, event-by-event
  float                JpsiK_e_genPV_vx       ;
  float                JpsiK_e_genPV_vy       ;
  float                JpsiK_e_genPV_vz       ;
  float                JpsiK_e_genSV_vx       ;
  float                JpsiK_e_genSV_vy       ;
  float                JpsiK_e_genSV_vz       ;

  int                  JpsiK_e_ngenmuons      ;
  int                  JpsiK_e_isgenmatched;
  //  int                  JpsiK_e_isgen3;
  //  int                  JpsiK_e_isgen3matched;


  int JpsiK_e_nch;
  //  int JpsiK_e_nch_after_dnn;
  int JpsiK_e_nch_before;
  //  int JpsiK_e_nch_qr;

//  int JpsiK_e_ngentau3;
//  int JpsiK_e_ngentau;
//  float JpsiK_e_gentaupt;
//  float JpsiK_e_gentaueta;
//  float JpsiK_e_gentauphi;
//  float JpsiK_e_gentaumass;
//  float JpsiK_e_gentaupt_bd;
//  float JpsiK_e_gentaueta_bd;
//  float JpsiK_e_gentauphi_bd;
//  float JpsiK_e_gentaumass_bd;
//  int JpsiK_e_gentaudm;
  float JpsiK_e_q2_gen;
  int JpsiK_e_nBc;
  float JpsiK_e_B_pt_gen;
  float JpsiK_e_B_eta_gen;
  float JpsiK_e_B_phi_gen;
  float JpsiK_e_B_mass_gen;




//  std::vector<int  >                BsTauTau_nCandidates ;
//
//  std::vector<float>                BsTauTau_mu1_pt      ;
//  std::vector<float>                BsTauTau_mu1_eta     ;
//  std::vector<float>                BsTauTau_mu1_phi     ;
//  std::vector<float>                BsTauTau_mu1_mass     ;
//  std::vector<float>                BsTauTau_mu1_unfit_pt      ;
//  std::vector<float>                BsTauTau_mu1_unfit_eta     ;
//  std::vector<float>                BsTauTau_mu1_unfit_phi     ;
//  std::vector<float>                BsTauTau_mu1_unfit_mass     ;
//  std::vector<int  >                BsTauTau_mu1_q      ;   
//  std::vector<int  >                BsTauTau_mu1_isLoose   ;
//  std::vector<int  >                BsTauTau_mu1_isTight   ;
//  std::vector<int  >                BsTauTau_mu1_isPF      ;
//  std::vector<int  >                BsTauTau_mu1_isGlobal  ;
//  std::vector<int  >                BsTauTau_mu1_isTracker ;
//  std::vector<int  >                BsTauTau_mu1_isSoft    ;
//  std::vector<float>                BsTauTau_mu1_vx       ;
//  std::vector<float>                BsTauTau_mu1_vy       ;
//  std::vector<float>                BsTauTau_mu1_vz       ;
//  std::vector<float>                BsTauTau_mu1_iso       ;
//  std::vector<float>                BsTauTau_mu1_dbiso       ;
//
//  std::vector<float>                BsTauTau_tau_pt      ;
//  std::vector<float>                BsTauTau_tau_eta     ;
//  std::vector<float>                BsTauTau_tau_phi     ;
//  std::vector<float>                BsTauTau_tau_mass     ;
//  std::vector<float>                BsTauTau_tau_rhomass1     ;
//  std::vector<float>                BsTauTau_tau_rhomass2     ;
//  std::vector<int  >                BsTauTau_tau_q   ;   
//  std::vector<float>                BsTauTau_tau_vx       ;
//  std::vector<float>                BsTauTau_tau_vy       ;
//  std::vector<float>                BsTauTau_tau_vz       ;
//
//  std::vector<float>       BsTauTau_tau_max_dr_3prong;
//  std::vector<float>       BsTauTau_tau_lip;
//  std::vector<float>       BsTauTau_tau_lips;
//  std::vector<float>       BsTauTau_tau_pvip;
//  std::vector<float>       BsTauTau_tau_pvips;
//  std::vector<float>       BsTauTau_tau_fl3d;
//  std::vector<float>       BsTauTau_tau_fls3d;
//  std::vector<float>       BsTauTau_tau_alpha;
//  std::vector<float>       BsTauTau_tau_vprob;
//  std::vector<bool>        BsTauTau_tau_isRight;
//  std::vector<bool>        BsTauTau_tau_isRight1;
//  std::vector<bool>        BsTauTau_tau_isRight2;
//  std::vector<bool>        BsTauTau_tau_isRight3;
////  std::vector<float>        BsTauTau_tau_dr1;
////  std::vector<float>        BsTauTau_tau_dr2;
////  std::vector<float>        BsTauTau_tau_dr3;
////  std::vector<float>        BsTauTau_tau_ptres1;
////  std::vector<float>        BsTauTau_tau_ptres2;
////  std::vector<float>        BsTauTau_tau_ptres3;
//  std::vector<int>         BsTauTau_tau_matched_ppdgId;
//  std::vector<float>       BsTauTau_tau_matched_gentaupt;
//  std::vector<float>       BsTauTau_tau_sumofdnn; 
//  std::vector<int>       BsTauTau_tau_pfidx1;
//  std::vector<int>       BsTauTau_tau_pfidx2;
//  std::vector<int>       BsTauTau_tau_pfidx3;
//  std::vector<float>       BsTauTau_tau_pi1_dnn;
//  std::vector<float>       BsTauTau_tau_pi2_dnn;
//  std::vector<float>       BsTauTau_tau_pi3_dnn;
//
//
//  std::vector<float>       BsTauTau_tau_pi1_pt;
//  std::vector<float>       BsTauTau_tau_pi1_eta;
//  std::vector<float>       BsTauTau_tau_pi1_phi;
//  std::vector<float>       BsTauTau_tau_pi1_mass;
//  std::vector<float>       BsTauTau_tau_pi2_pt;
//  std::vector<float>       BsTauTau_tau_pi2_eta;
//  std::vector<float>       BsTauTau_tau_pi2_phi;
//  std::vector<float>       BsTauTau_tau_pi2_mass;
//  std::vector<float>       BsTauTau_tau_pi3_pt;
//  std::vector<float>       BsTauTau_tau_pi3_eta;
//  std::vector<float>       BsTauTau_tau_pi3_phi;
//  std::vector<float>       BsTauTau_tau_pi3_mass;
//
//
//  std::vector<float>                BsTauTau_B_pt      ;
//  std::vector<float>                BsTauTau_B_eta     ;
//  std::vector<float>                BsTauTau_B_phi     ;
//  std::vector<float>                BsTauTau_B_mass    ;
//  std::vector<float>                BsTauTau_B_vprob ;
//  std::vector<float>                BsTauTau_B_lip;
//  std::vector<float>                BsTauTau_B_lips;
//  std::vector<float>                BsTauTau_B_pvip;
//  std::vector<float>                BsTauTau_B_pvips;
//  std::vector<float>                BsTauTau_B_fl3d;
//  std::vector<float>                BsTauTau_B_fls3d;
//  std::vector<float>                BsTauTau_B_alpha;
//  std::vector<float>                BsTauTau_B_maxdoca;
//  std::vector<float>                BsTauTau_B_mindoca;
//  std::vector<float>                BsTauTau_B_vx      ;
//  std::vector<float>                BsTauTau_B_vy      ;
//  std::vector<float>                BsTauTau_B_vz      ;
//  std::vector<float>                BsTauTau_B_iso;
//  std::vector<int  >                BsTauTau_B_iso_ntracks;
//  std::vector<float>                BsTauTau_B_iso_mindoca;
//  std::vector<float>                BsTauTau_B_unfit_pt      ;
//  std::vector<float>                BsTauTau_B_unfit_mass       ;
//  std::vector<float>                BsTauTau_B_unfit_vprob    ;
//  std::vector<float>                BsTauTau_B_unfit_vx;
//  std::vector<float>                BsTauTau_B_unfit_vy;
//  std::vector<float>                BsTauTau_B_unfit_vz;
//
//  std::vector<float>                BsTauTau_PV_vx       ;
//  std::vector<float>                BsTauTau_PV_vy       ;
//  std::vector<float>                BsTauTau_PV_vz       ;
//
//  std::vector<float>                BsTauTau_bbPV_vx       ;
//  std::vector<float>                BsTauTau_bbPV_vy       ;
//  std::vector<float>                BsTauTau_bbPV_vz       ;
//
//  std::vector<float>                BsTauTau_bbPV_refit_vx       ;
//  std::vector<float>                BsTauTau_bbPV_refit_vy       ;
//  std::vector<float>                BsTauTau_bbPV_refit_vz       ;
//
//  std::vector<float>                BsTauTau_genPV_vx       ;
//  std::vector<float>                BsTauTau_genPV_vy       ;
//  std::vector<float>                BsTauTau_genPV_vz       ;
//
//  std::vector<int  >                BsTauTau_ngenmuons      ;
//  std::vector<int  >                BsTauTau_isgen3;
//  std::vector<int  >                BsTauTau_isgen3matched;
//  std::vector<int> BsTauTau_nch;
//  std::vector<int> BsTauTau_nch_after_dnn;
//  std::vector<int> BsTauTau_nch_before_dnn;
//  std::vector<int> BsTauTau_nch_qr;
//  std::vector<int> BsTauTau_ngentau3;
//  std::vector<int> BsTauTau_ngentau;
//  std::vector<float> BsTauTau_gentaupt;
//  std::vector<int> BsTauTau_gentaudm;
//
//
//  //////////////////////
//
//
//
//  std::vector<int  >                BsTauTauFH_nCandidates ;
//  std::vector<int  >                BsTauTauFH_ntaus ;
//
//  std::vector<float>                BsTauTauFH_mu1_pt      ;
//  std::vector<float>                BsTauTauFH_mu1_eta     ;
//  std::vector<float>                BsTauTauFH_mu1_phi     ;
//  std::vector<float>                BsTauTauFH_mu1_mass     ;
//  std::vector<float>                BsTauTauFH_mu1_unfit_pt      ;
//  std::vector<float>                BsTauTauFH_mu1_unfit_eta     ;
//  std::vector<float>                BsTauTauFH_mu1_unfit_phi     ;
//  std::vector<float>                BsTauTauFH_mu1_unfit_mass     ;
//  std::vector<int  >                BsTauTauFH_mu1_q      ;   
//  std::vector<int  >                BsTauTauFH_mu1_isLoose   ;
//  std::vector<int  >                BsTauTauFH_mu1_isTight   ;
//  std::vector<int  >                BsTauTauFH_mu1_isPF      ;
//  std::vector<int  >                BsTauTauFH_mu1_isGlobal  ;
//  std::vector<int  >                BsTauTauFH_mu1_isTracker ;
//  std::vector<int  >                BsTauTauFH_mu1_isSoft    ;
//  std::vector<float>                BsTauTauFH_mu1_vx       ;
//  std::vector<float>                BsTauTauFH_mu1_vy       ;
//  std::vector<float>                BsTauTauFH_mu1_vz       ;
//  std::vector<float>                BsTauTauFH_mu1_iso       ;
//  std::vector<float>                BsTauTauFH_mu1_dbiso       ;
//
//  std::vector<float>                BsTauTauFH_tau1_pt      ;
//  std::vector<float>                BsTauTauFH_tau1_eta     ;
//  std::vector<float>                BsTauTauFH_tau1_phi     ;
//  std::vector<float>                BsTauTauFH_tau1_mass     ;
//  std::vector<float>                BsTauTauFH_tau1_rhomass1     ;
//  std::vector<float>                BsTauTauFH_tau1_rhomass2     ;
//  std::vector<int  >                BsTauTauFH_tau1_q   ;   
//  std::vector<float>                BsTauTauFH_tau1_vx       ;
//  std::vector<float>                BsTauTauFH_tau1_vy       ;
//  std::vector<float>                BsTauTauFH_tau1_vz       ;
//
//  std::vector<float>       BsTauTauFH_tau1_max_dr_3prong;
//  std::vector<float>       BsTauTauFH_tau1_lip;
//  std::vector<float>       BsTauTauFH_tau1_lips;
//  std::vector<float>       BsTauTauFH_tau1_pvip;
//  std::vector<float>       BsTauTauFH_tau1_pvips;
//  std::vector<float>       BsTauTauFH_tau1_fl3d;
//  std::vector<float>       BsTauTauFH_tau1_fls3d;
//  std::vector<float>       BsTauTauFH_tau1_alpha;
//  std::vector<float>       BsTauTauFH_tau1_vprob;
//  std::vector<bool>        BsTauTauFH_tau1_isRight;
//  std::vector<int>         BsTauTauFH_tau1_matched_ppdgId;
//  std::vector<float>       BsTauTauFH_tau1_matched_gentaupt;
//  std::vector<float>       BsTauTauFH_tau1_sumofdnn; 
//  std::vector<int>       BsTauTauFH_tau1_pfidx1;
//  std::vector<int>       BsTauTauFH_tau1_pfidx2;
//  std::vector<int>       BsTauTauFH_tau1_pfidx3;
//  std::vector<float>       BsTauTauFH_tau1_pi1_dnn;
//  std::vector<float>       BsTauTauFH_tau1_pi2_dnn;
//  std::vector<float>       BsTauTauFH_tau1_pi3_dnn;
//
//  std::vector<float>       BsTauTauFH_tau1_pi1_pt;
//  std::vector<float>       BsTauTauFH_tau1_pi1_eta;
//  std::vector<float>       BsTauTauFH_tau1_pi1_phi;
//  std::vector<float>       BsTauTauFH_tau1_pi1_mass;
//  std::vector<float>       BsTauTauFH_tau1_pi2_pt;
//  std::vector<float>       BsTauTauFH_tau1_pi2_eta;
//  std::vector<float>       BsTauTauFH_tau1_pi2_phi;
//  std::vector<float>       BsTauTauFH_tau1_pi2_mass;
//  std::vector<float>       BsTauTauFH_tau1_pi3_pt;
//  std::vector<float>       BsTauTauFH_tau1_pi3_eta;
//  std::vector<float>       BsTauTauFH_tau1_pi3_phi;
//  std::vector<float>       BsTauTauFH_tau1_pi3_mass;
//
//  std::vector<float>                BsTauTauFH_tau2_pt      ;
//  std::vector<float>                BsTauTauFH_tau2_eta     ;
//  std::vector<float>                BsTauTauFH_tau2_phi     ;
//  std::vector<float>                BsTauTauFH_tau2_mass     ;
//  std::vector<float>                BsTauTauFH_tau2_rhomass1     ;
//  std::vector<float>                BsTauTauFH_tau2_rhomass2     ;
//  std::vector<int  >                BsTauTauFH_tau2_q   ;   
//  std::vector<float>                BsTauTauFH_tau2_vx       ;
//  std::vector<float>                BsTauTauFH_tau2_vy       ;
//  std::vector<float>                BsTauTauFH_tau2_vz       ;
//
//  std::vector<float>       BsTauTauFH_tau2_max_dr_3prong;
//  std::vector<float>       BsTauTauFH_tau2_lip;
//  std::vector<float>       BsTauTauFH_tau2_lips;
//  std::vector<float>       BsTauTauFH_tau2_pvip;
//  std::vector<float>       BsTauTauFH_tau2_pvips;
//  std::vector<float>       BsTauTauFH_tau2_fl3d;
//  std::vector<float>       BsTauTauFH_tau2_fls3d;
//  std::vector<float>       BsTauTauFH_tau2_alpha;
//  std::vector<float>       BsTauTauFH_tau2_vprob;
//  std::vector<bool>        BsTauTauFH_tau2_isRight;
//  std::vector<int>         BsTauTauFH_tau2_matched_ppdgId;
//  std::vector<float>       BsTauTauFH_tau2_matched_gentaupt;
//  std::vector<float>       BsTauTauFH_tau2_sumofdnn; 
//  std::vector<int>       BsTauTauFH_tau2_pfidx1;
//  std::vector<int>       BsTauTauFH_tau2_pfidx2;
//  std::vector<int>       BsTauTauFH_tau2_pfidx3;
//  std::vector<float>       BsTauTauFH_tau2_pi1_dnn;
//  std::vector<float>       BsTauTauFH_tau2_pi2_dnn;
//  std::vector<float>       BsTauTauFH_tau2_pi3_dnn;
//
//  std::vector<float>       BsTauTauFH_tau2_pi1_pt;
//  std::vector<float>       BsTauTauFH_tau2_pi1_eta;
//  std::vector<float>       BsTauTauFH_tau2_pi1_phi;
//  std::vector<float>       BsTauTauFH_tau2_pi1_mass;
//  std::vector<float>       BsTauTauFH_tau2_pi2_pt;
//  std::vector<float>       BsTauTauFH_tau2_pi2_eta;
//  std::vector<float>       BsTauTauFH_tau2_pi2_phi;
//  std::vector<float>       BsTauTauFH_tau2_pi2_mass;
//  std::vector<float>       BsTauTauFH_tau2_pi3_pt;
//  std::vector<float>       BsTauTauFH_tau2_pi3_eta;
//  std::vector<float>       BsTauTauFH_tau2_pi3_phi;
//  std::vector<float>       BsTauTauFH_tau2_pi3_mass;
//
//  std::vector<float>                BsTauTauFH_B_pt      ;
//  std::vector<float>                BsTauTauFH_B_eta     ;
//  std::vector<float>                BsTauTauFH_B_phi     ;
//  std::vector<float>                BsTauTauFH_B_mass    ;
//  std::vector<float>                BsTauTauFH_B_vprob ;
//  std::vector<float>                BsTauTauFH_B_lip;
//  std::vector<float>                BsTauTauFH_B_lips;
//  std::vector<float>                BsTauTauFH_B_pvip;
//  std::vector<float>                BsTauTauFH_B_pvips;
//  std::vector<float>                BsTauTauFH_B_fl3d;
//  std::vector<float>                BsTauTauFH_B_fls3d;
//  std::vector<float>                BsTauTauFH_B_alpha;
//  std::vector<float>                BsTauTauFH_B_maxdoca;
//  std::vector<float>                BsTauTauFH_B_mindoca;
//  std::vector<float>                BsTauTauFH_B_vx      ;
//  std::vector<float>                BsTauTauFH_B_vy      ;
//  std::vector<float>                BsTauTauFH_B_vz      ;
//  std::vector<float>                BsTauTauFH_B_iso;
//  std::vector<int  >                BsTauTauFH_B_iso_ntracks;
//  std::vector<float>                BsTauTauFH_B_iso_mindoca;
//  std::vector<float>                BsTauTauFH_B_unfit_pt      ;
//  std::vector<float>                BsTauTauFH_B_unfit_mass       ;
//  std::vector<float>                BsTauTauFH_B_unfit_vprob    ;
//  std::vector<float>                BsTauTauFH_B_unfit_vx;
//  std::vector<float>                BsTauTauFH_B_unfit_vy;
//  std::vector<float>                BsTauTauFH_B_unfit_vz;
//
//  std::vector<float>                BsTauTauFH_PV_vx       ;
//  std::vector<float>                BsTauTauFH_PV_vy       ;
//  std::vector<float>                BsTauTauFH_PV_vz       ;
//
//  std::vector<float>                BsTauTauFH_bbPV_vx       ;
//  std::vector<float>                BsTauTauFH_bbPV_vy       ;
//  std::vector<float>                BsTauTauFH_bbPV_vz       ;
//
//  std::vector<float>                BsTauTauFH_bbPV_refit_vx       ;
//  std::vector<float>                BsTauTauFH_bbPV_refit_vy       ;
//  std::vector<float>                BsTauTauFH_bbPV_refit_vz       ;
//
//  std::vector<float>                BsTauTauFH_genPV_vx       ;
//  std::vector<float>                BsTauTauFH_genPV_vy       ;
//  std::vector<float>                BsTauTauFH_genPV_vz       ;
//
//  std::vector<int  >                BsTauTauFH_ngenmuons      ;
//  std::vector<int  >                BsTauTauFH_isgen3;
//  std::vector<int  >                BsTauTauFH_isgen3matched;
//  std::vector<int> BsTauTauFH_nch;
//  std::vector<int> BsTauTauFH_nch_after_dnn;
//  std::vector<int> BsTauTauFH_nch_before_dnn;
//  std::vector<int> BsTauTauFH_nch_qr;
//  std::vector<int> BsTauTauFH_ngentau3;
//  std::vector<int> BsTauTauFH_ngentau;
//  std::vector<float> BsTauTauFH_gentaupt;
//  std::vector<int> BsTauTauFH_gentaudm;
//
//
//  //////////////////////
//
//  std::vector<float>       BsTauTauFH_mr_tau_pi1_pt;
//  std::vector<float>       BsTauTauFH_mr_tau_pi1_eta;
//  std::vector<float>       BsTauTauFH_mr_tau_pi1_phi;
//  std::vector<float>       BsTauTauFH_mr_tau_pi1_mass;
//  std::vector<float>       BsTauTauFH_mr_tau_pi2_pt;
//  std::vector<float>       BsTauTauFH_mr_tau_pi2_eta;
//  std::vector<float>       BsTauTauFH_mr_tau_pi2_phi;
//  std::vector<float>       BsTauTauFH_mr_tau_pi2_mass;
//  std::vector<float>       BsTauTauFH_mr_tau_pi3_pt;
//  std::vector<float>       BsTauTauFH_mr_tau_pi3_eta;
//  std::vector<float>       BsTauTauFH_mr_tau_pi3_phi;
//  std::vector<float>       BsTauTauFH_mr_tau_pi3_mass;
//  std::vector<float>       BsTauTauFH_mr_tau_genpt;
//  std::vector<float>       BsTauTauFH_mr_tau_geneta;
//  std::vector<float>       BsTauTauFH_mr_tau_genphi;
//  std::vector<float>       BsTauTauFH_mr_tau_genmass;
//  std::vector<float>       BsTauTauFH_mr_tau_genpt_bd;
//  std::vector<float>       BsTauTauFH_mr_tau_geneta_bd;
//  std::vector<float>       BsTauTauFH_mr_tau_genphi_bd;
//  std::vector<float>       BsTauTauFH_mr_tau_genmass_bd;
//
//
//
//  std::vector<int  >                BsDstarTauNu_nCandidates ;
//
//  std::vector<float>                BsDstarTauNu_mu1_pt      ;
//  std::vector<float>                BsDstarTauNu_mu1_eta     ;
//  std::vector<float>                BsDstarTauNu_mu1_phi     ;
//  std::vector<float>                BsDstarTauNu_mu1_mass     ;
//  std::vector<float>                BsDstarTauNu_mu1_unfit_pt      ;
//  std::vector<float>                BsDstarTauNu_mu1_unfit_eta     ;
//  std::vector<float>                BsDstarTauNu_mu1_unfit_phi     ;
//  std::vector<float>                BsDstarTauNu_mu1_unfit_mass     ;
//  std::vector<int  >                BsDstarTauNu_mu1_q      ;   
//  std::vector<int  >                BsDstarTauNu_mu1_isLoose   ;
//  std::vector<int  >                BsDstarTauNu_mu1_isTight   ;
//  std::vector<int  >                BsDstarTauNu_mu1_isPF      ;
//  std::vector<int  >                BsDstarTauNu_mu1_isGlobal  ;
//  std::vector<int  >                BsDstarTauNu_mu1_isTracker ;
//  std::vector<int  >                BsDstarTauNu_mu1_isSoft    ;
//  std::vector<float>                BsDstarTauNu_mu1_vx       ;
//  std::vector<float>                BsDstarTauNu_mu1_vy       ;
//  std::vector<float>                BsDstarTauNu_mu1_vz       ;
//  std::vector<float>                BsDstarTauNu_mu1_iso       ;
//  std::vector<float>                BsDstarTauNu_mu1_dbiso       ;
//
//  std::vector<float>                BsDstarTauNu_tau_fullfit_pt      ;
//  std::vector<float>                BsDstarTauNu_tau_fullfit_eta     ;
//  std::vector<float>                BsDstarTauNu_tau_fullfit_phi     ;
//  std::vector<float>                BsDstarTauNu_tau_fullfit_mass     ;
//  std::vector<float>                BsDstarTauNu_tau_pt      ;
//  std::vector<float>                BsDstarTauNu_tau_eta     ;
//  std::vector<float>                BsDstarTauNu_tau_phi     ;
//  std::vector<float>                BsDstarTauNu_tau_mass     ;
//  std::vector<float>                BsDstarTauNu_tau_rhomass1    ;
//  std::vector<float>                BsDstarTauNu_tau_rhomass2     ;
//  std::vector<int  >                BsDstarTauNu_tau_q   ;   
//  std::vector<float>                BsDstarTauNu_tau_vx       ;
//  std::vector<float>                BsDstarTauNu_tau_vy       ;
//  std::vector<float>                BsDstarTauNu_tau_vz       ;
//
//  std::vector<float>       BsDstarTauNu_tau_pi1_pt;
//  std::vector<float>       BsDstarTauNu_tau_pi1_eta;
//  std::vector<float>       BsDstarTauNu_tau_pi1_phi;
//  std::vector<float>       BsDstarTauNu_tau_pi1_mass;
//  std::vector<float>       BsDstarTauNu_tau_pi2_pt;
//  std::vector<float>       BsDstarTauNu_tau_pi2_eta;
//  std::vector<float>       BsDstarTauNu_tau_pi2_phi;
//  std::vector<float>       BsDstarTauNu_tau_pi2_mass;
//  std::vector<float>       BsDstarTauNu_tau_pi3_pt;
//  std::vector<float>       BsDstarTauNu_tau_pi3_eta;
//  std::vector<float>       BsDstarTauNu_tau_pi3_phi;
//  std::vector<float>       BsDstarTauNu_tau_pi3_mass;
//
//  std::vector<float>       BsDstarTauNu_tau_max_dr_3prong;
//  std::vector<float>       BsDstarTauNu_tau_lip;
//  std::vector<float>       BsDstarTauNu_tau_lips;
//  std::vector<float>       BsDstarTauNu_tau_pvip;
//  std::vector<float>       BsDstarTauNu_tau_pvips;
//  std::vector<float>       BsDstarTauNu_tau_fl3d;
//  std::vector<float>       BsDstarTauNu_tau_fls3d;
//  std::vector<float>       BsDstarTauNu_tau_alpha;
//  std::vector<float>       BsDstarTauNu_tau_vprob;
//  std::vector<bool>        BsDstarTauNu_tau_isRight;
//  std::vector<int>         BsDstarTauNu_tau_matched_ppdgId;
//  std::vector<float>       BsDstarTauNu_tau_matched_gentaupt;
//  std::vector<float>       BsDstarTauNu_tau_sumofdnn; 
//  std::vector<int>       BsDstarTauNu_tau_pfidx1;
//  std::vector<int>       BsDstarTauNu_tau_pfidx2;
//  std::vector<int>       BsDstarTauNu_tau_pfidx3;
//
//
//  std::vector<float>                BsDstarTauNu_B_pt      ;
//  std::vector<float>                BsDstarTauNu_B_eta     ;
//  std::vector<float>                BsDstarTauNu_B_phi     ;
//  std::vector<float>                BsDstarTauNu_B_mass    ;
//  std::vector<float>                BsDstarTauNu_B_vprob ;
//  std::vector<float>                BsDstarTauNu_B_lip;
//  std::vector<float>                BsDstarTauNu_B_lips;
//  std::vector<float>                BsDstarTauNu_B_pvip;
//  std::vector<float>                BsDstarTauNu_B_pvips;
//  std::vector<float>                BsDstarTauNu_B_fl3d;
//  std::vector<float>                BsDstarTauNu_B_fls3d;
//  std::vector<float>                BsDstarTauNu_B_alpha;
//  std::vector<float>                BsDstarTauNu_B_maxdoca;
//  std::vector<float>                BsDstarTauNu_B_mindoca;
//  std::vector<float>                BsDstarTauNu_B_vx      ;
//  std::vector<float>                BsDstarTauNu_B_vy      ;
//  std::vector<float>                BsDstarTauNu_B_vz      ;
//  std::vector<float>                BsDstarTauNu_B_iso;
//  std::vector<int  >                BsDstarTauNu_B_iso_ntracks;
//  std::vector<float>                BsDstarTauNu_B_iso_mindoca;
//  std::vector<float>                BsDstarTauNu_B_unfit_pt      ;
//  std::vector<float>                BsDstarTauNu_B_unfit_mass       ;
//  std::vector<float>                BsDstarTauNu_B_unfit_vprob    ;
//  std::vector<float>                BsDstarTauNu_B_unfit_vx;
//  std::vector<float>                BsDstarTauNu_B_unfit_vy;
//  std::vector<float>                BsDstarTauNu_B_unfit_vz;
//  std::vector<float>                BsDstarTauNu_B_mm2;
//  std::vector<float>                BsDstarTauNu_B_q2;
//  std::vector<float>                BsDstarTauNu_B_Es;
//  std::vector<float>                BsDstarTauNu_B_ptback;
//
//  std::vector<float>                BsDstarTauNu_PV_vx       ;
//  std::vector<float>                BsDstarTauNu_PV_vy       ;
//  std::vector<float>                BsDstarTauNu_PV_vz       ;
//
//  std::vector<float>                BsDstarTauNu_bbPV_vx       ;
//  std::vector<float>                BsDstarTauNu_bbPV_vy       ;
//  std::vector<float>                BsDstarTauNu_bbPV_vz       ;
//
//  std::vector<float>                BsDstarTauNu_bbPV_refit_vx       ;
//  std::vector<float>                BsDstarTauNu_bbPV_refit_vy       ;
//  std::vector<float>                BsDstarTauNu_bbPV_refit_vz       ;
//
//
//  std::vector<float>                BsDstarTauNu_Ds_pt      ;
//  std::vector<float>                BsDstarTauNu_Ds_eta     ;
//  std::vector<float>                BsDstarTauNu_Ds_phi     ;
//  std::vector<float>                BsDstarTauNu_Ds_mass       ;
//  std::vector<float>                BsDstarTauNu_Ds_vprob    ;
//  std::vector<float>                BsDstarTauNu_Ds_lip;
//  std::vector<float>                BsDstarTauNu_Ds_lips;
//  std::vector<float>                BsDstarTauNu_Ds_pvip;
//  std::vector<float>                BsDstarTauNu_Ds_pvips;
//  std::vector<float>                BsDstarTauNu_Ds_fl3d;
//  std::vector<float>                BsDstarTauNu_Ds_fls3d;
//  std::vector<float>                BsDstarTauNu_Ds_alpha;
//  std::vector<float>                BsDstarTauNu_Ds_vx      ;
//  std::vector<float>                BsDstarTauNu_Ds_vy      ;
//  std::vector<float>                BsDstarTauNu_Ds_vz      ;
//  std::vector<float>                BsDstarTauNu_Ds_unfit_pt      ;
//  std::vector<float>                BsDstarTauNu_Ds_unfit_mass       ;
//  std::vector<float>                BsDstarTauNu_Ds_ptfrac       ;
//
//  std::vector<float>                BsDstarTauNu_k_charge;
//  std::vector<float>                BsDstarTauNu_pi_charge;
//  std::vector<float>                BsDstarTauNu_spi_charge;
//
//  std::vector<float>                BsDstarTauNu_D0_pt      ;
//  std::vector<float>                BsDstarTauNu_D0_eta     ;
//  std::vector<float>                BsDstarTauNu_D0_phi     ;
//  std::vector<float>                BsDstarTauNu_D0_mass       ;
//  std::vector<float>                BsDstarTauNu_D0_vprob    ;
//  std::vector<float>                BsDstarTauNu_D0_lip;
//  std::vector<float>                BsDstarTauNu_D0_lips;
//  std::vector<float>                BsDstarTauNu_D0_pvip;
//  std::vector<float>                BsDstarTauNu_D0_pvips;
//  std::vector<float>                BsDstarTauNu_D0_fl3d;
//  std::vector<float>                BsDstarTauNu_D0_fls3d;
//  std::vector<float>                BsDstarTauNu_D0_alpha;
//  std::vector<float>                BsDstarTauNu_D0_vx      ;
//  std::vector<float>                BsDstarTauNu_D0_vy      ;
//  std::vector<float>                BsDstarTauNu_D0_vz      ;
//  std::vector<float>                BsDstarTauNu_D0_unfit_pt      ;
//  std::vector<float>                BsDstarTauNu_D0_unfit_mass       ;
//  std::vector<float>                BsDstarTauNu_D0_ptfrac       ;
//
//  std::vector<float>                BsDstarTauNu_genPV_vx       ;
//  std::vector<float>                BsDstarTauNu_genPV_vy       ;
//  std::vector<float>                BsDstarTauNu_genPV_vz       ;
//
//  std::vector<int  >                BsDstarTauNu_ngenmuons      ;
//  std::vector<bool  >                BsDstarTauNu_isgen3;
//  std::vector<bool  >                BsDstarTauNu_isgen3matched;
//  std::vector<int> BsDstarTauNu_nch;
//  std::vector<int> BsDstarTauNu_nch_after_dnn;
//  std::vector<int> BsDstarTauNu_nch_before_dnn;
//  std::vector<int> BsDstarTauNu_nch_qr;
//  std::vector<int> BsDstarTauNu_ngentau3;
//  std::vector<int> BsDstarTauNu_ngentau;
//  std::vector<float> BsDstarTauNu_gentaupt;
//  std::vector<int> BsDstarTauNu_gentaudm;





  //////////////////////////




  /** HLT trigger objects */
  /* std::vector<float>  		    triggerObject_pt	      ; */
  /* std::vector<float>  		    triggerObject_eta	      ; */
  /* std::vector<float>  		    triggerObject_phi	      ; */
  /* std::vector<float>  		    triggerObject_mass	      ; */
  /* std::vector<std::string>  		    triggerObject_lastname    ; */
  /* std::vector< std::vector<float> > triggerObject_filterIDs   ; // as defined in http://cmslxr.fnal.gov/lxr/source/DataFormats/HLTReco/interface/TriggerTypeDefs.h */
  /* //  std::vector< std::vector<std::string> > triggerObject_filterLabels; */
  /* std::map<std::string, std::vector<std::string> > triggerObject_filterLabels; */
  /* std::vector< std::vector<int> >   triggerObject_firedTrigger; // as defined in plugins/TriggersNtuplizer.cc */


  /*--------------------------PV infos--------------------------*/
  int                               PV_N		     ;
  bool                              PV_filter		 ;
  std::vector<float>                PV_chi2          ;
  std::vector<float>                PV_ndof          ;
  std::vector<float>                PV_rho           ;
  std::vector<float>                PV_z             ;
  std::vector<float>                BeamSpot_x0;
  std::vector<float>                BeamSpot_y0;
  std::vector<float>                BeamSpot_z0;
  /*--------------------------PU infos--------------------------*/  			       
  std::vector<float  >                nPuVtxTrue             ;// the *true* mean number of the poisson distribution for this event from which the number of interactions each bunch crossing has been sampled // In MC it can be float 
  std::vector<int  >                nPuVtx                 ;// the number of pileup interactions that have been added to the event in the current bunch crossing
  std::vector<int  >                bX                     ;// to which bunch crossing do these interaction belong?  
  
private:
  TTree* tree_;

};

#endif // NtupleBranches_H
