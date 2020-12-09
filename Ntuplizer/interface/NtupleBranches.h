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
  TH1F* cutflow_perevt ;
  TH1F* q2_nocut;
  TH1F* hammer_width;
  
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
  std::vector<int  >                IsJpsiTau   ;
  std::vector<int  >                IsBsTauTau   ;
  std::vector<int  >                IsBsTauTauFH   ;
  std::vector<int  >                IsBsDstarTauNu   ;

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

  std::vector<float> JpsiMu_hammer_ebe_a0_up;
  std::vector<float> JpsiMu_hammer_ebe_a0_down;
  std::vector<float> JpsiMu_hammer_ebe_a1_up;
  std::vector<float> JpsiMu_hammer_ebe_a1_down;
  std::vector<float> JpsiMu_hammer_ebe_a2_up;
  std::vector<float> JpsiMu_hammer_ebe_a2_down;

  std::vector<float> JpsiMu_hammer_ebe_b0_up;
  std::vector<float> JpsiMu_hammer_ebe_b0_down;
  std::vector<float> JpsiMu_hammer_ebe_b1_up;
  std::vector<float> JpsiMu_hammer_ebe_b1_down;
  std::vector<float> JpsiMu_hammer_ebe_b2_up;
  std::vector<float> JpsiMu_hammer_ebe_b2_down;

  std::vector<float> JpsiMu_hammer_ebe_c1_up;
  std::vector<float> JpsiMu_hammer_ebe_c1_down;
  std::vector<float> JpsiMu_hammer_ebe_c2_up;
  std::vector<float> JpsiMu_hammer_ebe_c2_down;

  std::vector<float> JpsiMu_hammer_ebe_d0_up;
  std::vector<float> JpsiMu_hammer_ebe_d0_down;
  std::vector<float> JpsiMu_hammer_ebe_d1_up;
  std::vector<float> JpsiMu_hammer_ebe_d1_down;
  std::vector<float> JpsiMu_hammer_ebe_d2_up;
  std::vector<float> JpsiMu_hammer_ebe_d2_down;



  /** HLT trigger decisions for Jpsi */
  std::map<std::string,bool> HLT_BPH_isFired;

  ////////////////////////////////////////

  std::vector<int  >                JpsiTau_nCandidates ;

  std::vector<float>                JpsiTau_mu1_pt      ;
  std::vector<float>                JpsiTau_mu1_eta     ;
  std::vector<float>                JpsiTau_mu1_phi     ;
  std::vector<float>                JpsiTau_mu1_mass     ;
  std::vector<int  >                JpsiTau_mu1_q      ;   
  std::vector<int  >                JpsiTau_mu1_isLoose   ;
  std::vector<int  >                JpsiTau_mu1_isTight   ;
  std::vector<int  >                JpsiTau_mu1_isPF      ;
  std::vector<int  >                JpsiTau_mu1_isGlobal  ;
  std::vector<int  >                JpsiTau_mu1_isTracker ;
  std::vector<int  >                JpsiTau_mu1_isSoft    ;
  std::vector<float>                JpsiTau_mu1_vx       ;
  std::vector<float>                JpsiTau_mu1_vy       ;
  std::vector<float>                JpsiTau_mu1_vz       ;
  std::vector<float>                JpsiTau_mu1_iso       ;
  std::vector<float>                JpsiTau_mu1_dbiso       ;

  std::vector<float>                JpsiTau_mu2_pt      ;
  std::vector<float>                JpsiTau_mu2_eta     ;
  std::vector<float>                JpsiTau_mu2_phi     ;
  std::vector<float>                JpsiTau_mu2_mass     ;
  std::vector<int  >                JpsiTau_mu2_q   ;   
  std::vector<int  >                JpsiTau_mu2_isLoose   ;
  std::vector<int  >                JpsiTau_mu2_isTight   ;
  std::vector<int  >                JpsiTau_mu2_isPF      ;
  std::vector<int  >                JpsiTau_mu2_isGlobal  ;
  std::vector<int  >                JpsiTau_mu2_isTracker ;
  std::vector<int  >                JpsiTau_mu2_isSoft    ;
  std::vector<float>                JpsiTau_mu2_vx       ;
  std::vector<float>                JpsiTau_mu2_vy       ;
  std::vector<float>                JpsiTau_mu2_vz       ;
  std::vector<float>                JpsiTau_mu2_iso       ;
  std::vector<float>                JpsiTau_mu2_dbiso       ;

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
  std::vector<bool>        JpsiTau_tau_isRight;
  std::vector<bool>        JpsiTau_tau_isRight1;
  std::vector<bool>        JpsiTau_tau_isRight2;
  std::vector<bool>        JpsiTau_tau_isRight3;
//  std::vector<float>        JpsiTau_tau_dr1;
//  std::vector<float>        JpsiTau_tau_dr2;
//  std::vector<float>        JpsiTau_tau_dr3;
//  std::vector<float>        JpsiTau_tau_ptres1;
//  std::vector<float>        JpsiTau_tau_ptres2;
//  std::vector<float>        JpsiTau_tau_ptres3;
  std::vector<int>         JpsiTau_tau_matched_ppdgId;
  std::vector<float>       JpsiTau_tau_matched_gentaupt;
  std::vector<float>       JpsiTau_tau_sumofdnn; 
  std::vector<float>       JpsiTau_tau_sumofdnn_others; 
  std::vector<int>       JpsiTau_tau_pfidx1;
  std::vector<int>       JpsiTau_tau_pfidx2;
  std::vector<int>       JpsiTau_tau_pfidx3;

  std::vector<float>       JpsiTau_tau_pi1_dnn;
  std::vector<float>       JpsiTau_tau_pi2_dnn;
  std::vector<float>       JpsiTau_tau_pi3_dnn;

  std::vector<float>       JpsiTau_tau_pi1_doca;
  std::vector<float>       JpsiTau_tau_pi2_doca;
  std::vector<float>       JpsiTau_tau_pi3_doca;

  std::vector<int>       JpsiTau_tau_pi1_pv;
  std::vector<int>       JpsiTau_tau_pi2_pv;
  std::vector<int>       JpsiTau_tau_pi3_pv;


  std::vector<float>                JpsiTau_tau_rhomass1     ;
  std::vector<float>                JpsiTau_tau_rhomass2     ;

  std::vector<float>       JpsiTau_tau_pi1_pt;
  std::vector<float>       JpsiTau_tau_pi1_eta;
  std::vector<float>       JpsiTau_tau_pi1_phi;
  std::vector<float>       JpsiTau_tau_pi1_mass;
  std::vector<int>       JpsiTau_tau_pi1_q;
  std::vector<float>       JpsiTau_tau_pi2_pt;
  std::vector<float>       JpsiTau_tau_pi2_eta;
  std::vector<float>       JpsiTau_tau_pi2_phi;
  std::vector<float>       JpsiTau_tau_pi2_mass;
  std::vector<int>       JpsiTau_tau_pi2_q;
  std::vector<float>       JpsiTau_tau_pi3_pt;
  std::vector<float>       JpsiTau_tau_pi3_eta;
  std::vector<float>       JpsiTau_tau_pi3_phi;
  std::vector<float>       JpsiTau_tau_pi3_mass;
  std::vector<int>       JpsiTau_tau_pi3_q;
  std::vector<float>       JpsiTau_ptbal;
  std::vector<float>       JpsiTau_jpsi_tau_alpha;

  std::vector<float>                JpsiTau_Jpsi_pt      ;
  std::vector<float>                JpsiTau_Jpsi_eta     ;
  std::vector<float>                JpsiTau_Jpsi_phi     ;
  std::vector<float>                JpsiTau_Jpsi_mass       ;
  std::vector<float>                JpsiTau_Jpsi_vprob    ;
  std::vector<float>                JpsiTau_Jpsi_lip;
  std::vector<float>                JpsiTau_Jpsi_lips;
  std::vector<float>                JpsiTau_Jpsi_pvip;
  std::vector<float>                JpsiTau_Jpsi_pvips;
  std::vector<float>                JpsiTau_Jpsi_fl3d;
  std::vector<float>                JpsiTau_Jpsi_fls3d;
  std::vector<float>                JpsiTau_Jpsi_alpha;
  std::vector<float>                JpsiTau_Jpsi_maxdoca;
  std::vector<float>                JpsiTau_Jpsi_mindoca;
  std::vector<float>                JpsiTau_Jpsi_vx      ;
  std::vector<float>                JpsiTau_Jpsi_vy      ;
  std::vector<float>                JpsiTau_Jpsi_vz      ;
  std::vector<float>                JpsiTau_Jpsi_unfit_pt      ;
  std::vector<float>                JpsiTau_Jpsi_unfit_mass       ;
  std::vector<float>                JpsiTau_Jpsi_unfit_vprob    ;
  std::vector<float>                JpsiTau_Jpsi_unfit_vx;
  std::vector<float>                JpsiTau_Jpsi_unfit_vy;
  std::vector<float>                JpsiTau_Jpsi_unfit_vz;


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
  std::vector<float>                JpsiTau_B_iso;
  std::vector<int  >                JpsiTau_B_iso_ntracks;
  std::vector<float>                JpsiTau_B_iso_mindoca;

  std::vector<float>                JpsiTau_B_iso_nocut;
  std::vector<int  >                JpsiTau_B_iso_ntracks_nocut;
  std::vector<float>                JpsiTau_B_iso_mindoca_nocut;

  std::vector<float>                JpsiTau_B_unfit_pt      ;
  std::vector<float>                JpsiTau_B_unfit_mass       ;
  std::vector<float>                JpsiTau_B_unfit_vprob    ;
  std::vector<float>                JpsiTau_B_unfit_vx;
  std::vector<float>                JpsiTau_B_unfit_vy;
  std::vector<float>                JpsiTau_B_unfit_vz;
  std::vector<float>                JpsiTau_B_q2;
  std::vector<float>                JpsiTau_B_mm2;
  std::vector<float>                JpsiTau_B_ptmiss;
  std::vector<float>                JpsiTau_B_Es;
  std::vector<float>                JpsiTau_B_ptback;

  std::vector<float>                JpsiTau_PV_vx       ;
  std::vector<float>                JpsiTau_PV_vy       ;
  std::vector<float>                JpsiTau_PV_vz       ;

  std::vector<float>                JpsiTau_bbPV_vx       ;
  std::vector<float>                JpsiTau_bbPV_vy       ;
  std::vector<float>                JpsiTau_bbPV_vz       ;

  //  std::vector<float>                JpsiTau_bbPV_refit_vx       ;
  //  std::vector<float>                JpsiTau_bbPV_refit_vy       ;
  //  std::vector<float>                JpsiTau_bbPV_refit_vz       ;

  std::vector<float>                JpsiTau_genPV_vx       ;
  std::vector<float>                JpsiTau_genPV_vy       ;
  std::vector<float>                JpsiTau_genPV_vz       ;

  std::vector<int  >                JpsiTau_ngenmuons      ;
  std::vector<int  >                JpsiTau_isgenmatched;
  std::vector<int  >                JpsiTau_isgen3;
  std::vector<int  >                JpsiTau_isgen3matched;
  std::vector<int> JpsiTau_nch;
  std::vector<int> JpsiTau_nch_after_dnn;
  std::vector<int> JpsiTau_nch_before_dnn;
  //  std::vector<int> JpsiTau_nch_pvipcut;
  //  std::vector<int> JpsiTau_nch_pvipcut_dz0;
  //  std::vector<int> JpsiTau_nch_pvipcut_dz2p5;
  std::vector<int> JpsiTau_nch_qr;
  std::vector<int> JpsiTau_ngentau3;
  std::vector<int> JpsiTau_ngentau;
  std::vector<float> JpsiTau_gentaupt;
  std::vector<float> JpsiTau_gentaueta;
  std::vector<float> JpsiTau_gentauphi;
  std::vector<float> JpsiTau_gentaumass;
  std::vector<float> JpsiTau_gentaupt_bd;
  std::vector<float> JpsiTau_gentaueta_bd;
  std::vector<float> JpsiTau_gentauphi_bd;
  std::vector<float> JpsiTau_gentaumass_bd;
  std::vector<int> JpsiTau_gentaudm;
  std::vector<float> JpsiTau_q2_gen;
  std::vector<float> JpsiTau_B_pt_gen;
  std::vector<float> JpsiTau_B_eta_gen;
  std::vector<float> JpsiTau_B_phi_gen;
  std::vector<float> JpsiTau_B_mass_gen;

  std::vector<float> JpsiTau_st_doca3d;
  std::vector<float> JpsiTau_st_doca2d;
  std::vector<float> JpsiTau_st_doca1d;
  std::vector<float> JpsiTau_st_doca3ds;
  std::vector<float> JpsiTau_st_doca2ds;
  std::vector<float> JpsiTau_st_doca1ds;
  std::vector<float> JpsiTau_st_doca3de;
  std::vector<float> JpsiTau_st_doca2de;
  std::vector<float> JpsiTau_st_doca1de;

  std::vector<float> JpsiTau_st_doca3d_max;
  std::vector<float> JpsiTau_st_doca2d_max;
  std::vector<float> JpsiTau_st_doca1d_max;
  std::vector<float> JpsiTau_st_doca3ds_max;
  std::vector<float> JpsiTau_st_doca2ds_max;
  std::vector<float> JpsiTau_st_doca1ds_max;
  std::vector<float> JpsiTau_st_dz_max;

  std::vector<bool> JpsiTau_st_isRight;
  std::vector<bool> JpsiTau_st_isBdecay;
  std::vector<bool> JpsiTau_st_isSignal;
  std::vector<int> JpsiTau_st_isBdecaypdg;
  std::vector<int> JpsiTau_st_isBdecayppdg;
  std::vector<int> JpsiTau_st_nprong;
  std::vector<float> JpsiTau_st_dz;
  std::vector<bool> JpsiTau_st_isAssociate;
  std::vector<bool> JpsiTau_st_isBelong;
  std::vector<int> JpsiTau_st_pvAssociationQuality;
  std::vector<float> JpsiTau_st_pt;
  std::vector<float> JpsiTau_st_eta;
  std::vector<float> JpsiTau_st_phi;

//  std::vector<float> JpsiTau_ed_pfeta;
//  std::vector<float> JpsiTau_ed_pfphi;
//  std::vector<int> JpsiTau_ed_isRight;
//  std::vector<int> JpsiTau_ed_id;
//  std::vector<float> JpsiTau_ed_pfdnn;
//  std::vector<float> JpsiTau_ed_genpt;
    
  std::vector<std::vector<float>> JpsiTau_hammer_ff;
  std::vector<float> JpsiTau_hammer_ebe;
  std::vector<std::vector<float>> JpsiTau_hammer_ebe_toy;
  //  std::vector<float> JpsiTau_hammer_ebe_up;
  //  std::vector<float> JpsiTau_hammer_ebe_down;
  //  std::vector<float> JpsiTau_hammer_ebe_rate_up;
//  std::vector<float> JpsiTau_hammer_ebe_rate_down;

  std::vector<float> JpsiTau_hammer_ebe_a0_up;
  std::vector<float> JpsiTau_hammer_ebe_a0_down;
  std::vector<float> JpsiTau_hammer_ebe_a1_up;
  std::vector<float> JpsiTau_hammer_ebe_a1_down;
  std::vector<float> JpsiTau_hammer_ebe_a2_up;
  std::vector<float> JpsiTau_hammer_ebe_a2_down;

  std::vector<float> JpsiTau_hammer_ebe_b0_up;
  std::vector<float> JpsiTau_hammer_ebe_b0_down;
  std::vector<float> JpsiTau_hammer_ebe_b1_up;
  std::vector<float> JpsiTau_hammer_ebe_b1_down;
  std::vector<float> JpsiTau_hammer_ebe_b2_up;
  std::vector<float> JpsiTau_hammer_ebe_b2_down;

  std::vector<float> JpsiTau_hammer_ebe_c1_up;
  std::vector<float> JpsiTau_hammer_ebe_c1_down;
  std::vector<float> JpsiTau_hammer_ebe_c2_up;
  std::vector<float> JpsiTau_hammer_ebe_c2_down;

  std::vector<float> JpsiTau_hammer_ebe_d0_up;
  std::vector<float> JpsiTau_hammer_ebe_d0_down;
  std::vector<float> JpsiTau_hammer_ebe_d1_up;
  std::vector<float> JpsiTau_hammer_ebe_d1_down;
  std::vector<float> JpsiTau_hammer_ebe_d2_up;
  std::vector<float> JpsiTau_hammer_ebe_d2_down;

  ////////////////////////



  std::vector<int  >                BsTauTau_nCandidates ;

  std::vector<float>                BsTauTau_mu1_pt      ;
  std::vector<float>                BsTauTau_mu1_eta     ;
  std::vector<float>                BsTauTau_mu1_phi     ;
  std::vector<float>                BsTauTau_mu1_mass     ;
  std::vector<float>                BsTauTau_mu1_unfit_pt      ;
  std::vector<float>                BsTauTau_mu1_unfit_eta     ;
  std::vector<float>                BsTauTau_mu1_unfit_phi     ;
  std::vector<float>                BsTauTau_mu1_unfit_mass     ;
  std::vector<int  >                BsTauTau_mu1_q      ;   
  std::vector<int  >                BsTauTau_mu1_isLoose   ;
  std::vector<int  >                BsTauTau_mu1_isTight   ;
  std::vector<int  >                BsTauTau_mu1_isPF      ;
  std::vector<int  >                BsTauTau_mu1_isGlobal  ;
  std::vector<int  >                BsTauTau_mu1_isTracker ;
  std::vector<int  >                BsTauTau_mu1_isSoft    ;
  std::vector<float>                BsTauTau_mu1_vx       ;
  std::vector<float>                BsTauTau_mu1_vy       ;
  std::vector<float>                BsTauTau_mu1_vz       ;
  std::vector<float>                BsTauTau_mu1_iso       ;
  std::vector<float>                BsTauTau_mu1_dbiso       ;

  std::vector<float>                BsTauTau_tau_pt      ;
  std::vector<float>                BsTauTau_tau_eta     ;
  std::vector<float>                BsTauTau_tau_phi     ;
  std::vector<float>                BsTauTau_tau_mass     ;
  std::vector<float>                BsTauTau_tau_rhomass1     ;
  std::vector<float>                BsTauTau_tau_rhomass2     ;
  std::vector<int  >                BsTauTau_tau_q   ;   
  std::vector<float>                BsTauTau_tau_vx       ;
  std::vector<float>                BsTauTau_tau_vy       ;
  std::vector<float>                BsTauTau_tau_vz       ;

  std::vector<float>       BsTauTau_tau_max_dr_3prong;
  std::vector<float>       BsTauTau_tau_lip;
  std::vector<float>       BsTauTau_tau_lips;
  std::vector<float>       BsTauTau_tau_pvip;
  std::vector<float>       BsTauTau_tau_pvips;
  std::vector<float>       BsTauTau_tau_fl3d;
  std::vector<float>       BsTauTau_tau_fls3d;
  std::vector<float>       BsTauTau_tau_alpha;
  std::vector<float>       BsTauTau_tau_vprob;
  std::vector<bool>        BsTauTau_tau_isRight;
  std::vector<bool>        BsTauTau_tau_isRight1;
  std::vector<bool>        BsTauTau_tau_isRight2;
  std::vector<bool>        BsTauTau_tau_isRight3;
//  std::vector<float>        BsTauTau_tau_dr1;
//  std::vector<float>        BsTauTau_tau_dr2;
//  std::vector<float>        BsTauTau_tau_dr3;
//  std::vector<float>        BsTauTau_tau_ptres1;
//  std::vector<float>        BsTauTau_tau_ptres2;
//  std::vector<float>        BsTauTau_tau_ptres3;
  std::vector<int>         BsTauTau_tau_matched_ppdgId;
  std::vector<float>       BsTauTau_tau_matched_gentaupt;
  std::vector<float>       BsTauTau_tau_sumofdnn; 
  std::vector<int>       BsTauTau_tau_pfidx1;
  std::vector<int>       BsTauTau_tau_pfidx2;
  std::vector<int>       BsTauTau_tau_pfidx3;
  std::vector<float>       BsTauTau_tau_pi1_dnn;
  std::vector<float>       BsTauTau_tau_pi2_dnn;
  std::vector<float>       BsTauTau_tau_pi3_dnn;


  std::vector<float>       BsTauTau_tau_pi1_pt;
  std::vector<float>       BsTauTau_tau_pi1_eta;
  std::vector<float>       BsTauTau_tau_pi1_phi;
  std::vector<float>       BsTauTau_tau_pi1_mass;
  std::vector<float>       BsTauTau_tau_pi2_pt;
  std::vector<float>       BsTauTau_tau_pi2_eta;
  std::vector<float>       BsTauTau_tau_pi2_phi;
  std::vector<float>       BsTauTau_tau_pi2_mass;
  std::vector<float>       BsTauTau_tau_pi3_pt;
  std::vector<float>       BsTauTau_tau_pi3_eta;
  std::vector<float>       BsTauTau_tau_pi3_phi;
  std::vector<float>       BsTauTau_tau_pi3_mass;


  std::vector<float>                BsTauTau_B_pt      ;
  std::vector<float>                BsTauTau_B_eta     ;
  std::vector<float>                BsTauTau_B_phi     ;
  std::vector<float>                BsTauTau_B_mass    ;
  std::vector<float>                BsTauTau_B_vprob ;
  std::vector<float>                BsTauTau_B_lip;
  std::vector<float>                BsTauTau_B_lips;
  std::vector<float>                BsTauTau_B_pvip;
  std::vector<float>                BsTauTau_B_pvips;
  std::vector<float>                BsTauTau_B_fl3d;
  std::vector<float>                BsTauTau_B_fls3d;
  std::vector<float>                BsTauTau_B_alpha;
  std::vector<float>                BsTauTau_B_maxdoca;
  std::vector<float>                BsTauTau_B_mindoca;
  std::vector<float>                BsTauTau_B_vx      ;
  std::vector<float>                BsTauTau_B_vy      ;
  std::vector<float>                BsTauTau_B_vz      ;
  std::vector<float>                BsTauTau_B_iso;
  std::vector<int  >                BsTauTau_B_iso_ntracks;
  std::vector<float>                BsTauTau_B_iso_mindoca;
  std::vector<float>                BsTauTau_B_unfit_pt      ;
  std::vector<float>                BsTauTau_B_unfit_mass       ;
  std::vector<float>                BsTauTau_B_unfit_vprob    ;
  std::vector<float>                BsTauTau_B_unfit_vx;
  std::vector<float>                BsTauTau_B_unfit_vy;
  std::vector<float>                BsTauTau_B_unfit_vz;

  std::vector<float>                BsTauTau_PV_vx       ;
  std::vector<float>                BsTauTau_PV_vy       ;
  std::vector<float>                BsTauTau_PV_vz       ;

  std::vector<float>                BsTauTau_bbPV_vx       ;
  std::vector<float>                BsTauTau_bbPV_vy       ;
  std::vector<float>                BsTauTau_bbPV_vz       ;

  std::vector<float>                BsTauTau_bbPV_refit_vx       ;
  std::vector<float>                BsTauTau_bbPV_refit_vy       ;
  std::vector<float>                BsTauTau_bbPV_refit_vz       ;

  std::vector<float>                BsTauTau_genPV_vx       ;
  std::vector<float>                BsTauTau_genPV_vy       ;
  std::vector<float>                BsTauTau_genPV_vz       ;

  std::vector<int  >                BsTauTau_ngenmuons      ;
  std::vector<int  >                BsTauTau_isgen3;
  std::vector<int  >                BsTauTau_isgen3matched;
  std::vector<int> BsTauTau_nch;
  std::vector<int> BsTauTau_nch_after_dnn;
  std::vector<int> BsTauTau_nch_before_dnn;
  std::vector<int> BsTauTau_nch_qr;
  std::vector<int> BsTauTau_ngentau3;
  std::vector<int> BsTauTau_ngentau;
  std::vector<float> BsTauTau_gentaupt;
  std::vector<int> BsTauTau_gentaudm;


  //////////////////////



  std::vector<int  >                BsTauTauFH_nCandidates ;
  std::vector<int  >                BsTauTauFH_ntaus ;

  std::vector<float>                BsTauTauFH_mu1_pt      ;
  std::vector<float>                BsTauTauFH_mu1_eta     ;
  std::vector<float>                BsTauTauFH_mu1_phi     ;
  std::vector<float>                BsTauTauFH_mu1_mass     ;
  std::vector<float>                BsTauTauFH_mu1_unfit_pt      ;
  std::vector<float>                BsTauTauFH_mu1_unfit_eta     ;
  std::vector<float>                BsTauTauFH_mu1_unfit_phi     ;
  std::vector<float>                BsTauTauFH_mu1_unfit_mass     ;
  std::vector<int  >                BsTauTauFH_mu1_q      ;   
  std::vector<int  >                BsTauTauFH_mu1_isLoose   ;
  std::vector<int  >                BsTauTauFH_mu1_isTight   ;
  std::vector<int  >                BsTauTauFH_mu1_isPF      ;
  std::vector<int  >                BsTauTauFH_mu1_isGlobal  ;
  std::vector<int  >                BsTauTauFH_mu1_isTracker ;
  std::vector<int  >                BsTauTauFH_mu1_isSoft    ;
  std::vector<float>                BsTauTauFH_mu1_vx       ;
  std::vector<float>                BsTauTauFH_mu1_vy       ;
  std::vector<float>                BsTauTauFH_mu1_vz       ;
  std::vector<float>                BsTauTauFH_mu1_iso       ;
  std::vector<float>                BsTauTauFH_mu1_dbiso       ;

  std::vector<float>                BsTauTauFH_tau1_pt      ;
  std::vector<float>                BsTauTauFH_tau1_eta     ;
  std::vector<float>                BsTauTauFH_tau1_phi     ;
  std::vector<float>                BsTauTauFH_tau1_mass     ;
  std::vector<float>                BsTauTauFH_tau1_rhomass1     ;
  std::vector<float>                BsTauTauFH_tau1_rhomass2     ;
  std::vector<int  >                BsTauTauFH_tau1_q   ;   
  std::vector<float>                BsTauTauFH_tau1_vx       ;
  std::vector<float>                BsTauTauFH_tau1_vy       ;
  std::vector<float>                BsTauTauFH_tau1_vz       ;

  std::vector<float>       BsTauTauFH_tau1_max_dr_3prong;
  std::vector<float>       BsTauTauFH_tau1_lip;
  std::vector<float>       BsTauTauFH_tau1_lips;
  std::vector<float>       BsTauTauFH_tau1_pvip;
  std::vector<float>       BsTauTauFH_tau1_pvips;
  std::vector<float>       BsTauTauFH_tau1_fl3d;
  std::vector<float>       BsTauTauFH_tau1_fls3d;
  std::vector<float>       BsTauTauFH_tau1_alpha;
  std::vector<float>       BsTauTauFH_tau1_vprob;
  std::vector<bool>        BsTauTauFH_tau1_isRight;
  std::vector<int>         BsTauTauFH_tau1_matched_ppdgId;
  std::vector<float>       BsTauTauFH_tau1_matched_gentaupt;
  std::vector<float>       BsTauTauFH_tau1_sumofdnn; 
  std::vector<int>       BsTauTauFH_tau1_pfidx1;
  std::vector<int>       BsTauTauFH_tau1_pfidx2;
  std::vector<int>       BsTauTauFH_tau1_pfidx3;
  std::vector<float>       BsTauTauFH_tau1_pi1_dnn;
  std::vector<float>       BsTauTauFH_tau1_pi2_dnn;
  std::vector<float>       BsTauTauFH_tau1_pi3_dnn;

  std::vector<float>       BsTauTauFH_tau1_pi1_pt;
  std::vector<float>       BsTauTauFH_tau1_pi1_eta;
  std::vector<float>       BsTauTauFH_tau1_pi1_phi;
  std::vector<float>       BsTauTauFH_tau1_pi1_mass;
  std::vector<float>       BsTauTauFH_tau1_pi2_pt;
  std::vector<float>       BsTauTauFH_tau1_pi2_eta;
  std::vector<float>       BsTauTauFH_tau1_pi2_phi;
  std::vector<float>       BsTauTauFH_tau1_pi2_mass;
  std::vector<float>       BsTauTauFH_tau1_pi3_pt;
  std::vector<float>       BsTauTauFH_tau1_pi3_eta;
  std::vector<float>       BsTauTauFH_tau1_pi3_phi;
  std::vector<float>       BsTauTauFH_tau1_pi3_mass;

  std::vector<float>                BsTauTauFH_tau2_pt      ;
  std::vector<float>                BsTauTauFH_tau2_eta     ;
  std::vector<float>                BsTauTauFH_tau2_phi     ;
  std::vector<float>                BsTauTauFH_tau2_mass     ;
  std::vector<float>                BsTauTauFH_tau2_rhomass1     ;
  std::vector<float>                BsTauTauFH_tau2_rhomass2     ;
  std::vector<int  >                BsTauTauFH_tau2_q   ;   
  std::vector<float>                BsTauTauFH_tau2_vx       ;
  std::vector<float>                BsTauTauFH_tau2_vy       ;
  std::vector<float>                BsTauTauFH_tau2_vz       ;

  std::vector<float>       BsTauTauFH_tau2_max_dr_3prong;
  std::vector<float>       BsTauTauFH_tau2_lip;
  std::vector<float>       BsTauTauFH_tau2_lips;
  std::vector<float>       BsTauTauFH_tau2_pvip;
  std::vector<float>       BsTauTauFH_tau2_pvips;
  std::vector<float>       BsTauTauFH_tau2_fl3d;
  std::vector<float>       BsTauTauFH_tau2_fls3d;
  std::vector<float>       BsTauTauFH_tau2_alpha;
  std::vector<float>       BsTauTauFH_tau2_vprob;
  std::vector<bool>        BsTauTauFH_tau2_isRight;
  std::vector<int>         BsTauTauFH_tau2_matched_ppdgId;
  std::vector<float>       BsTauTauFH_tau2_matched_gentaupt;
  std::vector<float>       BsTauTauFH_tau2_sumofdnn; 
  std::vector<int>       BsTauTauFH_tau2_pfidx1;
  std::vector<int>       BsTauTauFH_tau2_pfidx2;
  std::vector<int>       BsTauTauFH_tau2_pfidx3;
  std::vector<float>       BsTauTauFH_tau2_pi1_dnn;
  std::vector<float>       BsTauTauFH_tau2_pi2_dnn;
  std::vector<float>       BsTauTauFH_tau2_pi3_dnn;

  std::vector<float>       BsTauTauFH_tau2_pi1_pt;
  std::vector<float>       BsTauTauFH_tau2_pi1_eta;
  std::vector<float>       BsTauTauFH_tau2_pi1_phi;
  std::vector<float>       BsTauTauFH_tau2_pi1_mass;
  std::vector<float>       BsTauTauFH_tau2_pi2_pt;
  std::vector<float>       BsTauTauFH_tau2_pi2_eta;
  std::vector<float>       BsTauTauFH_tau2_pi2_phi;
  std::vector<float>       BsTauTauFH_tau2_pi2_mass;
  std::vector<float>       BsTauTauFH_tau2_pi3_pt;
  std::vector<float>       BsTauTauFH_tau2_pi3_eta;
  std::vector<float>       BsTauTauFH_tau2_pi3_phi;
  std::vector<float>       BsTauTauFH_tau2_pi3_mass;

  std::vector<float>                BsTauTauFH_B_pt      ;
  std::vector<float>                BsTauTauFH_B_eta     ;
  std::vector<float>                BsTauTauFH_B_phi     ;
  std::vector<float>                BsTauTauFH_B_mass    ;
  std::vector<float>                BsTauTauFH_B_vprob ;
  std::vector<float>                BsTauTauFH_B_lip;
  std::vector<float>                BsTauTauFH_B_lips;
  std::vector<float>                BsTauTauFH_B_pvip;
  std::vector<float>                BsTauTauFH_B_pvips;
  std::vector<float>                BsTauTauFH_B_fl3d;
  std::vector<float>                BsTauTauFH_B_fls3d;
  std::vector<float>                BsTauTauFH_B_alpha;
  std::vector<float>                BsTauTauFH_B_maxdoca;
  std::vector<float>                BsTauTauFH_B_mindoca;
  std::vector<float>                BsTauTauFH_B_vx      ;
  std::vector<float>                BsTauTauFH_B_vy      ;
  std::vector<float>                BsTauTauFH_B_vz      ;
  std::vector<float>                BsTauTauFH_B_iso;
  std::vector<int  >                BsTauTauFH_B_iso_ntracks;
  std::vector<float>                BsTauTauFH_B_iso_mindoca;
  std::vector<float>                BsTauTauFH_B_unfit_pt      ;
  std::vector<float>                BsTauTauFH_B_unfit_mass       ;
  std::vector<float>                BsTauTauFH_B_unfit_vprob    ;
  std::vector<float>                BsTauTauFH_B_unfit_vx;
  std::vector<float>                BsTauTauFH_B_unfit_vy;
  std::vector<float>                BsTauTauFH_B_unfit_vz;

  std::vector<float>                BsTauTauFH_PV_vx       ;
  std::vector<float>                BsTauTauFH_PV_vy       ;
  std::vector<float>                BsTauTauFH_PV_vz       ;

  std::vector<float>                BsTauTauFH_bbPV_vx       ;
  std::vector<float>                BsTauTauFH_bbPV_vy       ;
  std::vector<float>                BsTauTauFH_bbPV_vz       ;

  std::vector<float>                BsTauTauFH_bbPV_refit_vx       ;
  std::vector<float>                BsTauTauFH_bbPV_refit_vy       ;
  std::vector<float>                BsTauTauFH_bbPV_refit_vz       ;

  std::vector<float>                BsTauTauFH_genPV_vx       ;
  std::vector<float>                BsTauTauFH_genPV_vy       ;
  std::vector<float>                BsTauTauFH_genPV_vz       ;

  std::vector<int  >                BsTauTauFH_ngenmuons      ;
  std::vector<int  >                BsTauTauFH_isgen3;
  std::vector<int  >                BsTauTauFH_isgen3matched;
  std::vector<int> BsTauTauFH_nch;
  std::vector<int> BsTauTauFH_nch_after_dnn;
  std::vector<int> BsTauTauFH_nch_before_dnn;
  std::vector<int> BsTauTauFH_nch_qr;
  std::vector<int> BsTauTauFH_ngentau3;
  std::vector<int> BsTauTauFH_ngentau;
  std::vector<float> BsTauTauFH_gentaupt;
  std::vector<int> BsTauTauFH_gentaudm;


  //////////////////////

  std::vector<float>       BsTauTauFH_mr_tau_pi1_pt;
  std::vector<float>       BsTauTauFH_mr_tau_pi1_eta;
  std::vector<float>       BsTauTauFH_mr_tau_pi1_phi;
  std::vector<float>       BsTauTauFH_mr_tau_pi1_mass;
  std::vector<float>       BsTauTauFH_mr_tau_pi2_pt;
  std::vector<float>       BsTauTauFH_mr_tau_pi2_eta;
  std::vector<float>       BsTauTauFH_mr_tau_pi2_phi;
  std::vector<float>       BsTauTauFH_mr_tau_pi2_mass;
  std::vector<float>       BsTauTauFH_mr_tau_pi3_pt;
  std::vector<float>       BsTauTauFH_mr_tau_pi3_eta;
  std::vector<float>       BsTauTauFH_mr_tau_pi3_phi;
  std::vector<float>       BsTauTauFH_mr_tau_pi3_mass;
  std::vector<float>       BsTauTauFH_mr_tau_genpt;
  std::vector<float>       BsTauTauFH_mr_tau_geneta;
  std::vector<float>       BsTauTauFH_mr_tau_genphi;
  std::vector<float>       BsTauTauFH_mr_tau_genmass;
  std::vector<float>       BsTauTauFH_mr_tau_genpt_bd;
  std::vector<float>       BsTauTauFH_mr_tau_geneta_bd;
  std::vector<float>       BsTauTauFH_mr_tau_genphi_bd;
  std::vector<float>       BsTauTauFH_mr_tau_genmass_bd;



  std::vector<int  >                BsDstarTauNu_nCandidates ;

  std::vector<float>                BsDstarTauNu_mu1_pt      ;
  std::vector<float>                BsDstarTauNu_mu1_eta     ;
  std::vector<float>                BsDstarTauNu_mu1_phi     ;
  std::vector<float>                BsDstarTauNu_mu1_mass     ;
  std::vector<float>                BsDstarTauNu_mu1_unfit_pt      ;
  std::vector<float>                BsDstarTauNu_mu1_unfit_eta     ;
  std::vector<float>                BsDstarTauNu_mu1_unfit_phi     ;
  std::vector<float>                BsDstarTauNu_mu1_unfit_mass     ;
  std::vector<int  >                BsDstarTauNu_mu1_q      ;   
  std::vector<int  >                BsDstarTauNu_mu1_isLoose   ;
  std::vector<int  >                BsDstarTauNu_mu1_isTight   ;
  std::vector<int  >                BsDstarTauNu_mu1_isPF      ;
  std::vector<int  >                BsDstarTauNu_mu1_isGlobal  ;
  std::vector<int  >                BsDstarTauNu_mu1_isTracker ;
  std::vector<int  >                BsDstarTauNu_mu1_isSoft    ;
  std::vector<float>                BsDstarTauNu_mu1_vx       ;
  std::vector<float>                BsDstarTauNu_mu1_vy       ;
  std::vector<float>                BsDstarTauNu_mu1_vz       ;
  std::vector<float>                BsDstarTauNu_mu1_iso       ;
  std::vector<float>                BsDstarTauNu_mu1_dbiso       ;

  std::vector<float>                BsDstarTauNu_tau_fullfit_pt      ;
  std::vector<float>                BsDstarTauNu_tau_fullfit_eta     ;
  std::vector<float>                BsDstarTauNu_tau_fullfit_phi     ;
  std::vector<float>                BsDstarTauNu_tau_fullfit_mass     ;
  std::vector<float>                BsDstarTauNu_tau_pt      ;
  std::vector<float>                BsDstarTauNu_tau_eta     ;
  std::vector<float>                BsDstarTauNu_tau_phi     ;
  std::vector<float>                BsDstarTauNu_tau_mass     ;
  std::vector<float>                BsDstarTauNu_tau_rhomass1    ;
  std::vector<float>                BsDstarTauNu_tau_rhomass2     ;
  std::vector<int  >                BsDstarTauNu_tau_q   ;   
  std::vector<float>                BsDstarTauNu_tau_vx       ;
  std::vector<float>                BsDstarTauNu_tau_vy       ;
  std::vector<float>                BsDstarTauNu_tau_vz       ;

  std::vector<float>       BsDstarTauNu_tau_pi1_pt;
  std::vector<float>       BsDstarTauNu_tau_pi1_eta;
  std::vector<float>       BsDstarTauNu_tau_pi1_phi;
  std::vector<float>       BsDstarTauNu_tau_pi1_mass;
  std::vector<float>       BsDstarTauNu_tau_pi2_pt;
  std::vector<float>       BsDstarTauNu_tau_pi2_eta;
  std::vector<float>       BsDstarTauNu_tau_pi2_phi;
  std::vector<float>       BsDstarTauNu_tau_pi2_mass;
  std::vector<float>       BsDstarTauNu_tau_pi3_pt;
  std::vector<float>       BsDstarTauNu_tau_pi3_eta;
  std::vector<float>       BsDstarTauNu_tau_pi3_phi;
  std::vector<float>       BsDstarTauNu_tau_pi3_mass;

  std::vector<float>       BsDstarTauNu_tau_max_dr_3prong;
  std::vector<float>       BsDstarTauNu_tau_lip;
  std::vector<float>       BsDstarTauNu_tau_lips;
  std::vector<float>       BsDstarTauNu_tau_pvip;
  std::vector<float>       BsDstarTauNu_tau_pvips;
  std::vector<float>       BsDstarTauNu_tau_fl3d;
  std::vector<float>       BsDstarTauNu_tau_fls3d;
  std::vector<float>       BsDstarTauNu_tau_alpha;
  std::vector<float>       BsDstarTauNu_tau_vprob;
  std::vector<bool>        BsDstarTauNu_tau_isRight;
  std::vector<int>         BsDstarTauNu_tau_matched_ppdgId;
  std::vector<float>       BsDstarTauNu_tau_matched_gentaupt;
  std::vector<float>       BsDstarTauNu_tau_sumofdnn; 
  std::vector<int>       BsDstarTauNu_tau_pfidx1;
  std::vector<int>       BsDstarTauNu_tau_pfidx2;
  std::vector<int>       BsDstarTauNu_tau_pfidx3;


  std::vector<float>                BsDstarTauNu_B_pt      ;
  std::vector<float>                BsDstarTauNu_B_eta     ;
  std::vector<float>                BsDstarTauNu_B_phi     ;
  std::vector<float>                BsDstarTauNu_B_mass    ;
  std::vector<float>                BsDstarTauNu_B_vprob ;
  std::vector<float>                BsDstarTauNu_B_lip;
  std::vector<float>                BsDstarTauNu_B_lips;
  std::vector<float>                BsDstarTauNu_B_pvip;
  std::vector<float>                BsDstarTauNu_B_pvips;
  std::vector<float>                BsDstarTauNu_B_fl3d;
  std::vector<float>                BsDstarTauNu_B_fls3d;
  std::vector<float>                BsDstarTauNu_B_alpha;
  std::vector<float>                BsDstarTauNu_B_maxdoca;
  std::vector<float>                BsDstarTauNu_B_mindoca;
  std::vector<float>                BsDstarTauNu_B_vx      ;
  std::vector<float>                BsDstarTauNu_B_vy      ;
  std::vector<float>                BsDstarTauNu_B_vz      ;
  std::vector<float>                BsDstarTauNu_B_iso;
  std::vector<int  >                BsDstarTauNu_B_iso_ntracks;
  std::vector<float>                BsDstarTauNu_B_iso_mindoca;
  std::vector<float>                BsDstarTauNu_B_unfit_pt      ;
  std::vector<float>                BsDstarTauNu_B_unfit_mass       ;
  std::vector<float>                BsDstarTauNu_B_unfit_vprob    ;
  std::vector<float>                BsDstarTauNu_B_unfit_vx;
  std::vector<float>                BsDstarTauNu_B_unfit_vy;
  std::vector<float>                BsDstarTauNu_B_unfit_vz;
  std::vector<float>                BsDstarTauNu_B_mm2;
  std::vector<float>                BsDstarTauNu_B_q2;
  std::vector<float>                BsDstarTauNu_B_Es;
  std::vector<float>                BsDstarTauNu_B_ptback;

  std::vector<float>                BsDstarTauNu_PV_vx       ;
  std::vector<float>                BsDstarTauNu_PV_vy       ;
  std::vector<float>                BsDstarTauNu_PV_vz       ;

  std::vector<float>                BsDstarTauNu_bbPV_vx       ;
  std::vector<float>                BsDstarTauNu_bbPV_vy       ;
  std::vector<float>                BsDstarTauNu_bbPV_vz       ;

  std::vector<float>                BsDstarTauNu_bbPV_refit_vx       ;
  std::vector<float>                BsDstarTauNu_bbPV_refit_vy       ;
  std::vector<float>                BsDstarTauNu_bbPV_refit_vz       ;


  std::vector<float>                BsDstarTauNu_Ds_pt      ;
  std::vector<float>                BsDstarTauNu_Ds_eta     ;
  std::vector<float>                BsDstarTauNu_Ds_phi     ;
  std::vector<float>                BsDstarTauNu_Ds_mass       ;
  std::vector<float>                BsDstarTauNu_Ds_vprob    ;
  std::vector<float>                BsDstarTauNu_Ds_lip;
  std::vector<float>                BsDstarTauNu_Ds_lips;
  std::vector<float>                BsDstarTauNu_Ds_pvip;
  std::vector<float>                BsDstarTauNu_Ds_pvips;
  std::vector<float>                BsDstarTauNu_Ds_fl3d;
  std::vector<float>                BsDstarTauNu_Ds_fls3d;
  std::vector<float>                BsDstarTauNu_Ds_alpha;
  std::vector<float>                BsDstarTauNu_Ds_vx      ;
  std::vector<float>                BsDstarTauNu_Ds_vy      ;
  std::vector<float>                BsDstarTauNu_Ds_vz      ;
  std::vector<float>                BsDstarTauNu_Ds_unfit_pt      ;
  std::vector<float>                BsDstarTauNu_Ds_unfit_mass       ;
  std::vector<float>                BsDstarTauNu_Ds_ptfrac       ;

  std::vector<float>                BsDstarTauNu_k_charge;
  std::vector<float>                BsDstarTauNu_pi_charge;
  std::vector<float>                BsDstarTauNu_spi_charge;

  std::vector<float>                BsDstarTauNu_D0_pt      ;
  std::vector<float>                BsDstarTauNu_D0_eta     ;
  std::vector<float>                BsDstarTauNu_D0_phi     ;
  std::vector<float>                BsDstarTauNu_D0_mass       ;
  std::vector<float>                BsDstarTauNu_D0_vprob    ;
  std::vector<float>                BsDstarTauNu_D0_lip;
  std::vector<float>                BsDstarTauNu_D0_lips;
  std::vector<float>                BsDstarTauNu_D0_pvip;
  std::vector<float>                BsDstarTauNu_D0_pvips;
  std::vector<float>                BsDstarTauNu_D0_fl3d;
  std::vector<float>                BsDstarTauNu_D0_fls3d;
  std::vector<float>                BsDstarTauNu_D0_alpha;
  std::vector<float>                BsDstarTauNu_D0_vx      ;
  std::vector<float>                BsDstarTauNu_D0_vy      ;
  std::vector<float>                BsDstarTauNu_D0_vz      ;
  std::vector<float>                BsDstarTauNu_D0_unfit_pt      ;
  std::vector<float>                BsDstarTauNu_D0_unfit_mass       ;
  std::vector<float>                BsDstarTauNu_D0_ptfrac       ;

  std::vector<float>                BsDstarTauNu_genPV_vx       ;
  std::vector<float>                BsDstarTauNu_genPV_vy       ;
  std::vector<float>                BsDstarTauNu_genPV_vz       ;

  std::vector<int  >                BsDstarTauNu_ngenmuons      ;
  std::vector<bool  >                BsDstarTauNu_isgen3;
  std::vector<bool  >                BsDstarTauNu_isgen3matched;
  std::vector<int> BsDstarTauNu_nch;
  std::vector<int> BsDstarTauNu_nch_after_dnn;
  std::vector<int> BsDstarTauNu_nch_before_dnn;
  std::vector<int> BsDstarTauNu_nch_qr;
  std::vector<int> BsDstarTauNu_ngentau3;
  std::vector<int> BsDstarTauNu_ngentau;
  std::vector<float> BsDstarTauNu_gentaupt;
  std::vector<int> BsDstarTauNu_gentaudm;





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
