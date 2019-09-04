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

  //=================================================================================================================== 
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
  std::vector<std::vector<int> >  genParticle_dau      ;
  std::vector<float >  genParticle_tauvispt      ;
  std::vector<float >  genParticle_tauviseta      ;
  std::vector<float >  genParticle_tauvisphi      ;
  std::vector<float >  genParticle_tauvismass      ;
  std::vector<int >  genParticle_taudecay      ;

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

          
  /** electrons */
  int 	                      el_N		         ;
  std::vector<int>  	      el_pdgId	                 ;
  std::vector<float>  	      el_charge	                 ;
  std::vector<float>  	      el_e 		         ;
  std::vector<float>  	      el_eta		         ;
  std::vector<float>  	      el_phi		         ;
  std::vector<float>  	      el_mass		         ;
  std::vector<float>  	      el_pt		         ;
  std::vector<float>  	      el_et		         ;
  std::vector<float>  	      el_superCluster_eta	 ;  
  std::vector<float>  	      el_pfRhoCorrRelIso03  	 ;
  std::vector<float>  	      el_pfRhoCorrRelIso04  	 ;
  std::vector<float>  	      el_pfDeltaCorrRelIso  	 ;
  /* std::vector<float>  	      el_pfRelIso  	     	 ; */
  std::vector<float>  	      el_photonIso 	     	 ;
  std::vector<float>  	      el_neutralHadIso	     	 ;
  std::vector<float>  	      el_chargedHadIso	     	 ;
  /* std::vector<float>  	      el_trackIso	         ; */      
  /* std::vector<float>  	      el_pfRhoCorrRelIso03Boost  ; */
  /* std::vector<float>  	      el_pfRhoCorrRelIso04Boost  ; */
  /* std::vector<float>  	      el_pfDeltaCorrRelIsoBoost  ; */
  /* std::vector<float>  	      el_pfRelIsoBoost  	 ; */
  /* std::vector<float>  	      el_photonIsoBoost 	 ; */
  /* std::vector<float>  	      el_neutralHadIsoBoost      ; */
  /* std::vector<float>  	      el_chargedHadIsoBoost      ; */
  std::vector<int>	      el_passConversionVeto	 ;
  std::vector<float>          el_full5x5_sigmaIetaIeta	 ;
  std::vector<float>          el_dEtaIn		         ;
  std::vector<float>          el_dPhiIn		         ;
  std::vector<float>          el_hOverE		         ;
  std::vector<float>          el_relIsoWithDBeta	 ;
  std::vector<float>          el_ooEmooP		 ;
  std::vector<int>	      el_expectedMissingInnerHits;
  std::vector<float>          el_d0			 ;
  std::vector<float>          el_dz			 ;
  std::vector<float>          el_d0_allvertices		 ;
  std::vector<float>          el_dz_allvertices		 ;
  std::vector<float>          el_dr03EcalRecHitSumEt      ;
  std::vector<float>          el_dr03HcalDepth1TowerSumEt ;
  std::vector<float>          el_rho                      ;
  std::vector<bool>           el_ecalDriven               ;
  std::vector<float>          el_dEtaInSeed               ;
  std::vector<float>          el_full5x5_e2x5Max          ;
  std::vector<float>          el_full5x5_e5x5             ;
  std::vector<float>          el_full5x5_e1x5             ;
  std::vector<float>          el_full5x5_r9               ;
  std::vector<float>          el_dr03TkSumPt              ;
  std::vector<float>          el_superCluster_e           ;
  std::vector<float>          el_hadronicOverEm           ;
  std::vector<float>          el_seedEnergy               ;
  std::vector<int>	      el_isVetoElectron	          ;
  std::vector<int>  	      el_isLooseElectron	  ; 
  std::vector<int>	      el_isMediumElectron	  ;
  std::vector<int>	      el_isTightElectron	  ;    
  std::vector<int  >  	      el_isHltElectron            ;
  std::vector<int  >  	      el_isHeepElectron	          ;
  std::vector<int>	      el_isMVAMediumElectron	  ;
  std::vector<int>	      el_isMVATightElectron	  ;    
  std::vector<float>	      el_MVAscore	  ;    
  std::vector<int>	      el_MVAcategory	  ;    
  std::vector<int>	      el_isVetoElectronWithoutIPandIsolation	  ;
  std::vector<int>	      el_isMediumElectronWithoutIPandIsolation  ;
  std::vector<int>	      el_isTightElectronWithoutIPandIsolation	  ;    
  std::vector<int  >  	      el_isHeepElectronWithoutIPandIsolation	  ;
  std::vector<int  >  	      el_isLooseElectronWithoutIPandIsolation	  ;  
//  std::vector<int>	      el_isVetoElectronBoosted	  ;
//  std::vector<int>	      el_isMediumElectronBoosted  ;
//  std::vector<int>	      el_isTightElectronBoosted	  ;    
//  std::vector<int  >  	      el_isHeepElectronBoosted	  ;
//  std::vector<int  >  	      el_isLooseElectronBoosted	  ;  
  std::vector<float>  	      el_SemileptonicPFIso 	 ;//  Isolations for semileptonic tau channel  
  //  std::vector<float>  	      el_SemileptonicCorrPFIso   ;// the simple PF one and the corrected one for the tau presence

  /** muons */
  int 	                      mu_N		         ;
  std::vector<int>  	      mu_pdgId	                 ;
  std::vector<float>  	      mu_charge	                 ;
  std::vector<float>  	      mu_e 		         ;
  std::vector<float>  	      mu_eta		         ;
  std::vector<float>  	      mu_phi		         ;
  std::vector<float>  	      mu_mass		         ;
  std::vector<float>  	      mu_pt		         ;    
  std::vector<int  >          mu_isHighPtMuon		 ;
  std::vector<int  >          mu_isTightMuon		 ;
  std::vector<int  >  	      mu_isMediumMuon	         ; // YT added
  std::vector<int  >  	      mu_isMediumMuonGH	         ; // YT added
  std::vector<int  >          mu_isLooseMuon		 ;
  std::vector<int  >          mu_isPFMuon		 ;
  std::vector<int  >          mu_isSoftMuon		 ;   
  std::vector<int  >  	      mu_isGlobalMuon	         ;
  std::vector<int  >  	      mu_isTrackerMuon	         ; // YT added
  std::vector<int  >  	      mu_isTrackerHighPtMuon	         ;
  /* std::vector<float>  	      mu_pfRhoCorrRelIso03  	 ; */
  /* std::vector<float>  	      mu_pfRhoCorrRelIso04  	 ; */
  std::vector<float>  	      mu_pfDeltaCorrRelIso  	 ;
  /* std::vector<float>  	      mu_pfRelIso  	     	 ; */
  std::vector<float>  	      mu_photonIso 	     	 ;
  std::vector<float>  	      mu_neutralHadIso	     	 ;
  std::vector<float>  	      mu_chargedHadIso	     	 ;
  std::vector<float>  	      mu_trackIso	         ;
  std::vector<float>  	      mu_trackCorrIso	         ;
  std::vector<float>          mu_d0			 ;
  std::vector<float>          mu_dz			 ;
  std::vector<float>          mu_d0_allvertices		 ;
  std::vector<float>          mu_dz_allvertices		 ;
  std::vector<float>  	      mu_innerTrack_pt 	         ;
  std::vector<float>  	      mu_bestTrack_pt  	         ;
  std::vector<float>  	      mu_bestTrack_ptErr  	 ;    
  std::vector<float>  	      mu_tunePTrack_pt  	       ;
  std::vector<float>  	      mu_tunePTrack_ptErr  	 ;   
  /* std::vector<float>  	      mu_pfRhoCorrRelIso03Boost  ; */
  /* std::vector<float>  	      mu_pfRhoCorrRelIso04Boost  ; */
  /* std::vector<float>  	      mu_pfDeltaCorrRelIsoBoost  ; */
  /* std::vector<float>  	      mu_pfRelIsoBoost  	 ;     */
  /* std::vector<float>  	      mu_photonIsoBoost 	 ; */
  /* std::vector<float>  	      mu_neutralHadIsoBoost      ; */
  /* std::vector<float>  	      mu_chargedHadIsoBoost      ;   */
  std::vector<float>  	      mu_normChi2  	         ;
  std::vector<int  >  	      mu_trackerHits	         ;
  std::vector<int  >  	      mu_matchedStations         ;
  std::vector<int  >  	      mu_pixelHits 	         ;
  std::vector<int  >  	      mu_globalHits	         ;
  std::vector<float>  	      mu_SemileptonicPFIso 	 ;//  Isolations for semileptonic tau channel  
  //  std::vector<float>  	      mu_SemileptonicCorrPFIso   ;// the simple PF one and the corrected one for the tau presence

  /** taus */
  int 	                      tau_N		         ;
  std::vector<int>  	      tau_pdgId	                 ;
  std::vector<float>  	      tau_charge	         ;
  std::vector<float>  	      tau_e 		         ;
  std::vector<float>  	      tau_eta		         ;
  std::vector<float>  	      tau_phi		         ;
  std::vector<float>  	      tau_mass		         ;
  std::vector<float>  	      tau_pt		         ;  
  std::vector<float>          tau_d0			 ;  
  std::vector<float>          tau_dz			 ;  // YT added

  // YT added
//  std::vector<std::vector<int>>	      tau_associated_pdgId       ;  
//  std::vector<std::vector<float>>	      tau_associated_pt       ;  
//  std::vector<std::vector<float>>	      tau_associated_eta       ;  
//  std::vector<std::vector<float>>	      tau_associated_phi       ;  
//  std::vector<std::vector<float>>	      tau_associated_dr       ;  

//  std::vector<int>  	      tau_n_total                 ;
  std::vector<int>  	      tau_n_ch                 ;
  std::vector<int>  	      tau_n_nh                 ;
  //  std::vector<int>  	      tau_n_h_f                 ;
  //  std::vector<int>  	      tau_n_em_f                 ;
  std::vector<int>  	      tau_n_gamma                 ;
  //  std::vector<int>  	      tau_n_e                 ;
  //  std::vector<int>  	      tau_n_mu                 ;
  /* // */
  /* std::vector<float>  	      tau_pfRhoCorrRelIso03  	 ; */
  /* std::vector<float>  	      tau_pfRhoCorrRelIso04  	 ; */
  /* std::vector<float>  	      tau_pfDeltaCorrRelIso  	 ; */
  /* std::vector<float>  	      tau_pfRelIso  	     	 ; */
  /* std::vector<float>  	      tau_photonIso 	     	 ; */
  /* std::vector<float>  	      tau_neutralHadIso	     	 ; */
  /* std::vector<float>  	      tau_chargedHadIso	     	 ; */
  /* // */
  std::vector<float>  	      tau_photonPtSumOutsideSignalCone;
  /* // */
  /* std::vector<float>  	      tau_trackIso	         ; */
  /*  std::vector<float>  	      tau_pfRhoCorrRelIso03Boost ; */
  /* std::vector<float>  	      tau_pfRhoCorrRelIso04Boost ; */
  /* std::vector<float>  	      tau_pfDeltaCorrRelIsoBoost ; */
  /* std::vector<float>  	      tau_pfRelIsoBoost  	 ;     */
  /* std::vector<float>  	      tau_photonIsoBoost 	 ; */
  /* std::vector<float>  	      tau_neutralHadIsoBoost     ; */
  /* std::vector<float>  	      tau_chargedHadIsoBoost     ;   */
  /* // */
  std::vector<int  >  	      tau_TauType	         ;  
  std::vector<int  >  	      tau_decayMode	         ;  // YT added
  std::vector<float>  	      tau_chargedPionPt	         ;  // YT added
  std::vector<float>  	      tau_neutralPionPt	         ;  // YT added
  // YT added : newly added for the MVA training
  std::vector<float>  	      tau_nPhoton		 ;  
//  std::vector<float>  	      tau_nPhoton_1		 ;  
//  std::vector<float>  	      tau_nPhoton_1p5		 ;  
//  std::vector<float>  	      tau_nPhoton_2		 ;  
//  std::vector<float>  	      tau_nPhoton_2p5		 ;  
//  std::vector<float>  	      tau_nPhoton_3		 ;  
//  std::vector<float>  	      tau_nPhoton_4		 ;  
//  std::vector<float>  	      tau_nPhoton_5		 ;  
  std::vector<float>          tau_ptWeightedDetaStrip    ;
  std::vector<float>          tau_ptWeightedDphiStrip    ;
  std::vector<float>          tau_ptWeightedDrSignal     ;
  std::vector<float>          tau_ptWeightedDrIsolation  ;
  std::vector<float>          tau_leadingTrackChi2       ;
  std::vector<float>          tau_leadingTrackPt         ;
  std::vector<float>          tau_eRatio                 ;
  std::vector<float>          tau_dxy_Sig                ;
  std::vector<float>          tau_ip3d                   ;
  std::vector<float>          tau_ip3d_Sig               ;
  std::vector<bool>            tau_hasSecondaryVertex     ;
  //  std::vector<float>          tau_decayDistMag_x         ;
  //  std::vector<float>          tau_decayDistMag_y         ;
  //  std::vector<float>          tau_decayDistMag_z         ;
  std::vector<float>          tau_decayDistMag           ;
  std::vector<float>          tau_flightLenthSig         ;
  // YT added end : newly added for the MVA training


  /** tau discriminants */
  std::vector<bool>  	      tau_decayModeFindingNewDMs	              ;
  std::vector<bool>           tau_decayModeFinding  			      ;
  std::vector<bool>  	      tau_byLooseCombinedIsolationDeltaBetaCorr3Hits  ;
  std::vector<bool>  	      tau_byMediumCombinedIsolationDeltaBetaCorr3Hits ;
  std::vector<bool>  	      tau_byTightCombinedIsolationDeltaBetaCorr3Hits  ;
  std::vector<float>  	      tau_byCombinedIsolationDeltaBetaCorrRaw3Hits    ;
  std::vector<float>  	      tau_chargedIsoPtSum			      ;
  std::vector<float>  	      tau_neutralIsoPtSum			      ;
//  std::vector<float>  	      tau_neutralIsoPtSum_0p5			      ;
//  std::vector<float>  	      tau_neutralIsoPtSum_1			      ;
//  std::vector<float>  	      tau_neutralIsoPtSum_1p5			      ;
//  std::vector<float>  	      tau_neutralIsoPtSum_2			      ;
//  std::vector<float>  	      tau_neutralIsoPtSum_2p5			      ;
//  std::vector<float>  	      tau_neutralIsoPtSum_3			      ;
//  std::vector<float>  	      tau_neutralIsoPtSum_4			      ;
//  std::vector<float>  	      tau_neutralIsoPtSum_5			      ;
  std::vector<float>  	      tau_puCorrPtSum				      ;
  std::vector<float>  	      tau_chargedIsoPtSumdR03			      ;
  //  std::vector<float>  	      tau_footprintCorrectiondR03			      ;
  std::vector<float>  	      tau_neutralIsoPtSumdR03			      ;
  //  std::vector<float>  	      tau_neutralIsoPtSumWeight			      ;
  //  std::vector<float>  	      tau_neutralIsoPtSumWeightdR03			      ;
  std::vector<float>  	      tau_photonPtSumOutsideSignalConedR03			      ;

  std::vector<float>  	      tau_byIsolationMVArun2v1DBdR03oldDMwLTraw			      ;
  std::vector<float>  	      tau_byIsolationMVArun2v1DBnewDMwLTraw			      ;
  std::vector<float>  	      tau_byIsolationMVArun2v1DBoldDMwLTraw			      ;
  //  std::vector<float>  	      tau_byIsolationMVArun2v1DBoldDMwoLTraw			      ;
  //  std::vector<float>  	      tau_byIsolationMVArun2v1PWdR03oldDMwLTraw			      ;
  std::vector<float>  	      tau_byIsolationMVArun2v1PWnewDMwLTraw			      ;
  std::vector<float>  	      tau_byIsolationMVArun2v1PWoldDMwLTraw			      ;
  std::vector<bool>  	      tau_byLooseIsolationMVArun2v1DBdR03oldDMwLT			      ;
  std::vector<bool>  	      tau_byLooseIsolationMVArun2v1DBnewDMwLT			      ;
  std::vector<bool>  	      tau_byLooseIsolationMVArun2v1DBoldDMwLT			      ;
  //  std::vector<bool>  	      tau_byLooseIsolationMVArun2v1PWdR03oldDMwLT			      ;
  std::vector<bool>  	      tau_byLooseIsolationMVArun2v1PWnewDMwLT			      ;
  std::vector<bool>  	      tau_byLooseIsolationMVArun2v1PWoldDMwLT			      ;

  std::vector<bool>  	      tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT			      ;
  std::vector<bool>  	      tau_byMediumIsolationMVArun2v1DBnewDMwLT			      ;
  std::vector<bool>  	      tau_byMediumIsolationMVArun2v1DBoldDMwLT			      ;
  //  std::vector<bool>  	      tau_byMediumIsolationMVArun2v1PWdR03oldDMwLT			      ;
  std::vector<bool>  	      tau_byMediumIsolationMVArun2v1PWnewDMwLT			      ;
  std::vector<bool>  	      tau_byMediumIsolationMVArun2v1PWoldDMwLT;
 
  std::vector<bool>  	      tau_byTightIsolationMVArun2v1DBdR03oldDMwLT			      ;
  std::vector<bool>  	      tau_byTightIsolationMVArun2v1DBnewDMwLT			      ;
  std::vector<bool>  	      tau_byTightIsolationMVArun2v1DBoldDMwLT			      ;
  //  std::vector<bool>  	      tau_byTightIsolationMVArun2v1PWdR03oldDMwLT			      ;
  std::vector<bool>  	      tau_byTightIsolationMVArun2v1PWnewDMwLT			      ;
  std::vector<bool>  	      tau_byTightIsolationMVArun2v1PWoldDMwLT			      ;
  std::vector<bool>  	      tau_byVLooseIsolationMVArun2v1DBdR03oldDMwLT			      ;

  std::vector<bool>  	      tau_byVLooseIsolationMVArun2v1DBnewDMwLT			      ;
  std::vector<bool>  	      tau_byVLooseIsolationMVArun2v1DBoldDMwLT			      ;
  std::vector<bool>           tau_byVVLooseIsolationMVArun2v1DBoldDMwLT                        ;
   
  //  std::vector<bool>  	      tau_byVLooseIsolationMVArun2v1PWdR03oldDMwLT			      ;
  std::vector<bool>  	      tau_byVLooseIsolationMVArun2v1PWnewDMwLT			      ;
  std::vector<bool>  	      tau_byVLooseIsolationMVArun2v1PWoldDMwLT			      ;
  std::vector<bool>  	      tau_byVTightIsolationMVArun2v1DBdR03oldDMwLT			      ;
  std::vector<bool>  	      tau_byVTightIsolationMVArun2v1DBnewDMwLT			      ;
  std::vector<bool>  	      tau_byVTightIsolationMVArun2v1DBoldDMwLT			      ;

  //  std::vector<bool>  	      tau_byVTightIsolationMVArun2v1PWdR03oldDMwLT			      ;
  std::vector<bool>  	      tau_byVTightIsolationMVArun2v1PWnewDMwLT			      ;
  std::vector<bool>  	      tau_byVTightIsolationMVArun2v1PWoldDMwLT			      ;
  std::vector<bool>  	      tau_byVVTightIsolationMVArun2v1DBdR03oldDMwLT			      ;
  std::vector<bool>  	      tau_byVVTightIsolationMVArun2v1DBnewDMwLT			      ;
  std::vector<bool>  	      tau_byVVTightIsolationMVArun2v1DBoldDMwLT			      ;
  //  std::vector<bool>  	      tau_byVVTightIsolationMVArun2v1PWdR03oldDMwLT			      ;
  std::vector<bool>  	      tau_byVVTightIsolationMVArun2v1PWnewDMwLT			      ;

  std::vector<bool>  	      tau_byVVTightIsolationMVArun2v1PWoldDMwLT			      ;

  // YT added : newly added for the MVA training
  std::vector<float>  	      tau_byIsolationMVArun2v2DBoldDMwLTraw;
  std::vector<bool>  	      tau_byVVLooseIsolationMVArun2v2DBoldDMwLT;
  std::vector<bool>  	      tau_byVLooseIsolationMVArun2v2DBoldDMwLT;
  std::vector<bool>  	      tau_byLooseIsolationMVArun2v2DBoldDMwLT;
  std::vector<bool>  	      tau_byMediumIsolationMVArun2v2DBoldDMwLT;
  std::vector<bool>  	      tau_byTightIsolationMVArun2v2DBoldDMwLT;
  std::vector<bool>  	      tau_byVTightIsolationMVArun2v2DBoldDMwLT;
  std::vector<bool>  	      tau_byVVTightIsolationMVArun2v2DBoldDMwLT;
  // YT added : newly added for the MVA training

  std::vector<float>  	      tau_againstElectronMVA6raw			      ;
  std::vector<float>  	      tau_againstElectronMVA6category			      ;
  std::vector<bool>  	      tau_againstElectronVLooseMVA6			      ;
  std::vector<bool>  	      tau_againstElectronLooseMVA6			      ;
  std::vector<bool>  	      tau_againstElectronMediumMVA6			      ;
  std::vector<bool>  	      tau_againstElectronTightMVA6			      ;
  std::vector<bool>  	      tau_againstElectronVTightMVA6			      ;
      

  std::vector<bool>  	      tau_againstMuonLoose3			      ;
  std::vector<bool>  	      tau_againstMuonTight3			      ; 
      
 
  std::vector<float>  	      tau_byPhotonPtSumOutsideSignalCone			      ;
  //  std::vector<float>  	      tau_footprintCorrection			      ;



  
/* /\*----------------------Tau tracks---------------------------*\/ */
  
/*     std::vector<int>    nCharCand              ; */
/*     std::vector<int>    nNeuCand               ; */
/*     std::vector<int>    nGamCand               ; */
/*     std::vector<int>    decayMode              ; */
/*     std::vector<float>  leadTrack_dxy          ; */
/*     std::vector<float>  leadTrack_dxySig       ; */
/*     std::vector<float>  secVtxX                ; */
/*     std::vector<float>  secVtxY                ; */
/*     std::vector<float>  secVtxZ                ; */
/*     std::vector<float>  flightLengthX          ; */
/*     std::vector<float>  flightLengthY          ; */
/*     std::vector<float>  flightLengthZ          ; */
/*     std::vector<float>  flightLength           ; */
/*     std::vector<float>  flightLengthSig        ; */



  
   
  
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
  std::vector<int  >                IsJpsiEle   ;   

  std::vector<float>                Jpsi_mu1_pt      ;
  std::vector<float>                Jpsi_mu1_eta     ;
  std::vector<float>                Jpsi_mu1_phi     ;
  std::vector<int  >                Jpsi_mu1_ch      ;   
  std::vector<int  >                Jpsi_mu1_isLoose   ;
  std::vector<int  >                Jpsi_mu1_isTight   ;
  std::vector<int  >                Jpsi_mu1_isPF      ;
  std::vector<int  >                Jpsi_mu1_isGlobal  ;
  std::vector<int  >                Jpsi_mu1_isTracker ;
  std::vector<int  >                Jpsi_mu1_isSoft    ;

  std::vector<float>                Jpsi_mu2_pt      ;
  std::vector<float>                Jpsi_mu2_eta     ;
  std::vector<float>                Jpsi_mu2_phi     ;
  std::vector<int  >                Jpsi_mu2_ch   ;   
  std::vector<float>                Jpsi_mu2_isLoose   ;
  std::vector<float>                Jpsi_mu2_isTight   ;
  std::vector<float>                Jpsi_mu2_isPF      ;
  std::vector<float>                Jpsi_mu2_isGlobal  ;
  std::vector<float>                Jpsi_mu2_isTracker ;
  std::vector<float>                Jpsi_mu2_isSoft    ;

  std::vector<float>                Jpsi_mu3_pt      ;
  std::vector<float>                Jpsi_mu3_eta     ;
  std::vector<float>                Jpsi_mu3_phi     ;
  std::vector<int  >                Jpsi_mu3_ch   ;   
  std::vector<float>                Jpsi_mu3_isLoose   ;
  std::vector<float>                Jpsi_mu3_isTight   ;
  std::vector<float>                Jpsi_mu3_isPF      ;
  std::vector<float>                Jpsi_mu3_isGlobal  ;
  std::vector<float>                Jpsi_mu3_isTracker ;
  std::vector<float>                Jpsi_mu3_isSoft    ;
  std::vector<float>                Jpsi_mu3_x       ;
  std::vector<float>                Jpsi_mu3_y       ;
  std::vector<float>                Jpsi_mu3_z       ;

  std::vector<float>                Jpsi_mu3_isopt03      ;
  std::vector<float>                Jpsi_mu3_isopt04      ;
  std::vector<float>                Jpsi_mu3_isopt05      ;
  std::vector<float>                Jpsi_mu3_isopt06      ;
  std::vector<float>                Jpsi_mu3_isopt07      ;
  std::vector<float>                Jpsi_dr_mu3pf      ;


  std::vector<float>                Jpsi_dx      ;
  std::vector<float>                Jpsi_dy      ;
  std::vector<float>                Jpsi_dz      ;
  std::vector<float>                Jpsi_pt      ;
  std::vector<float>                Jpsi_eta     ;
  std::vector<float>                Jpsi_phi     ;
  std::vector<float>                Jpsi_mass       ;
  std::vector<float>                Jpsi_vtxprob    ;
  std::vector<float>                Jpsi_vtxz       ;


  std::vector<float>                Jpsi_trimu_dx      ;
  std::vector<float>                Jpsi_trimu_dy      ;
  std::vector<float>                Jpsi_trimu_dz      ;
  std::vector<float>                Jpsi_trimu_pt      ;
  std::vector<float>                Jpsi_trimu_eta     ;
  std::vector<float>                Jpsi_trimu_phi     ;
  std::vector<float>                Jpsi_trimu_mass    ;
  std::vector<float>                Jpsi_trimu_vtxprob ;
  std::vector<float>                Jpsi_trimu_vtxz    ;

  std::vector<float>                Jpsi_PV_x       ;
  std::vector<float>                Jpsi_PV_y       ;
  std::vector<float>                Jpsi_PV_z       ;

  std::vector<float>                Jpsi_flightSig3D             ; 
  std::vector<float>                Jpsi_flightLength3D          ;
  std::vector<float>                Jpsi_flightLengthErr3D       ;
  std::vector<float>                Jpsi_flightSig2D             ; 
  std::vector<float>                Jpsi_flightLength2D          ;
  std::vector<float>                Jpsi_flightLengthErr2D       ;
  std::vector<float>                Jpsi_trimu_flightSig3D       ; 
  std::vector<float>                Jpsi_trimu_flightLength3D    ;
  std::vector<float>                Jpsi_trimu_flightLengthErr3D ;
  std::vector<float>                Jpsi_trimu_flightSig2D       ; 
  std::vector<float>                Jpsi_trimu_flightLength2D    ;
  std::vector<float>                Jpsi_trimu_flightLengthErr2D ;


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
