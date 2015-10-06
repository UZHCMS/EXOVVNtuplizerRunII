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
  std::vector<float>              genParticle_px       ;
  std::vector<float>              genParticle_py       ;
  std::vector<float>              genParticle_pz       ;
  std::vector<float>              genParticle_e        ;
  std::vector<float>              genParticle_eta      ;
  std::vector<float>              genParticle_phi      ;
  std::vector<float>              genParticle_mass     ;
  std::vector<int  >              genParticle_pdgId    ;
  std::vector<int  >              genParticle_status   ;
  std::vector<int  >              genParticle_nDau     ;
  std::vector<int  >              genParticle_nMoth    ;
  std::vector<std::vector<int> >  genParticle_mother   ; 
  std::vector<std::vector<int> >  genParticle_dau      ;

  /** generator info */
  float                           lheV_pt              ;
  float                           lheHT                ;
  float                           lheNj                ;
  float                           genWeight            ;
  float                           qScale               ;
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
  std::vector<float>  	      el_pfRelIso  	     	 ;
  std::vector<float>  	      el_photonIso 	     	 ;
  std::vector<float>  	      el_neutralHadIso	     	 ;
  std::vector<float>  	      el_chargedHadIso	     	 ;
  std::vector<float>  	      el_trackIso	         ;            


  std::vector<float>  	      el_pfRhoCorrRelIso03Boost  ;
  std::vector<float>  	      el_pfRhoCorrRelIso04Boost  ;
  std::vector<float>  	      el_pfDeltaCorrRelIsoBoost  ;
  std::vector<float>  	      el_pfRelIsoBoost  	 ;    
  std::vector<float>  	      el_photonIsoBoost 	 ;
  std::vector<float>  	      el_neutralHadIsoBoost      ;
  std::vector<float>  	      el_chargedHadIsoBoost      ;
  
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

  std::vector<float>          el_dr03EcalRecHitSumEt      ;
  std::vector<float>          el_dr03HcalDepth1TowerSumEt ;
  std::vector<float>          el_rho                      ;
  std::vector<bool>           el_ecalDriven               ;
  std::vector<float>          el_dEtaInSeed               ;
  std::vector<float>          el_full5x5_e2x5Max          ;
  std::vector<float>          el_full5x5_e5x5             ;
  std::vector<float>          el_full5x5_e1x5             ;
  std::vector<float>          el_dr03TkSumPt              ;
  std::vector<float>          el_superCluster_e           ;
  std::vector<float>          el_hadronicOverEm           ;

  std::vector<int>	      el_isVetoElectron	          ;
  std::vector<int>	      el_isMediumElectron	  ;
  std::vector<int>	      el_isTightElectron	  ;    
  std::vector<int  >  	      el_isHeepElectron	          ;
  std::vector<int  >  	      el_isHeep51Electron         ;
  std::vector<int  >  	      el_isLooseElectron	  ;
 
  std::vector<int>	      el_isVetoElectronBoosted	  ;
  std::vector<int>	      el_isMediumElectronBoosted  ;
  std::vector<int>	      el_isTightElectronBoosted	  ;    
  std::vector<int  >  	      el_isHeepElectronBoosted	  ;
  std::vector<int  >  	      el_isHeep51ElectronBoosted  ;
  std::vector<int  >  	      el_isLooseElectronBoosted	  ;
  
  std::vector<float>  	      el_SemileptonicPFIso 	 ;//  Isolations for semileptonic tau channel  
  std::vector<float>  	      el_SemileptonicCorrPFIso   ;// the simple PF one and the corrected one for the tau presence


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
  std::vector<int  >          mu_isLooseMuon		 ;
  std::vector<int  >          mu_isPFMuon		 ;
  std::vector<int  >          mu_isSoftMuon		 ; 
  
  std::vector<float>  	      mu_pfRhoCorrRelIso03  	 ;
  std::vector<float>  	      mu_pfRhoCorrRelIso04  	 ;
  std::vector<float>  	      mu_pfDeltaCorrRelIso  	 ;
  std::vector<float>  	      mu_pfRelIso  	     	 ;
  std::vector<float>  	      mu_photonIso 	     	 ;
  std::vector<float>  	      mu_neutralHadIso	     	 ;
  std::vector<float>  	      mu_chargedHadIso	     	 ;
  std::vector<float>  	      mu_trackIso	         ;
  std::vector<float>          mu_d0			 ;
  std::vector<float>          mu_dz			 ;
  std::vector<float>  	      mu_bestTrack_pt  	         ;
  std::vector<float>  	      mu_bestTrack_ptErr  	 ;
    
  std::vector<float>  	      mu_pfRhoCorrRelIso03Boost  ;
  std::vector<float>  	      mu_pfRhoCorrRelIso04Boost  ;
  std::vector<float>  	      mu_pfDeltaCorrRelIsoBoost  ;
  std::vector<float>  	      mu_pfRelIsoBoost  	 ;    
  std::vector<float>  	      mu_photonIsoBoost 	 ;
  std::vector<float>  	      mu_neutralHadIsoBoost      ;
  std::vector<float>  	      mu_chargedHadIsoBoost      ;
  
  std::vector<float>  	      mu_normChi2  	         ;
  std::vector<int  >  	      mu_isGlobalMuon	         ;
  std::vector<int  >  	      mu_trackerHits	         ;
  std::vector<int  >  	      mu_matchedStations         ;
  std::vector<int  >  	      mu_pixelHits 	         ;
  std::vector<int  >  	      mu_globalHits	         ;

  std::vector<float>  	      mu_SemileptonicPFIso 	 ;//  Isolations for semileptonic tau channel  
  std::vector<float>  	      mu_SemileptonicCorrPFIso   ;// the simple PF one and the corrected one for the tau presence

  /** taus */
  int 	                      tau_N		         ;
  std::vector<int>  	      tau_pdgId	                 ;
  std::vector<float>  	      tau_charge	         ;
  std::vector<float>  	      tau_e 		         ;
  std::vector<float>  	      tau_eta		         ;
  std::vector<float>  	      tau_phi		         ;
  std::vector<float>  	      tau_mass		         ;
  std::vector<float>  	      tau_pt		         ;
  
  std::vector<float>  	      tau_pfRhoCorrRelIso03  	 ;
  std::vector<float>  	      tau_pfRhoCorrRelIso04  	 ;
  std::vector<float>  	      tau_pfDeltaCorrRelIso  	 ;
  std::vector<float>  	      tau_pfRelIso  	     	 ;
  std::vector<float>  	      tau_photonIso 	     	 ;
  std::vector<float>  	      tau_neutralHadIso	     	 ;
  std::vector<float>  	      tau_chargedHadIso	     	 ;
  std::vector<float>  	      tau_trackIso	         ;
  std::vector<float>          tau_d0			 ;
  
  std::vector<float>  	      tau_pfRhoCorrRelIso03Boost ;
  std::vector<float>  	      tau_pfRhoCorrRelIso04Boost ;
  std::vector<float>  	      tau_pfDeltaCorrRelIsoBoost ;
  std::vector<float>  	      tau_pfRelIsoBoost  	 ;    
  std::vector<float>  	      tau_photonIsoBoost 	 ;
  std::vector<float>  	      tau_neutralHadIsoBoost     ;
  std::vector<float>  	      tau_chargedHadIsoBoost     ;
  
  std::vector<int  >  	      tau_TauType	         ;  
  
  /** tau discriminants */
  std::vector<float>  	      tau_decayModeFindingNewDMs	              ;
  std::vector<float>          tau_decayModeFinding  			      ;
  std::vector<float>  	      tau_byLooseCombinedIsolationDeltaBetaCorr3Hits  ;
  std::vector<float>  	      tau_byMediumCombinedIsolationDeltaBetaCorr3Hits ;
  std::vector<float>  	      tau_byTightCombinedIsolationDeltaBetaCorr3Hits  ;
  std::vector<float>  	      tau_byCombinedIsolationDeltaBetaCorrRaw3Hits    ;
  std::vector<float>  	      tau_chargedIsoPtSum			      ;
  std::vector<float>  	      tau_neutralIsoPtSum			      ;
  std::vector<float>  	      tau_puCorrPtSum				      ;
  std::vector<float>  	      tau_byIsolationMVA3oldDMwoLTraw		      ;
  std::vector<float>  	      tau_byVLooseIsolationMVA3oldDMwoLT	      ;
  std::vector<float>  	      tau_byLooseIsolationMVA3oldDMwoLT		      ;
  std::vector<float>  	      tau_byMediumIsolationMVA3oldDMwoLT	      ;
  std::vector<float>  	      tau_byTightIsolationMVA3oldDMwoLT		      ;
  std::vector<float>  	      tau_byVTightIsolationMVA3oldDMwoLT	      ;
  std::vector<float>  	      tau_byVVTightIsolationMVA3oldDMwoLT	      ;
  std::vector<float>  	      tau_byIsolationMVA3oldDMwLTraw		      ;
  std::vector<float>  	      tau_byVLooseIsolationMVA3oldDMwLT		      ;
  std::vector<float>  	      tau_byLooseIsolationMVA3oldDMwLT		      ;
  std::vector<float>  	      tau_byMediumIsolationMVA3oldDMwLT		      ;
  std::vector<float>  	      tau_byTightIsolationMVA3oldDMwLT		      ;
  std::vector<float>  	      tau_byVTightIsolationMVA3oldDMwLT		      ;
  std::vector<float>  	      tau_byVVTightIsolationMVA3oldDMwLT	      ;
  std::vector<float>  	      tau_byIsolationMVA3newDMwoLTraw		      ;
  std::vector<float>  	      tau_byVLooseIsolationMVA3newDMwoLT	      ;
  std::vector<float>  	      tau_byLooseIsolationMVA3newDMwoLT		      ;
  std::vector<float>  	      tau_byMediumIsolationMVA3newDMwoLT	      ;
  std::vector<float>  	      tau_byTightIsolationMVA3newDMwoLT		      ;
  std::vector<float>  	      tau_byVTightIsolationMVA3newDMwoLT	      ;
  std::vector<float>  	      tau_byVVTightIsolationMVA3newDMwoLT	      ;
  std::vector<float>  	      tau_byIsolationMVA3newDMwLTraw		      ;
  std::vector<float>  	      tau_byVLooseIsolationMVA3newDMwLT		      ;
  std::vector<float>  	      tau_byLooseIsolationMVA3newDMwLT		      ;
  std::vector<float>  	      tau_byMediumIsolationMVA3newDMwLT		      ;
  std::vector<float>  	      tau_byTightIsolationMVA3newDMwLT		      ;
  std::vector<float>  	      tau_byVTightIsolationMVA3newDMwLT		      ;
  std::vector<float>  	      tau_byVVTightIsolationMVA3newDMwLT	      ;
  std::vector<float>  	      tau_againstElectronLoose			      ;
  std::vector<float>  	      tau_againstElectronMedium			      ;
  std::vector<float>  	      tau_againstElectronTight			      ;
  std::vector<float>  	      tau_againstElectronMVA5raw		      ;
  std::vector<float>  	      tau_againstElectronMVA5category		      ;
  std::vector<float>  	      tau_againstElectronVLooseMVA5 		      ;
  std::vector<float>  	      tau_againstElectronLooseMVA5  		      ;
  std::vector<float>  	      tau_againstElectronMediumMVA5 		      ;
  std::vector<float>  	      tau_againstElectronTightMVA5  		      ;
  std::vector<float>  	      tau_againstElectronVTightMVA5 		      ;
  std::vector<float>  	      tau_againstMuonLoose  			      ;
  std::vector<float>  	      tau_againstMuonMedium 			      ;
  std::vector<float>  	      tau_againstMuonTight  			      ;
  std::vector<float>  	      tau_againstMuonLoose2 			      ;
  std::vector<float>  	      tau_againstMuonMedium2			      ;
  std::vector<float>  	      tau_againstMuonTight2 			      ;
  std::vector<float>  	      tau_againstMuonLoose3 			      ;
  std::vector<float>  	      tau_againstMuonTight3 			      ;
  std::vector<float>  	      tau_againstMuonMVAraw 			      ;
  std::vector<float>  	      tau_againstMuonLooseMVA			      ;
  std::vector<float>  	      tau_againstMuonMediumMVA			      ;
  std::vector<float>  	      tau_againstMuonTightMVA			      ;
    
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


  /** energy density */
  float                             rho;


    
  /** AK4 jets */
  int                             jetAK4_N             ;
  std::vector<float>              jetAK4_pt            ;
  std::vector<float>              jetAK4_eta           ;
  std::vector<float>              jetAK4_mass          ;
  std::vector<float>              jetAK4_phi           ;
  std::vector<float>              jetAK4_e             ;
  std::vector<float>              jetAK4_jec           ;
  //std::vector<float>	      jetAK4_jecUp	       ;
  //std::vector<float>	      jetAK4_jecDown	       ; 
  std::vector<bool>              jetAK4_IDLoose        ;
  std::vector<bool>              jetAK4_IDTight        ;

  std::vector<float>              jetAK4_muf           ;
  std::vector<float>              jetAK4_phf           ;
  std::vector<float>              jetAK4_emf           ;
  std::vector<float>              jetAK4_nhf           ;
  std::vector<float>              jetAK4_chf           ;
  std::vector<float>              jetAK4_area          ;
  std::vector<int>                jetAK4_cm            ;
  std::vector<int>                jetAK4_nm            ;
  std::vector<float>              jetAK4_che           ;
  std::vector<float>              jetAK4_ne            ;
  std::vector<float>              jetAK4_hf_hf         ;
  std::vector<float>              jetAK4_hf_emf        ;
  std::vector<float>              jetAK4_hof           ;
  std::vector<int>                jetAK4_chm           ;
  std::vector<int>                jetAK4_neHadMult     ;
  std::vector<int>                jetAK4_phoMult       ;
  std::vector<float>              jetAK4_nemf          ;
  std::vector<float>              jetAK4_cemf          ;
  std::vector<int>                jetAK4_charge	       ;
  std::vector<int>                jetAK4_partonFlavour ;
  std::vector<int>                jetAK4_hadronFlavour ;
  std::vector<int>                jetAK4_genParton_pdgID;
  std::vector<int>                jetAK4_nbHadrons     ;
  std::vector<int>                jetAK4_ncHadrons     ;
//std::vector<float>              jetAK4_ssv           ;
  std::vector<float>              jetAK4_cisv 	       ;         
//std::vector<float>  	          jetAK4_tchp 	       ;
//std::vector<float>  	          jetAK4_tche 	       ;
//std::vector<float>  	          jetAK4_jp	       ;
//std::vector<float>  	          jetAK4_jbp  	       ;
  std::vector<float>              jetAK4_vtxMass       ;
  std::vector<float>              jetAK4_vtxNtracks    ;
  std::vector<float>              jetAK4_vtx3DVal      ;
  std::vector<float>              jetAK4_vtx3DSig      ; 
//std::vector<int  >              jetAK4_nSVs          ;



  /** AK8 jets */
  int 	        	      jetAK8_N                 ;
  std::vector<float>  	      jetAK8_pt                ;
  std::vector<float>  	      jetAK8_eta               ;
  std::vector<float>  	      jetAK8_mass              ;
  std::vector<float>  	      jetAK8_phi               ;
  std::vector<float>  	      jetAK8_e                 ;
  std::vector<float>  	      jetAK8_jec               ;
  //std::vector<float>        jetAK8_jecUp             ;
  //std::vector<float>        jetAK8_jecDown           ;
  std::vector<bool >  	      jetAK8_IDLoose      ;
  std::vector<bool >  	      jetAK8_IDTight      ;
  std::vector<float>  	      jetAK8_muf          ;
  std::vector<float>  	      jetAK8_phf          ;

  std::vector<float>  	      jetAK8_emf	       ;
  std::vector<float>  	      jetAK8_nhf	       ;
  std::vector<float>  	      jetAK8_chf	       ;
  std::vector<float>  	      jetAK8_area	       ;
  std::vector<int  >          jetAK8_cm 	       ;
  std::vector<int  >          jetAK8_nm 	       ;
  std::vector<float>          jetAK8_che	       ;
  std::vector<float>          jetAK8_ne 	       ;
  std::vector<float>          jetAK8_hf_hf         ;
  std::vector<float>          jetAK8_hf_emf        ;
  std::vector<float>          jetAK8_hof           ;
  std::vector<int>            jetAK8_chm           ;
  std::vector<int>            jetAK8_neHadMult     ;
  std::vector<int>            jetAK8_phoMult       ;
  std::vector<float>          jetAK8_nemf          ;
  std::vector<float>          jetAK8_cemf          ;
  std::vector<int  >  	      jetAK8_charge	       ;
  std::vector<int>            jetAK8_partonFlavour ;
  std::vector<int>            jetAK8_hadronFlavour ;
  std::vector<int>            jetAK8_genParton_pdgID;
  std::vector<int>            jetAK8_nbHadrons     ;
  std::vector<int>            jetAK8_ncHadrons     ;
  std::vector<float>  	      jetAK8_Hbbtag 	       ;
//std::vector<float>          jetAK8_ssv               ;
  std::vector<float>  	      jetAK8_csv 	       ;    
//std::vector<float>          jetAK8_tchp              ;
//std::vector<float>          jetAK8_tche              ;
//std::vector<float>          jetAK8_jp                ;
//std::vector<float>          jetAK8_jbp               ;
  std::vector<float>  	      jetAK8_tau1              ;
  std::vector<float>  	      jetAK8_tau2              ;
  std::vector<float>  	      jetAK8_tau3              ;    
  std::vector<float>  	      jetAK8_pruned_mass       ;
  std::vector<float>  	      jetAK8_softdrop_mass     ;
  std::vector<float>	      jetAK10_trimmed_mass     ;
  std::vector<float>  	      jetAK8_pruned_massCorr   ;
  std::vector<float>  	      jetAK8_softdrop_massCorr ;
  std::vector<float>	      jetAK10_trimmed_massCorr ;
  std::vector<float>  	      jetAK8_pruned_jec        ;
  std::vector<float>  	      jetAK8_softdrop_jec      ;
  std::vector<float>	      jetAK10_trimmed_jec      ;
  //std::vector<float>	      jetAK8_filtered_mass     ;
  //std::vector<float>	      jetAK8_nSubJets	       ;
  std::vector<float>  	      jetAK10_ecf1              ;
  std::vector<float>  	      jetAK10_ecf2              ;
  std::vector<float>  	      jetAK10_ecf3              ;    
  std::vector<float>  	      jetAK8_puppi_tau1              ;
  std::vector<float>  	      jetAK8_puppi_tau2              ;
  std::vector<float>  	      jetAK8_puppi_tau3              ;    
  std::vector<float>  	      jetAK8_puppi_pruned_mass       ;
  std::vector<float>  	      jetAK8_puppi_softdrop_mass     ;
  std::vector<float>  	      jetAK8_puppi_pruned_massCorr   ;
  std::vector<float>  	      jetAK8_puppi_softdrop_massCorr ;
  std::vector<float>  	      jetAK8_puppi_pruned_jec        ;
  std::vector<float>  	      jetAK8_puppi_softdrop_jec      ;
    
  /** AK8 jets pruned */
  //int 		     njetsAK8_pruned	   ;
  //std::vector<float>       jetAK8_pruned_pt	   ;
  //std::vector<float>       jetAK8_pruned_eta     ;
  //std::vector<float>       jetAK8_pruned_mass    ;
  //std::vector<float>       jetAK8_pruned_phi     ;
  //std::vector<float>       jetAK8_pruned_e	   ;
  //std::vector<int  >       jetAK8_pruned_charge  ;
  //std::vector<int  >       jetAK8_pruned_flavour ;
  //std::vector<float>       jetAK8_pruned_ssv     ;
  //std::vector<float>       jetAK8_pruned_csv     ;
  //std::vector<float>       jetAK8_pruned_tchp    ;
  //std::vector<float>       jetAK8_pruned_tche    ;
  //std::vector<float>       jetAK8_pruned_jp	   ;
  //std::vector<float>       jetAK8_pruned_jbp     ;
  //std::vector<int  >       jetAK8_pruned_nSVs    ;
	
  /** AK8 jets softdrop */
  //int 	             njetsAK8_softdrop         ;
  //std::vector<float>	     jetAK8_softdrop_pt        ;
  //std::vector<float>	     jetAK8_softdrop_eta       ;
  //std::vector<float>	     jetAK8_softdrop_mass      ;
  //std::vector<float>	     jetAK8_softdrop_phi       ;
  //std::vector<float>	     jetAK8_softdrop_e	       ;
  //std::vector<int  >	     jetAK8_softdrop_charge    ;
  //std::vector<int  >	     jetAK8_softdrop_flavour   ;
  //std::vector<float>	     jetAK8_softdrop_ssv       ;
  //std::vector<float>	     jetAK8_softdrop_csv       ;
  //std::vector<float>	     jetAK8_softdrop_tchp      ;
  //std::vector<float>	     jetAK8_softdrop_tche      ;
  //std::vector<float>	     jetAK8_softdrop_jp        ;
  //std::vector<float>	     jetAK8_softdrop_jbp       ;
  //std::vector<int  >       jetAK8_softdrop_nSVs      ;
    
  /** pruned AK8 subjets  */
  std::vector<int>    	            subjetAK8_pruned_N      ;
  std::vector< std::vector<float> > subjetAK8_pruned_pt     ;
  std::vector< std::vector<float> > subjetAK8_pruned_eta    ;
  std::vector< std::vector<float> > subjetAK8_pruned_mass   ;
  std::vector< std::vector<float> > subjetAK8_pruned_phi    ;
  std::vector< std::vector<float> > subjetAK8_pruned_e      ;
  std::vector< std::vector<int  > > subjetAK8_pruned_charge ;
  std::vector< std::vector<int  > > subjetAK8_pruned_partonFlavour;
  std::vector< std::vector<int  > > subjetAK8_pruned_hadronFlavour;
//std::vector< std::vector<float> > subjetAK8_pruned_ssv    ;
  std::vector< std::vector<float> > subjetAK8_pruned_csv    ;    
//std::vector< std::vector<float> > subjetAK8_pruned_tchp   ;
//std::vector< std::vector<float> > subjetAK8_pruned_tche   ;
//std::vector< std::vector<float> > subjetAK8_pruned_jp     ;
//std::vector< std::vector<float> > subjetAK8_pruned_jbp    ;
	
  /** softdrop AK8 subjets */
  std::vector<int>                  subjetAK8_softdrop_N      ;
  std::vector< std::vector<float> > subjetAK8_softdrop_pt     ;
  std::vector< std::vector<float> > subjetAK8_softdrop_eta    ;
  std::vector< std::vector<float> > subjetAK8_softdrop_mass   ;
  std::vector< std::vector<float> > subjetAK8_softdrop_phi    ;
  std::vector< std::vector<float> > subjetAK8_softdrop_e      ;
  std::vector< std::vector<int  > > subjetAK8_softdrop_charge ;
  std::vector< std::vector<int  > > subjetAK8_softdrop_partonFlavour;
  std::vector< std::vector<int  > > subjetAK8_softdrop_hadronFlavour;
//std::vector< std::vector<float> > subjetAK8_softdrop_ssv    ;
  std::vector< std::vector<float> > subjetAK8_softdrop_csv    ;
//std::vector< std::vector<float> > subjetAK8_softdrop_tchp   ;
//std::vector< std::vector<float> > subjetAK8_softdrop_tche   ;
//std::vector< std::vector<float> > subjetAK8_softdrop_jp     ;
//std::vector< std::vector<float> > subjetAK8_softdrop_jbp    ;

  /** AK4 genJets*/
  int			      genJetAK4_N               ;
  std::vector<float>  	      genJetAK4_pt              ;
  std::vector<float>  	      genJetAK4_eta             ;
  std::vector<float>  	      genJetAK4_mass            ;
  std::vector<float>  	      genJetAK4_phi             ;
  std::vector<float>  	      genJetAK4_e               ;
  std::vector<float>  	      genJetNoNuAK4_pt          ;
  std::vector<float>  	      genJetNoNuAK4_mass        ;
  std::vector<float>  	      genJetNoNuAK4_e           ;
   
  /*-------------------------AK8 genJets---------------------------*/   
  int			      genJetAK8_N               ;

  std::vector<float>  	      genJetAK8_pt              ;
  std::vector<float>  	      genJetAK8_eta             ;
  std::vector<float>  	      genJetAK8_mass            ;
  std::vector<float>  	      genJetAK8_phi             ;
  std::vector<float>  	      genJetAK8_e               ;
  std::vector<float>  	      genJetAK8_prunedmass      ;
  std::vector<float>  	      genJetAK8_softdropmass    ;
  
  /** HLT trigger decisions */
  std::map<std::string,bool> HLT_isFired;
	 
  /** HLT trigger objects */
  std::vector<float>  		    triggerObject_pt	      ;
  std::vector<float>  		    triggerObject_eta	      ;
  std::vector<float>  		    triggerObject_phi	      ;
  std::vector<float>  		    triggerObject_mass	      ;
  std::vector< std::vector<float> > triggerObject_filterIDs   ; // as defined in http://cmslxr.fnal.gov/lxr/source/DataFormats/HLTReco/interface/TriggerTypeDefs.h
  std::vector< std::vector<int> >   triggerObject_firedTrigger; // as defined in plugins/TriggersNtuplizer.cc


  /** HLT filter decisions */
  bool passFilter_HBHE_;
  bool passFilter_CSCHalo_;
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
  
  /** MET */
  std::vector<float>                METraw_et		      ;	 
  std::vector<float>                METraw_phi		      ;
  std::vector<float>  	            METraw_sumEt	      ;
  std::vector<float>  	      	    MET_corrPx  	      ;
  std::vector<float>  	            MET_corrPy  	      ;
  std::vector<float>  	      	    MET_et		      ;
  std::vector<float>  	      	    MET_phi		      ;
  std::vector<float>  	      	    MET_sumEt		      ;
  std::vector<float>  	      	    MET_T1Uncertainty	      ;


  /*------------------------EVENT infos-------------------------*/    
  int                               EVENT_event            ;
  int                               EVENT_run              ;
  int                               EVENT_lumiBlock        ;

  /*--------------------------PV infos--------------------------*/
  int                               PV_N		     ;
  bool                              PV_filter		 ;
  std::vector<float>                PV_chi2          ;
  std::vector<float>                PV_ndof          ;
  std::vector<float>                PV_rho           ;
  std::vector<float>                PV_z             ;
    
  /*--------------------------PU infos--------------------------*/  			       
  std::vector<int  >                nPuVtxTrue             ;// the *true* mean number of the poisson distribution for this event from which the number of interactions each bunch crossing has been sampled 
  std::vector<int  >                nPuVtx                 ;// the number of pileup interactions that have been added to the event in the current bunch crossing
  std::vector<int  >                bX                     ;// to which bunch crossing do these interaction belong?  
  
private:
  TTree* tree_;

};

#endif // NtupleBranches_H
