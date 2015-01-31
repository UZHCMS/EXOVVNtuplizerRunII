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
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
//#include "TauAnalysis/CandidateTools/interface/NSVfitStandaloneAlgorithm.h"
#include "FWCore/Framework/interface/Event.h"
#include "TMatrixD.h"
#include "TMatrix.h"

#include "TTree.h"
#include "TLorentzVector.h"

/*here we declare the input and output variables*/

class NtupleBranches {

public:
   NtupleBranches( TTree* tree = 0 );
   ~NtupleBranches( void );
   
   void branch( void );
   void getEventByLabels( edm::EventBase const & event );
   void reset( void );
   void fillTree( void ){ tree_->Fill(); };

//=================================================================================================================== 
  /* output tree variables*/
    
    /*----------------------gen particles-------------------------*/
    float                               lheV_pt              ;
    float                               lheHT                ;
    float                               lheNj                ;
    std::vector<float>  	        genParticle_pt	     ;
    std::vector<float>  	        genParticle_px	     ;
    std::vector<float>  	        genParticle_py	     ;
    std::vector<float>  	        genParticle_pz	     ;
    std::vector<float>  	        genParticle_e	     ;
    std::vector<float>  	        genParticle_eta	     ;
    std::vector<float>  	        genParticle_phi	     ;
    std::vector<float>  	        genParticle_mass     ;
    std::vector<int  >  	        genParticle_pdgId    ;
    std::vector<int  >  	        genParticle_status   ;
    std::vector<int  >  	        genParticle_nDau     ;
    std::vector<int  >  	        genParticle_nMoth    ;
    std::vector<std::vector<int> >      genParticle_mother   ; 
    std::vector<std::vector<int> >      genParticle_dau      ;
        
    /*-------------------------leptons----------------------------*/
    int 	                      nlep		             ;
    std::vector<int  >  	      lep_type		         ;
    std::vector<float>  	      lep_charge	         ;
    std::vector<float>  	      lep_e 		         ;
    std::vector<float>  	      lep_eta		         ;
    std::vector<float>  	      lep_mass		         ;
    std::vector<float>  	      lep_pt		         ;
    std::vector<float>  	      lep_phi		         ;
    std::vector<int  >  	      lep_isHEEP	         ;
    std::vector<int  >                lep_isHighPtMuon       	 ;
    std::vector<float>  	      lep_pfRhoCorrRelIso03  	 ;
    std::vector<float>  	      lep_pfRhoCorrRelIso04  	 ;
    std::vector<float>  	      lep_pfDeltaCorrRelIso  	 ;
    std::vector<float>  	      lep_pfRelIso  	     	 ;
    std::vector<float>  	      lep_photonIso 	     	 ;
    std::vector<float>  	      lep_neutralHadIso	     	 ;
    std::vector<float>  	      lep_chargedHadIso	     	 ;
    std::vector<float>  	      lep_trackIso	         ;    
        
    /*-------------------------AK5 jets---------------------------*/   
    int 	        	      njetsAK4               ;
    std::vector<float>  	      jetAK4_pt              ;
    std::vector<float>  	      jetAK4_eta             ;
    std::vector<float>  	      jetAK4_mass            ;
    std::vector<float>  	      jetAK4_phi             ;
    std::vector<float>  	      jetAK4_e               ;
    // std::vector<float>  	      jetAK4_jec             ;
    // std::vector<float>  	      jetAK4_jecUp           ;
    // std::vector<float>  	      jetAK4_jecDown         ; 
    std::vector<bool  >  	      jetAK4_IDLoose	     ;
    std::vector<int  >                jetAK4_cm              ;
    std::vector<int  >                jetAK4_nm              ;
    std::vector<float>                jetAK4_che             ;
    std::vector<float>                jetAK4_ne              ;
    std::vector<int  >  	      jetAK4_charge	     ;
    std::vector<int  >  	      jetAK4_flavour	     ;
    std::vector<float>  	      jetAK4_ssv 	     ;
    std::vector<float>  	      jetAK4_csv 	     ;         
    std::vector<float>  	      jetAK4_cisv 	     ;         
    std::vector<float>  	      jetAK4_tchp            ;
    std::vector<float>  	      jetAK4_tche            ;
    std::vector<float>  	      jetAK4_jp              ;
    std::vector<float>  	      jetAK4_jbp             ;
    std::vector<float  >  	      jetAK4_vtxMass            ;
	std::vector<float  >  	      jetAK4_vtxNtracks            ;
	std::vector<float  >  	      jetAK4_vtx3DVal            ;
	std::vector<float  >  	      jetAK4_vtx3DSig            ; 
    //std::vector<int  >                jetAK4_nSVs            ;
    std::vector<int  >  	      jetAK4_nconstituents   ;     

    /*-------------------------CA8 jets---------------------------*/
    int 	        	      	  njetsAK8               ;
    std::vector<float>  	      jetAK8_pt              ;
    std::vector<float>  	      jetAK8_eta             ;
    std::vector<float>  	      jetAK8_mass            ;
    std::vector<float>  	      jetAK8_phi             ;
    std::vector<float>  	      jetAK8_e               ;
    // std::vector<float>  	      jetAK8_jec             ;
    // std::vector<float>  	      jetAK8_jecUp           ;
    // std::vector<float>  	      jetAK8_jecDown         ;
    std::vector<bool >  	      jetAK8_IDLoose	     ;
    std::vector<int  >            jetAK8_cm         	 ;
    std::vector<int  >            jetAK8_nm         	 ;
    std::vector<float>            jetAK8_che       		 ;
    std::vector<float>            jetAK8_ne         	 ;
    std::vector<int  >  	      jetAK8_charge	     	 ;
    std::vector<int  >  	      jetAK8_flavour	     ;
    std::vector<float>  	      jetAK8_ssv 	     	 ;
    std::vector<float>  	      jetAK8_csv 	     	 ;    
    std::vector<float>  	      jetAK8_tchp            ;
    std::vector<float>  	      jetAK8_tche            ;
    std::vector<float>  	      jetAK8_jp              ;
    std::vector<float>  	      jetAK8_jbp             ;
    std::vector<float>  	      jetAK8_tau1            ;
    std::vector<float>  	      jetAK8_tau2            ; 
    std::vector<float>  	      jetAK8_tau3            ;    
    std::vector<int  >  	      jetAK8_nconstituents   ;   
	std::vector<float>  	      jetAK8pruned_mass      ;
	std::vector<float>  	      jetAK8trimmed_mass     ;
	std::vector<float>  	      jetAK8filtered_mass    ;  
	std::vector<float>  	      jetAK8_nSubJets        ;
    // /*----------------------CA8 jets pruned-----------------------*/
//     int 	        	      njetsAK8pruned         ;
//     std::vector<float>  	      jetAK8pruned_pt        ;
//     std::vector<float>  	      jetAK8pruned_eta       ;
//     std::vector<float>  	      jetAK8pruned_mass      ;
//     std::vector<float>  	      jetAK8pruned_phi       ;
//     std::vector<float>  	      jetAK8pruned_e         ;
//     std::vector<float>  	      jetAK8pruned_jec       ;
//     std::vector<float>  	      jetAK8pruned_jecUp     ;
//     std::vector<float>  	      jetAK8pruned_jecDown   ;
//     std::vector<int  >  	      jetAK8pruned_charge    ;
//     std::vector<int  >  	      jetAK8pruned_flavour   ;
//     std::vector<float>  	      jetAK8pruned_ssv 	     ;
//     std::vector<float>  	      jetAK8pruned_csv 	     ;
//     std::vector<float>  	      jetAK8pruned_tchp      ;
//     std::vector<float>  	      jetAK8pruned_tche      ;
//     std::vector<float>  	      jetAK8pruned_jp        ;
//     std::vector<float>  	      jetAK8pruned_jbp       ;
//     std::vector<int  >  	      jetAK8pruned_nSVs      ;
    
    /*----------------------CA8 subjets---------------------------*/
    //std::vector<int>    	      nsubjets               ;
    std::vector< std::vector<float> > subjetAK8pruned_pt     ;
    std::vector< std::vector<float> > subjetAK8pruned_eta    ;
    std::vector< std::vector<float> > subjetAK8pruned_mass   ;
    std::vector< std::vector<float> > subjetAK8pruned_phi    ;
    std::vector< std::vector<float> > subjetAK8pruned_e      ;
    std::vector< std::vector<int  > > subjetAK8pruned_charge ;
    std::vector< std::vector<int  > > subjetAK8pruned_flavour;
    std::vector< std::vector<float> > subjetAK8pruned_ssv    ;
    std::vector< std::vector<float> > subjetAK8pruned_csv    ;    
    std::vector< std::vector<float> > subjetAK8pruned_tchp   ;
    std::vector< std::vector<float> > subjetAK8pruned_tche   ;
    std::vector< std::vector<float> > subjetAK8pruned_jp     ;
    std::vector< std::vector<float> > subjetAK8pruned_jbp    ;
   
    /*---------------------HLT triggers---------------------------*/    
    bool                              isFired_HLT_Mu40_eta2p1_v9 ;
    bool                              isFired_HLT_Mu40_eta2p1_v10;
    bool                              isFired_HLT_Mu40_eta2p1_v11;    
    bool                              isFired_HLT_Ele80_CaloIdVT_TrkIdT_v8 ;
    bool                              isFired_HLT_Ele80_CaloIdVT_TrkIdT_v9 ;
    bool                              isFired_HLT_Ele80_CaloIdVT_TrkIdT_v10;
    bool                              isFired_HLT_Ele80_CaloIdVT_GsfTrkIdT_v1 ;
    bool                              isFired_HLT_Ele80_CaloIdVT_GsfTrkIdT_v2 ;
    bool			      isFired_HLT_PFHT650_v5;
    bool			      isFired_HLT_PFHT650_v6;
    bool			      isFired_HLT_PFHT650_v7;
    bool			      isFired_HLT_PFHT650_v8;
    bool			      isFired_HLT_PFHT650_v9;
    bool			      isFired_HLT_PFNoPUHT650_v1;
    bool			      isFired_HLT_PFNoPUHT650_v3;
    bool			      isFired_HLT_PFNoPUHT650_v4;
    bool			      isFired_HLT_PFJet320_v3;
    bool			      isFired_HLT_PFJet320_v4;
    bool			      isFired_HLT_PFJet320_v5;
    bool			      isFired_HLT_PFJet320_v6;
    bool			      isFired_HLT_PFJet320_v8;
    bool			      isFired_HLT_PFJet320_v9;			      
    
    /*-------------------------MET--------------------------------*/
    std::vector<float>                METraw_et              ;	 
    std::vector<float>                METraw_phi             ;     
    std::vector<float>  	      MET_et                 ;
    std::vector<float>  	      MET_phi                ;

    /*------------------------EVENT infos-------------------------*/    
    int                               EVENT_event            ;
    int                               EVENT_run              ;
    int                               EVENT_lumiBlock        ;

    /*--------------------------PV infos--------------------------*/
    int                               nPVs		             ;
    
    /*--------------------------PU infos--------------------------*/  			       
    std::vector<int  >                nPuVtxTrue             ;// the *true* mean number of the poisson distribution for this event from which the number of interactions each bunch crossing has been sampled 
    std::vector<int  >                nPuVtx                 ;// the number of pileup interactions that have been added to the event in the current bunch crossing
    std::vector<int  >                bX                     ;// to which bunch crossing do these interaction belong?  
  
private:
   TTree* tree_;

};

#endif // NtupleBranches_H
