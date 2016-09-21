#ifndef METsNtuplizer_H
#define METsNtuplizer_H

#include "../interface/CandidateNtuplizer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
class TFormula;

class METsNtuplizer : public CandidateNtuplizer {

public:
   METsNtuplizer( edm::EDGetTokenT<pat::METCollection>     mettoken    , 
		  edm::EDGetTokenT<pat::METCollection>     metpuppitoken    , 
		  edm::EDGetTokenT<pat::METCollection>     metmvatoken    , 
 	 	  edm::EDGetTokenT<pat::JetCollection>	   jettoken    ,
		  edm::EDGetTokenT<pat::MuonCollection>    muontoken   ,
		  edm::EDGetTokenT<double>		   rhotoken    ,
		  edm::EDGetTokenT<reco::VertexCollection> vtxtoken    ,
		  edm::EDGetTokenT<double>		   metSigtoken    ,
		  edm::EDGetTokenT<math::Error<2>::type>    metCovtoken ,
		  std::vector<std::string>		   jecAK4labels,
		  std::vector<std::string>		   corrformulas,
		  NtupleBranches*			   nBranches   ,
		  std::map< std::string, bool >&                       runFlags );
   ~METsNtuplizer( void );
   
   void 	fillBranches      ( edm::Event const & event, const edm::EventSetup& iSetup );
   void   	addTypeICorr      ( edm::Event const & event );
   double 	getJEC            ( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ );
   double 	getJECOffset      ( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ );
   void   	initJetCorrFactors( void );
   
   private:
    edm::EDGetTokenT<pat::METCollection>	 metInputToken_  ;
    edm::EDGetTokenT<pat::METCollection>	 metpuppiInputToken_  ;
    edm::EDGetTokenT<pat::METCollection>	 metmvaInputToken_  ;
    edm::EDGetTokenT<pat::JetCollection>	 jetInputToken_  ; 
    edm::EDGetTokenT<pat::MuonCollection>	 muonInputToken_ ;
    edm::EDGetTokenT<double>			 rhoToken_	 ;  
    edm::EDGetTokenT<reco::VertexCollection>	 verticeToken_   ;   
    edm::EDGetTokenT<double>			 metSigToken_	 ;  
    edm::EDGetTokenT<math::Error<2>::type>	 metCovToken_   ;   

    std::vector<std::string>			 jetCorrLabel_   ;
    std::vector<std::string>			 corrFormulas_   ;

    std::vector<std::string>			 offsetCorrLabel_;   
    boost::shared_ptr<FactorizedJetCorrector>	 jecAK4_	 ;
    boost::shared_ptr<FactorizedJetCorrector>	 jecOffset_	 ;
    
    edm::Handle<pat::METCollection>		 METs_  	 ;
    edm::Handle<pat::METCollection>		 METspuppi_  	 ;
    edm::Handle<pat::METCollection>		 METsmva_  	 ;
    edm::Handle<pat::MuonCollection>		 muons_ 	 ;
    edm::Handle<pat::JetCollection>		 jets_  	 ;
    edm::Handle< double >			 rho_		 ;
    edm::Handle<reco::VertexCollection> 	 vertices_	 ;

    edm::Handle<double>                          significanceHandle;
    edm::Handle<math::Error<2>::type>            covHandle;
 
   
    std::map<std::string,double>		 TypeICorrMap_   ;
    
    TFormula* corrPx_;
    TFormula* corrPy_;

    bool doCorrOnTheFly_;
    bool doMETSVFIT_;
    bool doMVAMET_;

};

#endif // METsNtuplizer_H
