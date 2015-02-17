#ifndef METsNtuplizer_H
#define METsNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class TFormula;

class METsNtuplizer : public CandidateNtuplizer {

public:
   METsNtuplizer( std::vector<edm::EDGetTokenT<pat::METCollection>> 				mettokens,
   	 							edm::EDGetTokenT<pat::JetCollection> 				jettoken,
  								edm::EDGetTokenT<pat::MuonCollection> 				muontoken,
  								std::vector<std::string> 							jecAK4labels,
  								std::vector<std::string>						 	corrformulas,
  								edm::EDGetTokenT<double> 							rhotoken,
  								edm::EDGetTokenT<reco::VertexCollection> 			vtxtoken,
  								NtupleBranches* 									nBranches );
   ~METsNtuplizer( void );
   
   void 	fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
   void   	addTypeICorr      ( edm::Event const & event );
   double 	getJEC            ( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ );
   double 	getJECOffset      ( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ );
   void   	initJetCorrFactors( void );
   
private:      
   //edm::InputTag METsRawLabel_;
   edm::EDGetTokenT<pat::METCollection> 	metInputToken_   	;
   edm::EDGetTokenT<pat::JetCollection> 	jetInputToken_   	; 
   edm::EDGetTokenT<pat::MuonCollection>	muonInputToken_  	;
   std::vector<std::string>	                jetCorrLabel_		;
   std::vector<std::string> 				corrFormulas_   	;
   edm::EDGetTokenT<double>  		    	rhoToken_   		;  
   edm::EDGetTokenT<reco::VertexCollection> verticeToken_    	;
   
   std::vector<std::string> 				offsetCorrLabel_	;
   
   boost::shared_ptr<FactorizedJetCorrector>   jecAK4_         	;
   boost::shared_ptr<FactorizedJetCorrector>   jecOffset_      	;
   
     
   //edm::Handle< std::vector<pat::MET> > METsRaw_;
   edm::Handle<pat::METCollection> 			METs_				;
   edm::Handle<pat::MuonCollection> 		muons_				;
   edm::Handle<pat::JetCollection>      	jets_      			;
   edm::Handle< double >                    rho_        		;
   edm::Handle<reco::VertexCollection> 	    vertices_     		;
   
   std::map<std::string,double> 			TypeICorrMap_   	;
  
  
   TFormula* corrPx_;
   TFormula* corrPy_;
};

#endif // METsNtuplizer_H
