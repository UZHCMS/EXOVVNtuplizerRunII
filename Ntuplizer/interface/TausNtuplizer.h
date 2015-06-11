#ifndef TausNtuplizer_H
#define TausNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class TausNtuplizer : public CandidateNtuplizer {

public:
  TausNtuplizer( edm::EDGetTokenT<pat::TauCollection> tauToken, edm::EDGetTokenT<pat::TauCollection> tauEleTauToken , edm::EDGetTokenT<pat::TauCollection> tauMuTauToken , edm::EDGetTokenT<double> rhoToken, edm::EDGetTokenT<reco::VertexCollection> verticeToken,  NtupleBranches* nBranches );
   ~TausNtuplizer( void );
   
   void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
   
 private:
   
   
   edm::EDGetTokenT<pat::TauCollection> 			tauInputToken_     		;
   edm::EDGetTokenT<pat::TauCollection> 			tauEleTauInputToken_  		;
   edm::EDGetTokenT<pat::TauCollection> 			tauMuTauInputToken_    		;
   edm::EDGetTokenT<double> 		   	      		rhoToken_               	;
   edm::EDGetTokenT<reco::VertexCollection>  			verticeToken_   	 	;
   
  
   edm::Handle< std::vector<pat::Tau> > taus_   ;
   edm::Handle< std::vector<pat::Tau> > eleTaus_;
   edm::Handle< std::vector<pat::Tau> > muTaus_ ;
   edm::Handle< double >                rho_     ;
   edm::Handle<reco::VertexCollection> 	     			vertices_     			;

};

#endif // TausNtuplizer_H
