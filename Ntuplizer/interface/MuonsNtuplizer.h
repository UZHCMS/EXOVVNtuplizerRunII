#ifndef MuonsNtuplizer_H
#define MuonsNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class MuonsNtuplizer : public CandidateNtuplizer {

public:
   MuonsNtuplizer( edm::EDGetTokenT<pat::MuonCollection>    muonToken   , 
                   edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
		   edm::EDGetTokenT<double>		    rhoToken	, 
		   edm::EDGetTokenT<pat::TauCollection>     boostedtauToken  ,
		   NtupleBranches*                          nBranches   ,
		   std::map< std::string, bool >&           runFlags
		   );
   ~MuonsNtuplizer( void );

   void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );

private:
   edm::EDGetTokenT<pat::MuonCollection>    muonToken_   ;
   edm::EDGetTokenT<reco::VertexCollection> verticeToken_;
   edm::EDGetTokenT<double> 		    rhoToken_    ;
   edm::EDGetTokenT<pat::TauCollection>     boostedtauToken_  ; 
  
   edm::Handle<pat::MuonCollection> 	    muons_       ;
   edm::Handle<reco::VertexCollection> 	    vertices_    ;
   edm::Handle<double> 			    rho_         ;
   edm::Handle<pat::TauCollection> 	    taus_        ;
   bool doBoostedTaus_;
   bool isJpsiMu_;
   bool isJpsiEle_;
};

#endif // MuonsNtuplizer_H
