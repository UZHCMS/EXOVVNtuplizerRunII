#ifndef MuonsNtuplizer_H
#define MuonsNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class MuonsNtuplizer : public CandidateNtuplizer {

public:
   MuonsNtuplizer( edm::EDGetTokenT<pat::MuonCollection> muonToken, edm::EDGetTokenT<reco::VertexCollection> verticeToken,  NtupleBranches* nBranches );
   ~MuonsNtuplizer( void );

   void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );

private:
   edm::EDGetTokenT<pat::MuonCollection> 	muonsToken_   ;
   edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
      
   edm::Handle<pat::MuonCollection> 		muons_;
   edm::Handle<reco::VertexCollection> 		vertices_;

};

#endif // MuonsNtuplizer_H
