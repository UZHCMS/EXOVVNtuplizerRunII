#ifndef VerticesNtuplizer_H
#define VerticesNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class VerticesNtuplizer : public CandidateNtuplizer {

public:
  VerticesNtuplizer( std::vector<edm::EDGetTokenT<reco::VertexCollection>> tokens,
		     edm::EDGetTokenT<reco::BeamSpot>             beamToken , 
		     NtupleBranches* nBranches,
		     std::map< std::string, bool >&  runFlags  );
  ~VerticesNtuplizer( void );
  
  bool fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
  
private:
   edm::EDGetTokenT<reco::VertexCollection> vtxToken_   ;
   edm::EDGetTokenT<reco::BeamSpot> beamToken_   ;

   edm::Handle< reco::VertexCollection >  vertices_;
   edm::Handle< reco::BeamSpot >  beamSpot_;

};

#endif // VerticesNtuplizer_H
