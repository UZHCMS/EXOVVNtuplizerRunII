#ifndef ElectronsNtuplizer_H
#define ElectronsNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class ElectronsNtuplizer : public CandidateNtuplizer {

public:
  ElectronsNtuplizer( edm::EDGetTokenT<pat::ElectronCollection> electronToken, edm::EDGetTokenT<reco::VertexCollection> verticesToken, NtupleBranches* nBranches );
  ~ElectronsNtuplizer( void );
  
  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
  
private:
    edm::EDGetTokenT<pat::ElectronCollection> 	electronToken_   ;
    edm::EDGetTokenT<reco::VertexCollection> verticesToken_;
   
   edm::Handle<pat::ElectronCollection> 	electrons_;
   edm::Handle<reco::VertexCollection> 		vertices_;

      
};

#endif // ElectronsNtuplizer_H
