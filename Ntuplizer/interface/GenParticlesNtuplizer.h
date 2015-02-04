#ifndef GenParticlesNtuplizer_H
#define GenParticlesNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class GenParticlesNtuplizer : public CandidateNtuplizer {

public:
  GenParticlesNtuplizer( std::vector<edm::EDGetTokenT<reco::GenParticleCollection>> tokens, NtupleBranches* nBranches );
  ~GenParticlesNtuplizer( void ); 

  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );

private:
   edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
   edm::Handle< reco::GenParticleCollection >  genParticles_;
      
};

#endif // GenParticlesNtuplizer_H
