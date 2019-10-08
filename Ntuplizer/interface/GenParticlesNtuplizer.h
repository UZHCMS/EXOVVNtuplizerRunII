#ifndef GenParticlesNtuplizer_H
#define GenParticlesNtuplizer_H

#include "../interface/CandidateNtuplizer.h"
#include <algorithm>
#include <vector>


class GenParticlesNtuplizer : public CandidateNtuplizer {

public:
  GenParticlesNtuplizer( std::vector<edm::EDGetTokenT<reco::GenParticleCollection>> tokens, 
			 NtupleBranches* nBranches, std::map< std::string, bool >&  runFlags );

  ~GenParticlesNtuplizer( void ); 

  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup  );

private:
   edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
   edm::Handle< reco::GenParticleCollection >  genParticles_;
   bool isJpsiMu_;
   bool isJpsiEle_;
   bool doGenHist_;
};

#endif // GenParticlesNtuplizer_H
