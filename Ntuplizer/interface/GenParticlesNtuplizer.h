#ifndef HELPER_H
#include "../interface/helper.h"
#endif

#ifndef GenParticlesNtuplizer_H
#define GenParticlesNtuplizer_H

#include "../interface/CandidateNtuplizer.h"
#include "PhysicsTools/HepMCCandAlgos/interface/GenParticlesHelper.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <algorithm>
#include <vector>


class GenParticlesNtuplizer : public CandidateNtuplizer {

public:
  GenParticlesNtuplizer( std::vector<edm::EDGetTokenT<reco::GenParticleCollection>> tokens, 
			 NtupleBranches* nBranches, std::map< std::string, bool >&  runFlags );

  ~GenParticlesNtuplizer( void ); 

  bool fillBranches( edm::Event const & event, const edm::EventSetup& iSetup  );

  void recursiveDaughters(size_t index, int rank, const reco::GenParticleCollection &src, std::vector<size_t> &allIndices, std::vector<int> &pdgs, std::vector<int> &layers, std::vector<float> &ppt, std::vector<float> &peta, std::vector<float> &pphi);

  

private:
   edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
   edm::Handle< reco::GenParticleCollection >  genParticles_;
   bool isJpsiMu_;
   bool isJpsiEle_;
   bool isJpsiTau_;
   bool doGenHist_;
   bool verbose_;

   helper aux;

};

#endif // GenParticlesNtuplizer_H
