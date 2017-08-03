#ifndef GenParticlesNtuplizer_H
#define GenParticlesNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class GenParticlesNtuplizer : public CandidateNtuplizer {

public:
  GenParticlesNtuplizer( std::vector<edm::EDGetTokenT<reco::GenParticleCollection>> tokens, 
			 std::vector<edm::EDGetTokenT<bool>> tauSpinnerTokens_bool,
			 std::vector<edm::EDGetTokenT<double>> tauSpinnerTokens_double,
			 NtupleBranches* nBranches );

  ~GenParticlesNtuplizer( void ); 

  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );

private:
   edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
   edm::EDGetTokenT<bool> tauSpinnerWTisValidToken_;
   edm::EDGetTokenT<double> tauSpinnerWTToken_;
   edm::EDGetTokenT<double> tauSpinnerWThminusToken_;
   edm::EDGetTokenT<double> tauSpinnerWThplusToken_;
   edm::EDGetTokenT<double> tauSpinnerTauPolFromZToken_;
   edm::EDGetTokenT<double> tauSpinnerWRightToken_;
   edm::EDGetTokenT<double> tauSpinnerWLeftToken_;
   edm::EDGetTokenT<double> tauSpinnerIsRightLeftToken_;
   


   edm::Handle< reco::GenParticleCollection >  genParticles_;
   edm::Handle< bool > tauSpinnerWTisValid_;
   edm::Handle< double > tauSpinnerWT_;
   edm::Handle< double > tauSpinnerWThminus_;
   edm::Handle< double > tauSpinnerWThplus_;
   edm::Handle< double > tauSpinnerTauPolFromZ_;
   edm::Handle< double > tauSpinnerWRight_;
   edm::Handle< double > tauSpinnerWLeft_;
   edm::Handle< double > tauSpinnerIsRightLeft_;

};

#endif // GenParticlesNtuplizer_H
