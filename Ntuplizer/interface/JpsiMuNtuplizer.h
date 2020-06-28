#ifndef HELPER_H
#include "../interface/helper.h"
#endif

#ifndef JpsiMuNtuplizer_H
#define JpsiMuNtuplizer_H

class JpsiMuNtuplizer : public CandidateNtuplizer {


 public:
  JpsiMuNtuplizer( edm::EDGetTokenT<pat::MuonCollection>    muonToken   , 
		   edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
		   edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken,
		   edm::EDGetTokenT<edm::TriggerResults> triggertoken,
		   edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobject,
		   edm::EDGetTokenT<reco::GenParticleCollection> genptoken, 
		   std::map< std::string, bool >& runFlags,
		   NtupleBranches* nBranches );

  ~JpsiMuNtuplizer( void );


  bool fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
  
private:
   edm::EDGetTokenT<pat::MuonCollection>    muonToken_   ;
   edm::EDGetTokenT<reco::VertexCollection> verticeToken_   ;
   edm::EDGetTokenT<reco::BeamSpot> bsToken_   ;
   edm::EDGetTokenT<pat::PackedCandidateCollection>   		packedpfcandidatesToken_;
   edm::EDGetTokenT<edm::TriggerResults> 		     HLTtriggersToken_;
   edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection>  triggerObjects_;
   edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;

   edm::Handle<pat::MuonCollection>      		       muons_		       ;
   edm::Handle< reco::VertexCollection >  vertices_;
   edm::Handle< reco::BeamSpot >  beamspot_;
   edm::Handle< std::vector<pat::PackedCandidate> > packedpfcandidates_   ;
   edm::Handle< edm::TriggerResults> 			     HLTtriggers_;
   edm::Handle<pat::TriggerObjectStandAloneCollection>	     triggerObjects;
   edm::Handle< reco::GenParticleCollection >  genParticles_;

   edm::ESHandle<TransientTrackBuilder> builder;

   const MagneticField                 *fMagneticField;

   float chi = 0.;
   float ndf = 0.;

   helper aux;

   bool runOnMC_;   
   bool useHammer_;   
   bool verbose_;

   Hammer::Hammer hammer;

   std::vector<std::string> _FFErrNames = {"delta_a0","delta_a1","delta_a2","delta_b0","delta_b1","delta_b2","delta_c1","delta_c2","delta_d0","delta_d1","delta_d2"};

   std::vector<double> _FFErr = {0.0009, 0.03, 0.6, 0.0005, 0.02, 0.6, 0.003, 0.08, 0.1, 0.16, 0.009};
   
};




#endif // JpsiMuNtuplizer_H

