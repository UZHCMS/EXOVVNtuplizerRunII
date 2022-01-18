#ifndef HELPER_H
#include "../interface/helper.h"
#endif

#ifndef JpsiKNtuplizerE_H
#define JpsiKNtuplizerE_H

#include <TRandom3.h>


class JpsiKNtuplizerE : public CandidateNtuplizer {


 public:
  JpsiKNtuplizerE( edm::EDGetTokenT<pat::MuonCollection>    muonToken   , 
		    edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
		    edm::EDGetTokenT<reco::BeamSpot>             beamToken , 
		    edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken,
		    edm::EDGetTokenT<edm::TriggerResults> triggertoken,
		    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobject,
		    edm::EDGetTokenT<reco::GenParticleCollection> genptoken, 
		    edm::EDGetTokenT<pat::PackedGenParticleCollection> packedgenptoken, 
		    std::map< std::string, bool >& runFlags,
		    std::map< std::string, double >& runValues,
		  std::map< std::string, std::string >& runStrings, NtupleBranches* nBranches , TH1F* fileweight );


  ~JpsiKNtuplizerE( void );
  
  bool fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
  const reco::Candidate*  checkMom(const reco::Candidate * candMom);

private:
   edm::EDGetTokenT<pat::MuonCollection>    muonToken_   ;
   edm::EDGetTokenT<reco::VertexCollection> verticeToken_   ;
   edm::EDGetTokenT<reco::BeamSpot> bsToken_   ;
   edm::EDGetTokenT<pat::PackedCandidateCollection>   		packedpfcandidatesToken_;
   edm::EDGetTokenT<edm::TriggerResults> 		     HLTtriggersToken_;
   edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection>  triggerObjects_;
   edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
   edm::EDGetTokenT<pat::PackedGenParticleCollection> packedgenParticlesToken_;
   //   edm::EDGetTokenT<std::vector<reco::GenJet>> genTauToken_;

   edm::Handle<pat::MuonCollection>      		       muons_		       ;
   edm::Handle< reco::VertexCollection >  vertices_;
   edm::Handle< reco::BeamSpot >  beamspot_;
   edm::Handle< std::vector<pat::PackedCandidate> > packedpfcandidates_   ;
   edm::Handle< edm::TriggerResults> 			     HLTtriggers_;
   edm::Handle<pat::TriggerObjectStandAloneCollection>	     triggerObjects;
   edm::Handle< reco::GenParticleCollection >  genParticles_;
   edm::Handle< std::vector<pat::PackedGenParticle> >  packedgenParticles_;
   //   edm::Handle< std::vector<reco::GenJet> >  genTaus_;

   edm::ESHandle<TransientTrackBuilder> builder;

   const MagneticField                 *fMagneticField;


   bool runOnMC_;   
   bool isTruth_;
   bool verbose_;

   float c_dz;

   float chi = 0.;
   float ndf = 0.;

   helper aux;

   bool isBkgBSample_;
   TH1F* histGenWeights_ ;
   float genWeightBkgB_ =1;
   int motherID_ =0;

};




#endif // JpsiKNtuplizerE_H

