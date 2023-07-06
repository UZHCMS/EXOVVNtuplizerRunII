#ifndef HELPER_H
#include "../interface/helper.h"
#endif

#ifndef ZTauTauNtuplizer_H
#define ZTauTauNtuplizer_H

#include <TRandom3.h>


class ZTauTauNtuplizer : public CandidateNtuplizer {


 public:
  ZTauTauNtuplizer( edm::EDGetTokenT<pat::MuonCollection>    muonToken   , 
		    edm::EDGetTokenT<pat::ElectronCollection>    electronToken   , 
		    edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
		    edm::EDGetTokenT<reco::BeamSpot>             beamToken , 
		    edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken,
		    edm::EDGetTokenT<edm::TriggerResults> triggertoken,
		    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobject,
		    edm::EDGetTokenT<reco::GenParticleCollection> genptoken, 
		    std::map< std::string, bool >& runFlags,
		    std::map< std::string, double >& runValues,
		    std::map< std::string, std::string >& runStrings,
                    NtupleBranches* nBranches);

  ~ZTauTauNtuplizer( void );
  

  bool fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
  const reco::Candidate*  checkMom(const reco::Candidate * candMom);
  FreeTrajectoryState initialFreeState(const reco::Track& tk, const MagneticField *field);

private:
   edm::EDGetTokenT<pat::MuonCollection>    muonToken_   ;
   edm::EDGetTokenT<pat::ElectronCollection>    electronToken_   ;
   edm::EDGetTokenT<reco::VertexCollection> verticeToken_   ;
   edm::EDGetTokenT<reco::BeamSpot> bsToken_   ;
   edm::EDGetTokenT<pat::PackedCandidateCollection>   		packedpfcandidatesToken_;
   edm::EDGetTokenT<edm::TriggerResults> 		     HLTtriggersToken_;
   edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection>  triggerObjects_;
   edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
   //   edm::EDGetTokenT<std::vector<reco::GenJet>> genTauToken_;

   edm::Handle<pat::MuonCollection>      		       muons_		       ;
   edm::Handle<pat::ElectronCollection>      		       electrons_		       ;
   edm::Handle< reco::VertexCollection >  vertices_;
   edm::Handle< reco::BeamSpot >  beamspot_;
   edm::Handle< std::vector<pat::PackedCandidate> > packedpfcandidates_   ;
   edm::Handle< edm::TriggerResults> 			     HLTtriggers_;
   edm::Handle<pat::TriggerObjectStandAloneCollection>	     triggerObjects;
   edm::Handle< reco::GenParticleCollection >  genParticles_;
   //   edm::Handle< std::vector<reco::GenJet> >  genTaus_;

   edm::ESHandle<TransientTrackBuilder> builder;

   const MagneticField                 *fMagneticField;


   bool runOnMC_;   
   bool isTruth_;
   bool verbose_;

  
   bool flag_fill = false;

   float c_dz;
   float c_fsig;
   float c_vprob;
   unsigned int c_charge;

   float chi = 0.;
   float ndf = 0.;

   int globalCounter = 0;

   const int numberofDNN = 14;

   const int numberofToys = 1000;
   std::vector<map<string, double>> FFdict;

   helper aux;




};




#endif // ZTauTauNtuplizer_H

