#ifndef HELPER_H
#include "../interface/helper.h"
#endif

#ifndef JpsiTauNtuplizer_H
#define JpsiTauNtuplizer_H


class JpsiTauNtuplizer : public CandidateNtuplizer {


 public:
  JpsiTauNtuplizer( edm::EDGetTokenT<pat::MuonCollection>    muonToken   , 
		    edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
		    edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken,
		    edm::EDGetTokenT<edm::TriggerResults> triggertoken,
		    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobject,
		    edm::EDGetTokenT<reco::GenParticleCollection> genptoken, 
		    edm::EDGetTokenT<std::vector<reco::GenJet>> genttoken,
		    std::map< std::string, bool >& runFlags,
		    std::map< std::string, double >& runValues,
		    std::map< std::string, std::string >& runStrings,
		    NtupleBranches* nBranches );

  ~JpsiTauNtuplizer( void );
  

  bool fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );

  
private:
   edm::EDGetTokenT<pat::MuonCollection>    muonToken_   ;
   edm::EDGetTokenT<reco::VertexCollection> verticeToken_   ;
   edm::EDGetTokenT<reco::BeamSpot> bsToken_   ;
   edm::EDGetTokenT<pat::PackedCandidateCollection>   		packedpfcandidatesToken_;
   edm::EDGetTokenT<edm::TriggerResults> 		     HLTtriggersToken_;
   edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection>  triggerObjects_;
   edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
   edm::EDGetTokenT<std::vector<reco::GenJet>> genTauToken_;

   edm::Handle<pat::MuonCollection>      		       muons_		       ;
   edm::Handle< reco::VertexCollection >  vertices_;
   edm::Handle< reco::BeamSpot >  beamspot_;
   edm::Handle< std::vector<pat::PackedCandidate> > packedpfcandidates_   ;
   edm::Handle< edm::TriggerResults> 			     HLTtriggers_;
   edm::Handle<pat::TriggerObjectStandAloneCollection>	     triggerObjects;
   edm::Handle< reco::GenParticleCollection >  genParticles_;
   edm::Handle< std::vector<reco::GenJet> >  genTaus_;

   edm::ESHandle<TransientTrackBuilder> builder;

   const MagneticField                 *fMagneticField;


   bool runOnMC_;   
   bool useDNN_;
   bool useHammer_;
   bool isTruth_;
   bool verbose_;

   float c_dz;
   float c_fsig;
   float c_vprob;
   float c_dnn;
   unsigned int c_charge;

   float chi = 0.;
   float ndf = 0.;

   ParticleMass muon_mass = 0.1056583;
   ParticleMass jpsi_mass = 3.09687;
   ParticleMass pion_mass = 0.139571;
   ParticleMass kaon_mass = 0.493677;
   
   float muon_sigma = 0.0000001;
   float jp_m_sigma = 0.00004;
   float pion_sigma = 0.000016;
   float kaon_sigma = 0.000016;

   Float_t mass_B0 = 5.27963;
   helper aux;
   
   tensorflow::MetaGraphDef* graphDef;
   tensorflow::Session* session;
   tensorflow::Tensor data; // (tensorflow::DT_FLOAT, { 1, 50, 8 }); // single batch of dimension 10
   tensorflow::Tensor label; // (tensorflow::DT_FLOAT, { 1,50}); 
   tensorflow::Tensor add_global; //(tensorflow::DT_FLOAT, { 1, 2 }); 
   tensorflow::Tensor isTraining; //(tensorflow::DT_FLOAT, { 1, 2 }); 
   tensorflow::Tensor norm; //(tensorflow::DT_FLOAT, { 1, 2 }); 

   std::string dnnfile_;

   Hammer::Hammer hammer;


   Hammer::Process Bc2JpsiLNu;
   
   int idxB = -1;
   std::vector<size_t> Bvtx_idxs;
   int idxTau = -1;
   std::vector<size_t> Tauvtx_idxs;
   int idxJpsi = -1;
   std::vector<size_t> Jpsivtx_idxs;
   

   std::vector<std::string> _FFErrNames = {"delta_a0","delta_a1","delta_a2","delta_b0","delta_b1","delta_b2","delta_c1","delta_c2","delta_d0","delta_d1","delta_d2"};

   std::vector<double> _FFErr = {0.0009, 0.03, 0.6, 0.0005, 0.02, 0.6, 0.003, 0.08, 0.1, 0.16, 0.009};

};




#endif // JpsiTauNtuplizer_H

