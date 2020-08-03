#ifndef HELPER_H
#include "../interface/helper.h"
#endif

#ifndef JpsiTauNtuplizer_H
#define JpsiTauNtuplizer_H

#include <TRandom3.h>

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

   std::vector<std::string> _FFErrNames = {"delta_a0","delta_a1","delta_a2","delta_b0","delta_b1","delta_b2","delta_c1","delta_c2","delta_d0","delta_d1","delta_d2"};

   std::vector<double> _FFErr = {0.6850051095, 0.3674778905, 0.08740600666, 0.05468939568, 0.01128871118, 0.006722573912, 0.006537973692, 0.003156718233, 0.0004988336396, 0.0004454548237, 0.};

//   std::vector<double> _FFmean = {0.217179,
//				  0.173516,
//				  -0.0304393,
//				  0.0108183,
//				  -0.284005,
//				  0.00280893,
//				  -0.0107516,
//				  -0.0613763,
//				  0.202081,
//				  -0.127171,
//				  0.102528};
   std::vector<double> _FFmean = {
     0.418268,
     -0.181944,
     0.114036,
     0.090122,
     0.00745944,
     0.0127816,
     -0.00283544,
     0.0310275,
     -0.00132216,
     -0.00438412,
     0.
   };


   std::vector<double> Inv = {2.13114,
			      7.40517,
			      130.878,
			      334.174,
			      7846.72,
			      21999.4,
			      23177.8,
			      37494.1,
			      101982,
			      4217630,
			      10000000000};
			      
   TRandom3 *ran;			      



   //   std::vector<double> _FFErr = {1, 1,1,1,1,1,1,1,1,1,1};


};




#endif // JpsiTauNtuplizer_H

