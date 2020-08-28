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

   const vector<string> parName = {"a0", "a1", "a2", "b0", "b1", "b2", "c1", "c2", "d0", "d1", "d2"};


   std::vector<std::string> _FFErrNames = {"delta_a0","delta_a1","delta_a2","delta_b0","delta_b1","delta_b2","delta_c1","delta_c2","delta_d0","delta_d1","delta_d2"};



   //   std::vector<std::string> _FFErrNames = {"delta_a0","delta_a1","delta_a2","delta_b0","delta_b1","delta_b2","delta_c1","delta_c2","delta_d0","delta_d1","delta_d2"};

   //   std::vector<double> _FFErr = {0.0009, 0.03, 0.6, 0.0005, 0.02, 0.6, 0.003, 0.08, 0.1, 0.16, 0.009};


   const double eigVar[11][11][2] = {

     {{ 2.98497949e-04 -2.98497949e-04}, { 5.99984947e-03 -5.99984947e-03}, {-4.80751102e-01, 4.80751102e-01}, { 1.92613583e-04 -1.92613583e-04}, { 1.22449212e-02 -1.22449212e-02}, {-4.68739027e-01, 4.68739027e-01}, { 1.97549526e-03 -1.97549526e-03}, {-4.56398505e-02, 4.56398505e-02}, {-2.60006981e-03, 2.60006981e-03}, { 1.26908422e-01 -1.26908422e-01}, {-2.74245631e-03, 2.74245631e-03}},

     {{-4.71133920e-04, 4.71133920e-04}, { 1.32879279e-02 -1.32879279e-02}, { 2.58680658e-01 -2.58680658e-01}, { 2.34768517e-04 -2.34768517e-04}, {-4.84884101e-03, 4.84884101e-03}, {-2.59229003e-01, 2.59229003e-01}, { 5.89170830e-04 -5.89170830e-04}, { 8.31133427e-03 -8.31133427e-03}, { 3.51683085e-03 -3.51683085e-03}, { 2.53291393e-02 -2.53291393e-02}, {-1.03239034e-03, 1.03239034e-03}},

     {{ 3.21051860e-04 -3.21051860e-04}, {-1.00848686e-02, 1.00848686e-02}, {-6.07635408e-03, 6.07635408e-03}, { 6.96528654e-05 -6.96528654e-05}, {-1.02757874e-02, 1.02757874e-02}, {-1.50516666e-02, 1.50516666e-02}, {-1.76069353e-03, 1.76069353e-03}, {-1.64231134e-02, 1.64231134e-02}, {-8.83373368e-04, 8.83373368e-04}, {-8.29672082e-02, 8.29672082e-02}, { 3.43573034e-03 -3.43573034e-03}},

     {{-0.00024581, 0.00024581}, { 0.00835942 -0.00835942}, { 0.00406215 -0.00406215}, {-0.0001332,  0.0001332 }, { 0.00802829 -0.00802829}, { 0.00329927 -0.00329927}, {-0.00040278, 0.00040278}, {-0.0522203,  0.0522203 }, {-0.00674185, 0.00674185}, { 0.00747539 -0.00747539}, {-0.00082702, 0.00082702}},
     
     {{-7.17542818e-05, 7.17542818e-05}, { 9.22993678e-03 -9.22993678e-03}, {-3.65319756e-04, 3.65319756e-04}, {-1.65621772e-04, 1.65621772e-04}, { 4.92781842e-03 -4.92781842e-03}, {-1.53923981e-04, 1.53923981e-04}, {-1.74730474e-04, 1.74730474e-04}, { 2.26653218e-03 -2.26653218e-03}, {-2.80312643e-03, 2.80312643e-03}, {-2.11471342e-03, 2.11471342e-03}, {-5.22260235e-04, 5.22260235e-04}}, 

     {{-1.43017501e-05, 1.43017501e-05}, {-3.56462486e-03, 3.56462486e-03}, { 1.44878230e-04 -1.44878230e-04}, {-4.81204615e-05, 4.81204615e-05}, { 5.00971340e-03 -5.00971340e-03}, {-1.67957772e-04, 1.67957772e-04}, {-2.23214866e-04, 2.23214866e-04}, { 4.46518626e-04 -4.46518626e-04}, {-2.01605459e-03, 2.01605459e-03}, {-3.00531649e-04, 3.00531649e-04}, {-1.71281800e-03, 1.71281800e-03}}, 

     {{-1.77628609e-05, 1.77628609e-05}, { 6.10991400e-04 -6.10991400e-04}, {-4.75623338e-05, 4.75623338e-05}, { 8.33003554e-06 -8.33003554e-06}, {-1.47200379e-03, 1.47200379e-03}, { 2.73280209e-05 -2.73280209e-05}, { 8.01445000e-05 -8.01445000e-05}, {-1.28153462e-04, 1.28153462e-04}, { 5.99241226e-04 -5.99241226e-04}, {-1.37455485e-04, 1.37455485e-04}, {-6.30878784e-03, 6.30878784e-03}}, 

     {{ 6.63504512e-06 -6.63504512e-06}, { 1.83093502e-04 -1.83093502e-04}, {-6.38128750e-06, 6.38128750e-06}, { 2.59030443e-04 -2.59030443e-04}, { 1.29025452e-03 -1.29025452e-03}, {-4.99350399e-06, 4.99350399e-06}, { 2.09936384e-04 -2.09936384e-04}, {-1.68740682e-04, 1.68740682e-04}, { 2.84494019e-03 -2.84494019e-03}, {-1.81897111e-04, 1.81897111e-04}, {-2.68290941e-06, 2.68290941e-06}}, 

     {{-1.49569119e-04, 1.49569119e-04}, { 4.84503067e-06 -4.84503067e-06}, {-4.84246511e-07, 4.84246511e-07}, { 1.50934586e-04 -1.50934586e-04}, {-6.34205309e-06, 6.34205309e-06}, { 4.85034013e-07 -4.85034013e-07}, {-4.50559903e-04, 4.50559903e-04}, { 1.85705251e-06 -1.85705251e-06}, { 2.30770019e-05 -2.30770019e-05}, { 8.59132053e-06 -8.59132053e-06}, {-1.18152144e-06, 1.18152144e-06}}, 

     {{-4.24937792e-04, 4.24937792e-04}, {-3.86021231e-06, 3.86021231e-06}, {-6.21564394e-07, 6.21564394e-07}, {-4.22286808e-05, 4.22286808e-05}, {-6.72777148e-08, 6.72777148e-08}, {-2.16603057e-07, 2.16603057e-07}, { 1.26574482e-04 -1.26574482e-04}, { 4.43939925e-07 -4.43939925e-07}, {-4.43928160e-06, 4.43928160e-06}, {-3.75944699e-06, 3.75944699e-06}, { 2.04546059e-06 -2.04546059e-06}}, 

     {{-8.15047610e-13, 8.15047610e-13}, {-3.72784833e-14, 3.72784833e-14}, {-5.47371076e-16, 5.47371076e-16}, {-6.05781075e-07, 6.05781075e-07}, { 1.16932422e-08 -1.16932422e-08}, {-2.25716457e-10, 2.25716457e-10}, {-1.99761013e-07, 1.99761013e-07}, {-3.14201906e-09, 3.14201906e-09}, { 6.45517340e-08 -6.45517340e-08}, { 2.26130492e-09 -2.26130492e-09}, { 7.91426459e-11 -7.91426459e-11}}
   };



   
};




#endif // JpsiMuNtuplizer_H

