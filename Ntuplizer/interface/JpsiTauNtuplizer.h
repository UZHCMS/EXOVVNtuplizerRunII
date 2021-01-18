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
		    edm::EDGetTokenT<pat::PackedGenParticleCollection> packedgenptoken, 
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
   bool useDNN_;
   bool useHammer_;
   bool isTruth_;
   bool verbose_;

   bool flag_fill = false;

   float c_dz;
   float c_fsig;
   float c_vprob;
   float c_dnn;
   unsigned int c_charge;

   float chi = 0.;
   float ndf = 0.;

   const int numberofDNN = 80;

   const int numberofToys = 1000;
   std::vector<map<string, double>> FFdict;

   helper aux;
   
   tensorflow::MetaGraphDef* graphDef_old;
   tensorflow::MetaGraphDef* graphDef_perPF;
   tensorflow::MetaGraphDef* graphDef_perEVT;

   tensorflow::Session* session_old;
   tensorflow::Session* session_perPF;
   tensorflow::Session* session_perEVT;

   tensorflow::Tensor data; // (tensorflow::DT_FLOAT, { 1, 50, 8 }); // single batch of dimension 10

   tensorflow::Tensor label_perPF; // (tensorflow::DT_FLOAT, { 1,50}); 
   tensorflow::Tensor label_perEVT; // (tensorflow::DT_FLOAT, { 1,50}); 
   //   tensorflow::Tensor add_global; //(tensorflow::DT_FLOAT, { 1, 2 }); 
   tensorflow::Tensor isTraining; //(tensorflow::DT_FLOAT, { 1, 2 }); 
   //   tensorflow::Tensor norm; //(tensorflow::DT_FLOAT, { 1, 2 }); 

   tensorflow::Tensor data_old; // (tensorflow::DT_FLOAT, { 1, 50, 8 }); // single batch of dimension 10
   tensorflow::Tensor label_old; // (tensorflow::DT_FLOAT, { 1,50}); 
   tensorflow::Tensor add_global_old; //(tensorflow::DT_FLOAT, { 1, 2 }); 
   tensorflow::Tensor isTraining_old; //(tensorflow::DT_FLOAT, { 1, 2 }); 
   tensorflow::Tensor norm_old; //(tensorflow::DT_FLOAT, { 1, 2 }); 



   std::string dnnfile_old_;
   std::string dnnfile_perPF_;
   std::string dnnfile_perEVT_;

   //   Hammer::Hammer hammer;

   const vector<string> parName = {"a0", "a1", "a2", "b0", "b1", "b2", "c1", "c2", "d0", "d1", "d2"};
//   
   const vector<string> _FFErrNames = {"delta_a0","delta_a1","delta_a2","delta_b0","delta_b1","delta_b2","delta_c1","delta_c2","delta_d0","delta_d1","delta_d2"};
//
//   //   std::vector<double> _FFErr = {0.6850051095, 0.3674778905, 0.08740600666, 0.05468939568, 0.01128871118, 0.006722573912, 0.006537973692, 0.003156718233, 0.0004988336396, 0.0004454548237, 6.41244005e-07};
   std::vector<double> _FFErr = {0.6850042974, 0.3674780815, 0.0874061782, 0.0546894675, 0.0112887042, 0.0067225920, 0.0065379763, 0.0031567090, 0.0004988304, 0.0004454549, 0.0000006412};
   
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
//   const double avec_cent[3] = {-0.4182681183, -0.1819441298, 0.1140362192};
//   const double bvec_cent[3] = {-0.0901217995, 0.0074593144, -0.0127818946};
//   const double cvec_cent[2] = {0.0028349558, 0.0310274035};
//   const double dvec_cent[3] = {0.0013219585, -0.0043841356, -0.0000000258};



   std::vector<double> _FFmean = {-0.4182681183, -0.1819441298, 0.1140362192, -0.0901217995, 0.0074593144, -0.0127818946,  0.0028349558, 0.0310274035, 0.0013219585, -0.0043841356, -0.0000000258};



   const double eigVec[11][11] = {
     {4.35760696e-04, -1.28207353e-03,  3.67310259e-03, -4.49466512e-03, -6.35629036e-03, -2.12741605e-03, -2.71687448e-03,  2.10188686e-03, -2.99839624e-01, -9.53941269e-01, -1.27104129e-06},
     {8.75884938e-03,  3.61597837e-02, -1.15379356e-01,	1.52852441e-01,  8.17625885e-01, -5.30245605e-01, 9.34526790e-02,  5.80013879e-02,  9.71278152e-03, -8.66577627e-03, -5.81346305e-08},
     {-7.01822023e-01,  7.03934930e-01, -6.95185879e-02, 7.42766306e-02, -3.23615314e-02,  2.15509479e-02, -7.27477918e-03, -2.02150010e-03, -9.70763836e-04, -1.39534760e-03, -8.53608098e-10},
     {2.81185948e-04,  6.38864000e-04,  7.96887209e-04, -2.43558178e-03, -1.46714600e-02, -7.15802202e-03, 1.27409999e-03,  8.20571187e-02,  3.02576961e-01, -9.47990085e-02, -9.44696667e-01},
     {1.78756853e-02, -1.31949122e-02, -1.17563628e-01, 1.46797723e-01,  4.36526489e-01,  7.45205631e-01, -2.25146700e-01,  4.08734073e-01, -1.27138465e-02, -1.51031492e-04,  1.82352460e-02},
     {-6.84286257e-01, -7.05427114e-01, -1.72203692e-01, 6.03273549e-02, -1.36352214e-02, -2.49840794e-02, 4.17988987e-03, -1.58187024e-03,  9.72342534e-04, -4.86251398e-04, -3.51997765e-04},
     {2.88391659e-03,  1.60328155e-03, -2.01438109e-02, -7.36486847e-03, -1.54783464e-02, -3.32036908e-02, 1.22583039e-02,  6.65048271e-02, -9.03232652e-01, 2.84146584e-01, -3.11521062e-01},
     {-6.66271010e-02,  2.26172245e-02, -1.87894194e-01, -9.54851181e-01,  2.00778773e-01,  6.64206049e-02, -1.96013960e-02, -5.34546212e-02,  3.72281344e-03, 9.96599088e-04, -4.89988059e-03},
     {-3.79569853e-03,  9.57017854e-03, -1.01065324e-02, -1.23275108e-01, -2.48312506e-01, -2.99892451e-01, 9.16554601e-02,  9.01236134e-01,  4.62622205e-02, -9.96572677e-03,  1.00666413e-01},
     {1.85266608e-01,  6.89269388e-02, -9.49214460e-01,	1.36687843e-01, -1.87330041e-01, -4.47047285e-02, -2.10241638e-02, -5.76223886e-02,  1.72229290e-02, -8.43956857e-03,  3.52643440e-03},
     {-4.00356073e-03, -2.80939296e-03,  3.93076372e-02, -1.51221856e-02, -4.62639666e-02, -2.54785358e-01, -9.64945046e-01, -8.49907112e-04, -2.36858347e-03, 4.59184689e-03,  1.23420485e-04}
   };





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


//   const double eigVarInv[11][11][2] = {
//     
//     {{ 2.98497949e-04 -2.98497949e-04}, {-8.78225875e-04, 8.78225875e-04}, { 2.51609106e-03 -2.51609106e-03}, {-3.07886492e-03, 3.07886492e-03}, {-4.35408621e-03, 4.35408621e-03}, {-1.45728914e-03, 1.45728914e-03}, {-1.86107069e-03, 1.86107069e-03}, { 1.43980153e-03 -1.43980153e-03}, {-2.05391431e-01, 2.05391431e-01}, {-6.53453869e-01, 6.53453869e-01}, {-8.70668866e-07, 8.70668866e-07}},
//     {{ 3.21868517e-03 -3.21868517e-03}, { 1.32879279e-02 -1.32879279e-02}, {-4.23993846e-02, 4.23993846e-02}, { 5.61699219e-02 -5.61699219e-02}, { 3.00459592e-01 -3.00459592e-01}, {-1.94853638e-01, 1.94853638e-01}, { 3.43418112e-02 -3.43418112e-02}, { 2.13142387e-02 -2.13142387e-02}, { 3.56923432e-03 -3.56923432e-03}, {-3.18448284e-03, 3.18448284e-03}, {-2.13632009e-08, 2.13632009e-08}},
//     {{-6.13435808e-02, 6.13435808e-02}, { 6.15282619e-02 -6.15282619e-02}, {-6.07635408e-03, 6.07635408e-03}, { 6.49223641e-03 -6.49223641e-03}, {-2.82859778e-03, 2.82859778e-03}, { 1.88368600e-03 -1.88368600e-03}, {-6.35860645e-04, 6.35860645e-04}, {-1.76691598e-04, 1.76691598e-04}, {-8.48507568e-05, 8.48507568e-05}, {-1.21962001e-04, 1.21962001e-04}, {-7.46106760e-11, 7.46106760e-11}},
//     {{ 1.53779098e-05 -1.53779098e-05}, { 3.49391319e-05 -3.49391319e-05}, { 4.35813371e-05 -4.35813371e-05}, {-1.33200670e-04, 1.33200670e-04}, {-8.02374334e-04, 8.02374334e-04}, {-3.91468412e-04, 3.91468412e-04}, { 6.96798500e-05 -6.96798500e-05}, { 4.48766012e-03 -4.48766012e-03}, { 1.65477729e-02 -1.65477729e-02}, {-5.18450729e-03, 5.18450729e-03}, {-5.16649577e-02, 5.16649577e-02}},
//     {{ 2.01793324e-04 -2.01793324e-04}, {-1.48953461e-04, 1.48953461e-04}, {-1.32714102e-03, 1.32714102e-03}, { 1.65715608e-03 -1.65715608e-03}, { 4.92781842e-03 -4.92781842e-03}, { 8.41240595e-03 -8.41240595e-03}, {-2.54161450e-03, 2.54161450e-03}, { 4.61407806e-03 -4.61407806e-03}, {-1.43522852e-04, 1.43522852e-04}, {-1.70494984e-06, 1.70494984e-06}, { 2.05852299e-04 -2.05852299e-04}},
//     {{-4.60017731e-03, 4.60017731e-03}, {-4.74229867e-03, 4.74229867e-03}, {-1.15765516e-03, 1.15765516e-03}, { 4.05556193e-04 -4.05556193e-04}, {-9.16640301e-05, 9.16640301e-05}, {-1.67957772e-04, 1.67957772e-04}, { 2.80996941e-05 -2.80996941e-05}, {-1.06342682e-05, 1.06342682e-05}, { 6.53666213e-06 -6.53666213e-06}, {-3.26886975e-06, 3.26886975e-06}, {-2.36633735e-06, 2.36633735e-06}},
//     {{ 1.88549783e-05 -1.88549783e-05}, { 1.04822168e-05 -1.04822168e-05}, {-1.31699758e-04, 1.31699758e-04}, {-4.81513354e-05, 4.81513354e-05}, {-1.01197062e-04, 1.01197062e-04}, {-2.17084944e-04, 2.17084944e-04}, { 8.01445000e-05 -8.01445000e-05}, { 4.34806983e-04 -4.34806983e-04}, {-5.90531366e-03, 5.90531366e-03}, { 1.85774363e-03 -1.85774363e-03}, {-2.03671732e-03, 2.03671732e-03}},
//     {{-2.10322368e-04, 2.10322368e-04}, { 7.13959955e-05 -7.13959955e-05}, {-5.93127289e-04, 5.93127289e-04}, {-3.01418729e-03, 3.01418729e-03}, { 6.33800153e-04 -6.33800153e-04}, { 2.09670519e-04 -2.09670519e-04}, {-6.18759026e-05, 6.18759026e-05}, {-1.68740682e-04, 1.68740682e-04}, { 1.17518386e-05 -1.17518386e-05}, { 3.14597328e-06 -3.14597328e-06}, {-1.54674970e-05, 1.54674970e-05}},
//     {{-1.89340981e-06, 1.89340981e-06}, { 4.77389597e-06 -4.77389597e-06}, {-5.04144556e-06, 5.04144556e-06}, {-6.14933714e-05, 6.14933714e-05}, {-1.23865826e-04, 1.23865826e-04}, {-1.49595471e-04, 1.49595471e-04}, { 4.57205297e-05 -4.57205297e-05}, { 4.49563979e-04 -4.49563979e-04}, { 2.30770019e-05 -2.30770019e-05}, {-4.97120745e-06, 4.97120745e-06}, { 5.02154670e-05 -5.02154670e-05}},
//     {{ 8.25279141e-05 -8.25279141e-05}, { 3.07038411e-05 -3.07038411e-05}, {-4.22832212e-04, 4.22832212e-04}, { 6.08882665e-05 -6.08882665e-05}, {-8.34470804e-05, 8.34470804e-05}, {-1.99139394e-05, 1.99139394e-05}, {-9.36531634e-06, 9.36531634e-06}, {-2.56681741e-05, 2.56681741e-05}, { 7.67203774e-06 -7.67203774e-06}, {-3.75944699e-06, 3.75944699e-06}, { 1.57086741e-06 -1.57086741e-06}},
//     {{-2.56725932e-09, 2.56725932e-09}, {-1.80150640e-09, 1.80150640e-09}, { 2.52057867e-08 -2.52057867e-08}, {-9.69701084e-09, 9.69701084e-09}, {-2.96664913e-08, 2.96664913e-08}, {-1.63379583e-07, 1.63379583e-07}, {-6.18765226e-07, 6.18765226e-07}, {-5.44997841e-10, 5.44997841e-10}, {-1.51883995e-09, 1.51883995e-09}, { 2.94449429e-09 -2.94449429e-09}, { 7.91426459e-11 -7.91426459e-11}}
//     
//
//     
//   };
   



   std::vector<double> Inv = {2.13114700,
			      7.40520557,
			      1.30892794e+02,
			      3.34343289e+02,			      
			      7.84714743e+03,
			      2.21272016e+04,
			      2.33944757e+04,
			      1.00353127e+05,
			      4.01877950e+06,
			      5.03955932e+06,
			      2.43194284e+12};
			      
   TRandom3 *ran;			      



   //   std::vector<double> _FFErr = {1, 1,1,1,1,1,1,1,1,1,1};




};




#endif // JpsiTauNtuplizer_H

