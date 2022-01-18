#include "../interface/Ntuplizer.h"
#include "../interface/CandidateNtuplizer.h"
#include "../interface/METsNtuplizer.h"
#include "../interface/PileUpNtuplizer.h"
#include "../interface/GenEventNtuplizer.h"
#include "../interface/GenParticlesNtuplizer.h"
#include "../interface/VerticesNtuplizer.h"
#include "../interface/JpsiMuNtuplizer.h"
#include "../interface/JpsiTauNtuplizer.h"
#include "../interface/JpsiKNtuplizer.h"
#include "../interface/JpsiKNtuplizerE.h"
//#include "../interface/BsTauTauNtuplizer.h"
//#include "../interface/BsTauTauFHNtuplizer.h"
//#include "../interface/BsTauTauFHNtuplizer_mr.h"
//#include "../interface/BsDstarTauNuNtuplizer.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "FWCore/ParameterSet/interface/FileInPath.h"
//#include <Math/SMatrix.h>

// #include "DataFormats/METReco/interface/PFMET.h"


///////////////////////////////////////////////////////////////////////////////////
Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig):
 
  beamToken_                  (consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
	vtxToken_             	    (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
	rhoToken_             	    (consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  	packedpfcandidatesToken_    (consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("packedpfcandidates"))), 
//        svToken_                    (consumes<std::vector<reco::VertexCompositePtrCandidate>>(iConfig.getParameter<edm::InputTag>("SecondaryVertices"))), 
	puinfoToken_          	    (consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("PUInfo"))),
	geneventToken_        	    (consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo"))),     
	lheEventProductToken_       (consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("externallheProducer"))),     
	genparticleToken_     	    (consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genparticles"))),
  packedgenparticleToken_     (consumes<pat::PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("packedgenparticles"))),
//  gentauToken_     	    (consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("gentaus"))),

	muonToken_	      	    (consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
//	electronToken_	      	    (consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),

	//mvaValuesMapToken_          (consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap"))),
	//mvaCategoriesMapToken_      (consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap"))),
//	ebRecHitsToken_             (consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(iConfig.getParameter<edm::InputTag>("ebRecHits"))),

//	tauToken_	      	    (consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
				     //tauBoostedTauToken_	    (consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("tausBoostedTau"))),

//	metToken_	      	    (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
//	metpuppiToken_	      	    (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets_puppi"))),
//	metmvaToken_	      	    (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets_mva"))),
//	metSigToken_	      	    (consumes<double>(edm::InputTag("METSignificance","METSignificance"))),
//	metCovToken_	      	    (consumes<math::Error<2>::type>(edm::InputTag("METSignificance","METCovariance"))),

//	jetForMetCorrToken_   	    (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsForMetCorr"))),

	triggerToken_	      	    (consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLT"))),
	triggerObjects_	      	    (consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerobjects")))
//  isBkgBSample_ ( iConfig.getParameter<bool>("isBkgBSample")),
//  bweightfile_ (runStrings["bweightfile"])      
//	triggerPrescales_     	    (consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("triggerprescales"))),
//        noiseFilterToken_     	    (consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("noiseFilter"))),
//        HBHENoiseFilterLooseResultToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("noiseFilterSelection_HBHENoiseFilterLoose"))),
//        HBHENoiseFilterTightResultToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("noiseFilterSelection_HBHENoiseFilterTight"))),
//	HBHENoiseIsoFilterResultToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("noiseFilterSelection_HBHENoiseIsoFilter"))),
//	 ecalBadCalibFilterUpdateToken_(consumes< bool >(iConfig.getParameter<edm::InputTag>("noiseFilterSelection_ecalBadCalibReducedMINIAODFilter")))

{


  nevents = 0;
  nevents_all = 0;

  /*=======================================================================================*/
  edm::Service<TFileService> fs;
  TTree* tree = fs->make<TTree>( "tree", "tree" );
 
  //std::map< std::string, bool > runFlags;
  runFlags["runOnMC"] = iConfig.getParameter<bool>("runOnMC");

  runFlags["useHammer"] = iConfig.getParameter<bool>("useHammer");
  runFlags["doGenParticles"] = iConfig.getParameter<bool>("doGenParticles");
  runFlags["doGenEvent"] = iConfig.getParameter<bool>("doGenEvent");
  runFlags["doPileUp"] = iConfig.getParameter<bool>("doPileUp");
  runFlags["doVertices"] = iConfig.getParameter<bool>("doVertices");
  runFlags["doMissingEt"] = iConfig.getParameter<bool>("doMissingEt");
  runFlags["doJpsiMu"] = iConfig.getParameter<bool>("doJpsiMu");
  runFlags["doJpsiTau"] = iConfig.getParameter<bool>("doJpsiTau");
  runFlags["doJpsiK"] = iConfig.getParameter<bool>("doJpsiK");
  runFlags["doJpsiKE"] = iConfig.getParameter<bool>("doJpsiKE");
  runFlags["isBkgBSample"] = iConfig.getParameter<bool>("isBkgBSample");
  //  runFlags["doBsTauTau"] = iConfig.getParameter<bool>("doBsTauTau");
  //  runFlags["doBsTauTauFH"] = iConfig.getParameter<bool>("doBsTauTauFH");
  //  runFlags["doBsTauTauFH_mr"] = iConfig.getParameter<bool>("doBsTauTauFH_mr");
  //  runFlags["doBsDstarTauNu"] = iConfig.getParameter<bool>("doBsDstarTauNu");
  runFlags["doGenHist"] = iConfig.getParameter<bool>("doGenHist");
  runFlags["verbose"] = iConfig.getParameter<bool>("verbose");
  runValues["dzcut"] = iConfig.getParameter<double>("dzcut");
  runValues["fsigcut"] = iConfig.getParameter<double>("fsigcut");
  runValues["vprobcut"] = iConfig.getParameter<double>("vprobcut");
  runValues["tau_charge"] = iConfig.getParameter<unsigned int>("tau_charge");

  //  runStrings["dnnfile_old"] = iConfig.getParameter<std::string>("dnnfile_old");  
  runStrings["dnnfile_perPF"] = iConfig.getParameter<std::string>("dnnfile_perPF");  
  runStrings["dnnfile_perEVT_mc"] = iConfig.getParameter<std::string>("dnnfile_perEVT_mc");  
  runStrings["dnnfile_perEVT_data"] = iConfig.getParameter<std::string>("dnnfile_perEVT_data");  
  runStrings["bweightfile"] = iConfig.getParameter<std::string>("bweightfile");


  //  std::cout << "Ntuplizer: dnn_old file: " << runStrings["dnnfile_old"] << std::endl;
  std::cout << "Ntuplizer: dnn_perPF file: " << runStrings["dnnfile_perPF"] << std::endl;
  std::cout << "Ntuplizer: dnn_perEVT file (mc): " << runStrings["dnnfile_perEVT_mc"] << std::endl;
  std::cout << "Ntuplizer: dnn_perEVT file (data): " << runStrings["dnnfile_perEVT_data"] << std::endl;


  std::cout << "Ntuplizer: (dzcut, fsigcut, vprobcut, dnn cut, tau_charge) = " << runValues["dzcut"] << " " << runValues["fsigcut"] << " " << runValues["vprobcut"] << " " << " " << runValues["tau_charge"] << std::endl;
  
  std::cout << "useHammer = " << runFlags["useHammer"] << std::endl;

  std::string jecpath = iConfig.getParameter<std::string>("jecpath");
  jecpath = "EXOVVNtuplizerRunII/Ntuplizer/data/" + jecpath;
  //  std::cout << "jecpath  "<< jecpath  <<std::endl;
 
  nBranches_ = new NtupleBranches( runFlags, tree );
  
  /*=======================================================================================*/
  /* Histogram buildinng, definition in NtupleBrances */
  /* Histogram for cutflow */

  nBranches_->cutflow = fs->make<TH1F>("cutflow", "Per Event Ntuplizer Cutflow", 15, 0, 15);
  nBranches_->bweight = fs->make<TH1F>("bweight", "b weight", 10000, 0, 40);
  //nBranches_->hist_BkgB_weight = fs->make<TH2F>("hist_BkgB_weight", "Bkg B sample decay chain and weight", 11, 0, 11, 100,0 ,100);

  nBranches_->nmuon = fs->make<TH1F>("nmuon", "number of muon", 10, 0, 10);


  if(runFlags["runOnMC"]){
    nBranches_->q2_nocut = fs->make<TH1F>("q2_nocut", "q2 before any cut", 40, 0, 15);
  }

  if(runFlags["useHammer"]){
    nBranches_->hammer_width = fs->make<TH1F>("hammer_width", "Hammer width", 32, 0, 32);  
    nBranches_->hammer_width_lattice = fs->make<TH1F>("hammer_width_lattice", "Hammer width lattice", 32, 0, 32);
  }

  /* Histogram for genParticles */ 
  if (runFlags["doGenHist"]){
      nBranches_->genParticle_Bdau_X_id=fs->make<TH1F>("genParticle_Bdau_X_id", "Identity of X in B->J/#psi+X;id;Events;", 18, 0, 18);
      nBranches_->genParticle_Bdau_X_pt =fs->make<TH1F>("genParticle_Bdau_X_pt", "p_{T} of X in B->J/#psi+X;pt[GeV];Events;", 100, 0, 20);
      nBranches_->genParticle_Bdau_X_eta =fs->make<TH1F>("genParticle_Bdau_X_eta", "#eta of X in B->J/#psi+X", 48, -2.4, 2.4);
      nBranches_->genParticle_Bdau_X_phi =fs->make<TH1F>("genParticle_Bdau_X_phi", "#phi of X in B->J/#psi+X", 64, -3.2, 3.2);
      nBranches_->genParticle_Bdau_X_mass =fs->make<TH1F>("genParticle_Bdau_X_mass", "#phi of X in B->J/#psi+X", 64, -3.2, 3.2);
      nBranches_->genParticle_Bdau_mu1_pt =fs->make<TH1F>("genParticle_Bdau_mu1_pt", "p_{T} of #mu_{J/#psi,1} in B->J/#psi+X;pt[GeV];Events;", 50, 0, 10);
      nBranches_->genParticle_Bdau_mu1_eta =fs->make<TH1F>("genParticle_Bdau_mu1_eta", "#eta of #mu_{J/#psi,1} in B->J/#psi+X", 48, -2.4, 2.4);
      nBranches_->genParticle_Bdau_mu1_phi =fs->make<TH1F>("genParticle_Bdau_mu1_phi", "#phi of #mu_{J/#psi,1} in B->J/#psi+X", 64, -3.2, 3.2);
      nBranches_->genParticle_Bdau_mu2_pt =fs->make<TH1F>("genParticle_Bdau_mu2_pt", "p_{T} of #mu_{J/#psi,2} in B->J/#psi+X;pt[GeV];Events;", 50, 0, 10);
      nBranches_->genParticle_Bdau_mu2_eta =fs->make<TH1F>("genParticle_Bdau_mu2_eta", "#eta of #mu_{J/#psi,2} in B->J/#psi+X", 48, -2.4, 2.4);
      nBranches_->genParticle_Bdau_mu2_phi =fs->make<TH1F>("genParticle_Bdau_mu2_phi", "#phi of #mu_{J/#psi,2} in B->J/#psi+X", 64, -3.2, 3.2);
      nBranches_->genParticle_Bdau_Jpsi_pt =fs->make<TH1F>("genParticle_Bdau_Jpsi_pt", "p_{T} of J/#psi in B->J/#psi+X;pt[GeV];Events;", 100, 0, 20);
      nBranches_->genParticle_Bdau_Jpsi_eta =fs->make<TH1F>("genParticle_Bdau_Jpsi_eta", "#eta of J/#psi in B->J/#psi+X", 48, -2.4, 2.4);
      nBranches_->genParticle_Bdau_Jpsi_phi =fs->make<TH1F>("genParticle_Bdau_Jpsi_phi", "#phi of J/#psi in B->J/#psi+X", 64, -3.2, 3.2);
      nBranches_->genParticle_Bdau_Jpsi_mass =fs->make<TH1F>("genParticle_Bdau_Jpsi_mass", "mass of J/#psi in B->J/#psi+X", 100, 0, 20);
      nBranches_->genParticle_Bvis_pt =fs->make<TH1F>("genParticle_Bvis_pt", "Visible p_{T} of B in B->J/#psi+X;pt[GeV];Events;", 100, 0, 20);
      nBranches_->genParticle_Bvis_eta =fs->make<TH1F>("genParticle_Bvis_eta", "Visible #eta of B in B->J/#psi+X", 48, -2.4, 2.4);
      nBranches_->genParticle_Bvis_phi =fs->make<TH1F>("genParticle_Bvis_phi", "Visible #phi of B in B->J/#psi+X", 64, -3.2, 3.2);
      nBranches_->genParticle_Bvis_mass =fs->make<TH1F>("genParticle_Bvis_mass", "Visible mass of B in B->J/#psi+X", 100, 0, 20);
  } 
  //  if (runFlags["doGenHist"]) {
  nBranches_-> LabelHistograms( runFlags );
      //  }

  /*=======================================================================================*/


 
  /*=======================================================================================*/
  if (runFlags["doMissingEt"]) {
    std::vector<std::string> corrFormulas;
    corrFormulas.push_back(iConfig.getParameter<std::string>("corrMetPx"));
    corrFormulas.push_back(iConfig.getParameter<std::string>("corrMetPy"));

    std::vector<std::string> jecAK4Labels;
    std::vector<std::string> tmpVec = iConfig.getParameter<std::vector<std::string> >("jecAK4forMetCorr");
    std::string tmpString = "";
    for( unsigned int v = 0; v < tmpVec.size(); ++v ){
       tmpString = jecpath + tmpVec[v];
       std::cout << " tmpString "<< tmpString <<std::endl;
       jecAK4Labels.push_back(edm::FileInPath(tmpString).fullPath());
    }
    
    nTuplizers_["MET"] = new METsNtuplizer( metToken_          , 
					    metpuppiToken_     , 
					    metmvaToken_     , 
                        jetForMetCorrToken_, 
					    muonToken_         ,
					    rhoToken_	       ,
					    vtxToken_	       ,
					    metSigToken_       ,
					    metCovToken_       ,
					    jecAK4Labels       ,
                        corrFormulas       ,
					    nBranches_         ,
					    runFlags  );
  }
    
  
  /*=======================================================================================*/  
  std::vector<edm::EDGetTokenT<reco::VertexCollection>> vtxTokens;
  vtxTokens.push_back( vtxToken_  );  


  /*=======================================================================================*/  

						      
  if (runFlags["doVertices"]) {
    nTuplizers_["vertices"] = new VerticesNtuplizer( vtxTokens   , 
						     beamToken_,
                                                     nBranches_  ,
						     runFlags    );
  }

  if ( runFlags["runOnMC"] ){
      if (runFlags["isBkgBSample"] ){
          
          std::string bweightfilepath = edm::FileInPath("EXOVVNtuplizerRunII/Ntuplizer/" + runStrings["bweightfile"]).fullPath();
          //  fileWeights = new TFile("root://eoscms.cern.ch//eos/cms/store/user/fiorendi/p5prime/HBMC/decay_weight.root", "READ");                                                                                                                           
          
          std::cout << "b weight file = "<< bweightfilepath << std::endl;
          fileWeights = new TFile((TString)bweightfilepath);
          histGenWeights=(TH1F*)fileWeights->Get("weight");
          std::cout << "weight hist first bin is : "<<  histGenWeights->GetBinContent(1)<<std::endl ;
      }
  }

  if (runFlags["doJpsiMu"]) {
    std::cout<<"\n\n --->GETTING INSIDE doJpsiMu<---\n\n"<<std::endl;
    nTuplizers_["JpsiMu"] = new JpsiMuNtuplizer( muonToken_   , 
						 vtxToken_   , 
						 packedpfcandidatesToken_,
						 triggerToken_,
						 triggerObjects_,
						 genparticleToken_,
						 runFlags,
						 nBranches_ );
  }
  if (runFlags["doJpsiTau"]) {
    std::cout<<"\n\n --->GETTING INSIDE doJpsiTau<---\n\n"<<std::endl;
    nTuplizers_["JpsiTau"] = new JpsiTauNtuplizer( muonToken_   , 
						   vtxToken_   , 
						   beamToken_ ,
						   packedpfcandidatesToken_,
						   triggerToken_,
						   triggerObjects_,
						   genparticleToken_,
						   packedgenparticleToken_,
						   runFlags,
						   runValues,
						   runStrings,
                           nBranches_ , 
                           histGenWeights);
  }


  if (runFlags["doJpsiK"]) {
    std::cout<<"\n\n --->GETTING INSIDE doJpsiK<---\n\n"<<std::endl;
    nTuplizers_["JpsiK"] = new JpsiKNtuplizer( muonToken_   , 
					       vtxToken_   , 
					       beamToken_ ,
					       packedpfcandidatesToken_,
					       triggerToken_,
					       triggerObjects_,
					       genparticleToken_,
					       packedgenparticleToken_,
					       runFlags,
					       runValues,
					       runStrings,
					       nBranches_, histGenWeights );
  }

  if (runFlags["doJpsiKE"]) {
    std::cout<<"\n\n --->GETTING INSIDE doJpsiKE <---\n\n"<<std::endl;
    nTuplizers_["JpsiKE"] = new JpsiKNtuplizerE( muonToken_   , 
						   vtxToken_   , 
						   beamToken_ ,
						   packedpfcandidatesToken_,
						   triggerToken_,
						   triggerObjects_,
						   genparticleToken_,
						   packedgenparticleToken_,
						   runFlags,
						   runValues,
						   runStrings,
						   nBranches_, histGenWeights );
  }


//  if (runFlags["doBsTauTau"]) {
//    std::cout<<"\n\n --->GETTING INSIDE doBsTauTau<---\n\n"<<std::endl;
//    nTuplizers_["BsTauTau"] = new BsTauTauNtuplizer( muonToken_   , 
//						     vtxToken_   , 
//						     packedpfcandidatesToken_,
//						     triggerToken_,
//						     triggerObjects_,
//						     genparticleToken_,
//						     gentauToken_,
//						     runFlags,
//						     runValues,
//						     runStrings,
//						     nBranches_ );
//  }
//
//  if (runFlags["doBsTauTauFH"]) {
//    std::cout<<"\n\n --->GETTING INSIDE doBsTauTauFH<---\n\n"<<std::endl;
//    nTuplizers_["BsTauTauFH"] = new BsTauTauFHNtuplizer( muonToken_   , 
//							 vtxToken_   , 
//							 packedpfcandidatesToken_,
//							 triggerToken_,
//							 triggerObjects_,
//							 genparticleToken_,
//							 gentauToken_,
//							 runFlags,
//							 runValues,
//							 runStrings,
//							 nBranches_ );
//  }
//
//  if (runFlags["doBsTauTauFH_mr"]) {
//    std::cout<<"\n\n --->GETTING INSIDE doBsTauTauFH_mr<---\n\n"<<std::endl;
//    nTuplizers_["BsTauTauFH_mr"] = new BsTauTauFHNtuplizer_mr( muonToken_   , 
//							 vtxToken_   , 
//							 packedpfcandidatesToken_,
//							 triggerToken_,
//							 triggerObjects_,
//							 genparticleToken_,
//							 gentauToken_,
//							 runFlags,
//							 runValues,
//							 runStrings,
//							 nBranches_ );
//  }
//
//  if (runFlags["doBsDstarTauNu"]) {
//    std::cout<<"\n\n --->GETTING INSIDE doBsDstarTauNu<---\n\n"<<std::endl;
//    nTuplizers_["BsDstarTauNu"] = new BsDstarTauNuNtuplizer( muonToken_   , 
//							 vtxToken_   , 
//							 packedpfcandidatesToken_,
//							 triggerToken_,
//							 triggerObjects_,
//							 genparticleToken_,
//							 gentauToken_,
//							 runFlags,
//							 runValues,
//							 runStrings,
//							 nBranches_ );
//  }
  

 
  /*=======================================================================================*/    
  if ( runFlags["runOnMC"] ){
     
    if (runFlags["doGenParticles"]) {
      std::vector<edm::EDGetTokenT<reco::GenParticleCollection>> genpTokens;
      genpTokens.push_back( genparticleToken_ );

      nTuplizers_["genParticles"] = new GenParticlesNtuplizer( genpTokens, nBranches_, runFlags//, histGenWeights 
                                                               );
    }

    if (runFlags["doPileUp"]) {
      std::vector<edm::EDGetTokenT< std::vector<PileupSummaryInfo> > > puTokens;
      puTokens.push_back( puinfoToken_ );
      nTuplizers_["PU"] = new PileUpNtuplizer( puTokens, nBranches_, runFlags );
    }

    if (runFlags["doGenEvent"]) {
      std::vector<edm::EDGetTokenT< GenEventInfoProduct > > geneTokens;
      geneTokens.push_back( geneventToken_ );
      std::vector<edm::EDGetTokenT<  LHEEventProduct > > lheTokens;
      lheTokens.push_back( lheEventProductToken_);
      nTuplizers_["genEvent"] = new GenEventNtuplizer( geneTokens, nBranches_ , lheTokens, runFlags);
    }
  }


  // if (runFlags["doGenHist"] || runFlags["doGenHist"]){
  //     LabelHistogram(nBranches_, runFlags );


  // }
}

///////////////////////////////////////////////////////////////////////////////////
Ntuplizer::~Ntuplizer()
{
	  
  for( std::unordered_map<std::string,CandidateNtuplizer*>::iterator it = nTuplizers_.begin(); it != nTuplizers_.end(); ++it ){
    //    std::cout << "deconstructor: Branches: " << it->first << std::endl;
    delete it->second;
  }
   nTuplizers_.clear();
   
   delete nBranches_;
   
   std::cout << "Total number of filled events = " << nevents << "/" << nevents_all << ", eff = " << Float_t(nevents)/Float_t(nevents_all) << std::endl;

}

///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  //  std::cout << "check1" << std::endl;
  
  nBranches_->reset();

  //  std::cout << "check2" << std::endl;

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if( vertices->empty() ) return; // skip the event if no PV found
  
  //  std::cout << "check3" << std::endl;

  nBranches_->EVENT_event     = iEvent.id().event();
  nBranches_->EVENT_run       = iEvent.id().run();
  nBranches_->EVENT_lumiBlock = iEvent.id().luminosityBlock();  

  //  std::cout<<" ----------------- before the branches loop"<<std::endl; 

  bool isSave = true;
  for( std::unordered_map<std::string,CandidateNtuplizer*>::iterator it = nTuplizers_.begin(); it != nTuplizers_.end(); ++it ){

    isSave = (it->second)->fillBranches( iEvent, iSetup );
    //    std::cout << "Fill Branchines: " << it->first << ", result = " << isSave << std::endl;
    if(!isSave) break;
  }

  //  std::cout << "isSave =  " << isSave << std::endl;
  if(isSave){
    //    std::cout << "-------------- save -----------" << std::endl;
    nBranches_->fillTree();
    nevents++;
  }
  
  nBranches_->reset();    
 
  nevents_all++;
}






///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::beginJob(){
}

///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::endJob( ) {
    if (runFlags["isBkgBSample"]) fileWeights->Close();
}

///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::beginRun(edm::Run const&, edm::EventSetup const&){
}


///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::endRun(edm::Run const&, edm::EventSetup const&){
}

///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
}


///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&){
}


///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Ntuplizer);
