#include "../interface/Ntuplizer.h"
#include "../interface/CandidateNtuplizer.h"
#include "../interface/JetsNtuplizer.h"
#include "../interface/MuonsNtuplizer.h"
#include "../interface/ElectronsNtuplizer.h"
#include "../interface/METsNtuplizer.h"
#include "../interface/PileUpNtuplizer.h"
#include "../interface/GenEventNtuplizer.h"
#include "../interface/GenParticlesNtuplizer.h"
#include "../interface/TriggersNtuplizer.h"
#include "../interface/VerticesNtuplizer.h"


///////////////////////////////////////////////////////////////////////////////////
Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig):
	
	vtxToken_		(consumes<reco::VertexCollection>				(iConfig.getParameter<edm::InputTag>("vertices"))),
	rhoToken_		(consumes<double>						(iConfig.getParameter<edm::InputTag>("rho"))),
	puinfoToken_    	(consumes<std::vector<PileupSummaryInfo> >		        (iConfig.getParameter<edm::InputTag>("PUInfo"))),
	geneventToken_    	(consumes<GenEventInfoProduct>				        (iConfig.getParameter<edm::InputTag>("genEventInfo"))),
	
	genparticleToken_ 	(consumes<reco::GenParticleCollection>			        (iConfig.getParameter<edm::InputTag>("genparticles"))),
	
	jetToken_		(consumes<pat::JetCollection>					(iConfig.getParameter<edm::InputTag>("jets"))),
	fatjetToken_		(consumes<pat::JetCollection>					(iConfig.getParameter<edm::InputTag>("fatjets"))),
	prunedjetToken_		(consumes<pat::JetCollection>					(iConfig.getParameter<edm::InputTag>("prunedjets"))),
	softdropjetToken_	(consumes<pat::JetCollection>					(iConfig.getParameter<edm::InputTag>("softdropjets"))),
	
	flavourToken_		(consumes<reco::JetFlavourMatchingCollection>	                (iConfig.getParameter<edm::InputTag>("subjetflavour"))),

	muonToken_		(consumes<pat::MuonCollection>					(iConfig.getParameter<edm::InputTag>("muons"))),
	electronToken_		(consumes<pat::ElectronCollection>				(iConfig.getParameter<edm::InputTag>("electrons"))),
	tauToken_		(consumes<pat::TauCollection>					(iConfig.getParameter<edm::InputTag>("taus"))),
	metToken_		(consumes<pat::METCollection>					(iConfig.getParameter<edm::InputTag>("mets"))),
	
	triggerToken_		(consumes<edm::TriggerResults>					(iConfig.getParameter<edm::InputTag>("HLT")))

{
  
  /*=======================================================================================*/
  edm::Service<TFileService> fs;
  TTree* tree = fs->make<TTree>( "tree", "tree" );
  nBranches_ = new NtupleBranches( tree );
  runOnMC = iConfig.getParameter<bool>("runOnMC");
  
  /*=======================================================================================*/
  std::vector<edm::EDGetTokenT<pat::JetCollection>> jetTokens;
  jetTokens.push_back( jetToken_ 			);
  jetTokens.push_back( fatjetToken_ 		);
  jetTokens.push_back( prunedjetToken_ 		);
  jetTokens.push_back( softdropjetToken_	);
  // jetTokens.push_back( flavourToken_	 		);  

  /*=======================================================================================*/  
  std::vector<edm::EDGetTokenT<pat::METCollection>> metTokens;
  //METsLabels.push_back( iConfig.getParameter<edm::InputTag>("METrawColl") );
  metTokens.push_back( metToken_ );

  /*=======================================================================================*/  
  std::vector<edm::EDGetTokenT<reco::VertexCollection>> vtxTokens;
  vtxTokens.push_back( vtxToken_  );
  
  /*=======================================================================================*/
  std::vector<edm::EDGetTokenT<edm::TriggerResults>> triggerTokens;
  triggerTokens.push_back( triggerToken_ );
  
  /*=======================================================================================*/
  
  nTuplizers_["jets"]  	   = new JetsNtuplizer	    ( jetTokens		, flavourToken_	, nBranches_ );
  nTuplizers_["muons"] 	   = new MuonsNtuplizer     ( muonToken_	, vtxToken_	, rhoToken_, nBranches_ );
  nTuplizers_["electrons"] = new ElectronsNtuplizer ( electronToken_    , vtxToken_     , rhoToken_, nBranches_ );
  nTuplizers_["MET"]       = new METsNtuplizer      ( metTokens		, nBranches_ );
  nTuplizers_["vertices"]  = new VerticesNtuplizer  ( vtxTokens		, nBranches_ );
  nTuplizers_["triggers"]  = new TriggersNtuplizer  ( triggerTokens     , nBranches_ );

  /*=======================================================================================*/    
  if ( runOnMC ){
     std::vector<edm::EDGetTokenT<reco::GenParticleCollection>> genpTokens;
     genpTokens.push_back( genparticleToken_ );
     nTuplizers_["genParticles"] = new GenParticlesNtuplizer( genpTokens, nBranches_ );

     std::vector<edm::EDGetTokenT< std::vector<PileupSummaryInfo> > > puTokens;
     puTokens.push_back( puinfoToken_ );
     nTuplizers_["PU"] = new PileUpNtuplizer( puTokens, nBranches_ );

     std::vector<edm::EDGetTokenT< GenEventInfoProduct > > geneTokens;
     geneTokens.push_back( geneventToken_ );
     nTuplizers_["genEvent"] = new GenEventNtuplizer( geneTokens, nBranches_ );
  }
}

///////////////////////////////////////////////////////////////////////////////////
Ntuplizer::~Ntuplizer()
{
  
   for( std::map<std::string,CandidateNtuplizer*>::iterator it = nTuplizers_.begin(); it != nTuplizers_.end(); ++it )
      delete it->second;
   
   nTuplizers_.clear();
   
   delete nBranches_;
   
}

///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
	  
  if (vertices->empty()) return; // skip the event if no PV found
  

  nBranches_->EVENT_event     = iEvent.id().event();
  nBranches_->EVENT_run       = iEvent.id().run();
  nBranches_->EVENT_lumiBlock = iEvent.id().luminosityBlock();  
     
  for( std::map<std::string,CandidateNtuplizer*>::iterator it = nTuplizers_.begin(); it != nTuplizers_.end(); ++it )
    (it->second)->fillBranches( iEvent, iSetup );
  
  nBranches_->fillTree();
  
  nBranches_->reset();    
  
}


///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::beginJob(){
}

///////////////////////////////////////////////////////////////////////////////////
void Ntuplizer::endJob() {
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
