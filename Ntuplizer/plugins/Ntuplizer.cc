#include "../interface/Ntuplizer.h"
#include "../interface/CandidateNtuplizer.h"
#include "../interface/JetsNtuplizer.h"
#include "../interface/GenJetsNtuplizer.h"
#include "../interface/MuonsNtuplizer.h"
#include "../interface/ElectronsNtuplizer.h"
#include "../interface/TausNtuplizer.h"
#include "../interface/METsNtuplizer.h"
#include "../interface/PileUpNtuplizer.h"
#include "../interface/GenEventNtuplizer.h"
#include "../interface/GenParticlesNtuplizer.h"
#include "../interface/TriggersNtuplizer.h"
#include "../interface/VerticesNtuplizer.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

// #include "DataFormats/METReco/interface/PFMET.h"


///////////////////////////////////////////////////////////////////////////////////
Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig):
	

	vtxToken_             	    (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
	rhoToken_             	    (consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
	puinfoToken_          	    (consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("PUInfo"))),
	geneventToken_        	    (consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfo"))),     
	lheEventProductToken_       (consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("externallheProducer"))),     
	genparticleToken_     	    (consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genparticles"))),
	
	jetToken_             	    (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
	fatjetToken_          	    (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
	prunedjetToken_	      	    (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("prunedjets"))),
	softdropjetToken_     	    (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("softdropjets"))),
	trimmedjetToken_	    (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("trimmedjets"))),
	puppijetToken_	            (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("puppijets"))),
	genJetToken_	      	    (consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets"))),
        genJetAK8Token_	      	    (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("genJetsAK8"))),
	
	flavourToken_	      	    (consumes<reco::JetFlavourMatchingCollection>(iConfig.getParameter<edm::InputTag>("subjetflavour"))),

	muonToken_	      	    (consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
	electronToken_	      	    (consumes<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"))),
	eleHEEPIdMapToken_    	    (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHEEPIdMap"))),
        eleHEEPId51MapToken_  	    (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleHEEPId51Map"))),
	eleVetoIdMapToken_    	    (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"))),
	eleLooseIdMapToken_   	    (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
	eleMediumIdMapToken_  	    (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
	eleTightIdMapToken_   	    (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
	tauToken_	      	    (consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
	tauBoostedTauToken_	    (consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("tausBoostedTau"))),

	metToken_	      	    (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
	metpuppiToken_	      	    (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets_puppi"))),
	metmvaToken_	      	    (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets_mva"))),
	metSigToken_	      	    (consumes<double>(edm::InputTag("METSignificance","METSignificance"))),
	metCovToken_	      	    (consumes<math::Error<2>::type>(edm::InputTag("METSignificance","METCovariance"))),

	jetForMetCorrToken_   	    (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsForMetCorr"))),

	triggerToken_	      	    (consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("HLT"))),
	triggerObjects_	      	    (consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerobjects"))),
	triggerPrescales_     	    (consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("triggerprescales"))),
        noiseFilterToken_     	    (consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("noiseFilter"))),
        HBHENoiseFilterLooseResultToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("noiseFilterSelection_HBHENoiseFilterLoose"))),
        HBHENoiseFilterTightResultToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("noiseFilterSelection_HBHENoiseFilterTight"))),
        HBHENoiseIsoFilterResultToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("noiseFilterSelection_HBHENoiseIsoFilter")))

{

	
  /*=======================================================================================*/
  edm::Service<TFileService> fs;
  TTree* tree = fs->make<TTree>( "tree", "tree" );
  
  std::map< std::string, bool > runFlags;
  runFlags["runOnMC"] = iConfig.getParameter<bool>("runOnMC");
  runFlags["doGenParticles"] = iConfig.getParameter<bool>("doGenParticles");
  runFlags["doGenJets"] = iConfig.getParameter<bool>("doGenJets");
  runFlags["doGenEvent"] = iConfig.getParameter<bool>("doGenEvent");
  runFlags["doPileUp"] = iConfig.getParameter<bool>("doPileUp");
  runFlags["doElectrons"] = iConfig.getParameter<bool>("doElectrons");
  runFlags["doMuons"] = iConfig.getParameter<bool>("doMuons");
  runFlags["doTaus"] = iConfig.getParameter<bool>("doTaus");
  runFlags["doAK8Jets"] = iConfig.getParameter<bool>("doAK8Jets");
  runFlags["doAK4Jets"] = iConfig.getParameter<bool>("doAK4Jets");
  runFlags["doVertices"] = iConfig.getParameter<bool>("doVertices");
  runFlags["doTriggerDecisions"] = iConfig.getParameter<bool>("doTriggerDecisions");
  runFlags["doTriggerObjects"] = iConfig.getParameter<bool>("doTriggerObjects");
  runFlags["doHltFilters"] = iConfig.getParameter<bool>("doHltFilters");
  runFlags["doMissingEt"] = iConfig.getParameter<bool>("doMissingEt");
  runFlags["doBoostedTaus"] = iConfig.getParameter<bool>("doBoostedTaus");
  runFlags["doPrunedSubjets"] = iConfig.getParameter<bool>("doPrunedSubjets");
  runFlags["doTrimming"] = iConfig.getParameter<bool>("doTrimming");
  runFlags["doPuppi"] = iConfig.getParameter<bool>("doPuppi");
  runFlags["doHbbTag"] = iConfig.getParameter<bool>("doHbbTag");
  runFlags["doMETSVFIT"] = iConfig.getParameter<bool>("doMETSVFIT");
  runFlags["doMVAMET"] = iConfig.getParameter<bool>("doMVAMET");
  runFlags["doPuppiRecluster"] = iConfig.getParameter<edm::InputTag>("puppijets").label()!="";

  std::string jecpath = iConfig.getParameter<std::string>("jecpath");
  
  nBranches_ = new NtupleBranches( runFlags, tree );
  
  /*=======================================================================================*/
  if (runFlags["doAK4Jets"] || runFlags["doAK8Jets"]) {
  
    std::vector<edm::EDGetTokenT<pat::JetCollection>> jetTokens;
    jetTokens.push_back( jetToken_ 	   );
    jetTokens.push_back( fatjetToken_ 	   );
    jetTokens.push_back( prunedjetToken_   );
    jetTokens.push_back( softdropjetToken_ );
    jetTokens.push_back( trimmedjetToken_  );
    jetTokens.push_back( puppijetToken_    );
    //jetTokens.push_back( flavourToken_	 );  
  
    std::vector<std::string> jecAK8Labels;
    std::string tmpString = "";
    std::vector<std::string> tmpVec = iConfig.getParameter<std::vector<std::string> >("jecAK8chsPayloadNames");
    for( unsigned int v = 0; v < tmpVec.size(); ++v ){
       tmpString = jecpath + tmpVec[v];
       jecAK8Labels.push_back(tmpString);
    }    
    jecAK8Labels.push_back( iConfig.getParameter<std::string>("jecAK8chsUnc") );
    std::vector<std::string> jecAK8GroomedLabels;
    tmpVec.clear(); tmpVec = iConfig.getParameter<std::vector<std::string> >("jecAK8GroomedchsPayloadNames");
    for( unsigned int v = 0; v < tmpVec.size(); ++v ){
       tmpString = jecpath + tmpVec[v];
       jecAK8GroomedLabels.push_back(tmpString);
    }    
    
    std::vector<std::string> jecAK8PuppiLabels;
    tmpVec.clear(); tmpVec = iConfig.getParameter<std::vector<std::string> >("jecAK8PuppiPayloadNames");
    for( unsigned int v = 0; v < tmpVec.size(); ++v ){
       tmpString = jecpath + tmpVec[v];
       jecAK8PuppiLabels.push_back(tmpString);
    }    
    
    std::vector<std::string> jecAK4chsLabels;
    tmpVec.clear(); tmpVec = iConfig.getParameter<std::vector<std::string> >("jecAK4chsPayloadNames");
    for( unsigned int v = 0; v < tmpVec.size(); ++v ){
       tmpString = jecpath + tmpVec[v];
       jecAK4chsLabels.push_back(tmpString);
    }
    jecAK4chsLabels.push_back( iConfig.getParameter<std::string>("jecAK4chsUnc") );
    
  
    std::vector<std::string> jerAK8chsFileLabels;
    jerAK8chsFileLabels.push_back( iConfig.getParameter<std::string>("jerAK8chs_res_PayloadNames") );
    jerAK8chsFileLabels.push_back( iConfig.getParameter<std::string>("jerAK8chs_sf_PayloadNames") );
    std::vector<std::string> jerAK4chsFileLabels;
    jerAK4chsFileLabels.push_back( iConfig.getParameter<std::string>("jerAK4chs_res_PayloadNames") );
     jerAK4chsFileLabels.push_back( iConfig.getParameter<std::string>("jerAK4chs_sf_PayloadNames") );
     std::vector<std::string> jerAK8PuppiFileLabels;
    jerAK8PuppiFileLabels.push_back( iConfig.getParameter<std::string>("jerAK8Puppi_res_PayloadNames") );
    jerAK8PuppiFileLabels.push_back( iConfig.getParameter<std::string>("jerAK8Puppi_sf_PayloadNames") );

    std::vector<std::string> jerAK4PuppiFileLabels;
    jerAK4PuppiFileLabels.push_back( iConfig.getParameter<std::string>("jerAK4Puppi_res_PayloadNames") );
    jerAK4PuppiFileLabels.push_back( iConfig.getParameter<std::string>("jerAK4Puppi_sf_PayloadNames") );


    nTuplizers_["jets"] = new JetsNtuplizer( jetTokens      , 
                                             jecAK4chsLabels, 
					     jecAK8Labels   , 
					     jecAK8GroomedLabels   , 
					     jecAK8PuppiLabels   , 
					     flavourToken_  , 
					     rhoToken_      , 
					     vtxToken_      , 
					     nBranches_     ,
                                             runFlags	    ,
					     jerAK8chsFileLabels,
					     jerAK4chsFileLabels,
					     jerAK8PuppiFileLabels,
					     jerAK4PuppiFileLabels
					     ); 
  }

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
       jecAK4Labels.push_back(tmpString);
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

  if (runFlags["doMuons"]) {
    nTuplizers_["muons"]= new MuonsNtuplizer( muonToken_   , 
                                              vtxToken_    , 
					      rhoToken_    , 
					      tauBoostedTauToken_ ,
					      nBranches_  ,
					      runFlags     );
  }
						      
  if (runFlags["doElectrons"]) {
    
    std::vector<edm::EDGetTokenT<edm::ValueMap<bool> > > eleIdTokens;
    eleIdTokens.push_back(eleVetoIdMapToken_  );
    eleIdTokens.push_back(eleLooseIdMapToken_ );
    eleIdTokens.push_back(eleMediumIdMapToken_);
    eleIdTokens.push_back(eleTightIdMapToken_ );
    eleIdTokens.push_back(eleHEEPIdMapToken_  );
    eleIdTokens.push_back(eleHEEPId51MapToken_  );
    
    nTuplizers_["electrons"] = new ElectronsNtuplizer( electronToken_, 
                                                       vtxToken_     , 
						       rhoToken_     , 
						       eleIdTokens   , 
						       tauBoostedTauToken_ ,
						       nBranches_  ,
						       runFlags     );
  }    
						      
  if (runFlags["doVertices"]) {
    nTuplizers_["vertices"] = new VerticesNtuplizer( vtxTokens   , 
                                                     nBranches_ );
  }
  
  if (runFlags["doTriggerDecisions"] || runFlags["doTriggerObjects"] || runFlags["doTriggerDecisions"]) {
    nTuplizers_["triggers"] = new TriggersNtuplizer( triggerToken_, 
                                                     triggerObjects_, 
						     triggerPrescales_,
                                                     noiseFilterToken_,
						     HBHENoiseFilterLooseResultToken_,
						     HBHENoiseFilterTightResultToken_,
						     HBHENoiseIsoFilterResultToken_,
						     nBranches_,
                                                     iConfig,
                                                     runFlags );
  }

  if (runFlags["doTaus"]) {
     nTuplizers_["taus"] = new TausNtuplizer( tauToken_      ,  
                                              tauBoostedTauToken_,  
					      rhoToken_      ,  
					      vtxToken_      ,  
					      nBranches_     , 
                                              runFlags      ); 
  }
  /*=======================================================================================*/    
  if ( runFlags["runOnMC"] ){

    if (runFlags["doGenJets"]) 
      nTuplizers_["genJets"]   = new GenJetsNtuplizer   ( genJetToken_, genJetAK8Token_, nBranches_    );

    if (runFlags["doGenParticles"]) {
      std::vector<edm::EDGetTokenT<reco::GenParticleCollection>> genpTokens;
      genpTokens.push_back( genparticleToken_ );
      nTuplizers_["genParticles"] = new GenParticlesNtuplizer( genpTokens, nBranches_ );
    }

    if (runFlags["doPileUp"]) {
      std::vector<edm::EDGetTokenT< std::vector<PileupSummaryInfo> > > puTokens;
      puTokens.push_back( puinfoToken_ );
      nTuplizers_["PU"] = new PileUpNtuplizer( puTokens, nBranches_ );
    }

    if (runFlags["doGenEvent"]) {
      std::vector<edm::EDGetTokenT< GenEventInfoProduct > > geneTokens;
      geneTokens.push_back( geneventToken_ );
      std::vector<edm::EDGetTokenT<  LHEEventProduct > > lheTokens;
      lheTokens.push_back( lheEventProductToken_);
      nTuplizers_["genEvent"] = new GenEventNtuplizer( geneTokens, nBranches_ , lheTokens);
    }
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
  
  nBranches_->reset();

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if( vertices->empty() ) return; // skip the event if no PV found
           
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
