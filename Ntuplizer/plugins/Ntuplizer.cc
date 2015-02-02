
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "../interface/NtupleBranches.h"
#include "../interface/CandidateNtuplizer.h"
#include "../interface/JetsNtuplizer.h"
#include "../interface/MuonsNtuplizer.h"
#include "../interface/ElectronsNtuplizer.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"



// class declaration


class NtupleBranches;
class CandidateNtuplizer;

class Ntuplizer : public edm::EDAnalyzer {
public:
  explicit Ntuplizer(const edm::ParameterSet&);
  ~Ntuplizer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  //       virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //       virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //       virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //       virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  // ----------member data ---------------------------
  
  NtupleBranches* nBranches_;
  std::map<std::string,CandidateNtuplizer*> nTuplizers_;
  
  
  edm::EDGetTokenT<reco::VertexCollection>  vtxToken_;
  
  edm::EDGetTokenT<pat::JetCollection> 		jetToken_;
  edm::EDGetTokenT<pat::JetCollection> 		fatjetToken_;
  
  edm::EDGetTokenT<pat::MuonCollection> 	muonToken_;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<pat::TauCollection> 		tauToken_;
  
  edm::EDGetTokenT<pat::METCollection> 		metToken_;
  
  const reco::Vertex PV;
  
  
};


Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig):

	vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
	jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
	fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
	muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
	electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
	tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
	metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets")))
	

{
  
  edm::Service<TFileService> fs;
  TTree* tree = fs->make<TTree>( "tree", "tree" );
  nBranches_ = new NtupleBranches( tree );
  
//   edm::InputTag vertices_ = iConfig.getParameter<edm::InputTag>("vertices");
  
  // std::vector<edm::InputTag> jetLabels;
//   jetLabels.push_back( iConfig.getParameter<edm::InputTag>("jets") );
//   jetLabels.push_back( iConfig.getParameter<edm::InputTag>("fatjets") );
  
  std::vector<edm::EDGetTokenT<pat::JetCollection>> jetTokens;
  jetTokens.push_back( jetToken_ 	);
  jetTokens.push_back( fatjetToken_ );
  
  
  
  /*=======================================================================================*/
  
  nTuplizers_["jets"]  = new JetsNtuplizer	( jetTokens,  nBranches_ );
  nTuplizers_["muons"] = new MuonsNtuplizer ( muonToken_, vtxToken_ , nBranches_ );
  nTuplizers_["electrons"] = new ElectronsNtuplizer ( electronToken_, vtxToken_ , nBranches_ );
  
}


Ntuplizer::~Ntuplizer()
{
  
   //     for( std::map<std::string,CandidateNtuplizer*>::iterator it = nTuplizers_.begin(); it != nTuplizers_.end(); ++it )
   // *      delete it->second;
   // *
   // *   nTuplizers_.clear();
   // *
   // *   delete nBranches_;
}


//
// member functions
//

// ------------ method called for each event  ------------
void
Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;  
  
  edm::Handle<reco::VertexCollection> vertices;
      iEvent.getByToken(vtxToken_, vertices);
      if (vertices->empty()) return; // skip the event if no PV found
      const reco::Vertex &PV = vertices->front();
  
  for( std::map<std::string,CandidateNtuplizer*>::iterator it = nTuplizers_.begin(); it != nTuplizers_.end(); ++it )
    (it->second)->fillBranches( iEvent, iSetup );
  
  nBranches_->fillTree();
  
  nBranches_->reset();
  
  
  
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
Ntuplizer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Ntuplizer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
 * void 
 * Ntuplizer::beginRun(edm::Run const&, edm::EventSetup const&)
 * {
 * }
 */

// ------------ method called when ending the processing of a run  ------------
/*
 * void 
 * Ntuplizer::endRun(edm::Run const&, edm::EventSetup const&)
 * {
 * }
 */

// ------------ method called when starting to processes a luminosity block  ------------
/*
 * void 
 * Ntuplizer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
 * {
 * }
 */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
 * void 
 * Ntuplizer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
 * {
 * }
 */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Ntuplizer);
