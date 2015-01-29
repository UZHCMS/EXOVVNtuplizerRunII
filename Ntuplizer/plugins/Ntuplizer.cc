
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

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TLorentzVector.h>

// class declaration
//

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
      bool runOnMC;
      std::vector<unsigned long long> runRange_;

      std::string JetCorrectionService;
     
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig)

{
  
  edm::Service<TFileService> fs;
  TTree* tree = fs->make<TTree>( "tree", "tree" );
  nBranches_ = new NtupleBranches( tree );
  
  
  /*=======================================================================================*/
  
  nTuplizers_["jets"] = new JetsNtuplizer(  "selectedPatJetsCA8CHSwithQJetsForBoostedTaus", jecLevelsCA8, jecAK5Labels, nBranches_ );


}


Ntuplizer::~Ntuplizer()
{
 
/*     for( std::map<std::string,CandidateNtuplizer*>::iterator it = nTuplizers_.begin(); it != nTuplizers_.end(); ++it )
      delete it->second;
      
   nTuplizers_.clear();
   
   delete nBranches_; */  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   for( std::map<std::string,CandidateNtuplizer*>::iterator it = nTuplizers_.begin(); it != nTuplizers_.end(); ++it )
      (it->second)->fillBranches( iEvent, iSetup );
   
   nBranches_->fillTree();
   
   nBranches_->reset();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
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
void 
Ntuplizer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
Ntuplizer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
Ntuplizer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
Ntuplizer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
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
