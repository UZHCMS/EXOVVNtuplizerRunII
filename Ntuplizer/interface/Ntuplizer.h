#include <memory>
#include "DataFormats/METReco/interface/METCollection.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"

#include "../interface/NtupleBranches.h"

class NtupleBranches;
class CandidateNtuplizer;

class Ntuplizer : public edm::EDAnalyzer {
public:
  explicit Ntuplizer(const edm::ParameterSet&);
  ~Ntuplizer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() 																	override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) 							override;
  virtual void endJob() 																	override;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) 							override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) 								override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) 	override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const  &, edm::EventSetup const&)    override;
  
  // ----------member data ---------------------------
  
  NtupleBranches* nBranches_;
  std::map<std::string,CandidateNtuplizer*> 				nTuplizers_			;
  bool runOnMC;
  
  edm::EDGetTokenT<reco::VertexCollection>  				vtxToken_   		;
  edm::EDGetTokenT<double>  		    					rhoToken_   		;
  edm::EDGetTokenT< std::vector<PileupSummaryInfo> > 		puinfoToken_		;
  edm::EDGetTokenT< GenEventInfoProduct > 		geneventToken_		;
  
  edm::EDGetTokenT<reco::GenParticleCollection> 			genparticleToken_	;
  
  edm::EDGetTokenT<pat::JetCollection> 						jetToken_			;
  edm::EDGetTokenT<pat::JetCollection> 						fatjetToken_		;
  edm::EDGetTokenT<pat::JetCollection> 						prunedjetToken_		;
  edm::EDGetTokenT<pat::JetCollection> 						softdropjetToken_	;
  
  edm::EDGetTokenT<reco::JetFlavourMatchingCollection> 		flavourToken_		;

  edm::EDGetTokenT<pat::MuonCollection>     				muonToken_			;	
  edm::EDGetTokenT<pat::ElectronCollection> 				electronToken_		;
  edm::EDGetTokenT<pat::TauCollection> 	    				tauToken_			;
  
  edm::EDGetTokenT<pat::METCollection> 	    				metToken_			;
  edm::InputTag pfMETlabel;	     // input label for particle-flow MET
  
  edm::EDGetTokenT<edm::TriggerResults>						triggerToken_;
  

};
