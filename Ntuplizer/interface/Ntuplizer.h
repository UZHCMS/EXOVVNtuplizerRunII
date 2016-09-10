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
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "../interface/NtupleBranches.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"
class NtupleBranches;
class CandidateNtuplizer;

class Ntuplizer : public edm::EDAnalyzer {
public:
  explicit Ntuplizer(const edm::ParameterSet&);
  ~Ntuplizer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob()                                                                   override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&)                           override;
  virtual void endJob()                                                                     override;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&)                            override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&)                              override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)    override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const  &, edm::EventSetup const&)    override;
  
  // ----------member data ---------------------------
  
  NtupleBranches* nBranches_;
  std::map<std::string,CandidateNtuplizer*>                 nTuplizers_         ;
  
  edm::EDGetTokenT<reco::VertexCollection>                  vtxToken_           ;
  edm::EDGetTokenT<double>                                  rhoToken_           ;
  edm::EDGetTokenT< std::vector<PileupSummaryInfo> >        puinfoToken_        ;
  edm::EDGetTokenT< GenEventInfoProduct >                   geneventToken_      ;
  edm::EDGetTokenT<LHEEventProduct>	                    lheEventProductToken_;
  
  edm::EDGetTokenT<reco::GenParticleCollection>             genparticleToken_   ;
  
  edm::EDGetTokenT<pat::JetCollection>                      jetToken_           ;
  edm::EDGetTokenT<pat::JetCollection>                      fatjetToken_        ;
  edm::EDGetTokenT<pat::JetCollection>                      prunedjetToken_     ;
  edm::EDGetTokenT<pat::JetCollection>                      softdropjetToken_   ;
  edm::EDGetTokenT<pat::JetCollection>                      trimmedjetToken_     ;
  edm::EDGetTokenT<pat::JetCollection>                      puppijetToken_     ;
  edm::EDGetTokenT<pat::JetCollection>                      puppiprunedjetToken_     ;
  edm::EDGetTokenT<pat::JetCollection>                      puppisoftdropjetToken_     ;
  edm::EDGetTokenT<reco::GenJetCollection>                  genJetToken_	;
  edm::EDGetTokenT<pat::JetCollection>                      genJetAK8Token_	;
  

  edm::EDGetTokenT<reco::JetFlavourMatchingCollection> 	    flavourToken_	;

  edm::EDGetTokenT<pat::MuonCollection>     		    muonToken_  	;	
  edm::EDGetTokenT<edm::View<pat::Electron> >		    electronToken_	;
  edm::EDGetTokenT<edm::ValueMap<bool> >                    eleHEEPIdMapToken_  ;
  edm::EDGetTokenT<edm::ValueMap<bool> >                    eleHEEPId51MapToken_  ;
  edm::EDGetTokenT<edm::ValueMap<bool> >                    eleVetoIdMapToken_  ;
  edm::EDGetTokenT<edm::ValueMap<bool> >                    eleLooseIdMapToken_ ;
  edm::EDGetTokenT<edm::ValueMap<bool> >                    eleMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool> >                    eleTightIdMapToken_ ;
  edm::EDGetTokenT<pat::TauCollection> 	    		    tauToken_		;
  edm::EDGetTokenT<pat::TauCollection> 	    		    tauBoostedTauToken_	;


  edm::EDGetTokenT<pat::METCollection> 	    		    metToken_		;
  edm::EDGetTokenT<pat::METCollection> 	    		    metpuppiToken_		;
  edm::EDGetTokenT<pat::METCollection> 	    		    metmvaToken_		;
  edm::EDGetTokenT<double> 	    		            metSigToken_		;
  edm::EDGetTokenT<math::Error<2>::type> 	     	    metCovToken_		;
  edm::EDGetTokenT<pat::JetCollection>                      jetForMetCorrToken_ ;
  
  edm::EDGetTokenT<edm::TriggerResults>                     triggerToken_       ;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection>  triggerObjects_     ;
  edm::EDGetTokenT<pat::PackedTriggerPrescales>             triggerPrescales_   ;
  
  edm::EDGetTokenT<edm::TriggerResults>                     noiseFilterToken_;
  edm::EDGetTokenT<bool>                                    HBHENoiseFilterLooseResultToken_;
  edm::EDGetTokenT<bool>                                    HBHENoiseFilterTightResultToken_;
  edm::EDGetTokenT<bool>                                    HBHENoiseIsoFilterResultToken_;  

};
