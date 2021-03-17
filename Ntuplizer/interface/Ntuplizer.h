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
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"



class NtupleBranches;
class CandidateNtuplizer;

class Ntuplizer : public edm::EDAnalyzer {

public:
  explicit Ntuplizer(const edm::ParameterSet&);
  ~Ntuplizer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
 public:
  std::map< std::string, bool > runFlags;
  std::map< std::string, double > runValues;
  std::map< std::string, std::string > runStrings;
 

private:
  virtual void beginJob()                                                                   override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&)                           override;
  virtual void endJob()                                                                     override;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&)                            override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&)                              override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)    override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const  &, edm::EventSetup const&)    override;
  
  // ----------member data ---------------------------

  Int_t nevents;
  Int_t nevents_all;
  NtupleBranches* nBranches_;

  std::unordered_map<std::string,CandidateNtuplizer*>                 nTuplizers_         ;
  
  edm::EDGetTokenT<reco::BeamSpot>                          beamToken_          ;
  edm::EDGetTokenT<reco::VertexCollection>                  vtxToken_           ;
  edm::EDGetTokenT<double>                                  rhoToken_           ;
  edm::EDGetTokenT<pat::PackedCandidateCollection>          packedpfcandidatesToken_;
  edm::EDGetTokenT<std::vector<reco::VertexCompositePtrCandidate>>          svToken_;

  edm::EDGetTokenT< std::vector<PileupSummaryInfo> >        puinfoToken_        ;
  edm::EDGetTokenT< GenEventInfoProduct >                   geneventToken_      ;
  edm::EDGetTokenT<LHEEventProduct>	                    lheEventProductToken_;
  
  edm::EDGetTokenT<reco::GenParticleCollection>             genparticleToken_   ;
  edm::EDGetTokenT<pat::PackedGenParticleCollection>             packedgenparticleToken_   ;

  //  edm::EDGetTokenT<std::vector<reco::GenJet>>               gentauToken_   ;
 
  edm::EDGetTokenT<pat::MuonCollection>     		    muonToken_  	;	
  edm::EDGetTokenT<pat::ElectronCollection>     		    electronToken_  	;	
  //  edm::EDGetTokenT<edm::View<pat::Electron> >		    electronToken_	;
  /* edm::EDGetTokenT<edm::ValueMap<bool> >                    eleVetoIdMapToken_  ; */
  /* edm::EDGetTokenT<edm::ValueMap<bool> >                    eleLooseIdMapToken_ ; */
  /* edm::EDGetTokenT<edm::ValueMap<bool> >                    eleMediumIdMapToken_; */
  /* edm::EDGetTokenT<edm::ValueMap<bool> >                    eleTightIdMapToken_ ; */
  /* edm::EDGetTokenT<edm::ValueMap<bool> >                    eleHLTIdMapToken_  ; */
  /* edm::EDGetTokenT<edm::ValueMap<bool> >                    eleHEEPIdMapToken_  ; */
  /* edm::EDGetTokenT<edm::ValueMap<bool> >                    eleMVAMediumIdMapToken_; */
  /* edm::EDGetTokenT<edm::ValueMap<bool> >                    eleMVATightIdMapToken_ ; */
  /* edm::EDGetTokenT<edm::ValueMap<float> >                   mvaValuesMapToken_; */
  /* edm::EDGetTokenT<edm::ValueMap<int> >                     mvaCategoriesMapToken_; */
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>> ebRecHitsToken_;
  edm::EDGetTokenT<pat::TauCollection> 	    		    tauToken_		;
  

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
  edm::EDGetTokenT<bool>                                    ecalBadCalibFilterUpdateToken_;

  bool                                    isBkgBSample_;
  TFile  *fileWeights;
  TH1F  *histGenWeights;
};
