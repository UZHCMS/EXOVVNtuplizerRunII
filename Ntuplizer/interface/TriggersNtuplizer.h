#ifndef TriggersNtuplizer_H
#define TriggersNtuplizer_H

#include "../interface/CandidateNtuplizer.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

class TriggersNtuplizer : public CandidateNtuplizer {

public:
   TriggersNtuplizer(edm::EDGetTokenT<edm::TriggerResults> tokens, 
                     edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> object, 
		     edm::EDGetTokenT<pat::PackedTriggerPrescales> prescale, 
		     edm::EDGetTokenT<edm::TriggerResults> noiseFilterToken,
		     edm::EDGetTokenT<bool> HBHENoiseFilterLooseResultToken, 
		     edm::EDGetTokenT<bool> HBHENoiseFilterTightResultToken, 
		     edm::EDGetTokenT<bool> HBHENoiseIsoFilterResultToken, 
		     NtupleBranches* nBranches, 
		     const edm::ParameterSet& iConfig, 
		     std::map< std::string, bool >& runFlags );
   ~TriggersNtuplizer( void );
   
  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
  bool findTrigger( std::string trigName );

private:
   edm::EDGetTokenT<edm::TriggerResults> 		     HLTtriggersToken_;
   edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection>  triggerObjects_;
   edm::EDGetTokenT<pat::PackedTriggerPrescales>             triggerPrescales_;
   edm::EDGetTokenT<edm::TriggerResults> 		     noiseFilterToken_;
	
   edm::Handle< edm::TriggerResults> 			     HLTtriggers_;
   edm::Handle<pat::TriggerObjectStandAloneCollection>	     triggerObjects;
   edm::Handle<pat::PackedTriggerPrescales>	      	     triggerPrescales;
   edm::Handle< edm::TriggerResults> 			     noiseFilterBits_;
   
   // HLT Noise Filter names
   std::string HBHENoiseFilter_Selector_;
   edm::EDGetTokenT<bool> HBHENoiseFilterLoose_Selector_;
   edm::EDGetTokenT<bool> HBHENoiseFilterTight_Selector_;
   edm::EDGetTokenT<bool> HBHENoiseIsoFilter_Selector_;
   std::string CSCHaloNoiseFilter_Selector_;
   std::string CSCTightHalo2015Filter_Selector_;
   std::string HCALlaserNoiseFilter_Selector_;
   std::string ECALDeadCellNoiseFilter_Selector_;
   std::string GoodVtxNoiseFilter_Selector_;
   std::string TrkFailureNoiseFilter_Selector_;
   std::string EEBadScNoiseFilter_Selector_;
   std::string ECALlaserNoiseFilter_Selector_;
   std::string TrkPOGNoiseFilter_Selector_;
   std::string TrkPOG_manystrip_NoiseFilter_Selector_;
   std::string TrkPOG_toomanystrip_NoiseFilter_Selector_;
   std::string TrkPOG_logError_NoiseFilter_Selector_;
   std::string METFilters_Selector_;
   //NEW FOR ICHEP
   std::string CSCTightHaloTrkMuUnvetoFilter_Selector_  ;
   std::string globalTightHalo2016Filter_Selector_  ;
   std::string HcalStripHaloFilter_Selector_  ;
   std::string chargedHadronTrackResolutionFilter_Selector_ ;
   std::string muonBadTrackFilter_Selector_ ;
   
   bool doTriggerDecisions_;
   bool doTriggerObjects_;
   bool doHltFilters_;
   bool runOnMC_;
        
};

#endif // TriggersNtuplizer_H



