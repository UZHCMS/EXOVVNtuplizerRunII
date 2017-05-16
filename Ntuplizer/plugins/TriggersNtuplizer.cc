#include "../interface/TriggersNtuplizer.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"


//===================================================================================================================
TriggersNtuplizer::TriggersNtuplizer( edm::EDGetTokenT<edm::TriggerResults> tokens,
                                      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> object,
				      edm::EDGetTokenT<pat::PackedTriggerPrescales> prescale,
				      edm::EDGetTokenT<edm::TriggerResults> noiseFilterToken,
				      edm::EDGetTokenT<bool> HBHENoiseFilterLooseResultToken,
				      edm::EDGetTokenT<bool> HBHENoiseFilterTightResultToken,
				      edm::EDGetTokenT<bool> HBHENoiseIsoFilterResultToken,
				      NtupleBranches* nBranches,
				      const edm::ParameterSet& iConfig,
				      std::map< std::string, bool >& runFlags)
   : CandidateNtuplizer	( nBranches )
   , HLTtriggersToken_	( tokens )
   , triggerObjects_	( object )
   , triggerPrescales_	( prescale )
   , noiseFilterToken_	( noiseFilterToken )
   , HBHENoiseFilterLoose_Selector_( HBHENoiseFilterLooseResultToken )
   , HBHENoiseFilterTight_Selector_( HBHENoiseFilterTightResultToken )
   , HBHENoiseIsoFilter_Selector_( HBHENoiseIsoFilterResultToken )
   , doTriggerDecisions_( runFlags["doTriggerDecisions"] )
   , doTriggerObjects_	( runFlags["doTriggerObjects"] )
   , doHltFilters_	( runFlags["doHltFilters"] )
   , runOnMC_           ( runFlags["runOnMC"] )
{

  HBHENoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_HBHENoiseFilter");
  CSCHaloNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_CSCTightHaloFilter");
  CSCTightHalo2015Filter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_CSCTightHalo2015Filter");
  HCALlaserNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_hcalLaserEventFilter");
  ECALDeadCellNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter");
  GoodVtxNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_goodVertices");
  TrkFailureNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_trackingFailureFilter");
  EEBadScNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_eeBadScFilter");
  ECALlaserNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_ecalLaserCorrFilter");
  TrkPOGNoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_trkPOGFilters");
  TrkPOG_manystrip_NoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_trkPOG_manystripclus53X");
  TrkPOG_toomanystrip_NoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_trkPOG_toomanystripclus53X");
  TrkPOG_logError_NoiseFilter_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_trkPOG_logErrorTooManyClusters");
  METFilters_Selector_ =  iConfig.getParameter<std::string> ("noiseFilterSelection_metFilters");

  //NEW FOR ICHEP
  CSCTightHaloTrkMuUnvetoFilter_Selector_       =  iConfig.getParameter<std::string> ("noiseFilterSelection_CSCTightHaloTrkMuUnvetoFilter");
  globalTightHalo2016Filter_Selector_           =  iConfig.getParameter<std::string> ("noiseFilterSelection_globalTightHalo2016Filter");
  HcalStripHaloFilter_Selector_                 =  iConfig.getParameter<std::string> ("noiseFilterSelection_HcalStripHaloFilter");
  chargedHadronTrackResolutionFilter_Selector_  =  iConfig.getParameter<std::string> ("noiseFilterSelection_chargedHadronTrackResolutionFilter");
  muonBadTrackFilter_Selector_                  =  iConfig.getParameter<std::string> ("noiseFilterSelection_muonBadTrackFilter");
  
  //NEW FOR MORIOND

  badMuons_Selector_  =  iConfig.getParameter<std::string> ("noiseFilterSelection_badMuonsFilter");
  duplicateMuons_Selector_  =  iConfig.getParameter<std::string> ("noiseFilterSelection_duplicateMuonsFilter");
  nobadMuons_Selector_  =  iConfig.getParameter<std::string> ("noiseFilterSelection_nobadMuonsFilter");

}

//===================================================================================================================
TriggersNtuplizer::~TriggersNtuplizer( void )
{

}

//===================================================================================================================
bool TriggersNtuplizer::findTrigger( std::string trigName ){

   if( trigName.find("AK8PFJet360_TrimMass30") != std::string::npos ||
       trigName.find("AK8PFHT700_TrimR0p1PT0p03Mass50") != std::string::npos ||
       trigName.find("AK8PFHT650_TrimR0p1PT0p03Mass50") != std::string::npos ||
       trigName.find("AK8PFHT660_TrimR0p1PT0p03Mass50_BTagCSV_p20") != std::string::npos ||
       trigName.find("AK8DiPFJet280_200_TrimMass30_BTagCSV_p20") != std::string::npos ||
       trigName.find("AK8DiPFJet250_200_TrimMass30_BTagCSV_p20") != std::string::npos ||
       trigName.find("AK8PFJet450") != std::string::npos ||
       trigName.find("PFHT650_WideJetMJJ950DEtaJJ1p5") != std::string::npos ||
       trigName.find("PFHT650_WideJetMJJ900DEtaJJ1p5") != std::string::npos ||
       trigName.find("PFHT400_v") != std::string::npos ||
       trigName.find("PFHT425_v") != std::string::npos ||
       trigName.find("PFHT475_v") != std::string::npos ||
       trigName.find("PFHT600_v") != std::string::npos ||
       trigName.find("PFHT650_v") != std::string::npos ||
       trigName.find("PFHT800_v") != std::string::npos ||
       trigName.find("PFHT900_v") != std::string::npos ||
       trigName.find("PFJet320_v") != std::string::npos ||
       trigName.find("PFJet400_v") != std::string::npos ||
       trigName.find("PFJet450_v") != std::string::npos ||
       trigName.find("PFJet500_v") != std::string::npos ||
       trigName.find("CaloJet500_NoJetID") != std::string::npos ||

       trigName.find("HLT_IsoMu20_eta2p1") != std::string::npos ||
       trigName.find("HLT_IsoMu24_eta2p1") != std::string::npos ||
       trigName.find("HLT_IsoTkMu24_eta2p1") != std::string::npos ||
       trigName.find("HLT_IsoMu24_v") != std::string::npos ||
       trigName.find("HLT_IsoTkMu24_v") != std::string::npos ||
       trigName.find("HLT_IsoMu27_v") != std::string::npos ||
       trigName.find("HLT_IsoTkMu27_v") != std::string::npos ||
       trigName.find("HLT_Mu45_eta2p1") != std::string::npos ||
       //trigName.find("HLT_Mu50_eta2p1") != std::string::npos ||

       trigName.find("HLT_Mu50_v") != std::string::npos ||
       trigName.find("HLT_TkMu50_v") != std::string::npos ||
       trigName.find("HLT_Ele27_WPTight_Gsf") != std::string::npos ||
       trigName.find("HLT_Ele27_WPLoose_Gsf") != std::string::npos ||
       trigName.find("HLT_Ele27_WPLoose_Gsf_WHbbBoost") != std::string::npos ||
       trigName.find("HLT_Ele27_eta2p1_WPLoose") != std::string::npos ||
       trigName.find("HLT_Ele27_eta2p1_WP75_Gsf") != std::string::npos ||
       trigName.find("HLT_Ele23_CaloIdL_TrackIdL_IsoVL") != std::string::npos ||
       trigName.find("HLT_Ele32_eta2p1_WP75_Gsf") != std::string::npos ||
       trigName.find("HLT_Ele45_WPLoose_Gsf_v") != std::string::npos ||
       trigName.find("HLT_Ele105_CaloIdVT_GsfTrkIdT") != std::string::npos ||
       trigName.find("HLT_Ele115_CaloIdVT_GsfTrkIdT") != std::string::npos ||
       trigName.find("HLT_Ele145_CaloIdVT_GsfTrkIdT") != std::string::npos ||
       trigName.find("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165") != std::string::npos ||
       trigName.find("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20") != std::string::npos ||
       trigName.find("HLT_Ele27_eta2p1_WPLoose_Gsf_DoubleMediumIsoPFTau35") != std::string::npos ||
       trigName.find("HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20") != std::string::npos ||
       trigName.find("HLT_IsoMu16_eta2p1_MET30_JetIDCleaned_LooseIsoPFTau50") != std::string::npos ||
       trigName.find("HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1_v1") != std::string::npos ||
       //H->tautau triggers
       trigName.find("HLT_LooseIsoPFTau50") != std::string::npos||
       trigName.find("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_v") != std::string::npos||
       trigName.find("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v") != std::string::npos||
       trigName.find("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v") != std::string::npos||
       trigName.find("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v") != std::string::npos||
       trigName.find("HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_v") != std::string::npos||
       trigName.find("HLT_IsoMu24_eta2p1") != std::string::npos||
       trigName.find("HLT_IsoMu17_eta2p1") != std::string::npos||
       trigName.find("HLT_IsoMu18_v") != std::string::npos||
       trigName.find("HLT_IsoMu22_v") != std::string::npos||
       trigName.find("HLT_IsoMu27_v") != std::string::npos||
       trigName.find("HLT_IsoMu20_v") != std::string::npos||
       trigName.find("HLT_IsoMu22_eta2p1_v") != std::string::npos||
       trigName.find("HLT_IsoMu24_v") != std::string::npos||
       trigName.find("HLT_IsoTkMu18_v") != std::string::npos||
       trigName.find("HLT_IsoTkMu20_v") != std::string::npos||
       trigName.find("HLT_IsoTkMu22_v") != std::string::npos||
       trigName.find("HLT_IsoTkMu22_eta2p1_v") != std::string::npos||
       trigName.find("HLT_IsoTkMu24_v") != std::string::npos||
       trigName.find("HLT_IsoTkMu27_v") != std::string::npos||
       trigName.find("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v") != std::string::npos||
       trigName.find("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v") != std::string::npos||
       trigName.find("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v") != std::string::npos||
       trigName.find("HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v") != std::string::npos||
       trigName.find("HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v") != std::string::npos||
       trigName.find("HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20") != std::string::npos||
       trigName.find("HLT_Ele25_eta2p1") != std::string::npos||
       trigName.find("HLT_Ele22_eta2p1_WP75_Gsf_v") != std::string::npos||
       trigName.find("HLT_Ele32_eta2p1_WP75_Gsf") != std::string::npos||
       trigName.find("HLT_Ele23_WPLoose_Gsf") != std::string::npos||
       trigName.find("HLT_Ele23_WPTight_Gsf") != std::string::npos||
       trigName.find("HLT_Ele24_eta2p1_WPLoose_Gsf_v") != std::string::npos||
       trigName.find("HLT_Ele25_WPTight_Gsf_v") != std::string::npos||
       trigName.find("HLT_Ele25_eta2p1_WPLoose_Gsf_v") != std::string::npos||
       trigName.find("HLT_Ele25_eta2p1_WPTight_Gsf_v") != std::string::npos||
       trigName.find("HLT_Ele27_WPLoose_Gsf") != std::string::npos||
       trigName.find("HLT_Ele27_WPTight_Gsf") != std::string::npos||
       trigName.find("HLT_Ele27_eta2p1_WPLoose_Gsf_v") != std::string::npos||
       trigName.find("HLT_Ele27_eta2p1_WPTight_Gsf_v") != std::string::npos||
       trigName.find("HLT_Ele32_eta2p1_WPTight_Gsf_v") != std::string::npos||
       trigName.find("HLT_Ele32_WPTight_Gsf_v") != std::string::npos||
       trigName.find("HLT_Ele45_WPLoose_Gsf_L1JetTauSeeded") != std::string::npos||
       trigName.find("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg") != std::string::npos||
       trigName.find("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg") != std::string::npos||
       // Double leptons
       trigName.find("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") != std::string::npos||
       trigName.find("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v") != std::string::npos||
       trigName.find("HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL_v") != std::string::npos||
       trigName.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") != std::string::npos||
       trigName.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") != std::string::npos||
       trigName.find("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v") != std::string::npos||
       trigName.find("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v") != std::string::npos||
       trigName.find("HLT_Mu30_TkMu11_v") != std::string::npos||
       trigName.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos||
       trigName.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos||
       trigName.find("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos||
       trigName.find("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos||
       trigName.find("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v") != std::string::npos||
       trigName.find("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") != std::string::npos||
       trigName.find("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") != std::string::npos||
       // MET triggers
       trigName.find("HLT_PFMET110_PFMHT110_IDTight_v") != std::string::npos||
       trigName.find("HLT_PFMET120_PFMHT120_IDTight_v") != std::string::npos||
       trigName.find("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v") != std::string::npos||
       trigName.find("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v") != std::string::npos||
       trigName.find("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v") != std::string::npos||
       trigName.find("HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight_v") != std::string::npos||
       trigName.find("HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v") != std::string::npos||
       trigName.find("HLT_PFMET120_BTagCSV_p067_v") != std::string::npos||
       trigName.find("HLT_PFMET170_NoiseCleaned_v") != std::string::npos||
       trigName.find("HLT_PFMET170_HBHECleaned_v") != std::string::npos||
       trigName.find("HLT_PFMET170_HBHE_BeamHaloCleaned_v") != std::string::npos||
       //Alternative triggers
       trigName.find("HLT_MonoCentralPFJet80_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight")!= std::string::npos||
       trigName.find("HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight")!= std::string::npos
   ) return true;
   else
     return false;
}


bool TriggersNtuplizer::findFilter( std::string filterName ){

   if( filterName.find("hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09") != std::string::npos ||
       filterName.find("hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09") != std::string::npos ||
       filterName.find("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09") != std::string::npos ||
       filterName.find("hltL3crIsoL1sSingleMu20erL1f0L2f10QL3f22QL3trkIsoFiltered0p09") != std::string::npos ||
       filterName.find("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09") != std::string::npos ||
       filterName.find("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p09") != std::string::npos ||
       filterName.find("hltL3fL1sMu16L1f0Tkf18QL3trkIsoFiltered0p09") != std::string::npos ||
       filterName.find("hltL3fL1sMu18L1f0Tkf20QL3trkIsoFiltered0p09") != std::string::npos ||
       filterName.find("hltL3fL1sMu20erL1f0Tkf22QL3trkIsoFiltered0p09") != std::string::npos ||
       filterName.find("hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09") != std::string::npos ||
       filterName.find("hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09") != std::string::npos ||
       filterName.find("hltL3fL1sMu22Or25L1f0Tkf27QL3trkIsoFiltered0p09") != std::string::npos ||
       filterName.find("hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09") != std::string::npos ||
       filterName.find("hltOverlapFilterSingleIsoMu17LooseIsoPFTau20") != std::string::npos ||
       filterName.find("hltL3crIsoL1sMu16erTauJet20erL1f0L2f10QL3f17QL3trkIsoFiltered0p09") != std::string::npos ||
       filterName.find("hltOverlapFilterIsoMu17LooseIsoPFTau20") != std::string::npos ||
       filterName.find("hltL3crIsoL1sSingleMu18erIorSingleMu20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09") != std::string::npos ||
       filterName.find("hltOverlapFilterSingleIsoMu19LooseIsoPFTau20") != std::string::npos ||
       filterName.find("hltL3crIsoL1sMu18erTauJet20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09") != std::string::npos ||
       filterName.find("hltOverlapFilterIsoMu19LooseIsoPFTau20") != std::string::npos ||
       filterName.find("hltL3crIsoL1sSingleMu20erIorSingleMu22erL1f0L2f10QL3f21QL3trkIsoFiltered0p09") != std::string::npos ||
       filterName.find("hltOverlapFilterSingleIsoMu21LooseIsoPFTau20") != std::string::npos ||
       filterName.find("hltPFTau20TrackLooseIsoAgainstMuon") != std::string::npos ||
       filterName.find("hltEle23WPLooseGsfTrackIsoFilter") != std::string::npos ||
       filterName.find("hltSingleEle24WPLooseGsfTrackIsoFilter") != std::string::npos ||
       filterName.find("hltEle25WPTightGsfTrackIsoFilter") != std::string::npos ||
       filterName.find("hltEle25erWPLooseGsfTrackIsoFilter") != std::string::npos ||
       filterName.find("hltEle25erWPTightGsfTrackIsoFilter") != std::string::npos ||
       filterName.find("hltEle27noerWPLooseGsfTrackIsoFilter") != std::string::npos ||
       filterName.find("hltEle27WPTightGsfTrackIsoFilter") != std::string::npos ||
       filterName.find("hltEle27erWPLooseGsfTrackIsoFilter") != std::string::npos ||
       filterName.find("hltEle27erWPTightGsfTrackIsoFilter") != std::string::npos ||
       filterName.find("hltEle32WPTightGsfTrackIsoFilter") != std::string::npos ||
       filterName.find("hltEle22WPLooseL1SingleIsoEG20erGsfTrackIsoFilter") != std::string::npos ||
       filterName.find("hltOverlapFilterSingleIsoEle22WPLooseGsfLooseIsoPFTau20") != std::string::npos ||
       filterName.find("hltEle24WPLooseL1SingleIsoEG22erGsfTrackIsoFilter") != std::string::npos ||
       filterName.find("hltOverlapFilterSingleIsoEle24WPLooseGsfLooseIsoPFTau20") != std::string::npos ||
       filterName.find("hltEle24WPLooseL1IsoEG22erTau20erGsfTrackIsoFilter") != std::string::npos ||
       filterName.find("hltOverlapFilterIsoEle24WPLooseGsfLooseIsoPFTau20") != std::string::npos ||
       filterName.find("hltOverlapFilterIsoEle27WPLooseGsfLooseIsoPFTau20") != std::string::npos ||
       filterName.find("hltEle32WPLooseGsfTrackIsoFilter") != std::string::npos ||
       filterName.find("hltOverlapFilterIsoEle32WPLooseGsfLooseIsoPFTau20") != std::string::npos ||
       filterName.find("hltPFTau20TrackLooseIso") != std::string::npos ||
       filterName.find("hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter") != std::string::npos ||
       filterName.find("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter") != std::string::npos ||
       filterName.find("hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter") != std::string::npos ||
       filterName.find("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter") != std::string::npos ||
       filterName.find("hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter") != std::string::npos ||
       filterName.find("hltEle27erWPLooseGsfTrackIsoFilter") != std::string::npos ||
       filterName.find("hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8") != std::string::npos ||
       filterName.find("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8") != std::string::npos ||
       filterName.find("hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered17") != std::string::npos ||
       filterName.find("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23") != std::string::npos ||
       filterName.find("hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23") != std::string::npos ||
       filterName.find("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter") != std::string::npos ||
       filterName.find("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLDZFilter") != std::string::npos || 
       filterName.find("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4") != std::string::npos || 
       filterName.find("hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4") != std::string::npos
       ) return true;
   else
     return false;   
}



//===================================================================================================================
void TriggersNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  event.getByToken(HLTtriggersToken_, HLTtriggers_);
  event.getByToken(triggerObjects_  , triggerObjects);
  event.getByToken(triggerPrescales_, triggerPrescales);

  const edm::TriggerNames& trigNames = event.triggerNames(*HLTtriggers_);


  if (doTriggerDecisions_) {
  	 for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {
  	  if( findTrigger(trigNames.triggerName(i)) ){
   	     nBranches_->HLT_isFired[trigNames.triggerName(i)] = HLTtriggers_->accept(i);
	     //	     std::cout << "Trigger " << trigNames.triggerName(i) << ": " << (HLTtriggers_->accept(i) ? "PASS" : "fail (or not run)") << std::endl;
   	  }
   	}

  } //doTriggerDecisions_

  ////////////////// Trigger objects ///////////////////////////////////
  if (doTriggerObjects_) {

     	std::vector<int> vfiredTrigger; vfiredTrigger.clear();
	std::vector<float> vfilterIDs; vfilterIDs.clear();

  	for (pat::TriggerObjectStandAlone obj : *triggerObjects) {

  		obj.unpackPathNames(trigNames);

  		std::vector<std::string> pathNamesAll  = obj.pathNames(false);
		//  		std::vector<std::string> pathNamesLast = obj.pathNames(true);

		//		std::cout << "Size of pathNames All = " << pathNamesAll.size() << ", Last = " << pathNamesLast.size()  << std::endl;
		//		for (unsigned h = 0, n = pathNamesLast.size(); h < n; ++h) {
		//		  std::cout << "\t Lastname = " << pathNamesLast[h] << std::endl;
		//		}
		for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {

		  //		  bool isBoth = obj.hasPathName( pathNamesAll[h], true , true );
		  //		  bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );
		  //		  bool isBoth = obj.hasPathName( pathNamesLast[h], true , true );
		  //		  bool isL3   = obj.hasPathName( pathNamesLast[h], false, true );

		  //		  std::cout << "trigger object =" << pathNamesLast[h] << "(isBoth, isL3) = " << isBoth << " " << isL3 << std::endl;
		  //		  std::cout << "\t Finalname = " << pathNamesAll[h] << "(isBoth, isL3) = " << isBoth << " " << isL3 << std::endl;
			
		//			if( isBoth || isL3 ){


		  for (unsigned hh = 0; hh < obj.filterIds().size(); ++hh) vfilterIDs.push_back( obj.filterIds()[hh]); // as defined in http://cmslxr.fnal.gov/lxr/source/DataFormats/HLTReco/interface/TriggerTypeDefs.h

			   std::vector<std::string> vfilterLabels; vfilterLabels.clear();

			   bool isFilterExist = false;

			   for (unsigned hh = 0; hh < obj.filterLabels().size(); ++hh){
			     if(findFilter(obj.filterLabels()[hh])){
			       vfilterLabels.push_back( obj.filterLabels()[hh]);
			       isFilterExist = true;
			     }
			   }

			   
			   if(isFilterExist){
			     nBranches_->triggerObject_pt  .push_back(obj.pt());
			     nBranches_->triggerObject_eta .push_back(obj.eta());
			     nBranches_->triggerObject_phi .push_back(obj.phi());
			     nBranches_->triggerObject_mass.push_back(obj.mass());
			     nBranches_->triggerObject_lastname.push_back(pathNamesAll[h]);
			     //			     nBranches_->triggerObject_lastname.push_back("test");

			     //			     if(nBranches_->triggerObject_filterLabels.find(pathNamesLast[h]) == nBranches_->triggerObject_filterLabels.end()){
			     if(nBranches_->triggerObject_filterLabels.find(pathNamesAll[h]) == nBranches_->triggerObject_filterLabels.end()){
			       //  std::cout << "index NOT found !!!" << std::endl;
			       nBranches_->triggerObject_filterLabels[pathNamesAll[h]] = vfilterLabels;
			       //			       nBranches_->triggerObject_filterLabels["test"] = vfilterLabels;
			     }else{
			       // add to the original
			       // std::cout << "index found !!!" << std::endl;
			       
			       //			       std::vector<std::string> vec1 = nBranches_->triggerObject_filterLabels["test"];
			       std::vector<std::string> vec1 = nBranches_->triggerObject_filterLabels[pathNamesAll[h]];
			       
//			       std::cout << "before : size of vec1 = " << vec1.size() << std::endl;
//			       for(int ii = 0; ii < (int)vec1.size(); ii++){
//				 std::cout << "\t\t " << vec1.at(ii) << std::endl;
//			       }

			       vec1.insert(vec1.end(), vfilterLabels.begin(), vfilterLabels.end());

//			       std::cout << "after : size of vec1 = " << vec1.size() << std::endl;
//			       for(int ii = 0; ii < (int)vec1.size(); ii++){
//				 std::cout << "\t\t " << vec1.at(ii) << std::endl;
//			       }

			       nBranches_->triggerObject_filterLabels[pathNamesAll[h]] = vec1;
			     }

			   }


  			   if( pathNamesAll[h] == "HLT_AK8PFJet360_TrimMass30_v1") vfiredTrigger.push_back( 0 );
  			   if( pathNamesAll[h] == "HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v1") vfiredTrigger.push_back( 1 );
  			   if( pathNamesAll[h] == "HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p41_v1") vfiredTrigger.push_back( 2 );
  			   if( pathNamesAll[h] == "HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v1") vfiredTrigger.push_back( 3 );
  			   if( pathNamesAll[h] == "HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v1") vfiredTrigger.push_back( 4 );
  			   if( pathNamesAll[h] == "HLT_PFHT900_v1") vfiredTrigger.push_back( 5 );
  			   if( pathNamesAll[h] == "HLT_IsoMu24_eta2p1_v1") vfiredTrigger.push_back( 6 );
  			   if( pathNamesAll[h] == "HLT_IsoMu24_eta2p1_v2") vfiredTrigger.push_back( 7 );
  			   if( pathNamesAll[h] == "HLT_Mu45_eta2p1_v1") vfiredTrigger.push_back( 8 );
  			   if( pathNamesAll[h] == "HLT_Mu50_eta2p1_v1") vfiredTrigger.push_back( 9 );
  			   if( pathNamesAll[h] == "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v1") vfiredTrigger.push_back( 10 );
  			   if( pathNamesAll[h] == "HLT_Ele32_eta2p1_WP75_Gsf_v1") vfiredTrigger.push_back( 11 );
  			   if( pathNamesAll[h] == "HLT_Ele105_CaloIdVT_GsfTrkIdT_v1") vfiredTrigger.push_back( 12 );
  			   if( pathNamesAll[h] == "HLT_Ele105_CaloIdVT_GsfTrkIdT_v2") vfiredTrigger.push_back( 13 );
  			   if( pathNamesAll[h] == "HLT_Ele115_CaloIdVT_GsfTrkIdT_v1") vfiredTrigger.push_back( 14 );


  			   // else vfiredTrigger.push_back( -99 );
			   //  			}
			
		}

  		nBranches_->triggerObject_firedTrigger.push_back(vfiredTrigger);
		nBranches_->triggerObject_filterIDs.push_back(vfilterIDs);
  	}
  } //doTriggerObjects_


  // HLT Noise Filters
  // for deprecation see https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
  if (doHltFilters_) {

    event.getByToken(noiseFilterToken_, noiseFilterBits_);
    const edm::TriggerNames &names = event.triggerNames(*noiseFilterBits_);

    for (unsigned int i = 0, n = noiseFilterBits_->size(); i < n; ++i) {
      if (names.triggerName(i) == HBHENoiseFilter_Selector_)
        nBranches_->passFilter_HBHE_ = noiseFilterBits_->accept(i);
      if (names.triggerName(i) == CSCHaloNoiseFilter_Selector_)
        nBranches_->passFilter_CSCHalo_ = noiseFilterBits_->accept(i); // DEPRECATED
      if (names.triggerName(i) == CSCTightHalo2015Filter_Selector_)
        nBranches_->passFilter_CSCTightHalo2015_ = noiseFilterBits_->accept(i); // TO BE USED
      if (names.triggerName(i) == HCALlaserNoiseFilter_Selector_)
        nBranches_->passFilter_HCALlaser_ = noiseFilterBits_->accept(i); // DEPRECATED
      if (names.triggerName(i) == ECALDeadCellNoiseFilter_Selector_)
        nBranches_->passFilter_ECALDeadCell_ = noiseFilterBits_->accept(i); // under scrutiny
      if (names.triggerName(i) == GoodVtxNoiseFilter_Selector_)
        nBranches_->passFilter_GoodVtx_ = noiseFilterBits_->accept(i); // TO BE USED
      if (names.triggerName(i) == TrkFailureNoiseFilter_Selector_)
        nBranches_->passFilter_TrkFailure_ = noiseFilterBits_->accept(i); // DEPRECATED
      if (names.triggerName(i) == EEBadScNoiseFilter_Selector_)
        nBranches_->passFilter_EEBadSc_ = noiseFilterBits_->accept(i); // under scrutiny
      if (names.triggerName(i) == ECALlaserNoiseFilter_Selector_)
        nBranches_->passFilter_ECALlaser_ = noiseFilterBits_->accept(i); // DEPRECATED
      if (names.triggerName(i) == TrkPOGNoiseFilter_Selector_)
        nBranches_->passFilter_TrkPOG_ = noiseFilterBits_->accept(i); // DEPRECATED
      if (names.triggerName(i) == TrkPOG_manystrip_NoiseFilter_Selector_)
        nBranches_->passFilter_TrkPOG_manystrip_ = noiseFilterBits_->accept(i); // DEPRECATED
      if (names.triggerName(i) == TrkPOG_toomanystrip_NoiseFilter_Selector_)
        nBranches_->passFilter_TrkPOG_toomanystrip_ = noiseFilterBits_->accept(i); // DEPRECATED
      if (names.triggerName(i) == TrkPOG_logError_NoiseFilter_Selector_)
        nBranches_->passFilter_TrkPOG_logError_ = noiseFilterBits_->accept(i); // DEPRECATED
      if (names.triggerName(i) == METFilters_Selector_)
        nBranches_->passFilter_METFilters_ = noiseFilterBits_->accept(i); // DEPRECATED
       //NEW FOR ICHEP
      if (names.triggerName(i) == CSCTightHaloTrkMuUnvetoFilter_Selector_)
        nBranches_->passFilter_CSCTightHaloTrkMuUnvetoFilter_ = noiseFilterBits_->accept(i);
      if (names.triggerName(i) == globalTightHalo2016Filter_Selector_           )
        nBranches_->passFilter_globalTightHalo2016_ = noiseFilterBits_->accept(i); // TO BE USED FOR ICHEP 2016
      if (names.triggerName(i) == HcalStripHaloFilter_Selector_                 )
        nBranches_->passFilter_HcalStripHalo_ = noiseFilterBits_->accept(i);
      if (names.triggerName(i) == chargedHadronTrackResolutionFilter_Selector_  )
        nBranches_->passFilter_chargedHadronTrackResolution_ = noiseFilterBits_->accept(i); // TO BE USED FOR ICHEP 2016
      if (names.triggerName(i) == muonBadTrackFilter_Selector_                  )
        nBranches_->passFilter_muonBadTrack_ = noiseFilterBits_->accept(i); // TO BE USED FOR ICHEP 2016
      //NEW FOR MORIOND
      if (names.triggerName(i) == badMuons_Selector_                  )
        nBranches_->flag_badMuons_ = noiseFilterBits_->accept(i); // TO BE USED FOR ICHEP 2016

      if (names.triggerName(i) == duplicateMuons_Selector_                  )
        nBranches_->flag_duplicateMuons_ = noiseFilterBits_->accept(i); // TO BE USED FOR ICHEP 2016

      if (names.triggerName(i) == nobadMuons_Selector_                  )
        nBranches_->flag_nobadMuons_ = noiseFilterBits_->accept(i); // TO BE USED FOR ICHEP 2016

    }

    //if( !runOnMC_ /*&& event.id().run() < 251585*/ ){

       edm::Handle<bool> HBHENoiseFilterLooseResultHandle;
       event.getByToken(HBHENoiseFilterLoose_Selector_, HBHENoiseFilterLooseResultHandle);
       bool HBHENoiseFilterLooseResult = *HBHENoiseFilterLooseResultHandle;
       if (!HBHENoiseFilterLooseResultHandle.isValid()) {
         LogDebug("") << "CaloTowerAnalyzer: Could not find HBHENoiseFilterResult" << std::endl;
       }

       nBranches_->passFilter_HBHELoose_ = HBHENoiseFilterLooseResult;

       edm::Handle<bool> HBHENoiseFilterTightResultHandle;
       event.getByToken(HBHENoiseFilterTight_Selector_, HBHENoiseFilterTightResultHandle);
       bool HBHENoiseFilterTightResult = *HBHENoiseFilterTightResultHandle;
       if (!HBHENoiseFilterTightResultHandle.isValid()) {
         LogDebug("") << "CaloTowerAnalyzer: Could not find HBHENoiseFilterResult" << std::endl;
       }

       nBranches_->passFilter_HBHETight_ = HBHENoiseFilterTightResult;

        edm::Handle<bool> HBHENoiseIsoFilterResultHandle;
       event.getByToken(HBHENoiseIsoFilter_Selector_, HBHENoiseIsoFilterResultHandle);
       bool HBHENoiseIsoFilterResult = *HBHENoiseIsoFilterResultHandle;
       if (!HBHENoiseIsoFilterResultHandle.isValid()) {
         LogDebug("") << "CaloTowerAnalyzer: Could not find HBHENoiseFilterResult" << std::endl;
       }

       nBranches_->passFilter_HBHEIso_ = HBHENoiseIsoFilterResult;

    //}

  } //doHltFilters_
}
