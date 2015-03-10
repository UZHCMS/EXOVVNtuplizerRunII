#include "../interface/TriggersNtuplizer.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

//===================================================================================================================        
TriggersNtuplizer::TriggersNtuplizer(edm::EDGetTokenT<edm::TriggerResults> tokens, edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> object, edm::EDGetTokenT<pat::PackedTriggerPrescales> prescale,  NtupleBranches* nBranches )
   : CandidateNtuplizer	( nBranches )
   , HLTtriggersToken_	( tokens )
   , triggerObjects_	( object )	
   , triggerPrescales_	( prescale )		
{
   
}

//===================================================================================================================
TriggersNtuplizer::~TriggersNtuplizer( void )
{

}

//===================================================================================================================
void TriggersNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
	
	event.getByToken(HLTtriggersToken_, HLTtriggers_);	
	event.getByToken(triggerObjects_  , triggerObjects);
	event.getByToken(triggerPrescales_, triggerPrescales);
	
	const edm::TriggerNames& trigNames = event.triggerNames(*HLTtriggers_);
	
	
	// for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {
// 	  std::cout << "Trigger " << trigNames.triggerName(i) << std::endl;//<< ": " << (HLTtriggers_->accept(i) ? "PASS" : "fail (or not run)") << std::endl;
// 	}
	
	////////////////// dijet triggers ///////////////////////////////////
	unsigned int TrggIndex_HLT_AK8PFJet360TrimMod_Mass30_v1						(trigNames.triggerIndex("HLT_AK8PFJet360TrimMod_Mass30_v1"						));
	unsigned int TrggIndex_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v1				(trigNames.triggerIndex("HT_AK8PFHT700_TrimR0p1PT0p03Mass50_v1"				));
	unsigned int TrggIndex_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p3_v1	(trigNames.triggerIndex("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p3_v1"	));
	unsigned int TrggIndex_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v1				(trigNames.triggerIndex("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v1"				));
	unsigned int TrggIndex_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v1				(trigNames.triggerIndex("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v1"				));
	unsigned int TrggIndex_HLT_PFHT900_v1												(trigNames.triggerIndex("HLT_PFHT900_v1"												));

	nBranches_->isFired_HLT_AK8PFJet360TrimMod_Mass30_v1               = false;
	nBranches_->isFired_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v1         = false;
	nBranches_->isFired_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p3_v1 = false;
	nBranches_->isFired_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v1          = false;
	nBranches_->isFired_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v1          = false;
	nBranches_->isFired_HLT_PFHT900_v1                                 = false;
	
	if( TrggIndex_HLT_AK8PFJet360TrimMod_Mass30_v1						< HLTtriggers_->size() ) nBranches_->isFired_HLT_AK8PFJet360TrimMod_Mass30_v1						= HLTtriggers_->accept(TrggIndex_HLT_AK8PFJet360TrimMod_Mass30_v1						);
	if( TrggIndex_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v1				< HLTtriggers_->size() ) nBranches_->isFired_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v1				= HLTtriggers_->accept(TrggIndex_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v1				);      
	if( TrggIndex_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p3_v1	< HLTtriggers_->size() ) nBranches_->isFired_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p3_v1	= HLTtriggers_->accept(TrggIndex_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p3_v1	);
	if( TrggIndex_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v1				< HLTtriggers_->size() ) nBranches_->isFired_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v1				= HLTtriggers_->accept(TrggIndex_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v1				);
	if( TrggIndex_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v1				< HLTtriggers_->size() ) nBranches_->isFired_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v1				= HLTtriggers_->accept(TrggIndex_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v1				); 
	if( TrggIndex_HLT_PFHT900_v1												< HLTtriggers_->size() ) nBranches_->isFired_HLT_PFHT900_v1												= HLTtriggers_->accept(TrggIndex_HLT_PFHT900_v1												);

	////////////////// single lepton triggers ///////////////////////////////////
	unsigned int TrggIndex_HLT_Ele95_CaloIdVT_GsfTrkIdT_v1							(trigNames.triggerIndex("HLT_Ele95_CaloIdVT_GsfTrkIdT_v1"							));
	unsigned int TrggIndex_HLT_Mu40_v1														(trigNames.triggerIndex("HLT_Mu40_v1"														));
	unsigned int TrggIndex_HLT_IsoMu20_eta2p1_IterTrk02_v1							(trigNames.triggerIndex("HLT_IsoMu20_eta2p1_IterTrk02_v1"							));
	unsigned int TrggIndex_HLT_IsoMu20_eta2p1_IterTrk02_TriCentralPFJet40_v1	(trigNames.triggerIndex("HLT_IsoMu20_eta2p1_IterTrk02_TriCentralPFJet40_v1"	));
	unsigned int TrggIndex_HLT_Ele27_eta2p1_WP85_Gsf_v1								(trigNames.triggerIndex("HLT_Ele27_eta2p1_WP85_Gsf_v1"								));
	unsigned int TrggIndex_HLT_Ele27_eta2p1_WP85_Gsf_TriCentralPFJet40_v1		(trigNames.triggerIndex("HLT_Ele27_eta2p1_WP85_Gsf_TriCentralPFJet40_v1"		));
	
	nBranches_->isFired_HLT_Ele95_CaloIdVT_GsfTrkIdT_v1							= false;
	nBranches_->isFired_HLT_Mu40_v1 														= false;
	nBranches_->isFired_HLT_IsoMu20_eta2p1_IterTrk02_v1							= false;
	nBranches_->isFired_HLT_IsoMu20_eta2p1_IterTrk02_TriCentralPFJet40_v1	= false;
	nBranches_->isFired_HLT_Ele27_eta2p1_WP85_Gsf_v1								= false;
	nBranches_->isFired_HLT_Ele27_eta2p1_WP85_Gsf_TriCentralPFJet40_v1		= false;

	if( TrggIndex_HLT_Ele95_CaloIdVT_GsfTrkIdT_v1							< HLTtriggers_->size() ) nBranches_->isFired_HLT_Ele95_CaloIdVT_GsfTrkIdT_v1							= HLTtriggers_->accept(TrggIndex_HLT_Ele95_CaloIdVT_GsfTrkIdT_v1	);
	if( TrggIndex_HLT_Mu40_v1														< HLTtriggers_->size() ) nBranches_->isFired_HLT_Mu40_v1														= HLTtriggers_->accept(TrggIndex_HLT_Mu40_v1								);
	if( TrggIndex_HLT_IsoMu20_eta2p1_IterTrk02_v1							< HLTtriggers_->size() ) nBranches_->isFired_HLT_IsoMu20_eta2p1_IterTrk02_v1							= HLTtriggers_->accept(TrggIndex_HLT_IsoMu20_eta2p1_IterTrk02_v1	);
	if( TrggIndex_HLT_IsoMu20_eta2p1_IterTrk02_TriCentralPFJet40_v1	< HLTtriggers_->size() ) nBranches_->isFired_HLT_IsoMu20_eta2p1_IterTrk02_TriCentralPFJet40_v1	= HLTtriggers_->accept(TrggIndex_HLT_IsoMu20_eta2p1_IterTrk02_TriCentralPFJet40_v1);
	if( TrggIndex_HLT_Ele27_eta2p1_WP85_Gsf_TriCentralPFJet40_v1		< HLTtriggers_->size() ) nBranches_->isFired_HLT_Ele27_eta2p1_WP85_Gsf_TriCentralPFJet40_v1		= HLTtriggers_->accept(TrggIndex_HLT_Ele27_eta2p1_WP85_Gsf_TriCentralPFJet40_v1);
	if( TrggIndex_HLT_Ele27_eta2p1_WP85_Gsf_v1								< HLTtriggers_->size() ) nBranches_->isFired_HLT_Ele27_eta2p1_WP85_Gsf_v1								= HLTtriggers_->accept(TrggIndex_HLT_Ele27_eta2p1_WP85_Gsf_v1		);
	
		
	////////////////// Trigger objects ///////////////////////////////////	
	

   std::vector<float> vfilterIDs			;
   vfilterIDs.clear();
   std::vector<int> vfiredTrigger		;
   vfiredTrigger.clear();

	
	for (pat::TriggerObjectStandAlone obj : *triggerObjects) { 
	
		obj.unpackPathNames(trigNames);
		
		std::vector<std::string> pathNamesAll  = obj.pathNames(false);
		std::vector<std::string> pathNamesLast = obj.pathNames(true);
		
		for (unsigned h = 0, n = pathNamesLast.size(); h < n; ++h) {
		
			bool isBoth = obj.hasPathName( pathNamesLast[h], true , true );
                        bool isL3   = obj.hasPathName( pathNamesLast[h], false, true );
			
			if( isBoth || isL3 ){

				nBranches_->triggerObject_pt				.push_back(obj.pt());
			   nBranches_->triggerObject_eta				.push_back(obj.eta());
			   nBranches_->triggerObject_phi				.push_back(obj.phi());
			   nBranches_->triggerObject_mass			.push_back(obj.mass());

				
			   for (unsigned h = 0; h < obj.filterIds().size(); ++h) vfilterIDs.push_back( obj.filterIds()[h]); // as defined in http://cmslxr.fnal.gov/lxr/source/DataFormats/HLTReco/interface/TriggerTypeDefs.h
				
			      if( pathNamesLast[h] == "HLT_AK8PFJet360TrimMod_Mass30_v1") vfiredTrigger.push_back( 0 );
			      if( pathNamesLast[h] == "HT_AK8PFHT700_TrimR0p1PT0p03Mass50_v1") vfiredTrigger.push_back( 1 );
			      if( pathNamesLast[h] == "HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p3_v1")	vfiredTrigger.push_back( 2 );
			      if( pathNamesLast[h] == "HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v1") vfiredTrigger.push_back( 3 );
			      if( pathNamesLast[h] == "HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v1") vfiredTrigger.push_back( 4 );
			      if( pathNamesLast[h] == "HLT_PFHT900_v1") vfiredTrigger.push_back( 5 );
			      if( pathNamesLast[h] == "HLT_Ele95_CaloIdVT_GsfTrkIdT_v1") vfiredTrigger.push_back( 6 );
			      if( pathNamesLast[h] == "HLT_Mu40_v1") vfiredTrigger.push_back( 7 );
			      // else 																								vfiredTrigger   	    .push_back( -99 );
			}
			
		}
		
		nBranches_->triggerObject_filterIDs.push_back(vfilterIDs);
		nBranches_->triggerObject_firedTrigger.push_back(vfiredTrigger);
		
	}
	
}
// ###### Available triggers #######
// Trigger generation_step
// Trigger digitisation_step
// Trigger L1simulation_step
// Trigger digi2raw_step
// Trigger HLTriggerFirstPath
// Trigger HLT_Mu40_v1
// Trigger HLT_IsoMu20_eta2p1_IterTrk02_v1
// Trigger HLT_IsoTkMu20_eta2p1_IterTrk02_v1
// Trigger HLT_IsoMu24_eta2p1_IterTrk02_v1
// Trigger HLT_IsoTkMu24_eta2p1_IterTrk02_v1
// Trigger HLT_IsoMu24_IterTrk02_v1
// Trigger HLT_IsoTkMu24_IterTrk02_v1
// Trigger HLT_Mu17_Mu8_v1
// Trigger HLT_Mu17_TkMu8_v1
// Trigger HLT_Mu30_TkMu11_v1
// Trigger HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v1
// Trigger HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v1
// Trigger HLT_DoubleMu4_3_Bs_v1
// Trigger HLT_DoubleMu4_3_Jpsi_Displaced_v1
// Trigger HLT_Dimuon20_Jpsi_v1
// Trigger HLT_Dimuon13_PsiPrime_v1
// Trigger HLT_Dimuon13_Upsilon_v1
// Trigger HLT_Mu25_TkMu0_dEta18_Onia_v1
// Trigger HLT_DoubleMu4_JpsiTrk_Displaced_v1
// Trigger HLT_DoubleMu4_PsiPrimeTrk_Displaced_v1
// Trigger HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v1
// Trigger HLT_DoubleMu33NoFiltersNoVtx_v1
// Trigger HLT_DoubleMu38NoFiltersNoVtx_v1
// Trigger HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL_v1
// Trigger HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL_v1
// Trigger HLT_Ele27_eta2p1_WP85_Gsf_v1
// Trigger HLT_Ele32_eta2p1_WP85_Gsf_v1
// Trigger HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v1
// Trigger HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon18_AND_HE10_R9Id65_Mass95_v1
// Trigger HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_v1
// Trigger HLT_PFJet260_v1
// Trigger HLT_AK8PFJet360TrimMod_Mass30_v1
// Trigger HLT_L2Mu10_NoVertex_NoBPTX_v1
// Trigger HLT_L2Mu10_NoVertex_NoBPTX3BX_NoHalo_v1
// Trigger HLT_L2Mu20_NoVertex_3Sta_NoBPTX3BX_NoHalo_v1
// Trigger HLT_L2Mu30_NoVertex_3Sta_NoBPTX3BX_NoHalo_v1
// Trigger HLT_JetE30_NoBPTX_v1
// Trigger HLT_JetE30_NoBPTX3BX_NoHalo_v1
// Trigger HLT_JetE50_NoBPTX3BX_NoHalo_v1
// Trigger HLT_JetE70_NoBPTX3BX_NoHalo_v1
// Trigger HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1
// Trigger HLT_Ele22_eta2p1_WP85_Gsf_LooseIsoPFTau20_v1
// Trigger HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg_v1
// Trigger HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120_v1
// Trigger HLT_ReducedIterativeTracking_v1
// Trigger HLT_IsoMu20_eta2p1_IterTrk02_TriCentralPFJet60_50_35_v1
// Trigger HLT_IsoMu20_eta2p1_IterTrk02_TriCentralPFJet40_v1
// Trigger HLT_IsoMu20_eta2p1_IterTrk02_CentralPFJet30_BTagCSV_v1
// Trigger HLT_IsoMu24_eta2p1_IterTrk02_TriCentralPFJet60_50_35_v1
// Trigger HLT_IsoMu24_eta2p1_IterTrk02_TriCentralPFJet40_v1
// Trigger HLT_IsoMu24_eta2p1_IterTrk02_CentralPFJet30_BTagCSV_v1
// Trigger HLT_Ele27_eta2p1_WP85_Gsf_TriCentralPFJet40_v1
// Trigger HLT_Ele27_eta2p1_WP85_Gsf_TriCentralPFJet60_50_35_v1
// Trigger HLT_Ele27_eta2p1_WP85_Gsf_CentralPFJet30_BTagCSV_v1
// Trigger HLT_Ele32_eta2p1_WP85_Gsf_TriCentralPFJet40_v1
// Trigger HLT_Ele32_eta2p1_WP85_Gsf_TriCentralPFJet60_50_35_v1
// Trigger HLT_Ele32_eta2p1_WP85_Gsf_CentralPFJet30_BTagCSV_v1
// Trigger HLT_Mu40_eta2p1_PFJet200_PFJet50_v1
// Trigger HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1
// Trigger HLT_Ele23_Ele12_CaloId_TrackId_Iso_v1
// Trigger HLT_Ele17_Ele12_Ele10_CaloId_TrackId_v1
// Trigger HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v1
// Trigger HLT_Mu8_TrkIsoVVL_Ele23_Gsf_CaloId_TrackId_Iso_MediumWP_v1
// Trigger HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40_v1
// Trigger HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40_v1
// Trigger HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40_v1
// Trigger HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40_v1
// Trigger HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40_v1
// Trigger HLT_Photon135_PFMET40_v1
// Trigger HLT_Photon150_PFMET40_v1
// Trigger HLT_Photon160_PFMET40_v1
// Trigger HLT_Photon250_NoHE_PFMET40_v1
// Trigger HLT_Photon300_NoHE_PFMET40_v1
// Trigger HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF_v1
// Trigger HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF_v1
// Trigger HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF_v1
// Trigger HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF_v1
// Trigger HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF_v1
// Trigger HLT_Photon135_VBF_v1
// Trigger HLT_Photon150_VBF_v1
// Trigger HLT_Photon160_VBF_v1
// Trigger HLT_Photon250_NoHE_VBF_v1
// Trigger HLT_Photon300_NoHE_VBF_v1
// Trigger HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v1
// Trigger HLT_Ele95_CaloIdVT_GsfTrkIdT_v1
// Trigger HLT_DoublePhoton85_v1
// Trigger HLT_Photon155_v1
// Trigger HLT_Ele20WP60_Ele8_Mass55_v1
// Trigger HLT_Ele25WP60_SC4_Mass55_v1
// Trigger HLT_L2DoubleMu23_NoVertex_v1
// Trigger HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10_v1
// Trigger HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10_v1
// Trigger HLT_PFMET170_NoiseCleaned_v1
// Trigger HLT_PFMET120_NoiseCleaned_BTagCSV07_v1
// Trigger HLT_PFHT350_PFMET120_NoiseCleaned_v1
// Trigger HLT_PFHT900_v1
// Trigger HLT_ZeroBias_v1
// Trigger HLT_Physics_v1
// Trigger HLTriggerFinalPath
//
