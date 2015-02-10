#include "../interface/TriggersNtuplizer.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

//===================================================================================================================        
TriggersNtuplizer::TriggersNtuplizer( std::vector<edm::EDGetTokenT<edm::TriggerResults>> tokens, NtupleBranches* nBranches )
   : CandidateNtuplizer( nBranches )
   , HLTtriggersToken_( tokens[0] )
{
   
}

//===================================================================================================================
TriggersNtuplizer::~TriggersNtuplizer( void )
{

}

//===================================================================================================================
void TriggersNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
	
	event.getByToken(HLTtriggersToken_, HLTtriggers_);
	const edm::TriggerNames& trigNames = event.triggerNames(*HLTtriggers_);
	
// for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {
//         std::cout << "Trigger " << trigNames.triggerName(i) <<
//                 ": " << (HLTtriggers_->accept(i) ? "PASS" : "fail (or not run)")
//                 << std::endl;
//     }
	

	

	unsigned int TrggIndex_HLT_AK8PFJet360TrimMod_Mass30_v1					( trigNames.triggerIndex("HLT_AK8PFJet360TrimMod_Mass30_v1"));
	unsigned int TrggIndex_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v1			( trigNames.triggerIndex("HT_AK8PFHT700_TrimR0p1PT0p03Mass50_v1"));
	unsigned int TrggIndex_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p3_v1	( trigNames.triggerIndex("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p3_v1"));
	unsigned int TrggIndex_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v1			( trigNames.triggerIndex("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v1"));
	unsigned int TrggIndex_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v1			( trigNames.triggerIndex("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v1"));
	unsigned int TrggIndex_HLT_PFHT900_v1									( trigNames.triggerIndex("HLT_PFHT900_v1"));

	nBranches_->isFired_HLT_AK8PFJet360TrimMod_Mass30_v1               = false;
	nBranches_->isFired_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v1         = false;
	nBranches_->isFired_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p3_v1 = false;
	nBranches_->isFired_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v1          = false;
	nBranches_->isFired_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v1          = false;
	nBranches_->isFired_HLT_PFHT900_v1                                 = false;
	
	if( TrggIndex_HLT_AK8PFJet360TrimMod_Mass30_v1 					< HLTtriggers_->size() ) nBranches_->isFired_HLT_AK8PFJet360TrimMod_Mass30_v1 				= HLTtriggers_->accept(TrggIndex_HLT_AK8PFJet360TrimMod_Mass30_v1				);
	if( TrggIndex_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v1 			< HLTtriggers_->size() ) nBranches_->isFired_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v1 		= HLTtriggers_->accept(TrggIndex_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v1			);	
	if( TrggIndex_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p3_v1 	< HLTtriggers_->size() ) nBranches_->isFired_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p3_v1 = HLTtriggers_->accept(TrggIndex_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p3_v1	);
	if( TrggIndex_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v1 			< HLTtriggers_->size() ) nBranches_->isFired_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v1 			= HLTtriggers_->accept(TrggIndex_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v1 			);
	if( TrggIndex_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v1 			< HLTtriggers_->size() ) nBranches_->isFired_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v1 			= HLTtriggers_->accept(TrggIndex_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v1			); 
	if( TrggIndex_HLT_PFHT900_v1 									< HLTtriggers_->size() ) nBranches_->isFired_HLT_PFHT900_v1 								= HLTtriggers_->accept(TrggIndex_HLT_PFHT900_v1 								);

	
}
