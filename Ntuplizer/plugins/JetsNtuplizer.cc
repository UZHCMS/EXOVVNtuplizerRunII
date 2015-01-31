#include "../interface/JetsNtuplizer.h"
#include "DataFormats/JetReco/interface/CATopJetTagInfo.h"

//===================================================================================================================        
JetsNtuplizer::JetsNtuplizer( std::vector<edm::InputTag> labels, NtupleBranches* nBranches )
   : CandidateNtuplizer( nBranches )
   , jetInputTag_( labels[0] )
   , fatjetInputTag_( labels[1] )

{

   
}

//===================================================================================================================
JetsNtuplizer::~JetsNtuplizer( void )
{
}

void JetsNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  event.getByLabel(jetInputTag_, jets_);
  event.getByLabel(fatjetInputTag_, fatjets_);
  
  
  int nsubjets = 0;
 
  nBranches_->njetsAK4 = 0;
  for (const pat::Jet &j : *jets_) {
        if (j.pt() < 20) continue;
    bool IDLoose = false;
      
    if( j.nConstituents() > 1 && 
        j.muonEnergyFraction() < 0.99 && 
        j.photonEnergyFraction() < 0.99 && 
        j.chargedEmEnergyFraction() < 0.99 &&
        j.neutralHadronEnergyFraction() < 0.99 && 
        j.chargedHadronEnergyFraction() > 0. ) IDLoose = true;
      
    if( !IDLoose ) continue;
    nBranches_->njetsAK4++;
                
    nBranches_->jetAK4_pt     	    .push_back(j.pt());
    nBranches_->jetAK4_eta    	    .push_back(j.eta());
    nBranches_->jetAK4_mass   	    .push_back(j.mass());
    nBranches_->jetAK4_phi    	    .push_back(j.phi());   
    nBranches_->jetAK4_e      	    .push_back(j.energy());
    nBranches_->jetAK4_IDLoose      .push_back(IDLoose);
    nBranches_->jetAK4_cm     	    .push_back(j.chargedMultiplicity());
    nBranches_->jetAK4_nm     	    .push_back(j.neutralMultiplicity());
    nBranches_->jetAK4_che     	    .push_back(j.chargedHadronEnergy()+j.electronEnergy()+j.muonEnergy());
    nBranches_->jetAK4_ne     	    .push_back(j.neutralHadronEnergy()+j.photonEnergy());
    nBranches_->jetAK4_charge 	    .push_back(j.charge());
    nBranches_->jetAK4_ssv    	    .push_back(j.bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));
	nBranches_->jetAK4_cisv    	    .push_back(std::max(0.f,j.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags")));
    nBranches_->jetAK4_tchp   	    .push_back(j.bDiscriminator("trackCountingHighPurBJetTags"));
    nBranches_->jetAK4_tche   	    .push_back(j.bDiscriminator("trackCountingHighEffBJetTags"));
    nBranches_->jetAK4_jp     	    .push_back(j.bDiscriminator("jetProbabilityBJetTags"));
    nBranches_->jetAK4_jbp    	    .push_back(j.bDiscriminator("jetBProbabilityBJetTags"));
	nBranches_->jetAK4_flavour		.push_back(abs(j.partonFlavour()));
	
    nBranches_->jetAK4_vtxMass		.push_back(j.userFloat("vtxMass")); 
	nBranches_->jetAK4_vtxNtracks	.push_back(j.userFloat("vtxNtracks")); 
	nBranches_->jetAK4_vtx3DVal		.push_back(j.userFloat("vtx3DVal")); 
	nBranches_->jetAK4_vtx3DSig		.push_back(j.userFloat("vtx3DSig")); 
  
  }		
  

  nBranches_->njetsAK8 = 0;
  //int nsubjets = 0;

  std::vector<float> vSubjetpt     ;
  std::vector<float> vSubjeteta    ;
  std::vector<float> vSubjetmass   ;
  std::vector<float> vSubjetphi    ;
  std::vector<float> vSubjete	   ;
  std::vector<int  > vSubjetcharge ;
  std::vector<int  > vSubjetflavour;
  std::vector<float> vSubjetssv    ;
  std::vector<float> vSubjetcsv    ;
  std::vector<float> vSubjettchp   ;
  std::vector<float> vSubjettche   ; 
  std::vector<float> vSubjetjp     ;
  std::vector<float> vSubjetjbp    ; 
  
  for (const pat::Jet &fj : *fatjets_) {
        if (fj.pt() < 80) continue;

    bool IDLoose = false;
      
    if( fj.nConstituents() > 1 && 
        fj.muonEnergyFraction() < 0.99 && 
        fj.photonEnergyFraction() < 0.99 && 
        fj.chargedEmEnergyFraction() < 0.99 &&
        fj.neutralHadronEnergyFraction() < 0.99 && 
        fj.chargedHadronEnergyFraction() > 0. ) IDLoose = true;
    
    if( !IDLoose ) continue;
    
    nBranches_->njetsAK8++;	       
                   
    nBranches_->jetAK8_pt     	    .push_back(fj.pt());
    nBranches_->jetAK8_eta    	    .push_back(fj.eta());
    nBranches_->jetAK8_mass   	    .push_back(fj.mass());
    nBranches_->jetAK8_phi    	    .push_back(fj.phi());
    nBranches_->jetAK8_e      	    .push_back(fj.energy());
    nBranches_->jetAK8_cm     	    .push_back(fj.chargedMultiplicity());
    nBranches_->jetAK8_nm     	    .push_back(fj.neutralMultiplicity());     				       
    nBranches_->jetAK8_che     	    .push_back(fj.chargedHadronEnergy()+fj.electronEnergy()+fj.muonEnergy());
    nBranches_->jetAK8_ne     	    .push_back(fj.neutralHadronEnergy()+fj.photonEnergy());
    nBranches_->jetAK8_charge 	    .push_back(fj.charge());					       
    nBranches_->jetAK8_ssv    	    .push_back(fj.bDiscriminator("simpleSecondaryVertexHighPurBJetTags")); 
    nBranches_->jetAK8_csv    	    .push_back(fj.bDiscriminator("combinedInclusiveSecondaryVertexV2BJetTags"));      
    nBranches_->jetAK8_tchp   	    .push_back(fj.bDiscriminator("trackCountingHighPurBJetTags"));	       
    nBranches_->jetAK8_tche   	    .push_back(fj.bDiscriminator("trackCountingHighEffBJetTags"));	       
    nBranches_->jetAK8_jp     	    .push_back(fj.bDiscriminator("jetProbabilityBJetTags"));	       
    nBranches_->jetAK8_jbp    	    .push_back(fj.bDiscriminator("jetBProbabilityBJetTags"));
  
    nBranches_->jetAK8_flavour		.push_back(abs(fj.partonFlavour()));
    nBranches_->jetAK8_tau1			.push_back(fj.userFloat("NjettinessAK8:tau1"));	       
    nBranches_->jetAK8_tau2			.push_back(fj.userFloat("NjettinessAK8:tau2"));
    nBranches_->jetAK8_tau3			.push_back(fj.userFloat("NjettinessAK8:tau3"));	
	
	nBranches_->jetAK8pruned_mass	.push_back(fj.userFloat("ak8PFJetsCHSPrunedLinks"));
    nBranches_->jetAK8trimmed_mass	.push_back(fj.userFloat("ak8PFJetsCHSTrimmedLinks"));
	nBranches_->jetAK8filtered_mass	.push_back(fj.userFloat("ak8PFJetsCHSFilteredLinks"));
	
	
	// edm::Handle<edm::View<reco::CATopJetTagInfo>> tagInfo;
// 	event.getByLabel("tagInfo", tagInfo);

	//     const reco::SecondaryVertexTagInfo * svTagInfos = fj.tagInfoSecondaryVertex("secondaryVertex");
	//     nBranches_->jetCA8_nSVs   .push_back(svTagInfos->nVertices());
	//
	//  reco::CATopTag const * tagInfo =  dynamic_cast<reco::CATopTagInfo const *>( fj.tagInfo("caTop"));
	// 	int nSubJets;
	// if ( tagInfo != 0 ) {
	//    	// double minMass = tagInfo->properties().minMass;
	// 	//double topMass = tagInfo->properties().topMass;
	//     nSubJets = tagInfo->properties().nSubJets;
	// }
	//
	// nBranches_->jetAK8_nSubJets		.push_back(nSubJets);
	
    nsubjets = 0;
    for( unsigned int sj=0; sj<fj.numberOfDaughters(); ++sj ){

      const pat::Jet* subjet = dynamic_cast<const pat::Jet*>(fj.daughter(sj));
      if( subjet->pt() < 0.01 ) continue;
      nsubjets++;
  }
	
}
}



