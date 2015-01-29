#include "../interface/JetsNtuplizer.h"

//===================================================================================================================        
JetsNtuplizer::JetsNtuplizer( std::vector<edm::InputTag> labels, std::vector<std::string> jecCA8Labels, std::vector<std::string> jecAK5Labels, NtupleBranches* nBranches )
   : CandidateNtuplizer( nBranches )
   , jetsAK5Label_( labels[0] )
   , jetsCA8Label_( labels[1] )
   , jetsCA8prunedLabel_( labels[2] )
   , verticesLabel_( labels[3] )
   , rhoLabel_( labels[4] )
   , flavLabel_( labels[5] )
{
   jecCA8PayloadNames_ = jecCA8Labels;
   jecCA8PayloadNames_.pop_back();
   jecCA8UncName_ = jecCA8Labels.back();

   jecAK5PayloadNames_ = jecAK5Labels;
   jecAK5PayloadNames_.pop_back();
   jecAK5UncName_ = jecAK5Labels.back();
      
   initJetCorrFactors();
   
   
}

//===================================================================================================================
JetsNtuplizer::~JetsNtuplizer( void )
{
}

//===================================================================================================================


//===================================================================================================================
void JetsNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

event.getByLabel(jetLabel_      , jets_      );

nBranches_->njetsAK5 = 0;
   for( unsigned j=0; j<jets_->size(); ++j ){
     
     reco::Candidate::LorentzVector uncorrJet = (*jetsAK5_)[j].correctedP4(0);

     jecAK5_->setJetEta( uncorrJet.eta()          );
     
     nBranches_->njetsAK5++;
     
//   event.getByLabel(jetsAK5Label_      , jetsAK5_      );
//  
//   if ( flavLabel_.label() != "" )
//      event.getByLabel(flavLabel_         , theTagByValue );
//   	    
//  /*here we want to save the jets info*/   
//   nBranches_->njetsAK5 = 0;
//   for( unsigned j=0; j<jetsAK5_->size(); ++j ){
//       
//     bool IDLoose = false;
//       
//     if( (*jetsAK5_)[j].nConstituents() > 1 && 
//         (*jetsAK5_)[j].muonEnergyFraction() < 0.99 && 
//         (*jetsAK5_)[j].photonEnergyFraction() < 0.99 && 
//         (*jetsAK5_)[j].chargedEmEnergyFraction() < 0.99 &&
//         (*jetsAK5_)[j].neutralHadronEnergyFraction() < 0.99 && 
//         (*jetsAK5_)[j].chargedHadronEnergyFraction() > 0. ) IDLoose = true;
//       
//     //if( !IDLoose ) continue;
//       
//     reco::Candidate::LorentzVector uncorrJet = (*jetsAK5_)[j].correctedP4(0);
// 
//     jecAK5_->setJetEta( uncorrJet.eta()          );
//     jecAK5_->setJetPt ( uncorrJet.pt()           );
//     jecAK5_->setJetE  ( uncorrJet.energy()       );
//     jecAK5_->setJetA  ( (*jetsAK5_)[j].jetArea() );
//     jecAK5_->setRho   ( *(rho_.product())        );
//     jecAK5_->setNPV   ( vertices_->size()        );
//     double corr = jecAK5_->getCorrection();
//     
//     jecAK5Unc_->setJetEta( uncorrJet.eta() );
//     jecAK5Unc_->setJetPt( corr * uncorrJet.pt() );
//     double corrUp = corr * (1 + fabs(jecAK5Unc_->getUncertainty(1)));
//     jecAK5Unc_->setJetEta( uncorrJet.eta() );
//     jecAK5Unc_->setJetPt( corr * uncorrJet.pt() );
//     double corrDown = corr * ( 1 - fabs(jecAK5Unc_->getUncertainty(-1)) );
//                 
//     //if( corr*uncorrJet.pt() < 20. ) continue;
// 
//     nBranches_->njetsAK5++;
//                 
//     nBranches_->jetAK5_pt     	    .push_back(corr*uncorrJet.pt());
// 
// 
//     const reco::GenParticle * genP = (*jetsAK5_)[j].genParton();
//     if( genP ) nBranches_->jetAK5_flavour.push_back(abs(genP->pdgId()));
//     else nBranches_->jetAK5_flavour.push_back(abs((*jetsAK5_)[j].partonFlavour()));
//                
//   }
  
		
}
