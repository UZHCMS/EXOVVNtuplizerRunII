#ifndef GenJetsNtuplizer_H
#define GenJetsNtuplizer_H

#include "../interface/CandidateNtuplizer.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

class GenJetsNtuplizer : public CandidateNtuplizer {

public:

 GenJetsNtuplizer( edm::EDGetTokenT<reco::GenJetCollection> token, edm::EDGetTokenT<pat::JetCollection> AK8token, NtupleBranches* nBranches );
   ~GenJetsNtuplizer( void );

  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
 
private:
  
  edm::EDGetTokenT<reco::GenJetCollection> 					genJetInputToken_       ;
  edm::EDGetTokenT<pat::JetCollection> 				    	genJetAK8InputToken_    ;

  edm::Handle<reco::GenJetCollection>      					genJets_         				;
  edm::Handle<pat::JetCollection>      		    			genJetsAK8_      				;

};

#endif // GenJetsNtuplizer_H
