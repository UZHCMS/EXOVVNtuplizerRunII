#ifndef TausNtuplizer_H
#define TausNtuplizer_H

#include "../interface/CandidateNtuplizer.h"
#include "RecoTauTag/RecoTau/interface/PFRecoTauClusterVariables.h"

class TausNtuplizer : public CandidateNtuplizer {

public:
  TausNtuplizer( edm::EDGetTokenT<pat::TauCollection> tauToken, edm::EDGetTokenT<pat::TauCollection> tauBoostedTauToken , edm::EDGetTokenT<double> rhoToken, edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken, edm::EDGetTokenT<reco::VertexCollection> verticeToken,  NtupleBranches* nBranches, std::map< std::string, bool >& runFlags );
   ~TausNtuplizer( void );
   
   void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup /* , std::map< std::string, bool >& runFlags  */  );
   
 private:
   
   
   edm::EDGetTokenT<pat::TauCollection> 			tauInputToken_     		;
   edm::EDGetTokenT<pat::TauCollection> 			tauBoostedTauInputToken_  	;
   edm::EDGetTokenT<double> 		   	      		rhoToken_               	;
   edm::EDGetTokenT<pat::PackedCandidateCollection>   		packedpfcandidatesToken_;
   edm::EDGetTokenT<reco::VertexCollection>  			verticeToken_   	 	;

   
  
   edm::Handle< std::vector<pat::Tau> > taus_   ;
   edm::Handle< std::vector<pat::Tau> > boostedTaus_;
   edm::Handle< double >                rho_     ;
   edm::Handle< std::vector<pat::PackedCandidate> > packedpfcandidates_   ;
   edm::Handle<reco::VertexCollection> 	     			vertices_     			;

   
   bool doBoostedTaus_;
   TauIdMVAAuxiliaries clusterVariables_;
};

#endif // TausNtuplizer_H
