#include "../interface/MuonsNtuplizer.h"

#include <cmath>

//===================================================================================================================        
MuonsNtuplizer::MuonsNtuplizer( edm::EDGetTokenT<pat::MuonCollection> muonToken, edm::EDGetTokenT<reco::VertexCollection> verticeToken, NtupleBranches* nBranches )
   : CandidateNtuplizer( nBranches )
   , muonsToken_( muonToken )
   , verticesToken_( verticeToken )
	   
{
}

//===================================================================================================================
MuonsNtuplizer::~MuonsNtuplizer( void )
{

}
  
//===================================================================================================================
void MuonsNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

    event.getByToken(muonsToken_, muons_); 
    event.getByToken(verticesToken_, vertices_);
 

    int nmus = 0;
    
	for (const pat::Muon &mu : *muons_) {
        bool isHighPt = mu.isHighPtMuon(vertices_->at(0));
      
        if( !isHighPt ) continue;
      
      nmus++;
                  
      nBranches_->lep_type   .push_back(mu.pdgId() );   	 
      nBranches_->lep_charge .push_back(mu.charge());   	 
      nBranches_->lep_e      .push_back(mu.energy());   	 
      nBranches_->lep_eta    .push_back(mu.eta()   );   	 
      nBranches_->lep_mass   .push_back(mu.mass()  );   	 
      nBranches_->lep_pt     .push_back(mu.pt()    );   	 
      nBranches_->lep_phi    .push_back(mu.phi()   );
	  nBranches_->lep_isHighPtMuon.push_back(isHighPt);
      nBranches_->lep_isHEEP      .push_back(-99);
	  
      /*===== ISO ====*/
      
      nBranches_->lep_pfDeltaCorrRelIso     .push_back((mu.chargedHadronIso() + std::max(0., mu.neutralHadronIso() + mu.photonIso() - 0.5*mu.puChargedHadronIso()))/mu.pt());
      nBranches_->lep_pfRelIso              .push_back((mu.chargedHadronIso() + mu.neutralHadronIso()+ mu.photonIso())/mu.pt()) ; 
      nBranches_->lep_photonIso             .push_back(mu.photonIso());
      nBranches_->lep_neutralHadIso         .push_back(mu.neutralHadronIso());
      nBranches_->lep_chargedHadIso         .push_back(mu.chargedHadronIso());
      nBranches_->lep_trackIso              .push_back(mu.trackIso());
    } 
    
    nBranches_->nlep +=  nmus;
    
}