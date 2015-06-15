#ifndef ElectronsNtuplizer_H
#define ElectronsNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class ElectronsNtuplizer : public CandidateNtuplizer {

public:
  ElectronsNtuplizer( edm::EDGetTokenT<edm::View<pat::Electron> >          electronToken, 
                      edm::EDGetTokenT<reco::VertexCollection>             verticeToken , 
		      edm::EDGetTokenT<double>  		           rhoToken     ,
                      std::vector<edm::EDGetTokenT<edm::ValueMap<bool> > > eleIDtokens  ,
		      NtupleBranches*				           nBranches    );
  ~ElectronsNtuplizer( void );
  
  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
  
  float dEtaInSeed( const pat::Electron &ele );
  bool  eleIDpassed( std::string id, const pat::Electron &ele );
  
private:
   edm::EDGetTokenT<edm::View<pat::Electron> > electronToken_;
   edm::EDGetTokenT<reco::VertexCollection>  verticeToken_ ;
   edm::EDGetTokenT<double> 		     rhoToken_     ;   
   edm::EDGetTokenT<edm::ValueMap<bool> >    electronVetoIdMapToken_;  
   edm::EDGetTokenT<edm::ValueMap<bool> >    electronLooseIdMapToken_; 
   edm::EDGetTokenT<edm::ValueMap<bool> >    electronMediumIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> >    electronTightIdMapToken_;
   edm::EDGetTokenT<edm::ValueMap<bool> >    electronHEEPIdMapToken_;
   
   edm::Handle<edm::View<pat::Electron> >    electrons_    ;
   edm::Handle<reco::VertexCollection> 	     vertices_     ;
   edm::Handle<double> 			     rho_	   ;
   edm::Handle<edm::ValueMap<bool> >         heep_id_decisions;
   edm::Handle<edm::ValueMap<bool> >         veto_id_decisions;  
   edm::Handle<edm::ValueMap<bool> >         loose_id_decisions; 
   edm::Handle<edm::ValueMap<bool> >         medium_id_decisions;
   edm::Handle<edm::ValueMap<bool> >         tight_id_decisions;
      
};

#endif // ElectronsNtuplizer_H
