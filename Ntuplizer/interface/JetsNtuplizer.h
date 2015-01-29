#ifndef JetsNtuplizer_H
#define JetsNtuplizer_H

#include "../interface/CandidateNtuplizer.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"

class JetsNtuplizer : public CandidateNtuplizer {

public:
//   JetsNtuplizer( std::vector<edm::InputTag> labels, std::vector<std::string> jecCA8Labels, std::vector<std::string> jecAK5Labels, NtupleBranches* nBranches );
//   ~JetsNtuplizer( void );
  
  JetsNtuplizer(edm::InputTag> labels, std::vector<std::string> jecCA8Labels, NtupleBranches* nBranches );
  ~JetsNtuplizer( void );  

  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
//  void initJetCorrFactors( void );

 
private:
  edm::InputTag jetLabel ;
  edm::Handle< std::vector<pat::Jet> >      jets_      ;
//   edm::InputTag jetsAK5Label_      ;
//   edm::InputTag jetsCA8Label_      ;
//   edm::InputTag jetsCA8prunedLabel_;
//   edm::InputTag verticesLabel_     ;
//   edm::InputTag rhoLabel_          ;
//   edm::InputTag flavLabel_         ;
//   
//   edm::Handle< std::vector<pat::Jet> >      jetsAK5_      ;
//   edm::Handle< std::vector<pat::Jet> >      jetsCA8_      ;
//   edm::Handle< std::vector<pat::Jet> >      jetsCA8pruned_;        
//   edm::Handle< std::vector<reco::Vertex> >  vertices_     ; 
//   edm::Handle< double >                     rho_          ;
//   edm::Handle<reco::JetFlavourMatchingCollection> theTagByValue;
//   
//   std::vector<std::string>                    jecCA8PayloadNames_;
//   std::string                                 jecCA8UncName_     ;  
//   boost::shared_ptr<JetCorrectionUncertainty> jecCA8Unc_         ;
//   boost::shared_ptr<FactorizedJetCorrector>   jecCA8_            ;      
// 
//   std::vector<std::string>                    jecAK5PayloadNames_;
//   std::string                                 jecAK5UncName_     ;  
//   boost::shared_ptr<JetCorrectionUncertainty> jecAK5Unc_         ;
//   boost::shared_ptr<FactorizedJetCorrector>   jecAK5_            ;
  
};

#endif // JetsNtuplizer_H
