#ifndef PileUpNtuplizer_H
#define PileUpNtuplizer_H

#include "../interface/CandidateNtuplizer.h"

class PileUpNtuplizer : public CandidateNtuplizer {

public:
  PileUpNtuplizer( std::vector< edm::EDGetTokenT< std::vector<PileupSummaryInfo> > > tokens, NtupleBranches* nBranches, std::map< std::string, bool >&  runFlags );
  ~PileUpNtuplizer( void );
  
  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
  
private:
   edm::EDGetTokenT< std::vector<PileupSummaryInfo> > pileUpToken_; 
     
   edm::Handle< std::vector<PileupSummaryInfo> >  pileUpInfo_;
   bool isJpsiMu_;
   bool isJpsiEle_;


};

#endif // PileUpNtuplizer_H
