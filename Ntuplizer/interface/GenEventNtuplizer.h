#ifndef GenEventNtuplizer_H
#define GenEventNtuplizer_H

#include <algorithm>
#include <numeric>
#include <functional>
#include <iostream>

#include "../interface/CandidateNtuplizer.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"

class GenEventNtuplizer : public CandidateNtuplizer {

public:
  GenEventNtuplizer( std::vector< edm::EDGetTokenT< GenEventInfoProduct > > tokens, NtupleBranches* nBranches ,std::vector< edm::EDGetTokenT< LHEEventProduct > > tokens_lhe,std::map< std::string, bool >&  runFlags );
  ~GenEventNtuplizer( void );
  
  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
  
private:
   edm::EDGetTokenT< GenEventInfoProduct > geneventToken_; 
     
   edm::Handle< GenEventInfoProduct >  geneventInfo_;
  
   edm::EDGetTokenT<LHEEventProduct > lheEventProductToken_; 
   edm::Handle<LHEEventProduct> lheEventProduct_;
   bool isJpsiMu_;
   bool isJpsiEle_;
   bool isJpsiTau_;

};

#endif // GenEventNtuplizer_H
