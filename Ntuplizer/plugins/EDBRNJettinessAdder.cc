#include "ExoDiBosonResonances/EDBRJets/interface/EDBRNJettinessAdder.h"
//#include "fastjet/contrib/Njettiness.hh"

#include "FWCore/Framework/interface/MakerMacros.h"

void EDBRNJettinessAdder::produce(edm::Event & iEvent, const edm::EventSetup & iSetup) {
  // read input collection                                                                                                                                                                                                            
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByToken(src_token_, jets);

  // prepare room for output                                                                                                                                                                                                          
  std::auto_ptr<std::vector<pat::Jet> >  out(new std::vector<pat::Jet>());
  out->reserve(jets->size());

  for ( typename edm::View<pat::Jet>::const_iterator jetIt = jets->begin() ; jetIt != jets->end() ; ++jetIt ) {
    //pat::Jet newCand(*jetIt);                                                                                                                                                                                                       
    edm::Ptr<pat::Jet> jetPtr = jets->ptrAt(jetIt - jets->begin());
    pat::Jet* newJet = jetPtr->clone();

    float t1=getTau(1, jetPtr );
    float t2=getTau(2, jetPtr );
    float t3=getTau(3, jetPtr );

    newJet->addUserFloat("tau1",t1);
    newJet->addUserFloat("tau2",t2);
    newJet->addUserFloat("tau3",t3);
    newJet->addUserFloat("tau21",t2/t1);
    //std::cout << "Tau21 = " << t2/t1 << std::endl;                                                                                                                                                                                  
    //std::cout << "Pruned mass = " << jetPtr->userFloat("ak8PFJetsCHSPrunedLinks") << std::endl;                                                                                                                                     

    out->push_back(*newJet);
  }

  iEvent.put(out);
}

float EDBRNJettinessAdder::getTau(int num, edm::Ptr<pat::Jet> object) const
{
  /*std::vector<reco::CandidatePtr> all_particles;
  for (unsigned k =0; k < object->numberOfSourceCandidatePtrs(); k++) {
    all_particles.push_back( object->sourceCandidatePtr(k) );
  }

  std::vector<fastjet::PseudoJet> FJparticles;
  for (unsigned particle = 0; particle < all_particles.size(); particle++){
    const reco::CandidatePtr thisParticle = all_particles.at(particle);
    FJparticles.push_back( fastjet::PseudoJet( thisParticle->px(), thisParticle->py(), thisParticle->pz(), thisParticle->energy() ) );
  }
  fastjet::contrib::NsubParameters paraNsub = fastjet::contrib::NsubParameters(1.0, cone_); //assume R=0.7 jet clusering used                                                                                                         
  fastjet::contrib::Njettiness routine(fastjet::contrib::Njettiness::onepass_kt_axes, paraNsub);
  return routine.getTau(num, FJparticles);*/
  return 1.;
  
}



DEFINE_FWK_MODULE(EDBRNJettinessAdder);
