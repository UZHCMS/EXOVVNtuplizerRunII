#ifndef JpsiEleNtuplizer_H
#define JpsiEleNtuplizer_H

#include "../interface/CandidateNtuplizer.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include <CondFormats/JetMETObjects/interface/JetResolutionObject.h>


// for vertexing
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"  // edm::Handle
#include "TrackingTools/IPTools/interface/IPTools.h"
     
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"                                                                                                       
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxData.h" 
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonQuality.h"
             
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"                                                                                                 
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"



#include "TrackingTools/TransientTrack/interface/TransientTrack.h"             
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/IPTools/interface/IPTools.h"


#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include <tuple>


class JpsiEleNtuplizer : public CandidateNtuplizer {

public:
  JpsiEleNtuplizer( edm::EDGetTokenT<pat::Electron>    electronToken   , edm::EDGetTokenT<reco::VertexCollection> verticeToken, NtupleBranches* nBranches );
  ~JpsiEleNtuplizer( void );
  std::tuple<Float_t, Float_t> vertexProb(const std::vector<reco::TransientTrack> &tracks);

  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
  
private:
   edm::EDGetTokenT<pat::Electron> electronToken_;
   edm::EDGetTokenT<reco::VertexCollection> verticeToken_   ;

   //   edm::Handle<pat::JetCollection>      		       jets_		       ;

   edm::Handle< reco::VertexCollection >  vertices_;
   edm::Handle<pat::Electron>     electrons_    ;
   edm::ESHandle<TransientTrackBuilder> builder;

};

#endif // JpsiEleNtuplizer_H
