#ifndef JpsiMuNtuplizer_H
#define JpsiMuNtuplizer_H

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
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Common/interface/Handle.h"  // edm::Handle
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"   
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxData.h" 
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonQuality.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"             
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"



#include <tuple>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <algorithm> 
#include <stdio.h>
#include <stdlib.h>
#include <map>
//IN CASE THERE'S ONE LIBRARY MISSING OR SOMETHING


class JpsiMuNtuplizer : public CandidateNtuplizer {
// ----------------------------------------------------------------------

struct track_entry_t {
  track_entry_t(int ix, int pid, bool mfit) : trackIx(ix), particleID(pid), massFit(mfit) {}
  int trackIx;
  int particleID;
  bool massFit;
};


 public:
 typedef std::set<track_entry_t>::iterator HFDecayTreeTrackIterator;
 JpsiMuNtuplizer( edm::EDGetTokenT<pat::MuonCollection>    muonToken   , 
		  edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
		  edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken,
		  NtupleBranches* nBranches );
  ~JpsiMuNtuplizer( void );
  std::tuple<Float_t, Float_t> vertexProb(const std::vector<reco::TransientTrack> &tracks);
  void getAllTracks(std::vector<track_entry_t> *out_vector, int onlyThisVertex);


  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
  
private:
   edm::EDGetTokenT<pat::MuonCollection>    muonToken_   ;
   edm::EDGetTokenT<reco::VertexCollection> verticeToken_   ;

   //   edm::Handle<pat::JetCollection>      		       jets_		       ;

   edm::Handle< reco::VertexCollection >  vertices_;
   edm::Handle<pat::MuonCollection>      		       muons_		       ;
   edm::EDGetTokenT<pat::PackedCandidateCollection>   		packedpfcandidatesToken_;

   edm::Handle< std::vector<pat::PackedCandidate> > packedpfcandidates_   ;

   edm::ESHandle<TransientTrackBuilder> builder;
   AdaptiveVertexFitter avf;

   //BS part
   std::set<track_entry_t> fTrackIndices; // added tracks
   std::vector<HFDecayTreeTrackIterator> fSubVertices;
};

#endif // JpsiMuNtuplizer_H
