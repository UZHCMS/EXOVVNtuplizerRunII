#ifndef JpsiTauNtuplizer_H
#define JpsiTauNtuplizer_H

#include "../interface/CandidateNtuplizer.h"
#include "../interface/MyStruct.h"
//#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
//#include "SimDataFormats/JetMatching/interface/JetFlavour.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "JetMETCorrections/Modules/interface/JetResolution.h"
//#include <CondFormats/JetMETObjects/interface/JetResolutionObject.h>

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

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
//#include "DataFormats/Common/interface/TriggerResults.h"
//#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
//#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"             
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"


// kinematic fit package
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"


#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"

#include <tuple>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <algorithm> 
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>

//struct particle_cand 
//{
//
//  // Longitudinal impact parameter and its significance
//  Float_t lip;
//  Float_t lips;
//  
//  // Impact parameter for the PV and its significance
//  Float_t pvip;
//  Float_t pvips;
//
//  // Flight length and its significance
//  Float_t fl3d;
//  Float_t fls3d;
//  
//  // opening angle
//  Float_t alpha;
//
//};



class JpsiTauNtuplizer : public CandidateNtuplizer {


 public:
  JpsiTauNtuplizer( edm::EDGetTokenT<pat::MuonCollection>    muonToken   , 
		   edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
		   edm::EDGetTokenT<reco::BeamSpot> bsToken, 
		   edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken,
		   edm::EDGetTokenT<pat::PackedCandidateCollection> losttrackToken,
		   edm::EDGetTokenT<edm::TriggerResults> triggertoken,
		   edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobject,
		   edm::EDGetTokenT<reco::GenParticleCollection> genptoken, 
		    edm::EDGetTokenT<std::vector<reco::GenJet>> genttoken,
		   std::map< std::string, bool >& runFlags,
		   NtupleBranches* nBranches );

  ~JpsiTauNtuplizer( void );

  std::tuple<Float_t, TransientVertex> vertexProb( const std::vector<reco::TransientTrack>& tracks);

  particle_cand calculateIPvariables(AnalyticalImpactPointExtrapolator extrapolator,
				     RefCountedKinematicParticle particle,
				     RefCountedKinematicVertex vertex,
				     reco::Vertex wrtVertex);
  
  std::pair<bool, Measurement1D> absoluteImpactParameter(const TrajectoryStateOnSurface& tsos,
							 RefCountedKinematicVertex vertex,
							 VertexDistance& distanceComputer);


  math::PtEtaPhiMLorentzVector daughter_p4(std::vector< RefCountedKinematicParticle > fitted_children, size_t i);

  float MuonPFIso(pat::Muon muon);
  void fillBranches( edm::Event const & event, const edm::EventSetup& iSetup );
  Float_t getMaxDoca(std::vector<RefCountedKinematicParticle> &kinParticles);
  Float_t getMinDoca(std::vector<RefCountedKinematicParticle> &kinParticles);
  TVector3 getVertex(const reco::GenParticle& part);

  void printout(const RefCountedKinematicVertex& myVertex);
  void printout(const RefCountedKinematicParticle& myParticle);
  void printout(const RefCountedKinematicTree& myTree);

  Int_t decaymode_id(std::string str);
  
private:
   edm::EDGetTokenT<pat::MuonCollection>    muonToken_   ;
   edm::EDGetTokenT<reco::VertexCollection> verticeToken_   ;
   edm::EDGetTokenT<reco::BeamSpot> bsToken_   ;
   edm::EDGetTokenT<pat::PackedCandidateCollection>   		packedpfcandidatesToken_;
   edm::EDGetTokenT<pat::PackedCandidateCollection>   		losttrackToken_;
   edm::EDGetTokenT<edm::TriggerResults> 		     HLTtriggersToken_;
   edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection>  triggerObjects_;
   edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
   edm::EDGetTokenT<std::vector<reco::GenJet>> genTauToken_;

   edm::Handle<pat::MuonCollection>      		       muons_		       ;
   edm::Handle< reco::VertexCollection >  vertices_;
   edm::Handle< reco::BeamSpot >  beamspot_;
   edm::Handle< std::vector<pat::PackedCandidate> > packedpfcandidates_   ;
   edm::Handle< std::vector<pat::PackedCandidate> > losttrack_   ;
   edm::Handle< edm::TriggerResults> 			     HLTtriggers_;
   edm::Handle<pat::TriggerObjectStandAloneCollection>	     triggerObjects;
   edm::Handle< reco::GenParticleCollection >  genParticles_;
   edm::Handle< std::vector<reco::GenJet> >  genTaus_;

   edm::ESHandle<TransientTrackBuilder> builder;

   const MagneticField                 *fMagneticField;

   ParticleMass muon_mass = 0.1056583;
   ParticleMass jpsi_mass = 3.09687;
   ParticleMass pion_mass = 0.139571;
   ParticleMass kaon_mass = 0.493677;
   
   float muon_sigma = 0.0000001;
   float jp_m_sigma = 0.00004;
   float pion_sigma = 0.000016;
   float kaon_sigma = 0.000016;
   float chi = 0.;
   float ndf = 0.;

   bool runOnMC_;   
   bool doCutFlow_;
   
};




#endif // JpsiTauNtuplizer_H

