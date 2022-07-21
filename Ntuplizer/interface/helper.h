#ifndef HELPER_H
#define HELPER_H

#include "../interface/CandidateNtuplizer.h"
#include "../interface/MyStruct.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
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

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"


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
#include "TMatrixDSym.h"
#include "TVectorD.h"


#include "Hammer/Hammer.hh"
#include "Hammer/Process.hh"
#include "Hammer/Particle.hh"


//#include "../interface/MyStruct.h"
//#include <iostream>
//#include <string.h>
//#include <TVector3.h>
//#include <TRandom.h>

using namespace std;

// mass


class helper{
  
 public:


  ParticleMass muon_mass = 0.1056583;
  ParticleMass jpsi_mass = 3.09687;
  ParticleMass pion_mass = 0.139571;
  ParticleMass kaon_mass = 0.493677;
  ParticleMass ds_mass = 2.01026;
  ParticleMass d0_mass = 1.86483;
  
  float muon_sigma = 0.0000001;
  float jp_m_sigma = 0.00004;
  float pion_sigma = 0.000016;
  float kaon_sigma = 0.000016;
  float phi_sigma = 0.000016;
  float ds_sigma = 0.00005;
  float d0_sigma = 0.00005;
  
  Float_t mass_kaon = 0.493677;
  Float_t mass_pion = 0.139571;
  Float_t mass_D0 = 1.86483;
  Float_t mass_Dstar = 2.01026;
  Float_t mass_B0 = 5.27963;
  Float_t mass_Bc = 6.2749;
  


  // return vertex x, y, and z
  TVector3 getVertex(const reco::GenParticle& part);

  // return tau decay mode
  Int_t decaymode_id(std::string str);

  // return muon isolation
  float MuonPFIso(pat::Muon muon);

  // max. or min. of distance of closest approach
  Float_t getMaxDoca(std::vector<RefCountedKinematicParticle> &kinParticles);
  Float_t getMinDoca(std::vector<RefCountedKinematicParticle> &kinParticles);
  
  // nagivator of the kinematic fit class
  void printout(const RefCountedKinematicVertex& myVertex);
  void printout(const RefCountedKinematicParticle& myParticle);
  void printout(const RefCountedKinematicTree& myTree);
  

  // calculate IP related variables
  particle_cand calculateIPvariables(AnalyticalImpactPointExtrapolator extrapolator,
				     RefCountedKinematicParticle particle,
				     RefCountedKinematicVertex vertex,
				     reco::Vertex wrtVertex);

  std::pair<float, float> calculateIPvariables(RefCountedKinematicVertex tauVertex,
					       RefCountedKinematicVertex jpsiVertex,
					       reco::Vertex refitVertex);

  // absolute impact parameter
  std::pair<bool, Measurement1D> absoluteImpactParameter(const TrajectoryStateOnSurface& tsos,
							 RefCountedKinematicVertex vertex,
							 VertexDistance& distanceComputer);

  std::pair<bool, Measurement1D> absoluteImpactParameter3D(const TrajectoryStateOnSurface& tsos,
							   RefCountedKinematicVertex vertex);


  std::pair<bool, Measurement1D> absoluteTransverseImpactParameter(const TrajectoryStateOnSurface& tsos,
								   RefCountedKinematicVertex vertex);
  

  std::pair<bool, Measurement1D> signedTransverseImpactParameter(const TrajectoryStateOnSurface& tsos,
								 RefCountedKinematicVertex vertex,
								 reco::Vertex wrtVertex);

  std::pair<bool, Measurement1D> signedImpactParameter3D(const TrajectoryStateOnSurface& tsos,
							 RefCountedKinematicVertex vertex,
							 reco::Vertex wrtVertex);


  math::PtEtaPhiMLorentzVector daughter_p4(std::vector< RefCountedKinematicParticle > fitted_children, size_t i);
  
  // Vertex probability calculater
  std::tuple<Float_t, TransientVertex> vertexProb( const std::vector<reco::TransientTrack>& tracks);

  void recursiveDaughters(size_t index,
			  int rank, 
			  const reco::GenParticleCollection &src,
			  std::vector<size_t> &allIndices,
			  std::vector<int> &pdgs,
			  std::vector<int> &layers,
			  std::vector<float> &ppt,
			  std::vector<float> &peta,
			  std::vector<float> &pphi,
			  std::vector<int> &isfinal,
			  bool verbose
			  );
  

  reco::Track fix_track(const reco::Track *tk, double delta = 1e-8);
  reco::Track fix_track(const reco::TrackRef& tk);


  bool isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle);
  
  std::tuple<Bool_t, RefCountedKinematicParticle, RefCountedKinematicVertex, RefCountedKinematicTree> KinematicFit(std::vector<RefCountedKinematicParticle> particles, Float_t constrain_mass, Float_t constrain_error);

  bool basicPFcut(pat::PackedCandidate pf);

};





#endif
