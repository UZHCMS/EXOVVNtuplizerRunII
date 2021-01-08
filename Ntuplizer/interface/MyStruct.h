#ifndef MyStruct_H
#define MyStruct_H

#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MultiTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicTree.h"

// kinematic fit package
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"


#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"


struct particle_cand 
{
  
  // Longitudinal impact parameter and its significance
  Float_t lip;
  Float_t lips;
  
  // Impact parameter for the PV and its significance
  Float_t pvip;
  Float_t pvips;
  
  // Flight length and its significance
  Float_t fl3d;
  Float_t fls3d;
  
  // opening angle
  Float_t alpha;
  
};

struct attribute
{

  Float_t doca3d;
  Float_t doca3de;
  Float_t doca3ds;

  Float_t doca2d;
  Float_t doca2de;
  Float_t doca2ds;

//  Float_t doca1d;
//  Float_t doca1de;
//  Float_t doca1ds;

  Bool_t isRight;

  Float_t dz;
  
  Bool_t isAssociate;
  Int_t pvAssociationQuality;
  Float_t pt;
  Float_t eta;
  Float_t phi;
  Int_t charge;
  Float_t mass;
  
  Bool_t isBdecay;
  Int_t isBdecaypdg;
  Int_t isBdecayppdg;
  Bool_t isSignal;
  Int_t nprong;

  Float_t near_dz; 

};







struct taucand2{
  
  Int_t cand_tau_id1;
  Int_t cand_tau_id2;
  Int_t cand_tau_id3;
  Float_t cand_tau_pt;
  Float_t cand_tau_eta;
  Float_t cand_tau_phi;
  Float_t cand_tau_mass;
//  Float_t cand_tau_lip;
//  Float_t cand_tau_lips;
//  Float_t cand_tau_pvip;
//  Float_t cand_tau_pvips;
//  Float_t cand_tau_fl3d;
//  Float_t cand_tau_fls3d;
//  Float_t cand_tau_alpha;
  particle_cand cand_tau;
  Float_t cand_tau_vprob;
  Float_t cand_tau_vx;
  Float_t cand_tau_vy;
  Float_t cand_tau_vz;
  Float_t cand_tau_max_dr_3prong;	
  Int_t cand_tau_charge;
  Bool_t cand_tau_isRight;
  Bool_t cand_tau_isRight1;
  Bool_t cand_tau_isRight2;
  Bool_t cand_tau_isRight3;
  //  Float_t cand_tau_dr1;
  //  Float_t cand_tau_dr2;
  //  Float_t cand_tau_dr3;
  //  Float_t cand_tau_ptres1;
  //  Float_t cand_tau_ptres2;
  //  Float_t cand_tau_ptres3;
  Int_t cand_tau_matched_ppdgId;
  Float_t cand_tau_matched_gentaupt;
  Float_t cand_tau_sumofdnn;
  Float_t cand_tau_sumofdnn_others;
  Float_t cand_tau_pi1_dnn;
  Float_t cand_tau_pi2_dnn;
  Float_t cand_tau_pi3_dnn;
  Float_t cand_tau_pi1_doca;
  Float_t cand_tau_pi2_doca;
  Float_t cand_tau_pi3_doca;
  Int_t cand_tau_pi1_pv;
  Int_t cand_tau_pi2_pv;
  Int_t cand_tau_pi3_pv;
  Float_t cand_b_vprob;
  Float_t cand_b_vx;
  Float_t cand_b_vy;
  Float_t cand_b_vz;
  Float_t cand_b_pt;
  Float_t cand_b_eta; 
  Float_t cand_b_phi;
  Float_t cand_b_mass;
  particle_cand cand_b;
  //  Float_t cand_b_lip;
  //  Float_t cand_b_lips;
//  Float_t cand_b_pvip;
//  Float_t cand_b_pvips;
//  Float_t cand_b_fl3d;
//  Float_t cand_b_fls3d;
//  Float_t cand_b_alpha;
  Float_t cand_b_iso;
  Float_t cand_b_iso_ntracks;
  Float_t cand_b_iso_mindoca;
  Float_t cand_b_iso_nocut;
  Float_t cand_b_iso_ntracks_nocut;
  Float_t cand_b_iso_mindoca_nocut;
    
  bool operator<(const taucand2& another) const { 
    return cand_tau_pt > another.cand_tau_pt;
  }
};



struct taucandsimple{
  
  Int_t cand_tau_id1;
  Int_t cand_tau_id2;
  Int_t cand_tau_id3;
  Float_t cand_tau_pt;
  Int_t cand_tau_charge;

  bool operator<(const taucandsimple& another) const { 
    return cand_tau_pt > another.cand_tau_pt;
  }
};


struct taucandgen{
  
  Int_t cand_tau_id1;
  Int_t cand_tau_id2;
  Int_t cand_tau_id3;
  Float_t cand_gentau_pt;
  Float_t cand_gentau_eta;
  Float_t cand_gentau_phi;
  Float_t cand_gentau_mass;
  Float_t cand_gentau_pt_bd;
  Float_t cand_gentau_eta_bd;
  Float_t cand_gentau_phi_bd;
  Float_t cand_gentau_mass_bd;

  bool operator<(const taucandgen& another) const { 
    return cand_gentau_pt > another.cand_gentau_pt;
  }
};



struct pfcand_struct{
  
  Int_t idx;
  //  Float_t cand_absdz;
  pat::PackedCandidate pfcand;
  reco::TransientTrack track;
  attribute pfaux;
    
  bool operator<(const pfcand_struct& another) const { 
    //    return cand_absdz < another.cand_absdz;
    return TMath::Abs(pfaux.dz) < TMath::Abs(another.pfaux.dz);
  }
};



struct taucand{
  
  //  Int_t cand_tau_id1;
  //  Int_t cand_tau_id2;
  //  Int_t cand_tau_id3;
  math::PtEtaPhiMLorentzVector cand_tlv_tau_fit;
  //  Float_t cand_tau_pt;
  //  Float_t cand_tau_eta;
  //  Float_t cand_tau_phi;
  //  Float_t cand_tau_mass;
//  Float_t cand_tau_lip;
//  Float_t cand_tau_lips;
//  Float_t cand_tau_pvip;
//  Float_t cand_tau_pvips;
//  Float_t cand_tau_fl3d;
//  Float_t cand_tau_fls3d;
//  Float_t cand_tau_alpha;
  particle_cand cand_tau;
  RefCountedKinematicVertex cand_tau_vertex;
  //  Float_t cand_tau_vprob;
  //  Float_t cand_tau_vx;
  //  Float_t cand_tau_vy;
  //  Float_t cand_tau_vz;
  Float_t cand_tau_max_dr_3prong;	
  Int_t cand_tau_charge;
  Bool_t cand_tau_isRight;
  Bool_t cand_tau_isRight1;
  Bool_t cand_tau_isRight2;
  Bool_t cand_tau_isRight3;
  //  Float_t cand_tau_dr1;
  //  Float_t cand_tau_dr2;
  //  Float_t cand_tau_dr3;
  //  Float_t cand_tau_ptres1;
  //  Float_t cand_tau_ptres2;
  //  Float_t cand_tau_ptres3;
  Int_t cand_tau_matched_ppdgId;
  Float_t cand_tau_matched_gentaupt;
  Float_t cand_tau_sumofdnn;
  Float_t cand_tau_sumofdnn_others;
  Float_t cand_tau_pi1_dnn;
  Float_t cand_tau_pi2_dnn;
  Float_t cand_tau_pi3_dnn;
//  Float_t cand_tau_pi1_doca;
//  Float_t cand_tau_pi2_doca;
//  Float_t cand_tau_pi3_doca;
//  Int_t cand_tau_pi1_pv;
//  Int_t cand_tau_pi2_pv;
//  Int_t cand_tau_pi3_pv;
  RefCountedKinematicVertex cand_b_vertex;
  //  Float_t cand_b_vprob;
  //  Float_t cand_b_vx;
  //  Float_t cand_b_vy;
  //  Float_t cand_b_vz;
  RefCountedKinematicParticle cand_b_part;
//  Float_t cand_b_pt;
//  Float_t cand_b_eta; 
//  Float_t cand_b_phi;
//  Float_t cand_b_mass;
  particle_cand cand_b;
  //  Float_t cand_b_lip;
  //  Float_t cand_b_lips;
//  Float_t cand_b_pvip;
//  Float_t cand_b_pvips;
//  Float_t cand_b_fl3d;
//  Float_t cand_b_fls3d;
//  Float_t cand_b_alpha;
  std::vector<float> cand_b_iso;
  std::vector<int> cand_b_iso_ntracks;
  std::vector<float> cand_b_iso_mindoca;
  //  Float_t cand_b_iso_nocut;
  //  Float_t cand_b_iso_ntracks_nocut;
  //  Float_t cand_b_iso_mindoca_nocut;
  //  attribute aux1;
  //  attribute aux2;
  //  attribute aux3;
  pfcand_struct cand_pf1;
  pfcand_struct cand_pf2;
  pfcand_struct cand_pf3;
    
  bool operator<(const taucand& another) const { 
    return cand_tlv_tau_fit.Pt() > another.cand_tlv_tau_fit.Pt();
  }
};




#endif
