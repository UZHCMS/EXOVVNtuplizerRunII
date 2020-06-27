#ifndef MyStruct_H
#define MyStruct_H

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



struct taucand{
  
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
  Float_t cand_tau_pi1_dnn;
  Float_t cand_tau_pi2_dnn;
  Float_t cand_tau_pi3_dnn;
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
    
  bool operator<(const taucand& another) const { 
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



struct pfcand{
  
  Int_t cand_idx;
  Float_t cand_absdz;
    
  bool operator<(const pfcand& another) const { 
    return cand_absdz < another.cand_absdz;
  }
};






#endif
