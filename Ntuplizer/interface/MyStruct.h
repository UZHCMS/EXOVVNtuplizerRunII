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
  Float_t cand_tau_lip;
  Float_t cand_tau_lips;
  Float_t cand_tau_pvip;
  Float_t cand_tau_pvips;
  Float_t cand_tau_fl3d;
  Float_t cand_tau_fls3d;
  Float_t cand_tau_alpha;
  Float_t cand_tau_vprob;
  Float_t cand_tau_vx;
  Float_t cand_tau_vy;
  Float_t cand_tau_vz;
  Float_t cand_tau_iso;
  Float_t cand_tau_iso_ntracks;
  Float_t cand_tau_iso_mindoca;
  Float_t cand_tau_max_dr_3prong;	
  Int_t cand_tau_charge;
  Bool_t cand_tau_isRight;
  Int_t cand_tau_matched_ppdgId;
  Float_t cand_tau_matched_gentaupt;
  Float_t cand_b_vprob;
  Float_t cand_b_vx;
  Float_t cand_b_vy;
  Float_t cand_b_vz;
  Float_t cand_b_pt;
  Float_t cand_b_eta; 
  Float_t cand_b_phi;
  Float_t cand_b_mass;
  Float_t cand_b_lip;
  Float_t cand_b_lips;
  Float_t cand_b_pvip;
  Float_t cand_b_pvips;
  Float_t cand_b_fl3d;
  Float_t cand_b_fls3d;
  Float_t cand_b_alpha;
    
  bool operator<(const taucand& another) const { 
    return cand_tau_pt > another.cand_tau_pt;
  }
};


#endif
