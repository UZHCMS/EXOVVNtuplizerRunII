#include "../interface/ZTauTauNtuplizer.h"


//===================================================================================================================
ZTauTauNtuplizer::ZTauTauNtuplizer( edm::EDGetTokenT<pat::MuonCollection>    muonToken   ,
				    edm::EDGetTokenT<pat::ElectronCollection>    electronToken   ,
				    edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
				    edm::EDGetTokenT<reco::BeamSpot>             beamToken, 
				    edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken,
				    edm::EDGetTokenT<edm::TriggerResults> triggertoken,
				    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobject,
				    edm::EDGetTokenT<reco::GenParticleCollection> genptoken,
				    std::map< std::string, bool >& runFlags,
				    std::map< std::string, double >& runValues,
				    std::map< std::string, std::string >& runStrings,
                                    NtupleBranches* nBranches) 
: CandidateNtuplizer ( nBranches )
  , muonToken_	        ( muonToken )
  , electronToken_	        ( electronToken )
  , verticeToken_          ( verticeToken )
  , bsToken_ (beamToken)
  , packedpfcandidatesToken_(packedpfcandidatesToken) 
  , HLTtriggersToken_	( triggertoken )
  , triggerObjects_	( triggerobject )
  , genParticlesToken_( genptoken )
  , runOnMC_   (runFlags["runOnMC"])  , verbose_   (runFlags["verbose"])
  
{

  if(verbose_){
    std::cout << "[ZTauTauNtuplizer] runOnMC    = " << runOnMC_ << std::endl;
  }

}




//===================================================================================================================
ZTauTauNtuplizer::~ZTauTauNtuplizer( void )
{
}



bool ZTauTauNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  if(verbose_) std::cout << "[ZTauTauNtuplizer] ---------------- event, run, lumi = " << event.id().event() << " " << event.id().run() << " " << event.id().luminosityBlock() << "----------------" << std::endl;
  

  std::vector<TLorentzVector> match_tlv;
  std::vector<int> match_dm;
  Int_t n_tauhad = 0;
  Int_t n_tauhad_1prong = 0;
  Int_t n_tauhad_3prong = 0;
  Int_t n_muon = 0;
  Int_t n_electron = 0;
  
  if(runOnMC_){ 
    
    event.getByToken(genParticlesToken_ , genParticles_);   

    int pion_counter = 0;
    TLorentzVector pion1; pion1.SetPtEtaPhiM(0., 0., 0., 0.); 
    TLorentzVector pion2; pion2.SetPtEtaPhiM(0., 0., 0., 0.);
    TLorentzVector pion3; pion3.SetPtEtaPhiM(0., 0., 0., 0.);
    double pion1charge = 0;
    double pion2charge = 0;
    double pion3charge = 0;

    for(unsigned int i = 0; i < genParticles_->size(); ++i) {
      const reco::GenParticle* gen_particle = &genParticles_->at(i);
      if ( TMath::Abs(gen_particle->pdgId()) == 15 && gen_particle->status() ==2 && TMath::Abs(gen_particle->mother(0)->pdgId())==23 && pion_counter==0) {//  && gen_particle->status() == 62) {
	
        for (size_t j = 0; j < gen_particle->numberOfDaughters(); ++j) {
          if (TMath::Abs(gen_particle->daughter(j)->pdgId())==211)
            pion_counter++;
        }
        if (pion_counter!=3) {continue;}
        std::vector<double> pt, eta, phi, mass, charge;
        for (size_t j = 0; j < gen_particle->numberOfDaughters(); ++j) {
          if (TMath::Abs(gen_particle->daughter(j)->pdgId())==211) {
            pt.push_back(gen_particle->daughter(j)->pt());
            eta.push_back(gen_particle->daughter(j)->eta());
            phi.push_back(gen_particle->daughter(j)->phi());
            mass.push_back(gen_particle->daughter(j)->mass());
            if (gen_particle->daughter(j)->pdgId() > 0) { charge.push_back(1); } else { charge.push_back(-1); }
          }
        }

        if ((pt.at(0) > pt.at(1)) && (pt.at(0) > pt.at(2))) {
          pion1charge = charge.at(0);
          pion1.SetPtEtaPhiM( pt.at(0), eta.at(0), phi.at(0), mass.at(0) );
          if (pt.at(1)>pt.at(2)) {
            pion2charge = charge.at(1);
            pion3charge = charge.at(2);
            pion2.SetPtEtaPhiM( pt.at(1), eta.at(1), phi.at(1), mass.at(1) );
            pion3.SetPtEtaPhiM( pt.at(2), eta.at(2), phi.at(2), mass.at(2) );
          } else {
            pion2charge = charge.at(2);
            pion3charge = charge.at(1);
            pion2.SetPtEtaPhiM( pt.at(2), eta.at(2), phi.at(2), mass.at(2) );
            pion3.SetPtEtaPhiM( pt.at(1), eta.at(1), phi.at(1), mass.at(1) );
          }
        }

        if ((pt.at(1) > pt.at(0)) && (pt.at(1) > pt.at(2))) {
          pion1charge = charge.at(1);
          pion1.SetPtEtaPhiM( pt.at(1), eta.at(1), phi.at(1), mass.at(1) );
          if (pt.at(0)>pt.at(2)) {
            pion2charge = charge.at(0);
            pion3charge = charge.at(2);
            pion2.SetPtEtaPhiM( pt.at(0), eta.at(0), phi.at(0), mass.at(0) );
            pion3.SetPtEtaPhiM( pt.at(2), eta.at(2), phi.at(2), mass.at(2) );
          } else {
            pion2charge = charge.at(2);
            pion3charge = charge.at(0);
            pion2.SetPtEtaPhiM( pt.at(2), eta.at(2), phi.at(2), mass.at(2) );
            pion3.SetPtEtaPhiM( pt.at(0), eta.at(0), phi.at(0), mass.at(0) );
          }
        }

        if ((pt.at(2) > pt.at(0)) && (pt.at(2) > pt.at(1))) {
          pion1charge = charge.at(2);
          pion1.SetPtEtaPhiM( pt.at(2), eta.at(2), phi.at(2), mass.at(2) );
          if (pt.at(0)>pt.at(1)) {
            pion2charge = charge.at(0);
            pion3charge = charge.at(1);
            pion2.SetPtEtaPhiM( pt.at(0), eta.at(0), phi.at(0), mass.at(0) );
            pion3.SetPtEtaPhiM( pt.at(1), eta.at(1), phi.at(1), mass.at(1) );
          } else {
            pion2charge = charge.at(1);
            pion3charge = charge.at(0);
            pion2.SetPtEtaPhiM( pt.at(1), eta.at(1), phi.at(1), mass.at(1) );
            pion3.SetPtEtaPhiM( pt.at(0), eta.at(0), phi.at(0), mass.at(0) );
          }
        }



        if (pion1charge+pion2charge==0 && pion1charge+pion3charge==0) {
          nBranches_->truth_tau_dipion1_mass.push_back((pion1+pion2).M());
          nBranches_->truth_tau_dipion1_pt.push_back  ((pion1+pion2).Pt());
          nBranches_->truth_tau_dipion1_eta.push_back ((pion1+pion2).Eta());
          nBranches_->truth_tau_dipion1_phi.push_back ((pion1+pion2).Phi());

          nBranches_->truth_tau_dipion2_mass.push_back((pion1+pion3).M());
          nBranches_->truth_tau_dipion2_pt.push_back  ((pion1+pion3).Pt());
          nBranches_->truth_tau_dipion2_eta.push_back ((pion1+pion3).Eta());
          nBranches_->truth_tau_dipion2_phi.push_back ((pion1+pion3).Phi());
        } else if (pion1charge+pion2charge==0 && pion2charge+pion3charge==0) {
          nBranches_->truth_tau_dipion1_mass.push_back((pion1+pion2).M());
          nBranches_->truth_tau_dipion1_pt.push_back  ((pion1+pion2).Pt());
          nBranches_->truth_tau_dipion1_eta.push_back ((pion1+pion2).Eta());
          nBranches_->truth_tau_dipion1_phi.push_back ((pion1+pion2).Phi());

          nBranches_->truth_tau_dipion2_mass.push_back((pion2+pion3).M());
          nBranches_->truth_tau_dipion2_pt.push_back  ((pion2+pion3).Pt());
          nBranches_->truth_tau_dipion2_eta.push_back ((pion2+pion3).Eta());
          nBranches_->truth_tau_dipion2_phi.push_back ((pion2+pion3).Phi());
        } else if (pion1charge+pion3charge==0 && pion2charge+pion3charge==0) {
          nBranches_->truth_tau_dipion1_mass.push_back((pion1+pion3).M());
          nBranches_->truth_tau_dipion1_pt.push_back  ((pion1+pion3).Pt());
          nBranches_->truth_tau_dipion1_eta.push_back ((pion1+pion3).Eta());
          nBranches_->truth_tau_dipion1_phi.push_back ((pion1+pion3).Phi());

          nBranches_->truth_tau_dipion2_mass.push_back((pion2+pion3).M());
          nBranches_->truth_tau_dipion2_pt.push_back  ((pion2+pion3).Pt());
          nBranches_->truth_tau_dipion2_eta.push_back ((pion2+pion3).Eta());
          nBranches_->truth_tau_dipion2_phi.push_back ((pion2+pion3).Phi());
        }
      }
    }


    /////////////////////////////////////////////////



    for(size_t jpgp=0; jpgp < genParticles_->size(); ++jpgp){

      if(TMath::Abs((*genParticles_)[jpgp].pdgId())==15 && 
	 TMath::Abs((*genParticles_)[jpgp].status())==2 && 
	 TMath::Abs((*genParticles_)[jpgp].mother(0)->pdgId())==23){

	TLorentzVector p_gentau;
	Int_t n_charged_pions = 0;
	Int_t n_neutral_pions = 0;
	Int_t n_mu_decay = 0;
	Int_t n_e_decay = 0;
	Int_t n_occurance = 0;

	for(unsigned int jdau=0; jdau < (*genParticles_)[jpgp].numberOfDaughters(); jdau++){

	  Int_t decay_pdgid = TMath::Abs((*genParticles_)[jpgp].daughter(jdau)->pdgId());
	  
	  if(decay_pdgid==211) n_charged_pions += 1;
	  else if(decay_pdgid==111) n_neutral_pions += 1;
	  else if(decay_pdgid==13) n_mu_decay += 1;
	  else if(decay_pdgid==11) n_e_decay += 1;

	  if(!(decay_pdgid==12 || decay_pdgid==14 || decay_pdgid==16)){
	    TLorentzVector _genvis_;
	    _genvis_.SetPtEtaPhiM((*genParticles_)[jpgp].daughter(jdau)->pt(),
				  (*genParticles_)[jpgp].daughter(jdau)->eta(),
				  (*genParticles_)[jpgp].daughter(jdau)->phi(),
				  (*genParticles_)[jpgp].daughter(jdau)->mass()
				  );

	    p_gentau += _genvis_;
	    
	  }
	}

	match_tlv.push_back(p_gentau);
	
	Int_t dm = -1;

	if(n_e_decay==1){
	  dm = -2;
	  n_electron += 1;
	}else if(n_mu_decay==1){
	  dm = -1;
	  n_muon += 1;
	}else if(n_charged_pions==1 && n_neutral_pions==0){
	  dm = 0;
	  n_tauhad += 1; 
	  n_tauhad_1prong += 1;
	}else if(n_charged_pions==1 && n_neutral_pions>=1){
	  dm = 1; 
	  n_tauhad += 1; 
	  n_tauhad_1prong += 1;
	}else if(n_charged_pions==3 && n_neutral_pions==0){
	  dm = 10;
	  n_tauhad += 1; 
	  n_tauhad_3prong += 1;
	}else if(n_charged_pions==3 && n_neutral_pions>=1){
	  dm = 11;
	  n_tauhad += 1; 
	  n_tauhad_3prong += 1;
	}else{
	  dm = -9;
	  n_tauhad += 1; 
	}

	match_dm.push_back(dm);

	n_occurance+=1;
      }
    }
  }






  nBranches_->ZTauTau_n_muon = n_muon;
  nBranches_->ZTauTau_n_electron = n_electron;
  nBranches_->ZTauTau_n_tauhad = n_tauhad;
  nBranches_->ZTauTau_n_tauhad_1prong = n_tauhad_1prong;
  nBranches_->ZTauTau_n_tauhad_3prong = n_tauhad_3prong;


  
  /********************************************************************
   *
   * Step1: Check if the J/psi trigger is fired.
   * Namely, HLT_DoubleMu4_JpsiTrk_Displaced_v
   *
   ********************************************************************/

  event.getByToken(HLTtriggersToken_, HLTtriggers_);
  nBranches_->cutflow->Fill(10); // just to save total events here ... 

  bool isTriggered = false;
  const edm::TriggerNames& trigNames = event.triggerNames(*HLTtriggers_);

  for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {

    //    string trig = trigNames.triggerName(i);
    //    string str2find = "HLT_IsoMu24_v*";

    //    bool test = trigNames.triggerName(i).find("HLT_IsoMu24_v*")!= std::string::npos;
    //    boaol test = trig.find(str2find)!= std::string::npos;
    //
    //    std::cout << trigNames.triggerName(i) <<  " " << test << " " << HLTtriggers_->accept(i) << std::endl;
    

    if(trigNames.triggerName(i).find("HLT_IsoMu24_v")!= std::string::npos){
      if(HLTtriggers_->accept(i)){
	isTriggered = true;
	//	std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      }
    }
  }

  if(!isTriggered) return false;
  nBranches_->cutflow->Fill(1);

  if(verbose_) std::cout << "[ZTauTauNtuplizer] Trigger fired:" << isTriggered << std::endl;

  /********************************************************************
   *
   * Step2: pre-select muons for building J/psi candidates ... 
   * For muons, no requirement applied (soft muon ID is required once ther vertex is chosen)
   *
   ********************************************************************/

  event.getByToken(verticeToken_   , vertices_     );
  event.getByToken(bsToken_, beamspot_);
  event.getByToken(muonToken_	, muons_    );
  event.getByToken(triggerObjects_  , triggerObjects);

  //  reco::BeamSpot& beamspot = beamspot_;

  //  std::cout << "------------------------------------------" << std::endl;

  std::vector<pat::Muon> muoncollection;
  muoncollection.clear();


  size_t idx_muon = -1;
  
  for(size_t imuon = 0; imuon < muons_->size(); ++ imuon){

    const pat::Muon & muon = (*muons_)[imuon];

    if(TMath::Abs(muon.pdgId())!=13) continue;
      //      std::cout << "TEST_not_good:" <<muon.pdgId() <<std::endl;
    //    }

    if(muon.pt() < 26.) continue;
    if(TMath::Abs(muon.eta()) > 2.4) continue;
    if(!(muon.track().isNonnull())) continue;

    if(aux.MuonPFIso(muon) > 0.15) continue;
    if(!muon.isPFMuon()) continue;

    muoncollection.push_back(muon);
    idx_muon = imuon;
  }

  nBranches_->nmuon->Fill( muoncollection.size() );


  if(!( muoncollection.size() == 1)) return false;



  Int_t n_vetomuon = 0;

  for(size_t imuon = 0; imuon < muons_->size(); ++ imuon){

    const pat::Muon & muon = (*muons_)[imuon];

    if(imuon == idx_muon){
      //      std::cout << "exclude" << std::endl;
      continue;
    }
    if(TMath::Abs(muon.pdgId())!=13) continue;
    
    if(muon.pt() < 5.) continue;
    if(TMath::Abs(muon.eta()) > 2.4) continue;
    if(!(muon.track().isNonnull())) continue;

    //    if(aux.MuonPFIso(muon) > 0.15) continue;
    if(!muon.isLooseMuon()) continue;


    n_vetomuon += 1;
  }

  event.getByToken(electronToken_	, electrons_    );
  
  Int_t n_vetoelectron = 0;

  for(size_t ielectron = 0; ielectron < electrons_->size(); ++ ielectron){

    const pat::Electron & electron = (*electrons_)[ielectron];

    if(electron.pt() < 5.) continue;
    if(TMath::Abs(electron.eta()) > 2.4) continue;
    if (!electron.passConversionVeto()) continue;
    const reco::GsfTrackRef gsfTrk = electron.gsfTrack();

    if(!gsfTrk.isNonnull()) continue;

    n_vetoelectron += 1;
  }

  nBranches_->ZTauTau_n_vetomuon = n_vetomuon;
  nBranches_->ZTauTau_n_vetoelectron = n_vetoelectron;





  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

  const reco::TrackRef track_muon = muoncollection[0].muonBestTrack();
  reco::TransientTrack tt_muon = (*builder).build(aux.fix_track(track_muon));

  KinematicParticleFactoryFromTransientTrack pFactory;

  nBranches_->cutflow->Fill(2);

  if(verbose_) std::cout << "[ZTauTauNtuplizer] At least 2 muons are found" << std::endl;

  /********************************************************************
   *
   * Step3: building J/psi candidates 
   *
   ********************************************************************/

  nBranches_->cutflow->Fill(3);

  if(verbose_) std::cout << "[ZTauTauNtuplizer] J/psi found" << std::endl;

  /********************************************************************
   *
   * Step4: Kinematic fit for the J/psi candidate
   *        Use kinematicFitPackage
   *
   ********************************************************************/

  nBranches_->cutflow->Fill(4);

  /********************************************************************
   *
   * Step5: determine bbbar-PV, extrapolated back from the J/psi candidate
   *
   * We tried several possibilities
   * 
   * 1) using minimum lip
   * 2) using minimum pvip
   * 3) using PV (first content of the PV collection)
   *
   * and found that 1) or 2) is the best. Here, we use method (1).
   * 
   * Note: we can refit the vertex each time, by excluding the muons from J/psi
   * but since the efficiency of selecting the right vertex is already quite high (w/o vertex refit)
   * I don't think it is necessary to complicate things at this stage 
   * (we can do the refit once we select B candidate)
   *
   ********************************************************************/


  // define extrapolator
  edm::ESHandle<MagneticField> fieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(fieldHandle);
  fMagneticField = fieldHandle.product();
  
  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);
  //  TransverseImpactPointExtrapolator extrapolatort(fMagneticField);

  //  TSCBLBuilderNoMaterial blsBuilder;

//  Float_t max_criteria = 999;
  reco::Vertex closestVertex = *(vertices_->begin());

  nBranches_->cutflow->Fill(5);

  if(verbose_) std::cout << "[ZTauTauNtuplizer] J/psi soft muon ID passed" << std::endl;

  /********************************************************************
   *
   * Adding more attributes for each PF candidate
   *
   ********************************************************************/

  event.getByToken( packedpfcandidatesToken_               , packedpfcandidates_      ); 

  std::vector<pfcand_struct> pfcands;

  
  for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){   
      
    pat::PackedCandidate pf = (*packedpfcandidates_)[ii];

    if(!aux.basicPFcut(pf)) continue;

    if(pf.pt() < 0.5) continue;

    Float_t precut_dz = pf.vz() - closestVertex.position().z();
    if(TMath::Abs(precut_dz) > 0.12) continue;

    reco::TransientTrack  tt_track = (*builder).build(aux.fix_track(&pf.pseudoTrack()));
    
    pfcand_struct _cand_ = {
      (Int_t)ii,
      pf,
      tt_track,
    };
    
    pfcands.push_back(_cand_);


  }

  // sorted by pt
  sort(pfcands.begin(), pfcands.end());



  Int_t numOfch = (size_t)pfcands.size();
  
  if(numOfch<3) return false;
  nBranches_->cutflow->Fill(6);


  if(verbose_) std::cout << "[ZTauTauNtuplizer] Starts to build tau candidate out of " << numOfch << " pion candidates" << std::endl;

  std::vector<taucand_ztt> cands;
    

  for(int iii = 0; iii < numOfch; iii ++){
      
    pat::PackedCandidate pf1 = pfcands[iii].pfcand;
    reco::TransientTrack track1 = pfcands[iii].track;

    for(int jjj = iii+1; jjj < numOfch; jjj ++){
	
      pat::PackedCandidate pf2 = pfcands[jjj].pfcand;
      reco::TransientTrack track2 = pfcands[jjj].track;


      for(int kkk = jjj+1; kkk < numOfch; kkk ++){

	pat::PackedCandidate pf3 = pfcands[kkk].pfcand;
	reco::TransientTrack track3 = pfcands[kkk].track;

	Int_t tau_charge = pf1.charge() + pf2.charge() + pf3.charge(); 

	if(TMath::Abs(tau_charge)!=(int)1) continue; 


	/* reconstruct taus*/

	std::vector<RefCountedKinematicParticle> tauParticles;

	tauParticles.push_back(pFactory.particle(track1, aux.pion_mass, chi, ndf, aux.pion_sigma));
	tauParticles.push_back(pFactory.particle(track2, aux.pion_mass, chi, ndf, aux.pion_sigma));
	tauParticles.push_back(pFactory.particle(track3, aux.pion_mass, chi, ndf, aux.pion_sigma));

	RefCountedKinematicParticle tau_part;
	RefCountedKinematicVertex tau_vertex;
	RefCountedKinematicTree tauTree;
	Bool_t taufit_flag;
	std::tie(taufit_flag, tau_part, tau_vertex, tauTree) = aux.KinematicFit(tauParticles, -1, -1);
	
	if(!taufit_flag) continue;

	//	if(TMath::Prob(tau_vertex->chiSquared(), tau_vertex->degreesOfFreedom()) <= c_vprob) continue;
	if(TMath::Prob(tau_vertex->chiSquared(), tau_vertex->degreesOfFreedom()) <= 0.1) continue;
	  
	std::vector< RefCountedKinematicParticle > tau_children = tauTree->finalStateParticles();
	  
	math::PtEtaPhiMLorentzVector pi1_fit = aux.daughter_p4(tau_children, 0);
	math::PtEtaPhiMLorentzVector pi2_fit = aux.daughter_p4(tau_children, 1);
	math::PtEtaPhiMLorentzVector pi3_fit = aux.daughter_p4(tau_children, 2);

	if(pi1_fit.Pt() < 10.) continue;

	//	math::PtEtaPhiMLorentzVector tlv_tau_fit = tt1_fit + tt2_fit + tt3_fit;
	math::PtEtaPhiMLorentzVector tlv_tau_fit = pi1_fit + pi2_fit + pi3_fit;

	//	if(tlv_tau_fit.Pt() < 3.){
	  //	  std::cout <<"remove pt" << std::endl;
	//	  continue;
	//	}


	if(tlv_tau_fit.M() > 1.7) continue;



	taucand_ztt _cand_ = {

	  tlv_tau_fit,
	  pi1_fit,
	  pi2_fit,
	  pi3_fit,
	  tau_part,
	  tau_vertex,
	  (Int_t) tau_charge,
	  pfcands[iii],
	  pfcands[jjj],
	  pfcands[kkk],
	};
	  
	cands.push_back(_cand_);
      }
    }
  }
    
  //  sort(cands.begin(), cands.end());


  if(cands.size()==0) return false;

  if(verbose_) std::cout << "[ZTauTauNtuplizer] " << cands.size() << " tau candidates were found" << std::endl;

  nBranches_->cutflow->Fill(7);

  //  Int_t ncomb = 0;

  for(int ic=0; ic < (int)cands.size(); ic++){

    /********************************************************************
     *
     * Step9: Filling normal branches
     *
     ********************************************************************/


    nBranches_->ZTauTau_tau_pt.push_back(cands[ic].cand_tlv_tau_fit.Pt());
    nBranches_->ZTauTau_tau_eta.push_back(cands[ic].cand_tlv_tau_fit.Eta());
    nBranches_->ZTauTau_tau_phi.push_back(cands[ic].cand_tlv_tau_fit.Phi());
    nBranches_->ZTauTau_tau_mass.push_back(cands[ic].cand_tlv_tau_fit.M());
    nBranches_->ZTauTau_tau_q.push_back(cands[ic].cand_tau_charge);

    nBranches_->ZTauTau_tau_vx.push_back(cands[ic].cand_tau_vertex->vertexState().position().x());
    nBranches_->ZTauTau_tau_vy.push_back(cands[ic].cand_tau_vertex->vertexState().position().y());
    nBranches_->ZTauTau_tau_vz.push_back(cands[ic].cand_tau_vertex->vertexState().position().z());


    Float_t max_dr_3prong = -1;

    Float_t dR_12 = reco::deltaR(cands[ic].cand_tlv_pi1_fit.Eta(), cands[ic].cand_tlv_pi1_fit.Phi(), 
				 cands[ic].cand_tlv_pi2_fit.Eta(), cands[ic].cand_tlv_pi2_fit.Phi());

    Float_t dR_13 = reco::deltaR(cands[ic].cand_tlv_pi1_fit.Eta(), cands[ic].cand_tlv_pi1_fit.Phi(), 
				 cands[ic].cand_tlv_pi3_fit.Eta(), cands[ic].cand_tlv_pi3_fit.Phi());

    Float_t dR_23 = reco::deltaR(cands[ic].cand_tlv_pi2_fit.Eta(), cands[ic].cand_tlv_pi2_fit.Phi(), 
				 cands[ic].cand_tlv_pi3_fit.Eta(), cands[ic].cand_tlv_pi3_fit.Phi());


    if(max_dr_3prong < dR_12) max_dr_3prong = dR_12;
    if(max_dr_3prong < dR_13) max_dr_3prong = dR_13;
    if(max_dr_3prong < dR_23) max_dr_3prong = dR_23;



    nBranches_->ZTauTau_tau_max_dr_3prong.push_back(max_dr_3prong);
    nBranches_->ZTauTau_tau_vprob.push_back(TMath::Prob(cands[ic].cand_tau_vertex->chiSquared(), cands[ic].cand_tau_vertex->degreesOfFreedom()) ); 

    bool flag_tau_isSignal = false;

    for(int im=0; im < (int)match_tlv.size(); im++){
      if(match_dm[im]==10 || match_dm[im]==11){
	Float_t dR_tau_gen = reco::deltaR(cands[ic].cand_tlv_tau_fit.Eta(),
					  cands[ic].cand_tlv_tau_fit.Phi(),
					  match_tlv[im].Eta(),
					  match_tlv[im].Phi());
	
	if(dR_tau_gen < 0.2) flag_tau_isSignal = true;
      }
    }
    
    nBranches_->ZTauTau_tau_isSignal.push_back(flag_tau_isSignal);



    /////////////////////////////////////////

    float iso = 0;
    int ntracks = 0;


    for( int itrk = 0; itrk < (int)packedpfcandidates_->size(); ++itrk ){   
      
      pat::PackedCandidate pf = (*packedpfcandidates_)[itrk];
      
      if(!aux.basicPFcut(pf)) continue;
      
      if(pf.pt() < 0.5) continue;
      
      if(itrk==cands[ic].cand_pf1.idx || 
	 itrk==cands[ic].cand_pf2.idx || 
	 itrk==cands[ic].cand_pf3.idx){
	//	std::cout << "overlapped ... removed!" << std::endl;
	continue;
      }
      
      Float_t iso_pt = pf.pt();
      Float_t iso_eta = pf.eta();
      Float_t iso_phi = pf.phi();
	
      Float_t iso_dr = reco::deltaR(iso_eta, iso_phi, cands[ic].cand_tlv_tau_fit.Eta(), cands[ic].cand_tlv_tau_fit.Phi());
				    
      
      if(iso_dr > 0.4) continue;
      
      if(iso_pt > 0.8){
	iso += iso_pt;
	ntracks += 1;
      }	  
    }


    nBranches_->ZTauTau_tau_iso.push_back(iso);
    nBranches_->ZTauTau_tau_ntracks.push_back(ntracks);


    ////////////////////////////////////






    std::vector<TLorentzVector> rhomass;
      
    TLorentzVector tlv_pion1;// = cands[ic].cand_tlv_pi1_fit;
    TLorentzVector tlv_pion2;// = cands[ic].cand_tlv_pi2_fit;
    TLorentzVector tlv_pion3;// = cands[ic].cand_tlv_pi3_fit;

    tlv_pion1.SetPtEtaPhiM(cands[ic].cand_tlv_pi1_fit.Pt(), cands[ic].cand_tlv_pi1_fit.Eta(), cands[ic].cand_tlv_pi1_fit.Phi(), cands[ic].cand_tlv_pi1_fit.M());
    tlv_pion2.SetPtEtaPhiM(cands[ic].cand_tlv_pi2_fit.Pt(), cands[ic].cand_tlv_pi2_fit.Eta(), cands[ic].cand_tlv_pi2_fit.Phi(), cands[ic].cand_tlv_pi2_fit.M());
    tlv_pion3.SetPtEtaPhiM(cands[ic].cand_tlv_pi3_fit.Pt(), cands[ic].cand_tlv_pi3_fit.Eta(), cands[ic].cand_tlv_pi3_fit.Phi(), cands[ic].cand_tlv_pi3_fit.M());

    Int_t q1 = cands[ic].cand_pf1.pfcand.charge();
    Int_t q2 = cands[ic].cand_pf2.pfcand.charge();
    Int_t q3 = cands[ic].cand_pf3.pfcand.charge();

    if(q1*q2 == -1){
      TLorentzVector tlv_rho = tlv_pion1 + tlv_pion2;
      rhomass.push_back(tlv_rho);
    }else if(q1*q2==1){
      TLorentzVector tlv_rho = tlv_pion1 + tlv_pion2;
    }
      
    if(q1*q3 == -1){
      TLorentzVector tlv_rho = tlv_pion1 + tlv_pion3;
      rhomass.push_back(tlv_rho);
    }else if(q1*q3==1){
      TLorentzVector tlv_rho = tlv_pion1 + tlv_pion3;
    }
      
    if(q2*q3 == -1){
      TLorentzVector tlv_rho = tlv_pion2 + tlv_pion3;
      rhomass.push_back(tlv_rho);
    }else if(q2*q3==1){
      TLorentzVector tlv_rho = tlv_pion2 + tlv_pion3;
    }
      

    if(rhomass.size()==2){

      if(rhomass.at(0).Pt() > rhomass.at(1).Pt()){      
	nBranches_->ZTauTau_tau_rhomass1.push_back(rhomass.at(0).M());
	nBranches_->ZTauTau_tau_rhomass2.push_back(rhomass.at(1).M());
	nBranches_->ZTauTau_tau_rhopt1.push_back(rhomass.at(0).Pt());
	nBranches_->ZTauTau_tau_rhopt2.push_back(rhomass.at(1).Pt());
      }else{
	nBranches_->ZTauTau_tau_rhomass1.push_back(rhomass.at(1).M());
	nBranches_->ZTauTau_tau_rhomass2.push_back(rhomass.at(0).M());
	nBranches_->ZTauTau_tau_rhopt1.push_back(rhomass.at(1).Pt());
	nBranches_->ZTauTau_tau_rhopt2.push_back(rhomass.at(0).Pt());
      }

    }else{
      nBranches_->ZTauTau_tau_rhomass1.push_back(-1);
      nBranches_->ZTauTau_tau_rhomass2.push_back(-1);
      nBranches_->ZTauTau_tau_rhopt1.push_back(-1);
      nBranches_->ZTauTau_tau_rhopt2.push_back(-1);
    }


    nBranches_->ZTauTau_tau_pi1_pt.push_back(tlv_pion1.Pt());
    nBranches_->ZTauTau_tau_pi1_eta.push_back(tlv_pion1.Eta());
    nBranches_->ZTauTau_tau_pi1_phi.push_back(tlv_pion1.Phi());
    nBranches_->ZTauTau_tau_pi1_mass.push_back(tlv_pion1.M());
    nBranches_->ZTauTau_tau_pi1_q.push_back(q1);
      
    nBranches_->ZTauTau_tau_pi2_pt.push_back(tlv_pion2.Pt());
    nBranches_->ZTauTau_tau_pi2_eta.push_back(tlv_pion2.Eta());
    nBranches_->ZTauTau_tau_pi2_phi.push_back(tlv_pion2.Phi());
    nBranches_->ZTauTau_tau_pi2_mass.push_back(tlv_pion2.M());
    nBranches_->ZTauTau_tau_pi2_q.push_back(q2);


    nBranches_->ZTauTau_tau_pi3_pt.push_back(tlv_pion3.Pt());
    nBranches_->ZTauTau_tau_pi3_eta.push_back(tlv_pion3.Eta());
    nBranches_->ZTauTau_tau_pi3_phi.push_back(tlv_pion3.Phi());
    nBranches_->ZTauTau_tau_pi3_mass.push_back(tlv_pion3.M());
    nBranches_->ZTauTau_tau_pi3_q.push_back(q3);


    // kinematic fit

    std::vector<RefCountedKinematicParticle> allParticles;
      
    allParticles.push_back(pFactory.particle(cands[ic].cand_pf1.track, aux.pion_mass, chi, ndf, aux.pion_sigma));
    allParticles.push_back(pFactory.particle(cands[ic].cand_pf2.track, aux.pion_mass, chi, ndf, aux.pion_sigma));
    allParticles.push_back(pFactory.particle(cands[ic].cand_pf3.track, aux.pion_mass, chi, ndf, aux.pion_sigma));
    allParticles.push_back(pFactory.particle(tt_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));


    RefCountedKinematicParticle bc_part;
    RefCountedKinematicVertex bc_vertex;
    RefCountedKinematicTree bcTree;
    Bool_t bcfit_flag;
    std::tie(bcfit_flag, bc_part, bc_vertex, bcTree) = aux.KinematicFit(allParticles, -1, -1);
	
    if(bcfit_flag){
 
      particle_cand Bcand; 
      Bcand = aux.calculateIPvariables(extrapolator, bc_part, bc_vertex, closestVertex);
	  
      nBranches_->ZTauTau_Z_pt.push_back(bc_part->currentState().globalMomentum().perp());
      nBranches_->ZTauTau_Z_eta.push_back(bc_part->currentState().globalMomentum().eta());
      nBranches_->ZTauTau_Z_phi.push_back(bc_part->currentState().globalMomentum().phi());
      nBranches_->ZTauTau_Z_mass.push_back(bc_part->currentState().mass());
      nBranches_->ZTauTau_Z_vprob.push_back( TMath::Prob(bc_vertex->chiSquared(), bc_vertex->degreesOfFreedom()) ); //bc_part->currentState().mass());
      nBranches_->ZTauTau_Z_lip.push_back(Bcand.lip);
      nBranches_->ZTauTau_Z_lips.push_back(Bcand.lips);
      nBranches_->ZTauTau_Z_pvip.push_back(Bcand.pvip);
      nBranches_->ZTauTau_Z_pvips.push_back(Bcand.pvips);
      nBranches_->ZTauTau_Z_fls3d.push_back(Bcand.fls3d);
      nBranches_->ZTauTau_Z_fl3d.push_back(Bcand.fl3d);
      nBranches_->ZTauTau_Z_alpha.push_back(Bcand.alpha);

      nBranches_->ZTauTau_Z_vx.push_back(bc_vertex->vertexState().position().x());
      nBranches_->ZTauTau_Z_vy.push_back(bc_vertex->vertexState().position().y());
      nBranches_->ZTauTau_Z_vz.push_back(bc_vertex->vertexState().position().z());

    }
  }

  nBranches_->ZTauTau_mu_pt = muoncollection[0].pt();
  nBranches_->ZTauTau_mu_eta = muoncollection[0].eta();
  nBranches_->ZTauTau_mu_phi = muoncollection[0].phi();
  nBranches_->ZTauTau_mu_mass = muoncollection[0].mass();
  nBranches_->ZTauTau_mu_q = muoncollection[0].charge();
  nBranches_->ZTauTau_mu_isLoose = muoncollection[0].isLooseMuon();
  nBranches_->ZTauTau_mu_isTight = muoncollection[0].isTightMuon(closestVertex);
  nBranches_->ZTauTau_mu_isPF = muoncollection[0].isPFMuon();
  nBranches_->ZTauTau_mu_isGlobal = muoncollection[0].isGlobalMuon();
  nBranches_->ZTauTau_mu_isTracker = muoncollection[0].isTrackerMuon();
  nBranches_->ZTauTau_mu_isSoft = muoncollection[0].isSoftMuon(closestVertex);
  nBranches_->ZTauTau_mu_vx = muoncollection[0].vx();
  nBranches_->ZTauTau_mu_vy = muoncollection[0].vy();
  nBranches_->ZTauTau_mu_vz = muoncollection[0].vz();
  nBranches_->ZTauTau_mu_dbiso = aux.MuonPFIso(muoncollection[0]);

  bool flag_mu_isSignal = false;

  for(int im=0; im < (int)match_tlv.size(); im++){
    if(match_dm[im]==-1){
      Float_t dR_mu_gen = reco::deltaR(
					muoncollection[0].eta(),
					muoncollection[0].phi(),
					match_tlv[im].Eta(),
					match_tlv[im].Phi());
      
      if(dR_mu_gen < 0.2) flag_mu_isSignal = true;
    }
  }
  
  nBranches_->ZTauTau_mu_isSignal = flag_mu_isSignal;

  

  nBranches_->ZTauTau_PV_vx = vertices_->begin()->position().x();
  nBranches_->ZTauTau_PV_vy = vertices_->begin()->position().y();
  nBranches_->ZTauTau_PV_vz = vertices_->begin()->position().z();

  nBranches_->ZTauTau_ncands = cands.size();

  //  if(!runOnMC_) return true;

  return true;

}

const reco::Candidate*  ZTauTauNtuplizer::checkMom(const reco::Candidate * candMom){
    int diquarks[] = { 1103,2101,2103,2203,3101,3103,3201,3203,3303,4101,4103,4201,4203,4301,4303,4403,5101,5103,5201,5203,5301,5303,5401, 5403,5503};
    if (candMom == nullptr) return nullptr;

    if (candMom->mother(0) == nullptr) {
        return candMom;
    }
    std::vector<int> B_hadron = {511,521,531,541,5112,5122,5132,5212,5232};
    std::vector<int>::iterator it0 = std::find(B_hadron.begin(), B_hadron.end(),  abs(candMom->pdgId()));
    if (it0!= B_hadron.end()){
        return candMom;    
    }
    
    int * p = std::find (diquarks, diquarks+25, candMom->mother(0)->pdgId());
    //std::cout << "  check mother  "<< abs(candMom->mother(0)->pdgId()) << std::endl; 
   if (abs(candMom->mother(0)->pdgId()) < 8  ||    \
        abs(candMom->mother(0)->pdgId())== 21 ||    \
        abs(candMom->mother(0)->pdgId())== 2212 ||  \
        (p != (diquarks+25))
        ){
       // std::cout << " ---- original mother   "<<  abs(candMom->pdgId()) << std::endl;
        return candMom;
    }
    else {
        candMom = checkMom(candMom->mother(0));
        //std::cout << "  intermediate mother   "<<  abs(candMom->pdgId())<< std::endl;
        return candMom;
    }
}



FreeTrajectoryState ZTauTauNtuplizer::initialFreeState(const reco::Track& tk, const MagneticField *field) {
  Basic3DVector<float> pos(tk.vertex());
  GlobalPoint gpos(pos);
  Basic3DVector<float> mom(tk.momentum());
  GlobalVector gmom(mom);
  GlobalTrajectoryParameters par(gpos, gmom, tk.charge(), field);
  CurvilinearTrajectoryError err(tk.covariance());
  return FreeTrajectoryState(par, err);
}

