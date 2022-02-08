#include "../interface/JpsiKNtuplizerE.h"


//===================================================================================================================
JpsiKNtuplizerE::JpsiKNtuplizerE( edm::EDGetTokenT<pat::MuonCollection>    muonToken   ,
				    edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
				    edm::EDGetTokenT<reco::BeamSpot>             beamToken, 
				    edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken,
				    edm::EDGetTokenT<edm::TriggerResults> triggertoken,
				    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobject,
				    edm::EDGetTokenT<reco::GenParticleCollection> genptoken,
				    edm::EDGetTokenT<pat::PackedGenParticleCollection> packedgenptoken,
				    std::map< std::string, bool >& runFlags,
				    std::map< std::string, double >& runValues,
				std::map< std::string, std::string >& runStrings, NtupleBranches* nBranches, TH1F* histGenWeights)

: CandidateNtuplizer ( nBranches )
  , muonToken_	        ( muonToken )
  , verticeToken_          ( verticeToken )
  , bsToken_ (beamToken)
  , packedpfcandidatesToken_(packedpfcandidatesToken) 
  , HLTtriggersToken_	( triggertoken )
  , triggerObjects_	( triggerobject )
  , genParticlesToken_( genptoken )
  , packedgenParticlesToken_( packedgenptoken )
  , runOnMC_   (runFlags["runOnMC"])
  , verbose_   (runFlags["verbose"])
  , c_dz (runValues["dzcut"])
  , isBkgBSample_ (runFlags["isBkgBSample"]) 
  , histGenWeights_ (histGenWeights)        
{

  if(verbose_){
    std::cout << "[JpsiKNtuplizerE] runOnMC    = " << runOnMC_ << std::endl;
    std::cout << "[JpsiKNtuplizerE] dzcut      = " << c_dz << std::endl;
  }
}




//===================================================================================================================
JpsiKNtuplizerE::~JpsiKNtuplizerE( void )
{
}



bool JpsiKNtuplizerE::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  if(verbose_) std::cout << "[JpsiKNtuplizerE] ---------------- event, run, lumi = " << event.id().event() << " " << event.id().run() << " " << event.id().luminosityBlock() << "----------------" << std::endl;
  

  TVector3 genvertex(-9.,-9.,-9.);
  TVector3 genvertex_sv(-9.,-9.,-9.);

  std::vector<TLorentzVector> gen_nr_mu;
  std::vector<TLorentzVector> gen_jpsi_mu;

  TLorentzVector pB_gen;
  TLorentzVector pJpsi_gen;
  Int_t nBc = 0;

  if(runOnMC_){ 
    
    event.getByToken(genParticlesToken_ , genParticles_);   
    event.getByToken(packedgenParticlesToken_ , packedgenParticles_);   

    for( unsigned p=0; p < genParticles_->size(); ++p){
      
      if(!(TMath::Abs((*genParticles_)[p].pdgId())==521 && (*genParticles_)[p].status()==2)) continue;
      
      auto _part_ = (*genParticles_)[p];
      
      bool isJpsi = false;
      
      for(auto d : _part_.daughterRefVector()) {
	if(TMath::Abs(d->pdgId()) == 443){
	  genvertex_sv = TVector3(d->vx(), d->vy(), d->vz());
	  isJpsi = true;
	}
      }
      
      pB_gen.SetPtEtaPhiM(_part_.pt(), _part_.eta(), _part_.phi(), _part_.mass());
      genvertex = aux.getVertex(_part_);
      nBc += 1;

      if(!isJpsi) continue;
                 
      for(auto d : _part_.daughterRefVector()) {
	if(TMath::Abs(d->pdgId()) == 443){
	  
	  pJpsi_gen.SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
	  
	  for (auto dd : d->daughterRefVector()) {
	    if(TMath::Abs(dd->pdgId())==13){
	      TLorentzVector jpsi_muon_tlv;
	      jpsi_muon_tlv.SetPtEtaPhiM(dd->pt(), dd->eta(), dd->phi(), dd->mass());
	      gen_jpsi_mu.push_back(jpsi_muon_tlv);
	    }
	  }
	  
	}
	if(TMath::Abs(d->pdgId())==13){
	  TLorentzVector nojpsi_muon_tlv;
	  nojpsi_muon_tlv.SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
	  
	  gen_nr_mu.push_back(nojpsi_muon_tlv);
	}
      }
    }
  }



  TLorentzVector q2_gen; 
  q2_gen = pB_gen - pJpsi_gen;
  
  if(runOnMC_){ 
    nBranches_->q2_nocut->Fill(q2_gen.M2());
  }



  /*****************
   * Save mother and weight if needed 
   ****************/
  // Implementation fo the weight for the B chain decay in the generic background B sample
  if (isBkgBSample_){ 
    //    std::cout << "---------------- New event-------------"<< std::endl;
      for(size_t p=0; p < genParticles_->size(); ++p){
          if (abs((*genParticles_)[p].pdgId())==443 and abs((*genParticles_)[p].daughter(0)->pdgId())==13 ){ 

              //             std:: cout<< "A) Jpsi direct mother is "    << abs((*genParticles_)[p].mother(0)->pdgId()) << std::endl;
              std::vector<int> B_hadron = {511,521,531,541,5112,5122,5132,5212,5232};
              motherID_ = abs((JpsiKNtuplizerE::checkMom(&(*genParticles_)[p]))->pdgId());
              // std:: cout<< " Jpsi status is "    <<  (*genParticles_)[p].status() <<std::endl;
              
              // std:: cout<< " Jpsi daugh is "    <<   (*genParticles_)[p].daughter(0)->pdgId() << std::endl;
              // std:: cout<< " Jpsi mother is "    << motherID_ << std::endl;
              
              std::vector<int>::iterator it = std::find(B_hadron.begin(), B_hadron.end(), motherID_); 
              int index; 
              if (motherID_ == 541){ genWeightBkgB_ =0; //set as default to 0 to Bc events in generic Bkg sample  because the Bc sample should be used
              } else if (it != B_hadron.end()) {  
                  index = std::distance(B_hadron.begin(), it);    
                  //                  std:: cout<< "index is " << index <<" weight is " << histGenWeights_->GetBinContent(index+2)<<std::endl;
                  genWeightBkgB_ = histGenWeights_->GetBinContent(index+2);
              } else {    
                  genWeightBkgB_ = histGenWeights_->GetBinContent(1); // in the first bin of the hist there is a generic 'other' for all the b decays not contained in the B_hadron vector;
              }
          }
      }
  }else{
      genWeightBkgB_ = 1; 
  } 





  
  /********************************************************************
   *
   * Step1: Check if the J/psi trigger is fired.
   * Namely, HLT_DoubleMu4_JpsiTrk_Displaced_v
   *
   ********************************************************************/

  event.getByToken(HLTtriggersToken_, HLTtriggers_);

  nBranches_->cutflow->Fill(10); // just to save total events here ... 
  nBranches_->cutflow->Fill(0.,genWeightBkgB_); // it should be 1 if the sample is not Bkg B generic sample, so 1 for data and 1 for Bc if flags properly set
  nBranches_->bweight->Fill(genWeightBkgB_);

  
  bool isTriggered = false;
  const edm::TriggerNames& trigNames = event.triggerNames(*HLTtriggers_);
  std::string finalTriggerName="";
  std::string finalTriggerFilterObjName="";

  for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {

    if(trigNames.triggerName(i).find("HLT_DoubleMu4_JpsiTrk_Displaced_v")!= std::string::npos){
      nBranches_->HLT_BPH_isFired[trigNames.triggerName(i)] = HLTtriggers_->accept(i);
      if(HLTtriggers_->accept(i)){
	isTriggered = true;
	finalTriggerName=trigNames.triggerName(i);  
	finalTriggerFilterObjName="hltJpsiTkVertexFilter";
      }
    }
  }


  if(!isTriggered) return false;
  nBranches_->cutflow->Fill(1);

  if(verbose_) std::cout << "[JpsiKNtuplizerE] Trigger fired" << std::endl;

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

  for(size_t imuon = 0; imuon < muons_->size(); ++ imuon){

    const pat::Muon & muon = (*muons_)[imuon];

    if(muon.pt() < 4) continue;
    if(TMath::Abs(muon.eta()) > 2.4) continue;
    if(!(muon.track().isNonnull())) continue;


    // Trigger matching

    bool trigMatch = false;

    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    
      obj.unpackPathNames(trigNames);
      obj.unpackFilterLabels(event, *HLTtriggers_);

      std::vector<std::string> pathNamesAll  = obj.pathNames(false);

      bool isPathExist = false;

      for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
	if(pathNamesAll[h]==finalTriggerName) isPathExist = true;
      }
      
      if(!isPathExist) continue;

      bool isFilterExist = false;
    
      for (unsigned hh = 0; hh < obj.filterLabels().size(); ++hh){
	
	//	if(obj.filterLabels()[hh].find(finalTriggerFilterObjName) != std::string::npos){
	if(obj.filterLabels()[hh].find("hltDisplacedmumuFilterDoubleMu4Jpsi") != std::string::npos){
	  isFilterExist = true;
	}
      }
      
      if(!isFilterExist) continue;
      
      if(TMath::Abs(obj.pdgId()) != 13) continue;
      
      Float_t trigger_dR = reco::deltaR(obj.eta(), obj.phi(),
					muon.eta(), muon.phi());
      
      if(trigger_dR < 0.015 &&
	 obj.pt()/muon.pt() > 0.85 &&
	 obj.pt()/muon.pt() < 1.15
	 ){
	trigMatch = true;
	//	std::cout << "Muon" << imuon << " matches to the trigger object dR = " << trigger_dR << ", obj pT = " << obj.pt() << ", muon pT = " << muon.pt() << " " << obj.pdgId() << std::endl;
      }
    }

    if(!trigMatch) continue;

    muoncollection.push_back(muon);
    //    muoncollection_id.push_back(imuon);
  }

  nBranches_->nmuon->Fill( muoncollection.size() );

  if(!( muoncollection.size() >= 2)) return false;

  nBranches_->cutflow->Fill(2);

  if(verbose_) std::cout << "[JpsiKNtuplizerE] At least 2 muons are found" << std::endl;

  /********************************************************************
   *
   * Step3: building J/psi candidates 
   *
   ********************************************************************/

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

  Float_t jpsi_max_pt = -1;
  unsigned int mcidx_mu1 = -1;
  unsigned int mcidx_mu2 = -1;
  TLorentzVector jpsi_tlv_highest;

  for(int imu = 0; imu < (int)muoncollection.size(); imu++){
    for(int jmu = imu+1; jmu < (int)muoncollection.size(); jmu++){

      const pat::Muon mu1 = muoncollection[imu];
      const pat::Muon mu2 = muoncollection[jmu];

      TLorentzVector tlv_mu1;
      TLorentzVector tlv_mu2;

      tlv_mu1.SetPtEtaPhiM(mu1.pt(), mu1.eta(), mu1.phi(), mu1.mass());
      tlv_mu2.SetPtEtaPhiM(mu2.pt(), mu2.eta(), mu2.phi(), mu2.mass());

      TLorentzVector tlv_jpsi = (tlv_mu1 + tlv_mu2);

      Float_t jpsi_mass = tlv_jpsi.M();
      Float_t jpsi_pt = tlv_jpsi.Pt();

      if(mu1.charge() + mu2.charge() !=0) continue;
      if(jpsi_mass < 2.95) continue; // a little bit broad winder to take into account FSR ...
      if(jpsi_mass > 3.25) continue;

      if(jpsi_max_pt < jpsi_pt){
	jpsi_max_pt = jpsi_pt;
	mcidx_mu1 = imu;
	mcidx_mu2 = jmu;
	jpsi_tlv_highest = tlv_jpsi;
      }
    }
  }

  if(jpsi_max_pt == -1) return false;
  nBranches_->cutflow->Fill(3);

  if(verbose_) std::cout << "[JpsiKNtuplizerE] J/psi found" << std::endl;

  /********************************************************************
   *
   * Step4: Kinematic fit for the J/psi candidate
   *        Use kinematicFitPackage
   *
   ********************************************************************/

  const reco::TrackRef track1_muon = muoncollection[mcidx_mu1].muonBestTrack();
  const reco::TrackRef track2_muon = muoncollection[mcidx_mu2].muonBestTrack();
  reco::TransientTrack tt1_muon = (*builder).build(track1_muon);
  reco::TransientTrack tt2_muon = (*builder).build(track2_muon);

  KinematicParticleFactoryFromTransientTrack pFactory;
  std::vector<RefCountedKinematicParticle> muonParticles;

  muonParticles.push_back(pFactory.particle(tt1_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));
  muonParticles.push_back(pFactory.particle(tt2_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));
  
  RefCountedKinematicParticle jpsi_part;
  RefCountedKinematicVertex jpsi_vertex;
  RefCountedKinematicTree jpTree;
  Bool_t jpsifit_flag;

  std::tie(jpsifit_flag, jpsi_part, jpsi_vertex, jpTree) = aux.KinematicFit(muonParticles, -1, -1);

  if(!jpsifit_flag) return false;

  nBranches_->cutflow->Fill(4);


  std::vector< RefCountedKinematicParticle > jpsi_children = jpTree->finalStateParticles();

  math::PtEtaPhiMLorentzVector mu1_fit = aux.daughter_p4(jpsi_children, 0);
  math::PtEtaPhiMLorentzVector mu2_fit = aux.daughter_p4(jpsi_children, 1);

  std::vector<pat::Muon> muoncollection_selected;
  muoncollection_selected.push_back(muoncollection[mcidx_mu1]);
  muoncollection_selected.push_back(muoncollection[mcidx_mu2]);

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
  TransverseImpactPointExtrapolator extrapolatort(fMagneticField);

  Float_t max_criteria = 999;
  reco::Vertex closestVertex; 
  int counter = 0;

  for( reco::VertexCollection::const_iterator vtx = vertices_->begin(); vtx != vertices_->end(); ++vtx){

    particle_cand cand = aux.calculateIPvariables(extrapolator, jpsi_part, jpsi_vertex, *vtx);
 
    if(TMath::Abs(cand.lip) < max_criteria){
      max_criteria = TMath::Abs(cand.lip);
      closestVertex = *vtx;
    }

    counter += 1;
  }


  

  if(!(muoncollection[mcidx_mu1].isSoftMuon(closestVertex) > 0.5 && muoncollection[mcidx_mu2].isSoftMuon(closestVertex) > 0.5)) return false;
  nBranches_->cutflow->Fill(5);

  if(verbose_) std::cout << "[JpsiKNtuplizerE] J/psi soft muon ID passed" << std::endl;



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
    Float_t precut_dz = pf.vz() - closestVertex.position().z();
    if(TMath::Abs(precut_dz) > c_dz) continue;


    reco::TransientTrack  _track = (*builder).build(pf.pseudoTrack());
  
    TrajectoryStateOnSurface _tsos_pf = extrapolator.extrapolate(_track.impactPointState(), jpsi_vertex->position());

    TrajectoryStateOnSurface _tsost_pf = extrapolatort.extrapolate(_track.impactPointState(), jpsi_vertex->position());
    
    std::pair<bool,Measurement1D> _cur3DIP_pf = aux.signedImpactParameter3D(_tsos_pf, jpsi_vertex, closestVertex);

    std::pair<bool,Measurement1D> _cur3DIPt_pf = aux.signedTransverseImpactParameter(_tsost_pf, jpsi_vertex, closestVertex);

    
    Float_t doca3d = _cur3DIP_pf.second.value();
    Float_t doca3de = _cur3DIP_pf.second.error();
    Float_t doca3ds = _cur3DIP_pf.second.significance();

    Float_t doca2d = _cur3DIPt_pf.second.value();
    Float_t doca2de = _cur3DIPt_pf.second.error();
    Float_t doca2ds = _cur3DIPt_pf.second.significance();

    Float_t near_dz = -999;
    for( reco::VertexCollection::const_iterator vtx = vertices_->begin(); vtx != vertices_->end(); ++vtx){
      
      if(pf.vertexRef()->z()==vtx->position().z()){
	near_dz = closestVertex.position().z() - vtx->position().z();
      }
    }

    Bool_t isAssociate = (bool)(pf.vertexRef()->z()==closestVertex.position().z());


    // Trigger matching 
    bool trigMatch = false;
    float trigMatch_dr = -1;

    for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
    
      obj.unpackPathNames(trigNames);
      obj.unpackFilterLabels(event, *HLTtriggers_);

      std::vector<std::string> pathNamesAll  = obj.pathNames(false);

      bool isPathExist = false;

      for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
	if(pathNamesAll[h]==finalTriggerName) isPathExist = true;
      }
      
      if(!isPathExist) continue;

      bool isFilterExist = false;
    
      for (unsigned hh = 0; hh < obj.filterLabels().size(); ++hh){
	
	if(obj.filterLabels()[hh].find(finalTriggerFilterObjName) != std::string::npos){
	  isFilterExist = true;
	}
      }
      
      if(!isFilterExist) continue;

      if(TMath::Abs(obj.pdgId()) != 321) continue;
      
      Float_t trigger_dR = reco::deltaR(obj.eta(), obj.phi(),
					pf.eta(), pf.phi());

      if(trigger_dR < 0.015 &&
	 obj.pt()/pf.pt() > 0.85 &&
	 obj.pt()/pf.pt() < 1.15){
	 
	trigMatch = true;
	trigMatch_dr = trigger_dR;
      }
    }


    float dr_jpsi = reco::deltaR(pf.eta(),
				 pf.phi(),
				 jpsi_part->currentState().globalMomentum().eta(),
				 jpsi_part->currentState().globalMomentum().phi());


    if(doca3d < -0.04 || doca3d > 0.06) continue;
    if(dr_jpsi > 1.) continue;


    attribute attr = {
      (Float_t) doca3d,
      (Float_t) doca3de,
      (Float_t) doca3ds,
      (Float_t) doca2d,
      (Float_t) doca2de,
      (Float_t) doca2ds,
      (Float_t) precut_dz,
      (Bool_t) isAssociate,
      (Int_t) pf.pvAssociationQuality(),
      (Float_t) pf.pt(),
      (Float_t) pf.eta(),
      (Float_t) pf.phi(),
      (Int_t) pf.charge(),
      (Float_t) pf.mass(),
      (Bool_t) false,
      (Int_t) -1,
      (Int_t) -1,
      (Bool_t) false,
      (Int_t) -1,
      (Int_t) -1,
      (Float_t) near_dz,
      (Bool_t) trigMatch,
      (Float_t) trigMatch_dr,
      (Float_t) dr_jpsi
    };

    
    reco::TransientTrack  tt_track = (*builder).build(pf.pseudoTrack());
    
    pfcand_struct _cand_ = {
      (Int_t)ii,
      pf,
      tt_track,
      attr
    };
    
    pfcands.push_back(_cand_);

  }

  // sorted by pt
  sort(pfcands.begin(), pfcands.end());



  nBranches_->cutflow->Fill(6);


  int numOfch = (int)pfcands.size();
  if(verbose_) std::cout << "[JpsiKNtuplizerE] Starts to build tau candidate out of " << numOfch << " pion candidates" << std::endl;

  Int_t nkaon = 0;

  for(int iii = 0; iii < (int)pfcands.size(); iii ++){
      
    pat::PackedCandidate pf = pfcands[iii].pfcand;
    reco::TransientTrack track = pfcands[iii].track;
    attribute attr = pfcands[iii].pfaux;

    std::vector<RefCountedKinematicParticle> allParticles;
      
    allParticles.push_back(pFactory.particle(track, aux.kaon_mass, chi, ndf, aux.kaon_sigma));
    allParticles.push_back(pFactory.particle(tt1_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));
    allParticles.push_back(pFactory.particle(tt2_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));

    RefCountedKinematicParticle b_part;
    RefCountedKinematicVertex b_vertex;
    RefCountedKinematicTree bTree;
    Bool_t bfit_flag;
    std::tie(bfit_flag, b_part, b_vertex, bTree) = aux.KinematicFit(allParticles, -1, -1);
    if(!bfit_flag) continue;

    if(TMath::Prob(b_vertex->chiSquared(), b_vertex->degreesOfFreedom()) <= 0.1) continue;
        

    particle_cand Bcand; 
    Bcand = aux.calculateIPvariables(extrapolator, b_part, b_vertex, closestVertex);
	  
    nBranches_->JpsiK_e_B_pt.push_back(b_part->currentState().globalMomentum().perp());
    nBranches_->JpsiK_e_B_eta.push_back(b_part->currentState().globalMomentum().eta());
    nBranches_->JpsiK_e_B_phi.push_back(b_part->currentState().globalMomentum().phi());
    nBranches_->JpsiK_e_B_mass.push_back(b_part->currentState().mass());
    nBranches_->JpsiK_e_B_vprob.push_back( TMath::Prob(b_vertex->chiSquared(), b_vertex->degreesOfFreedom()) ); //b_part->currentState().mass());
    nBranches_->JpsiK_e_B_lip.push_back(Bcand.lip);
    nBranches_->JpsiK_e_B_lips.push_back(Bcand.lips);
    nBranches_->JpsiK_e_B_pvip.push_back(Bcand.pvip);
    nBranches_->JpsiK_e_B_pvips.push_back(Bcand.pvips);
    nBranches_->JpsiK_e_B_fls3d.push_back(Bcand.fls3d);
    nBranches_->JpsiK_e_B_fl3d.push_back(Bcand.fl3d);
    nBranches_->JpsiK_e_B_alpha.push_back(Bcand.alpha);    
    nBranches_->JpsiK_e_B_vx.push_back(b_vertex->vertexState().position().x());
    nBranches_->JpsiK_e_B_vy.push_back(b_vertex->vertexState().position().y());
    nBranches_->JpsiK_e_B_vz.push_back(b_vertex->vertexState().position().z());

    nBranches_->JpsiK_e_pi_pt.push_back(pf.pt());
    nBranches_->JpsiK_e_pi_eta.push_back(pf.eta());
    nBranches_->JpsiK_e_pi_phi.push_back(pf.phi());
    nBranches_->JpsiK_e_pi_mass.push_back(pf.mass());
    nBranches_->JpsiK_e_pi_q.push_back(pf.charge());
    nBranches_->JpsiK_e_pi_doca3d.push_back(attr.doca3d);
    nBranches_->JpsiK_e_pi_doca3de.push_back(attr.doca3de);
    nBranches_->JpsiK_e_pi_doca2d.push_back(attr.doca2d);
    nBranches_->JpsiK_e_pi_doca2de.push_back(attr.doca2de);
    nBranches_->JpsiK_e_pi_dz.push_back(attr.dz);
    nBranches_->JpsiK_e_pi_near_dz.push_back(attr.near_dz);
    nBranches_->JpsiK_e_pi_isAssociate.push_back(attr.isAssociate);
    nBranches_->JpsiK_e_pi_pvAssociationQuality.push_back(attr.pvAssociationQuality);
    nBranches_->JpsiK_e_pi_trigMatch.push_back(attr.trigMatch);
    nBranches_->JpsiK_e_pi_trigMatch_dr.push_back(attr.trigMatch_dr);

    nkaon+=1;
  }
    

  nBranches_->JpsiK_e_mu1_pt = muoncollection[mcidx_mu1].pt();
  nBranches_->JpsiK_e_mu1_eta = muoncollection[mcidx_mu1].eta();
  nBranches_->JpsiK_e_mu1_phi = muoncollection[mcidx_mu1].phi();
  nBranches_->JpsiK_e_mu1_mass = muoncollection[mcidx_mu1].mass();
  nBranches_->JpsiK_e_mu1_q = muoncollection[mcidx_mu1].charge();
  nBranches_->JpsiK_e_mu1_isLoose = muoncollection[mcidx_mu1].isLooseMuon();
  nBranches_->JpsiK_e_mu1_isTight = muoncollection[mcidx_mu1].isTightMuon(closestVertex);
  nBranches_->JpsiK_e_mu1_isPF = muoncollection[mcidx_mu1].isPFMuon();
  nBranches_->JpsiK_e_mu1_isGlobal = muoncollection[mcidx_mu1].isGlobalMuon();
  nBranches_->JpsiK_e_mu1_isTracker = muoncollection[mcidx_mu1].isTrackerMuon();
  nBranches_->JpsiK_e_mu1_isSoft = muoncollection[mcidx_mu1].isSoftMuon(closestVertex);
  nBranches_->JpsiK_e_mu1_vx = muoncollection[mcidx_mu1].vx();
  nBranches_->JpsiK_e_mu1_vy = muoncollection[mcidx_mu1].vy();
  nBranches_->JpsiK_e_mu1_vz = muoncollection[mcidx_mu1].vz();
  nBranches_->JpsiK_e_mu1_dbiso = aux.MuonPFIso(muoncollection[mcidx_mu1]);
  
  nBranches_->JpsiK_e_mu2_pt = muoncollection[mcidx_mu2].pt();
  nBranches_->JpsiK_e_mu2_eta = muoncollection[mcidx_mu2].eta();
  nBranches_->JpsiK_e_mu2_phi = muoncollection[mcidx_mu2].phi();
  nBranches_->JpsiK_e_mu2_mass = muoncollection[mcidx_mu2].mass();
  nBranches_->JpsiK_e_mu2_q = muoncollection[mcidx_mu2].charge();
  nBranches_->JpsiK_e_mu2_isLoose = muoncollection[mcidx_mu2].isLooseMuon();
  nBranches_->JpsiK_e_mu2_isTight = muoncollection[mcidx_mu2].isTightMuon(closestVertex);
  nBranches_->JpsiK_e_mu2_isPF = muoncollection[mcidx_mu2].isPFMuon();
  nBranches_->JpsiK_e_mu2_isGlobal = muoncollection[mcidx_mu2].isGlobalMuon();
  nBranches_->JpsiK_e_mu2_isTracker = muoncollection[mcidx_mu2].isTrackerMuon();
  nBranches_->JpsiK_e_mu2_isSoft = muoncollection[mcidx_mu2].isSoftMuon(closestVertex);
  nBranches_->JpsiK_e_mu2_vx = muoncollection[mcidx_mu2].vx();
  nBranches_->JpsiK_e_mu2_vy = muoncollection[mcidx_mu2].vy();
  nBranches_->JpsiK_e_mu2_vz = muoncollection[mcidx_mu2].vz();
  nBranches_->JpsiK_e_mu2_dbiso = aux.MuonPFIso(muoncollection[mcidx_mu2]);

  nBranches_->JpsiK_e_PV_vx = vertices_->begin()->position().x();
  nBranches_->JpsiK_e_PV_vy = vertices_->begin()->position().y();
  nBranches_->JpsiK_e_PV_vz = vertices_->begin()->position().z();

  nBranches_->JpsiK_e_bbPV_vx = closestVertex.position().x();
  nBranches_->JpsiK_e_bbPV_vy = closestVertex.position().y();
  nBranches_->JpsiK_e_bbPV_vz = closestVertex.position().z();
  nBranches_->JpsiK_e_bbPV_chi2 = closestVertex.chi2();
  nBranches_->JpsiK_e_bbPV_ndof = closestVertex.ndof();
  nBranches_->JpsiK_e_bbPV_rho = closestVertex.position().Rho();

  nBranches_->JpsiK_e_Jpsi_pt = jpsi_part->currentState().globalMomentum().perp();
  nBranches_->JpsiK_e_Jpsi_eta = jpsi_part->currentState().globalMomentum().eta();
  nBranches_->JpsiK_e_Jpsi_phi = jpsi_part->currentState().globalMomentum().phi();
  nBranches_->JpsiK_e_Jpsi_mass = jpsi_part->currentState().mass();
  nBranches_->JpsiK_e_Jpsi_vprob = TMath::Prob(jpsi_part->chiSquared(), jpsi_part->degreesOfFreedom());

  particle_cand JPcand;
  JPcand = aux.calculateIPvariables(extrapolator, jpsi_part, jpsi_vertex, closestVertex);

  nBranches_->JpsiK_e_Jpsi_lip = JPcand.lip;
  nBranches_->JpsiK_e_Jpsi_lips = JPcand.lips;
  nBranches_->JpsiK_e_Jpsi_pvip = JPcand.pvip;
  nBranches_->JpsiK_e_Jpsi_pvips = JPcand.pvips;
  nBranches_->JpsiK_e_Jpsi_fl3d = JPcand.fl3d;
  nBranches_->JpsiK_e_Jpsi_fls3d = JPcand.fls3d;
  nBranches_->JpsiK_e_Jpsi_alpha = JPcand.alpha;
  nBranches_->JpsiK_e_Jpsi_maxdoca = aux.getMaxDoca(muonParticles);
  nBranches_->JpsiK_e_Jpsi_mindoca = aux.getMinDoca(muonParticles);
  nBranches_->JpsiK_e_Jpsi_vx = jpsi_vertex->vertexState().position().x();
  nBranches_->JpsiK_e_Jpsi_vy = jpsi_vertex->vertexState().position().y();
  nBranches_->JpsiK_e_Jpsi_vz = jpsi_vertex->vertexState().position().z();  



  /********************************************************************
   *
   * Step10: check gen-matching and fill them
   *
   ********************************************************************/

  nBranches_->JpsiK_e_nch = numOfch;
  nBranches_->JpsiK_e_nch_before = pfcands.size();
  nBranches_->JpsiK_e_nCandidates = nkaon;

  if(!runOnMC_) return true;
  nBranches_->genWeightBkgB=genWeightBkgB_;

  nBranches_->JpsiK_e_q2_gen = q2_gen.M2();
  nBranches_->JpsiK_e_nBc = nBc;
  nBranches_->JpsiK_e_B_pt_gen = pB_gen.Pt();
  nBranches_->JpsiK_e_B_eta_gen = pB_gen.Eta();
  nBranches_->JpsiK_e_B_phi_gen = pB_gen.Phi();
  nBranches_->JpsiK_e_B_mass_gen = pB_gen.M();
  

  bool flag_jpsi_match = false;
  if(gen_jpsi_mu.size()==2){

    Float_t _dR_11 = reco::deltaR(gen_jpsi_mu[0].Eta(), gen_jpsi_mu[0].Phi(), 
				  mu1_fit.eta(), mu1_fit.phi());
    
    Float_t _dR_22 = reco::deltaR(gen_jpsi_mu[1].Eta(), gen_jpsi_mu[1].Phi(), 
				  mu2_fit.eta(), mu2_fit.phi());

    Float_t _dR_21 = reco::deltaR(gen_jpsi_mu[1].Eta(), gen_jpsi_mu[1].Phi(), 
				  mu1_fit.eta(), mu1_fit.phi());
    
    Float_t _dR_12 = reco::deltaR(gen_jpsi_mu[0].Eta(), gen_jpsi_mu[0].Phi(), 
				  mu2_fit.eta(), mu2_fit.phi());

    if(_dR_11 < 0.1 && _dR_22 < 0.1) flag_jpsi_match = true;
    else if(_dR_21 < 0.1 && _dR_12 < 0.1) flag_jpsi_match = true;
  }
  
  nBranches_->JpsiK_e_genPV_vx = genvertex.x();
  nBranches_->JpsiK_e_genPV_vy = genvertex.y();
  nBranches_->JpsiK_e_genPV_vz = genvertex.z();

  nBranches_->JpsiK_e_genSV_vx = genvertex_sv.x();
  nBranches_->JpsiK_e_genSV_vy = genvertex_sv.y();
  nBranches_->JpsiK_e_genSV_vz = genvertex_sv.z();

  nBranches_->JpsiK_e_ngenmuons = gen_nr_mu.size() + gen_jpsi_mu.size();
  nBranches_->JpsiK_e_isgenmatched = (int)flag_jpsi_match;

  
  return true;

}


const reco::Candidate*  JpsiKNtuplizerE::checkMom(const reco::Candidate * candMom){
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
