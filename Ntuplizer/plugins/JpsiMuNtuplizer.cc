#include "../interface/JpsiMuNtuplizer.h"


//===================================================================================================================
JpsiMuNtuplizer::JpsiMuNtuplizer( edm::EDGetTokenT<pat::MuonCollection>    muonToken   ,
				  edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
				  edm::EDGetTokenT<reco::BeamSpot> bsToken,
				  edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken,
				  edm::EDGetTokenT<edm::TriggerResults> triggertoken,
				  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobject,
				  edm::EDGetTokenT<reco::GenParticleCollection> genptoken,
				  std::map< std::string, bool >& runFlags,
				  NtupleBranches* nBranches )
  : CandidateNtuplizer ( nBranches )
  , muonToken_	        ( muonToken )
  , verticeToken_          ( verticeToken )
  , bsToken_          ( bsToken )
  , packedpfcandidatesToken_(packedpfcandidatesToken) 
  , HLTtriggersToken_	( triggertoken )
  , triggerObjects_	( triggerobject )
  , genParticlesToken_( genptoken )
  , runOnMC_   (runFlags["runOnMC"])
   
{
}

//===================================================================================================================
JpsiMuNtuplizer::~JpsiMuNtuplizer( void )
{

}


TVector3 JpsiMuNtuplizer::getVertex(const reco::GenParticle& part){
  return TVector3(part.vx(),part.vy(),part.vz());
}

float JpsiMuNtuplizer::MuonPFIso(pat::Muon muon, bool highpt){

  float sumChargedHadronPt = muon.pfIsolationR04().sumChargedHadronPt;
  float sumNeutralHadronEt = muon.pfIsolationR04().sumNeutralHadronEt;
  float sumPhotonEt = muon.pfIsolationR04().sumPhotonEt;
  float sumPUPt = muon.pfIsolationR04().sumPUPt;
  float iso = (sumChargedHadronPt + std::max( 0. ,sumNeutralHadronEt + sumPhotonEt - 0.5 * sumPUPt));// / muon.pt()
 
  return iso;
}


double JpsiMuNtuplizer::isoTrack(double docaCut, double r, double pmin) {

///  const double pCut(pmin), coneSize(r);
///  const bool verbose(false);

///  double iso(-1.), p(0.), sumP(0.);
///  TSimpleTrack *ps;
///  vector<int> cIdx, pIdx;
///  int pvIdx = pC->fPvIdx;
///
///  double trackP = pTrack->fPlab.Mag();
///
///  getSigTracks(cIdx, pC);
///
///  // -- look at all tracks that are associated to the same vertex
///  for (int i = 0; i < fpEvt->nSimpleTracks(); ++i) {
///    ps = fpEvt->getSimpleTrack(i);
///    if (verbose) {
///      cout << "   track " << i
///           << " with p = " << ps->getP().Mag()
///           << " eta = " << ps->getP().Eta()
///           << " pointing at PV " << ps->getPvIndex();
///    }
///
///
///    if (ps->getPvIndex() != pvIdx) {
///      if (verbose) cout << " skipped because of PV index mismatch" << endl;          //FIXME
///      continue;
///    }
///
///    // -- despite the name use momentum (FIXME?)
///    p = ps->getP().Mag();
///    if (p < pCut) {
///      if (verbose) cout << " skipped because of p = " << p << endl;
///      continue;
///    }
///    if (cIdx.end() != find(cIdx.begin(), cIdx.end(), i))  {
///      if (verbose) cout << " skipped because it is a sig track " << endl;
///      continue;
///    }
///    if (ps->getP().DeltaR(pTrack->fPlab) < coneSize) {
///      pIdx.push_back(i);
///      sumP += p;
///      if (verbose) cout << " USED THIS ONE" << endl;
///    }
///    else {
///      if (verbose) cout << " skipped because of deltaR = " << ps->getP().DeltaR(pC->fPlab) << endl;
///    }
///  }
///
///  // -- Now consider the DOCA tracks
///  int nsize = pC->fNstTracks.size();
///  if (nsize>0) {
///    for(int i = 0; i<nsize; ++i) {
///      int trkId = pC->fNstTracks[i].first;
///      double doca = pC->fNstTracks[i].second.first;
///      // double docaE = pC->fNstTracks[i].second.second;
///
///      if(doca > docaCut) continue; // check the doca cut
///
///      ps = fpEvt->getSimpleTrack(trkId);
///      p = ps->getP().Mag();
///
///
///      if ((ps->getPvIndex() > -1) && (ps->getPvIndex() != pvIdx)) {
///        if (verbose) cout << " doca track " << trkId  << " doca = " << doca  << " p = " << p
///                          << " skipped because it is from a different PV " << ps->getPvIndex() <<endl;
///        continue;
///      }
///
///      if (p < pCut) {
///        if (verbose) cout << " doca track " << trkId  << " doca = " << doca << " p = " << p
///                          << " skipped because of p = " << p << endl;
///        continue;
///      }
///
///      if (ps->getP().DeltaR(pC->fPlab) > coneSize) {
///        if (verbose) cout << " doca track " << trkId << " skipped because of deltaR = " << ps->getP().DeltaR(pC->fPlab) << endl;
///        continue;
///      }
///
///      // -- Skip tracks already included above
///      if (pIdx.end() != find(pIdx.begin(), pIdx.end(), trkId))  continue;
///      if (cIdx.end() != find(cIdx.begin(), cIdx.end(), trkId))  continue;
///
///      //      cout << "doca trk " << trkId << " doca = " << doca << endl;
///
///      sumP += p;
///      if (verbose) cout << " doca track " << trkId << " included "<<doca<<" "<<p<<endl;
///
///    } // for loop over tracks
///  } // end if

  
  //iso = trackP/(trackP + sumP);
  double iso = 1.;

  return iso;
}


Float_t JpsiMuNtuplizer::getMaxDoca(std::vector<RefCountedKinematicParticle> &kinParticles){

  double maxDoca = -1.0;

  TwoTrackMinimumDistance md;
  std::vector<RefCountedKinematicParticle>::iterator in_it, out_it;

  for (out_it = kinParticles.begin(); out_it != kinParticles.end(); ++out_it) {
    for (in_it = out_it + 1; in_it != kinParticles.end(); ++in_it) {
      md.calculate((*out_it)->currentState().freeTrajectoryState(),(*in_it)->currentState().freeTrajectoryState());
      if (md.distance() > maxDoca)
	maxDoca = md.distance();
    }
  }

  return maxDoca;
}



Float_t JpsiMuNtuplizer::getMinDoca(std::vector<RefCountedKinematicParticle> &kinParticles) {

  double minDoca = 99999.9;

  TwoTrackMinimumDistance md;
  unsigned j,k,n;

  n = kinParticles.size();
  for (j = 0; j < n; j++) {
    for (k = j+1; k < n; k++) {
      md.calculate(kinParticles[j]->currentState().freeTrajectoryState(),kinParticles[k]->currentState().freeTrajectoryState());
      if (md.distance() < minDoca)
	minDoca = md.distance();
    }
  }

  return minDoca;
}




std::tuple<Float_t, TransientVertex> JpsiMuNtuplizer::vertexProb( const std::vector<reco::TransientTrack>& tracks){

  Float_t vprob = -1;
  
  KalmanVertexFitter kalman_fitter;
  TransientVertex vertex;

  try{
    vertex = kalman_fitter.vertex(tracks);
  }catch(std::exception e){
    std::cout << "No vertex found ... return" << std::endl;
    return std::forward_as_tuple(-9, vertex);
  }

  if(vertex.isValid()){

    vprob =  TMath::Prob(vertex.totalChiSquared(), vertex.degreesOfFreedom());

    //    vx = vertex.position().x();
    //    vy = vertex.position().y();
    //    vz = vertex.position().z();
    
    return std::forward_as_tuple(vprob, vertex);

  }else{

    return std::forward_as_tuple(-9, vertex);

  }
}




particle_cand JpsiMuNtuplizer::calculateIPvariables(
						    AnalyticalImpactPointExtrapolator extrapolator,
						    RefCountedKinematicParticle particle,
						    RefCountedKinematicVertex vertex,
						    reco::Vertex wrtVertex
						    ){

  TrajectoryStateOnSurface tsos = extrapolator.extrapolate(particle->currentState().freeTrajectoryState(),
							   RecoVertex::convertPos(wrtVertex.position()));


  VertexDistance3D a3d;  

  std::pair<bool,Measurement1D> currentIp = IPTools::signedDecayLength3D(tsos, GlobalVector(0,0,1), wrtVertex);
  std::pair<bool,Measurement1D> cur3DIP = IPTools::absoluteImpactParameter(tsos, wrtVertex, a3d);
  
  // flight length
  Float_t fl3d = a3d.distance(wrtVertex, vertex->vertexState()).value();
  Float_t fl3de = a3d.distance(wrtVertex, vertex->vertexState()).error();
  Float_t fls3d = -1;

  if(fl3de!=0) fls3d = fl3d/fl3de;

  // longitudinal impact parameters
  Float_t lip = currentIp.second.value();
  Float_t lipe = currentIp.second.error();
  Float_t lips = -1;
  
  if(lipe!=0) lips = lip/lipe;

  // impact parameter to the PV
  Float_t pvip = cur3DIP.second.value();
  Float_t pvipe = cur3DIP.second.error();
  Float_t pvips = -1;
  
  if(pvipe!=0) pvips = pvip/pvipe;

  // opening angle
  TVector3 plab = TVector3(particle->currentState().globalMomentum().x(),
			   particle->currentState().globalMomentum().y(),
			   particle->currentState().globalMomentum().z());

  const TVector3 tv3diff = TVector3(vertex->vertexState().position().x() - wrtVertex.position().x(),
				    vertex->vertexState().position().y() - wrtVertex.position().y(),
				    vertex->vertexState().position().z() - wrtVertex.position().z()
				    );

  Float_t alpha = -1;

  if(plab.Mag() != 0. && tv3diff.Mag()!=0){
    alpha = plab.Dot(tv3diff) / (plab.Mag() * tv3diff.Mag());
  }

  particle_cand cand = {
    lip,
    lips,
    pvip, 
    pvips,
    fl3d,
    fls3d,
    alpha
  };


  return cand;
}





void JpsiMuNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  
  std::cout << "---------------- event, run, lumi = " << event.id().event() << " " << event.id().run() << " " << event.id().luminosityBlock() << "----------------" << std::endl;
  
  /********************************************************************
   *
   * Step1: check if the J/psi trigger is fired.
   * Namely, HLT_DoubleMu4_JpsiTrk_Displaced_v
   *
   ********************************************************************/

  event.getByToken(HLTtriggersToken_, HLTtriggers_);

  bool isTriggered = false;
  const edm::TriggerNames& trigNames = event.triggerNames(*HLTtriggers_);
  std::string finalTriggerName = "";

  for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {

    if(trigNames.triggerName(i).find("HLT_DoubleMu4_JpsiTrk_Displaced_v")!= std::string::npos){
      if(HLTtriggers_->accept(i)){
	isTriggered = true;
	finalTriggerName = trigNames.triggerName(i);
      }
    }
  }

  if(!isTriggered) return;

  //  std::cout << "finalTriggerName = "  << finalTriggerName << std::endl;



  /********************************************************************
   *
   * Step2: pre-select muons for building J/psi candidates ... 
   * For muons, no requirement applied
   *
   ********************************************************************/

  event.getByToken(verticeToken_   , vertices_     );
  event.getByToken(bsToken_   , bs_     );
  event.getByToken(muonToken_	, muons_    );
  event.getByToken(triggerObjects_  , triggerObjects);

  std::vector<pat::Muon> muoncollection;
  std::vector<int> muoncollection_id;
  muoncollection.clear();
  muoncollection_id.clear();

  for(size_t imuon = 0; imuon < muons_->size(); ++ imuon){

    const pat::Muon & muon = (*muons_)[imuon];

    if(muon.pt() < 4) continue;
    if(TMath::Abs(muon.eta()) > 2.4) continue;
    if(!(muon.track().isNonnull())) continue;

    //    bool isSoft = muon.isSoftMuon(*firstGoodVertex);
    //    bool isGlobal = muon.isGlobalMuon();
    //    bool isTracker = muon.isTrackerMuon();
    //    bool isLoose = muon.isLooseMuon();
    //    bool isTight =  muon.isTightMuon(*firstGoodVertex);
    //    bool isPF = muon.isPFMuon();
    //    if(!(isSoft && isGlobal)) continue;
    //    if(TMath::Abs(muon.muonBestTrack()->dz(firstGoodVertex->position())) > 0.5) continue;

    
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
	
	if(obj.filterLabels()[hh].find("hltJpsiTkVertexFilter") != std::string::npos){
	  isFilterExist = true;
	}
      }
      
      if(!isFilterExist) continue;
      
      Float_t trigger_dR = reco::deltaR(obj.eta(), obj.phi(),
					muon.eta(), muon.phi());
      
      if(trigger_dR < 0.1) trigMatch = true;
    }

    if(!trigMatch) continue;

    muoncollection.push_back(muon);
    muoncollection_id.push_back(imuon);
  }

  if(!( muoncollection.size() >= 2)) return;




  /********************************************************************
   *
   * Step3: building J/psi candidates 
   *
   ********************************************************************/

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

  Float_t jpsi_max_pt = -1;
  unsigned int idx_mu1 = -1;
  unsigned int idx_mu2 = -1;
  unsigned int mcidx_mu1 = -1;
  unsigned int mcidx_mu2 = -1;
  TLorentzVector jpsi_tlv_highest;
  Float_t jpsi_vprob_highest = -9;
  TransientVertex jpsi_vertex_highest;

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

      
      std::vector<reco::TransientTrack> transient_tracks_dimuon;
      
      transient_tracks_dimuon.push_back((*builder).build(muoncollection[imu].muonBestTrack()));
      transient_tracks_dimuon.push_back((*builder).build(muoncollection[jmu].muonBestTrack()));
      
      Float_t vprob_jpsi = -9;
      TransientVertex vertex_jpsi;
      std::tie(vprob_jpsi, vertex_jpsi) = vertexProb(transient_tracks_dimuon);
      //      if(!(vprob_jpsi > 0)) continue;

      if(jpsi_max_pt < jpsi_pt){
	jpsi_max_pt = jpsi_pt;
	idx_mu1 = muoncollection_id[imu];
	idx_mu2 = muoncollection_id[jmu];
	mcidx_mu1 = imu;
	mcidx_mu2 = jmu;
	jpsi_tlv_highest = tlv_jpsi;
	jpsi_vprob_highest = vprob_jpsi;
	jpsi_vertex_highest = vertex_jpsi;
      }
    }
  }

  if(jpsi_max_pt == -1) return;




  /********************************************************************
   *
   * Step4: Kinematic fit for the J/psi candidate
   *
   ********************************************************************/

  const reco::TrackRef track1_muon = muoncollection[mcidx_mu1].muonBestTrack();
  const reco::TrackRef track2_muon = muoncollection[mcidx_mu2].muonBestTrack();
  reco::TransientTrack tt1_muon = (*builder).build(track1_muon);
  reco::TransientTrack tt2_muon = (*builder).build(track2_muon);

  KinematicParticleFactoryFromTransientTrack pFactory;
  std::vector<RefCountedKinematicParticle> muonParticles;

  muonParticles.push_back(pFactory.particle(tt1_muon, muon_mass, chi, ndf, muon_sigma));
  muonParticles.push_back(pFactory.particle(tt2_muon, muon_mass, chi, ndf, muon_sigma));
  
  //creating the vertex fitter
  KinematicParticleVertexFitter kpvFitter;

  //reconstructing a J/Psi decay
  RefCountedKinematicTree jpTree = kpvFitter.fit(muonParticles);

  if(jpTree->isEmpty() || !jpTree->isValid() || !jpTree->isConsistent()) return;

  //creating the particle fitter
  KinematicParticleFitter csFitter;

  // creating the constraint
  KinematicConstraint* jpsi_constraint = new MassKinematicConstraint(jpsi_mass, jp_m_sigma);

  //the constrained fit
  jpTree = csFitter.fit(jpsi_constraint, jpTree);

  //getting the J/Psi KinematicParticle
  jpTree->movePointerToTheTop();
  RefCountedKinematicParticle jpsi_part = jpTree->currentParticle();
  if(!jpsi_part->currentState().isValid()) return;

  RefCountedKinematicVertex jpsi_vertex = jpTree->currentDecayVertex();
  if(!jpsi_vertex->vertexIsValid()) return; 

  if(TMath::Prob(jpsi_vertex->chiSquared(), jpsi_vertex->degreesOfFreedom()) <=0) return;


  /********************************************************************
   *
   * Step5: determine bbbar-PV, extrapolated back from the J/psi candidate
   * There are several ways to do this ...
   * 1) using minimum lip
   * 2) using minimum pvip
   * 3) using PV (first content of the PV collection)
   *
   ********************************************************************/

  // define extrapolator
  edm::ESHandle<MagneticField> fieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(fieldHandle);
  fMagneticField = fieldHandle.product();

  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);

  Float_t max_criteria = 999;
  reco::Vertex closestVertex; 

  for( reco::VertexCollection::const_iterator vtx = vertices_->begin(); vtx != vertices_->end(); ++vtx){
    
    //    bool isFake = (vtx->chi2()==0 && vtx->ndof()==0);
    
//    if(
//       !(!isFake && vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0)
//       ) continue;
    

    particle_cand cand = calculateIPvariables(extrapolator, jpsi_part, jpsi_vertex, *vtx);

    
    if(TMath::Abs(cand.lip) < max_criteria){
      //    if(TMath::Abs(cand.pvip) < max_criteria){
      max_criteria = TMath::Abs(cand.lip);
      //      max_criteria = TMath::Abs(cand.pvip);
      closestVertex = *vtx;
    }
  }


  
  /********************************************************************
   *
   * Step6: 3rd muon selection
   *
   ********************************************************************/

  //  event.getByToken( packedpfcandidatesToken_               , packedpfcandidates_      ); 

  pat::Muon mu3;
  Float_t max_pt3 = -1;
  
  for(size_t imuon = 0; imuon < muons_->size(); ++ imuon){

    const pat::Muon & muon = (*muons_)[imuon];

    if(muon.pt() < 4) continue;
    if(TMath::Abs(muon.eta()) > 2.4) continue;
    if(!(muon.track().isNonnull())) continue;
    if(imuon==idx_mu1 || imuon==idx_mu2) continue;

    if(muon.pt() > max_pt3){
      max_pt3 = muon.pt();
      mu3 = muon;
    }
  }

  if(max_pt3==-1) return;


  std::vector<RefCountedKinematicParticle> allParticles;

  const reco::TrackRef track3_muon = mu3.muonBestTrack();
  reco::TransientTrack tt3_muon = (*builder).build(track3_muon);

  allParticles.push_back(pFactory.particle(tt3_muon, muon_mass, chi, ndf, muon_sigma));
  allParticles.push_back(jpsi_part);

  RefCountedKinematicTree bcTree = kpvFitter.fit(allParticles);

  if(bcTree->isEmpty() || !bcTree->isValid() || !bcTree->isConsistent()) return;

  RefCountedKinematicParticle bc_part = bcTree->currentParticle();
  if(!bc_part->currentState().isValid()) return;

  RefCountedKinematicVertex bc_vertex = bcTree->currentDecayVertex();
  if(!bc_vertex->vertexIsValid()) return; 

  particle_cand JPcand = calculateIPvariables(extrapolator, jpsi_part, jpsi_vertex, closestVertex);
  particle_cand Bcand = calculateIPvariables(extrapolator, bc_part, bc_vertex, closestVertex);


  std::cout << "J/psi candidate (lip, lips, pvip, pvips, fl3d, fls3d, alpha) = " << JPcand.lip << " " << JPcand.lips << " " << JPcand.pvip << " " << JPcand.pvips << " " << JPcand.fl3d << " " << JPcand.fls3d << " " << JPcand.alpha << std::endl;

  std::cout << "B candidate (lip, lips, pvip, pvips, fl3d, fls3d, alpha) = " << Bcand.lip << " " << Bcand.lips << " " << Bcand.pvip << " " << Bcand.pvips << " " << Bcand.fl3d << " " << Bcand.fls3d << " " << Bcand.alpha << std::endl;



  // for unfit variables from here 
  TLorentzVector tlv_mu3;
  tlv_mu3.SetPtEtaPhiM(mu3.pt(), mu3.eta(), mu3.phi(), mu3.mass());

  TLorentzVector tlv_B = jpsi_tlv_highest + tlv_mu3;

  std::vector<reco::TransientTrack> transient_tracks_trimuon;
  transient_tracks_trimuon.push_back(tt1_muon);
  transient_tracks_trimuon.push_back(tt2_muon);
  transient_tracks_trimuon.push_back(tt3_muon);

  Float_t vprob_bc = -9;
  TransientVertex vertex_bc;
  std::tie(vprob_bc, vertex_bc) = vertexProb(transient_tracks_trimuon);
  // to here


//	Float_t min_dR_pf = 999.;
//	Float_t iso_pt03 =0.;
//	Float_t iso_pt04 =0.;
//	Float_t iso_pt05 =0.;
//	Float_t iso_pt06 =0.;
//	Float_t iso_pt07 =0.;
//
//	//isolation part
//	for(int kkk = 0; kkk < numOfch; kkk ++){
//	  
//	  pat::PackedCandidate _pf = pfcollection[kkk];
//	  if(TMath::Abs(_pf.pdgId())!=211) continue;
//	  TLorentzVector tlv_iso;
//	  tlv_iso.SetPtEtaPhiM(_pf.pt(), _pf.eta(), _pf.phi(), _pf.mass());
//	  Float_t dR_iso = tlv_muon3.DeltaR(tlv_iso);
//	  if(min_dR_pf>dR_iso) min_dR_pf=dR_iso;
//	  if(dR_iso>0.7) continue; 
//	  if(dR_iso<0.7)  iso_pt07 += tlv_iso.Pt();
//	  if(dR_iso<0.6)  iso_pt06 += tlv_iso.Pt();
//	  if(dR_iso<0.5)  iso_pt05 += tlv_iso.Pt();
//	  if(dR_iso<0.4)  iso_pt04 += tlv_iso.Pt();
//	  if(dR_iso<0.3)  iso_pt03 += tlv_iso.Pt();
//	  //std::cout<<dR_iso<<"<-dr, pt03->"<<iso_pt07<<std::endl;
//	}
	//isolation part
//	nBranches_->Jpsi_mu3_isopt03.push_back(iso_pt03);
//	nBranches_->Jpsi_mu3_isopt04.push_back(iso_pt04);
//	nBranches_->Jpsi_mu3_isopt05.push_back(iso_pt05);
//	nBranches_->Jpsi_mu3_isopt06.push_back(iso_pt06);
//	nBranches_->Jpsi_mu3_isopt07.push_back(iso_pt07);
//	nBranches_->Jpsi_dr_mu3pf.push_back(min_dR_pf);  
	//flight part
	//jpsi part


  /********************************************************************
   *
   * Step7: Filling branches ...
   *
   ********************************************************************/

  bool isMC = runOnMC_;
  
  TVector3 genvertex(-9.,-9.,-9.);
  
  std::vector<const reco::Candidate*> gen_nr_mu;
  std::vector<const reco::Candidate*> gen_jpsi_mu;
  

  if(isMC){
    event.getByToken(genParticlesToken_ , genParticles_); 
    
    for( unsigned p=0; p < genParticles_->size(); ++p){

      //      std::cout << "gen: " << (*genParticles_)[p].pdgId() << " " << (*genParticles_)[p].status() << std::endl;

      // Bc daughters loop
      if(TMath::Abs((*genParticles_)[p].pdgId())==541 && (*genParticles_)[p].status()==2){

	// retrieve production vertex
	genvertex = getVertex((*genParticles_)[p]);

	for(int idd = 0; idd < (int)(*genParticles_)[p].numberOfDaughters(); idd++){
	  Int_t dpid = (*genParticles_)[p].daughter(idd)->pdgId();
	  //	  std::cout << "\t -> " <<  << " " << (*genParticles_)[p].daughter(idd)->status()<< std::endl;
	  if(TMath::Abs(dpid)==13) gen_nr_mu.push_back((*genParticles_)[p].daughter(idd));
	}
      }

      // J/psi loop
      if(TMath::Abs((*genParticles_)[p].pdgId())==443 && 
	 (*genParticles_)[p].status()==2 && 
	 TMath::Abs((*genParticles_)[p].mother(0)->pdgId())==541){

	//	std::cout << "nMon = " << (*genParticles_)[p].numberOfMothers() << std::endl;
	//	std::cout << "mother pdgId = " <<  << std::endl;

	for(int idd = 0; idd < (int)(*genParticles_)[p].numberOfDaughters(); idd++){
	  Int_t dpid = (*genParticles_)[p].daughter(idd)->pdgId();
	  if(TMath::Abs(dpid)==13) gen_jpsi_mu.push_back((*genParticles_)[p].daughter(idd));
	  //	  std::cout << "\t -> " << (*genParticles_)[p].daughter(idd)->pdgId() << " " << (*genParticles_)[p].daughter(idd)->status()<< std::endl;
	}
      }

      
    }
  }



  nBranches_->Jpsi_mu1_isLoose.push_back(muoncollection[mcidx_mu1].isLooseMuon());
  nBranches_->Jpsi_mu1_isTight.push_back(muoncollection[mcidx_mu1].isTightMuon(closestVertex));
  nBranches_->Jpsi_mu1_isPF.push_back(muoncollection[mcidx_mu1].isPFMuon());
  nBranches_->Jpsi_mu1_isGlobal.push_back(muoncollection[mcidx_mu1].isGlobalMuon());
  nBranches_->Jpsi_mu1_isTracker.push_back(muoncollection[mcidx_mu1].isTrackerMuon());
  nBranches_->Jpsi_mu1_isSoft.push_back(muoncollection[mcidx_mu1].isSoftMuon(closestVertex));
  nBranches_->Jpsi_mu1_pt.push_back(muoncollection[mcidx_mu1].pt());
  nBranches_->Jpsi_mu1_eta.push_back(muoncollection[mcidx_mu1].eta());
  nBranches_->Jpsi_mu1_phi.push_back(muoncollection[mcidx_mu1].phi());
  nBranches_->Jpsi_mu1_q.push_back(muoncollection[mcidx_mu1].charge());
  nBranches_->Jpsi_mu1_vx.push_back(muoncollection[mcidx_mu1].vx());
  nBranches_->Jpsi_mu1_vy.push_back(muoncollection[mcidx_mu1].vy());
  nBranches_->Jpsi_mu1_vz.push_back(muoncollection[mcidx_mu1].vz());
  
  nBranches_->Jpsi_mu2_isLoose.push_back(muoncollection[mcidx_mu2].isLooseMuon());
  nBranches_->Jpsi_mu2_isTight.push_back(muoncollection[mcidx_mu2].isTightMuon(closestVertex));
  nBranches_->Jpsi_mu2_isPF.push_back(muoncollection[mcidx_mu2].isPFMuon());
  nBranches_->Jpsi_mu2_isGlobal.push_back(muoncollection[mcidx_mu2].isGlobalMuon());
  nBranches_->Jpsi_mu2_isTracker.push_back(muoncollection[mcidx_mu2].isTrackerMuon());
  nBranches_->Jpsi_mu2_isSoft.push_back(muoncollection[mcidx_mu2].isSoftMuon(closestVertex));
  nBranches_->Jpsi_mu2_pt.push_back(muoncollection[mcidx_mu2].pt());
  nBranches_->Jpsi_mu2_eta.push_back(muoncollection[mcidx_mu2].eta());
  nBranches_->Jpsi_mu2_phi.push_back(muoncollection[mcidx_mu2].phi());
  nBranches_->Jpsi_mu2_q.push_back(muoncollection[mcidx_mu2].charge());
  nBranches_->Jpsi_mu2_vx.push_back(muoncollection[mcidx_mu2].vx());
  nBranches_->Jpsi_mu2_vy.push_back(muoncollection[mcidx_mu2].vy());
  nBranches_->Jpsi_mu2_vz.push_back(muoncollection[mcidx_mu2].vz());

  nBranches_->Jpsi_mu3_isLoose.push_back(mu3.isLooseMuon());
  nBranches_->Jpsi_mu3_isTight.push_back(mu3.isTightMuon(closestVertex));
  nBranches_->Jpsi_mu3_isPF.push_back(mu3.isPFMuon());
  nBranches_->Jpsi_mu3_isGlobal.push_back(mu3.isGlobalMuon());
  nBranches_->Jpsi_mu3_isTracker.push_back(mu3.isTrackerMuon());
  nBranches_->Jpsi_mu3_isSoft.push_back(mu3.isSoftMuon(closestVertex));
  nBranches_->Jpsi_mu3_pt.push_back(mu3.pt());
  nBranches_->Jpsi_mu3_eta.push_back(mu3.eta());
  nBranches_->Jpsi_mu3_phi.push_back(mu3.phi());
  nBranches_->Jpsi_mu3_q.push_back(mu3.charge());
  nBranches_->Jpsi_mu3_vx.push_back(mu3.vx());
  nBranches_->Jpsi_mu3_vy.push_back(mu3.vy());
  nBranches_->Jpsi_mu3_vz.push_back(mu3.vz());

  nBranches_->Jpsi_PV_vx.push_back(vertices_->begin()->position().x());
  nBranches_->Jpsi_PV_vy.push_back(vertices_->begin()->position().y());
  nBranches_->Jpsi_PV_vz.push_back(vertices_->begin()->position().z());

  nBranches_->Jpsi_bbPV_vx.push_back(closestVertex.position().x());
  nBranches_->Jpsi_bbPV_vy.push_back(closestVertex.position().y());
  nBranches_->Jpsi_bbPV_vz.push_back(closestVertex.position().z());

  // -9 if there is no Bc found 
  nBranches_->Jpsi_genPV_vx.push_back(genvertex.x());
  nBranches_->Jpsi_genPV_vy.push_back(genvertex.y());
  nBranches_->Jpsi_genPV_vz.push_back(genvertex.z());

  nBranches_->Jpsi_pt.push_back(jpsi_part->currentState().globalMomentum().perp());
  nBranches_->Jpsi_eta.push_back(jpsi_part->currentState().globalMomentum().eta());
  nBranches_->Jpsi_phi.push_back(jpsi_part->currentState().globalMomentum().phi());
  nBranches_->Jpsi_mass.push_back(jpsi_part->currentState().mass());
  nBranches_->Jpsi_vprob.push_back(TMath::Prob(jpsi_part->chiSquared(), jpsi_part->degreesOfFreedom()));
  nBranches_->Jpsi_lip.push_back(JPcand.lip);
  nBranches_->Jpsi_lips.push_back(JPcand.lips);
  nBranches_->Jpsi_pvip.push_back(JPcand.pvip);
  nBranches_->Jpsi_pvips.push_back(JPcand.pvips);
  nBranches_->Jpsi_fl3d.push_back(JPcand.fl3d);
  nBranches_->Jpsi_fls3d.push_back(JPcand.fls3d);
  nBranches_->Jpsi_alpha.push_back(JPcand.alpha);
  nBranches_->Jpsi_maxdoca.push_back(getMaxDoca(muonParticles));
  nBranches_->Jpsi_mindoca.push_back(getMinDoca(muonParticles));
  nBranches_->Jpsi_vx.push_back(jpsi_vertex->vertexState().position().x());
  nBranches_->Jpsi_vy.push_back(jpsi_vertex->vertexState().position().y());
  nBranches_->Jpsi_vz.push_back(jpsi_vertex->vertexState().position().z());  
  nBranches_->Jpsi_unfitpt.push_back(jpsi_tlv_highest.Pt());
  nBranches_->Jpsi_unfitmass.push_back(jpsi_tlv_highest.M());
  nBranches_->Jpsi_unfitvprob.push_back(jpsi_vprob_highest);

  if(jpsi_vprob_highest!=-9){
    nBranches_->Jpsi_unfit_vx.push_back(jpsi_vertex_highest.position().x());
    nBranches_->Jpsi_unfit_vy.push_back(jpsi_vertex_highest.position().y());
    nBranches_->Jpsi_unfit_vz.push_back(jpsi_vertex_highest.position().z());
  }


  nBranches_->Jpsi_trimu_pt.push_back(bc_part->currentState().globalMomentum().perp());
  nBranches_->Jpsi_trimu_eta.push_back(bc_part->currentState().globalMomentum().eta());
  nBranches_->Jpsi_trimu_phi.push_back(bc_part->currentState().globalMomentum().phi());
  nBranches_->Jpsi_trimu_mass.push_back(bc_part->currentState().mass());
  nBranches_->Jpsi_trimu_vprob.push_back(TMath::Prob(bc_part->chiSquared(), bc_part->degreesOfFreedom()));
  nBranches_->Jpsi_trimu_lip.push_back(Bcand.lip);
  nBranches_->Jpsi_trimu_lips.push_back(Bcand.lips);
  nBranches_->Jpsi_trimu_pvip.push_back(Bcand.pvip);
  nBranches_->Jpsi_trimu_pvips.push_back(Bcand.pvips);
  nBranches_->Jpsi_trimu_fl3d.push_back(Bcand.fl3d);
  nBranches_->Jpsi_trimu_fls3d.push_back(Bcand.fls3d);
  nBranches_->Jpsi_trimu_alpha.push_back(Bcand.alpha);

  std::vector<RefCountedKinematicParticle> allParticles4doc;

  allParticles4doc.push_back(pFactory.particle(tt1_muon, muon_mass, chi, ndf, muon_sigma));
  allParticles4doc.push_back(pFactory.particle(tt2_muon, muon_mass, chi, ndf, muon_sigma));
  allParticles4doc.push_back(pFactory.particle(tt3_muon, muon_mass, chi, ndf, muon_sigma));

  nBranches_->Jpsi_trimu_maxdoca.push_back(getMaxDoca(allParticles4doc));
  nBranches_->Jpsi_trimu_mindoca.push_back(getMinDoca(allParticles4doc));
  nBranches_->Jpsi_trimu_vx.push_back(bc_vertex->vertexState().position().x());
  nBranches_->Jpsi_trimu_vy.push_back(bc_vertex->vertexState().position().y());
  nBranches_->Jpsi_trimu_vz.push_back(bc_vertex->vertexState().position().z());  

  nBranches_->Jpsi_trimu_unfitpt.push_back(tlv_B.Pt());
  nBranches_->Jpsi_trimu_unfitmass.push_back(tlv_B.M());
  nBranches_->Jpsi_trimu_unfitvprob.push_back(vprob_bc);
  
  if(vprob_bc!=-9){
    nBranches_->Jpsi_trimu_unfit_vx.push_back(vertex_bc.position().x());
    nBranches_->Jpsi_trimu_unfit_vy.push_back(vertex_bc.position().y());
    nBranches_->Jpsi_trimu_unfit_vz.push_back(vertex_bc.position().z());
  }


  
  nBranches_->Jpsi_ngenmuons.push_back(gen_nr_mu.size() + gen_jpsi_mu.size());

  bool flag_nr_match = false;
  if(gen_nr_mu.size()==1){
    Float_t _dR = reco::deltaR(
			       gen_nr_mu[0]->eta(), gen_nr_mu[0]->phi(), 
			       mu3.eta(), mu3.phi());

    if(_dR < 0.1) flag_nr_match = true;
  }

  bool flag_jpsi_match = false;
  if(gen_jpsi_mu.size()==2){

    Float_t _dR_11 = reco::deltaR(gen_jpsi_mu[0]->eta(), gen_jpsi_mu[0]->phi(), 
				muoncollection[mcidx_mu1].eta(), muoncollection[mcidx_mu1].phi());
    
    Float_t _dR_22 = reco::deltaR(gen_jpsi_mu[1]->eta(), gen_jpsi_mu[1]->phi(), 
				muoncollection[mcidx_mu2].eta(), muoncollection[mcidx_mu2].phi());

    Float_t _dR_21 = reco::deltaR(gen_jpsi_mu[1]->eta(), gen_jpsi_mu[1]->phi(), 
				muoncollection[mcidx_mu1].eta(), muoncollection[mcidx_mu1].phi());
    
    Float_t _dR_12 = reco::deltaR(gen_jpsi_mu[0]->eta(), gen_jpsi_mu[0]->phi(), 
				muoncollection[mcidx_mu2].eta(), muoncollection[mcidx_mu2].phi());
      

    if(_dR_11 < 0.1 && _dR_22 < 0.1) flag_jpsi_match = true;
    if(_dR_21 < 0.1 && _dR_12 < 0.1) flag_jpsi_match = true;
  }

  
  nBranches_->Jpsi_isgenmatched.push_back((int)flag_jpsi_match);
  nBranches_->Jpsi_mu3_isgenmatched.push_back((int)flag_nr_match);

  nBranches_->IsJpsiMu.push_back(1.);


}


void JpsiMuNtuplizer::printout(const RefCountedKinematicVertex& myVertex){
  std::cout << "Vertex:" << std::endl;
  if (myVertex->vertexIsValid()) {
    std::cout << "\t Decay vertex: " << myVertex->position() << myVertex->chiSquared() << " " << myVertex->degreesOfFreedom()
	      << std::endl;
  } else
    std::cout << "\t Decay vertex Not valid\n";
}

void JpsiMuNtuplizer::printout(const RefCountedKinematicParticle& myParticle){
  std::cout << "Particle:" << std::endl;
  //accessing the reconstructed Bs meson parameters:
  //SK: uncomment if needed  AlgebraicVector7 bs_par = myParticle->currentState().kinematicParameters().vector();

  //and their joint covariance matrix:
  //SK:uncomment if needed  AlgebraicSymMatrix77 bs_er = myParticle->currentState().kinematicParametersError().matrix();
  std::cout << "\t Momentum at vertex: " << myParticle->currentState().globalMomentum() << std::endl;
  std::cout << "\t Parameters at vertex: " << myParticle->currentState().kinematicParameters().vector() << std::endl;
}

void JpsiMuNtuplizer::printout(const RefCountedKinematicTree& myTree){
  if (!myTree->isValid()) {
    std::cout << "Tree is invalid. Fit failed.\n";
    return;
  }

  //accessing the tree components, move pointer to top
  myTree->movePointerToTheTop();

  //We are now at the top of the decay tree getting the B_s reconstructed KinematicPartlcle
  RefCountedKinematicParticle b_s = myTree->currentParticle();
  printout(b_s);

  // The B_s decay vertex
  RefCountedKinematicVertex b_dec_vertex = myTree->currentDecayVertex();
  printout(b_dec_vertex);

  // Get all the children of Bs:
  //In this way, the pointer is not moved
  std::vector<RefCountedKinematicParticle> bs_children = myTree->finalStateParticles();

  for (unsigned int i = 0; i < bs_children.size(); ++i) {
    printout(bs_children[i]);
  }

  std::cout << "\t ------------------------------------------" << std::endl;

  //Now navigating down the tree , pointer is moved:
  bool child = myTree->movePointerToTheFirstChild();

  if (child)
    while (myTree->movePointerToTheNextChild()) {
      RefCountedKinematicParticle aChild = myTree->currentParticle();
      printout(aChild);
    }
}
