#include "../interface/BsDstarTauNuNtuplizer.h"


//===================================================================================================================
BsDstarTauNuNtuplizer::BsDstarTauNuNtuplizer( edm::EDGetTokenT<pat::MuonCollection>    muonToken   ,
					      edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
					      edm::EDGetTokenT<reco::BeamSpot> bsToken,
					      edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken,
					      edm::EDGetTokenT<pat::PackedCandidateCollection> losttrackToken,
					      edm::EDGetTokenT<edm::TriggerResults> triggertoken,
					      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobject,
					      edm::EDGetTokenT<reco::GenParticleCollection> genptoken,
					      edm::EDGetTokenT<std::vector<reco::GenJet>> genttoken,
					      std::map< std::string, bool >& runFlags,
					      std::map< std::string, double >& runValues,
					      std::map< std::string, std::string >& runStrings,
					      NtupleBranches* nBranches )
: CandidateNtuplizer ( nBranches )
  , muonToken_	        ( muonToken )
  , verticeToken_          ( verticeToken )
  , bsToken_          ( bsToken )
  , packedpfcandidatesToken_(packedpfcandidatesToken) 
  , losttrackToken_(losttrackToken) 
  , HLTtriggersToken_	( triggertoken )
  , triggerObjects_	( triggerobject )
  , genParticlesToken_( genptoken )
  , genTauToken_( genttoken )
  , runOnMC_   (runFlags["runOnMC"])
  , useDNN_   (runFlags["useDNN"])
  , c_dz (runValues["dzcut"])
  , c_fsig (runValues["fsigcut"])
  , c_vprob (runValues["vprobcut"])
  , c_dnn (runValues["dnncut"])
  , dnnfile_ (runStrings["dnnfile"])   
{

  std::cout << "UseDNN = " << useDNN_ << std::endl;
  std::cout << "DNN file =" << dnnfile_ << std::endl;
  std::cout << "Confirm (dzcut, fsigcut, vprobcut) = " << c_dz << " " << c_fsig << " " << c_vprob << " " << c_dnn<< std::endl;
  
  if(useDNN_){
    //    std::string dnnfilepath = "EXOVVNtuplizerRunII/Ntuplizer/" +  dnnfile_;
    //    //    TFile * dnnfile = new TFile((TString)dnnfile_);
    //    std::cout << "dnn file input = " << edm::FileInPath(dnnfilepath).fullPath() << std::endl;
    //    TFile * dnnfile = new TFile((TString)edm::FileInPath(dnnfilepath).fullPath());
    //    TTree *tree = (TTree*) dnnfile->Get("tree");
    //    
    //    
    //    Float_t index[50];
    //    Float_t dnn[50];
    //    Float_t match[7];
    //    
    //    tree->SetBranchAddress("idx",&index);
    //    tree->SetBranchAddress("DNN",&dnn);
    //    tree->SetBranchAddress("match",&match);
    //    
    //    
    //    for(int ii=0; ii < tree->GetEntries(); ii++){
    //      
    //      Long64_t tentry = tree->LoadTree(ii);
    //      tree->GetEntry(tentry);
    //      
    //      std::vector<Int_t> vecidx;
    //      std::vector<Float_t> vecval;
    //      
    //      for(int jj=0; jj < 50; jj++){
    //	
    //	if(index[jj]!=-99){
    //
    //	  if(dnn[jj] > c_dnn){
    //
    //	    vecidx.push_back((int)index[jj]);
    //	    vecval.push_back(dnn[jj]);
    //
    //	  }
    //
    //	}
    //      }
    //      
    //      //      std::cout << "registered: " <<(Int_t)match[0] << " " << (Float_t)match[1] << " " << (Float_t)match[2] << " " << (Float_t)match[3] << " " << vecidx.size() << " " << vecval.size() << std::endl;
    //      DNNidx[std::make_tuple((Int_t)match[0], (Float_t)match[1], (Float_t)match[2], (Float_t)match[3])] = vecidx;
    //      DNNval[std::make_tuple((Int_t)match[0], (Float_t)match[1], (Float_t)match[2], (Float_t)match[3])] = vecval;
    //    }
    //    
    //    std::cout << tree->GetEntries() << " events read from DNN files" << std::endl;
    //  }
    std::string dnnfilepath = edm::FileInPath("EXOVVNtuplizerRunII/Ntuplizer/" +  dnnfile_).fullPath();

    std::cout << "dnn file path = " << dnnfilepath << std::endl;
    //    std::replace(dnnfilepath.begin(), dnnfilepath.end(), "D", "");

    std::string tbr = "DUMMY";  // to be replaced
    auto pos = dnnfilepath.find(tbr);
    auto len = tbr.length();
    if (pos != std::string::npos) {
      dnnfilepath.replace(pos, len, ""); // s == "a|b"
    }
    std::cout << "dnn file dir. = " << dnnfilepath << std::endl;    

    
    graphDef = tensorflow::loadMetaGraph(dnnfilepath);
    session = tensorflow::createSession(graphDef, dnnfilepath);
    //    graphDef = tensorflow::loadMetaGraph(edm::FileInPath(dnnfilepath).fullPath());
    //    session = tensorflow::createSession(graphDef, edm::FileInPath(dnnfilepath).fullPath());
    
    data = tensorflow::Tensor(tensorflow::DT_FLOAT, { 1, 50, 8 }); // single batch of dimension 10
    label = tensorflow::Tensor(tensorflow::DT_INT32, { 1,50}); 
    add_global = tensorflow::Tensor(tensorflow::DT_FLOAT, { 1, 2 }); 
    isTraining = tensorflow::Tensor(tensorflow::DT_BOOL, tensorflow::TensorShape()); 
    norm = tensorflow::Tensor(tensorflow::DT_FLOAT, { 1, 50 }); 
    
    //  for (size_t i = 0; i < 10; i++) input.matrix<float>()(0, i) = float(i);


  }

}

//===================================================================================================================
BsDstarTauNuNtuplizer::~BsDstarTauNuNtuplizer( void )
{

}


Int_t BsDstarTauNuNtuplizer::decaymode_id(std::string str){
  if(str=="electron") return -2;
  else if(str=="muon") return -1;
  else if(str=="oneProng0Pi0") return 0;
  else if(str=="oneProng1Pi0") return 1;
  else if(str=="oneProng2Pi0") return 2;
  else if(str=="oneProng3Pi0") return 3;
  else if(str=="oneProngOther") return 4;  
  else if(str=="threeProng0Pi0") return 10;
  else if(str=="threeProng1Pi0") return 11;
  else if(str=="threeProngOther") return 14;
  else if(str=="rare") return 15;
  else return -9;
}


TVector3 BsDstarTauNuNtuplizer::getVertex(const reco::GenParticle& part){
  return TVector3(part.vx(),part.vy(),part.vz());
}

float BsDstarTauNuNtuplizer::MuonPFIso(pat::Muon muon){

  float sumChargedHadronPt = muon.pfIsolationR04().sumChargedHadronPt;
  float sumNeutralHadronEt = muon.pfIsolationR04().sumNeutralHadronEt;
  float sumPhotonEt = muon.pfIsolationR04().sumPhotonEt;
  float sumPUPt = muon.pfIsolationR04().sumPUPt;
  float iso = (sumChargedHadronPt + std::max( 0. ,sumNeutralHadronEt + sumPhotonEt - 0.5 * sumPUPt));// / muon.pt()
 
  return iso;
}




Float_t BsDstarTauNuNtuplizer::getMaxDoca(std::vector<RefCountedKinematicParticle> &kinParticles){

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



Float_t BsDstarTauNuNtuplizer::getMinDoca(std::vector<RefCountedKinematicParticle> &kinParticles) {

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




std::tuple<Float_t, TransientVertex> BsDstarTauNuNtuplizer::vertexProb( const std::vector<reco::TransientTrack>& tracks){

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


//adapt absoluteImpactParameter functionality for RefCountedKinematicVertex
std::pair<bool, Measurement1D> BsDstarTauNuNtuplizer::absoluteImpactParameter(const TrajectoryStateOnSurface& tsos,
									      RefCountedKinematicVertex vertex,
									      VertexDistance& distanceComputer){
  if (!tsos.isValid()) {
    return std::pair<bool, Measurement1D>(false, Measurement1D(0., 0.));
  }
  GlobalPoint refPoint = tsos.globalPosition();
  GlobalError refPointErr = tsos.cartesianError().position();
  GlobalPoint vertexPosition = vertex->vertexState().position();
  GlobalError vertexPositionErr = RecoVertex::convertError(vertex->vertexState().error());
  return std::pair<bool, Measurement1D>(
					true,
					distanceComputer.distance(VertexState(vertexPosition, vertexPositionErr), VertexState(refPoint, refPointErr)));
}




particle_cand BsDstarTauNuNtuplizer::calculateIPvariables(
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


math::PtEtaPhiMLorentzVector BsDstarTauNuNtuplizer::daughter_p4(std::vector< RefCountedKinematicParticle > fitted_children, size_t i){
  const auto& state = fitted_children.at(i)->currentState();

  return math::PtEtaPhiMLorentzVector(
				      state.globalMomentum().perp(), 
				      state.globalMomentum().eta() ,
				      state.globalMomentum().phi() ,
				      state.mass()
				      );
}


bool BsDstarTauNuNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  
  // std::cout << "---------------- event, run, lumi = " << event.id().event() << " " << event.id().run() << " " << event.id().luminosityBlock() << "----------------" << std::endl;
  
  /********************************************************************
   *
   * Step1: check if the J/psi trigger is fired.
   * Namely, HLT_DoubleMu4_JpsiTrk_Displaced_v
   * and  HLT_Dimuon0_Jpsi3p5_Muon2_v
   ********************************************************************/

  event.getByToken(HLTtriggersToken_, HLTtriggers_);
  nBranches_->cutflow_perevt->Fill(0);

  bool isTriggered = false;
  const edm::TriggerNames& trigNames = event.triggerNames(*HLTtriggers_);
  std::vector<std::string> finalTriggerName;
  //    std::string finalTriggerFilterObjName="";

  for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {

    // if(trigNames.triggerName(i).find("HLT_DoubleMu4_JpsiTrk_Displaced_v")!= std::string::npos || trigNames.triggerName(i).find("HLT_Dimuon0_Jpsi3p5_Muon2_v")!= std::string::npos ){
           
    //     nBranches_->HLT_BPH_isFired[trigNames.triggerName(i)] = HLTtriggers_->accept(i);
    if(trigNames.triggerName(i).find("HLT_Mu8_IP3")!= std::string::npos || 
       trigNames.triggerName(i).find("HLT_Mu8_IP5")!= std::string::npos || 
       trigNames.triggerName(i).find("HLT_Mu8_IP6")!= std::string::npos || 
       trigNames.triggerName(i).find("HLT_Mu8p5_IP3p5")!= std::string::npos || 
       trigNames.triggerName(i).find("HLT_Mu9_IP4")!= std::string::npos || 
       trigNames.triggerName(i).find("HLT_Mu9_IP5")!= std::string::npos || 
       trigNames.triggerName(i).find("HLT_Mu9_IP6")!= std::string::npos || 
       trigNames.triggerName(i).find("HLT_Mu10p5_IP3p5")!= std::string::npos || 
       trigNames.triggerName(i).find("HLT_Mu12_IP6")!= std::string::npos || 
       trigNames.triggerName(i).find("HLT_Mu7_IP4")!= std::string::npos
       ){
      nBranches_->HLT_BPH_isFired[trigNames.triggerName(i)] = HLTtriggers_->accept(i);
      if(HLTtriggers_->accept(i)){
	isTriggered = true;
	//	std::cout << "This trigger is fired:" << trigNames.triggerName(i) << std::endl;
	finalTriggerName.push_back(trigNames.triggerName(i));
	//                finalTriggerFilterObjName="hltJpsiTkVertexFilter";
	// std::cout << "finalTriggerName = "  << finalTriggerName << std::endl;
          
      }
    }
  }


  /********************************************************************
   *
   * Step2: pre-select muons for building J/psi candidates ... 
   * For muons, no requirement applied
   *
   ********************************************************************/

  event.getByToken(verticeToken_   , vertices_     );
  event.getByToken(bsToken_   , beamspot_     );
  event.getByToken(muonToken_	, muons_    );
  event.getByToken(triggerObjects_  , triggerObjects);

  std::vector<pat::Muon> muoncollection;
  muoncollection.clear();

  if(!isTriggered) return false;
  nBranches_->cutflow_perevt->Fill(1);

  // evt Triggered

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
	      
	for(int iname=0; iname < (int)finalTriggerName.size(); iname ++){
	  if(pathNamesAll[h]==finalTriggerName[iname]) isPathExist = true;
	}
      }
      
      if(!isPathExist) continue;

      //            bool isFilterExist = false;
      //    
      //            for (unsigned hh = 0; hh < obj.filterLabels().size(); ++hh){
      //	
      //                if(obj.filterLabels()[hh].find(finalTriggerFilterObjName) != std::string::npos){
      //                    isFilterExist = true;
      //                }
      //            }
      //      
      //            if(!isFilterExist) continue;
      
      Float_t trigger_dR = reco::deltaR(obj.eta(), obj.phi(),
					muon.eta(), muon.phi());
      
      if(trigger_dR < 0.1) trigMatch = true;
    }

    if(!trigMatch) continue;

    //	std::cout << "imuon = " << imuon << " "  << muon.pt() << std::endl;

    muoncollection.push_back(muon);
  }


  //  std::cout << "number of matched muon = " << muoncollection.size() << std::endl;
  if(!( muoncollection.size() >= 1)) return false;
  nBranches_->cutflow_perevt->Fill(2);
    
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
    
  const reco::TrackRef track_muon = muoncollection[0].muonBestTrack();
  reco::TransientTrack tt_muon = (*builder).build(track_muon);

  KinematicParticleFactoryFromTransientTrack pFactory;
  KinematicParticleVertexFitter kpvFitter;

  edm::ESHandle<MagneticField> fieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(fieldHandle);
  fMagneticField = fieldHandle.product();

  AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);


  reco::Vertex closestVertex; 
  closestVertex = *(vertices_->begin());
  Float_t max_dz = 999;

  for( reco::VertexCollection::const_iterator vtx = vertices_->begin(); vtx != vertices_->end(); ++vtx){
      
    Float_t _dz_ = TMath::Abs(vtx->position().z() - muoncollection[0].vz());
      
    //      bool isFake = (vtx->chi2()==0 && vtx->ndof()==0);
    //      if( !isFake && vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0) {
    //	firstGoodVertex = vtx;
    //	break;
    //      }


      
    if(_dz_ < max_dz){
      max_dz = _dz_;
      closestVertex = *vtx;
    }
  }



  /********************************************************************
   *
   * Step6: Tau selection
   *        Just select highest in pT but there might be better selection ... 
   *
   ********************************************************************/
    
  event.getByToken( packedpfcandidatesToken_               , packedpfcandidates_      ); 
  event.getByToken( losttrackToken_               , losttrack_      ); 

  // This is for constructing D* meson
  std::vector<pat::PackedCandidate> pfcollection_pre; 
  std::vector<reco::TransientTrack> mytracks_pre;
    
  // This is for constructing tau
  std::vector<pat::PackedCandidate> pfcollection; 
  std::vector<reco::TransientTrack> mytracks;
  std::vector<Float_t> mydnn;

  Float_t ptsum = 0.;

  for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){   
    
    pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
    
    if(pf.pt() < 0.5) continue;
    if(!pf.hasTrackDetails()) continue;
    
    // use the PF candidates that come from closestVertex
    //      if(pf.vertexRef()->z()!=closestVertex.position().z()) continue;
    
    //      Float_t precut_dz = pf.vertexRef()->z() - closestVertex.position().z();
    Float_t precut_dz = pf.vz() - closestVertex.position().z();
    if(TMath::Abs(precut_dz) > c_dz) continue;
    
    Bool_t hpflag = pf.trackHighPurity();
    if(!hpflag) continue;
    if(pf.pseudoTrack().hitPattern().numberOfValidPixelHits() < 0) continue;
    if(pf.pseudoTrack().hitPattern().numberOfValidHits() < 3) continue;
    if(pf.pseudoTrack().normalizedChi2() > 100) continue;
    
    if(TMath::Abs(pf.pdgId())!=211) continue; 
    if(TMath::Abs(pf.eta()) > 2.5) continue; 
    
    pfcollection_pre.push_back(pf);
    reco::TransientTrack  tt_track = (*builder).build(pf.pseudoTrack());
    mytracks_pre.push_back(tt_track);
    
    ptsum += pf.pt();

  }
    
  //    for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){   
  //	
  //      pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
  //	
  //      pfcollection_pre.push_back(pf);
  //    }

  Int_t npf_before_dnn = 0;
  Int_t npf_qr = 0;


  if(useDNN_){


    ///////////////////////////////
    std::vector<pfcand> pfcands;
    //      std::vector<pat::PackedCandidate> pfmuon;


    Int_t count_dnn = 0;
    Int_t count_dnn_muon = 0;

    for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){   
	
      pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
	
      // first for the muon ... 

      //	if(TMath::Abs(pf.pdgId())==13){
      //	  std::cout << pf.pt() << " " << pf.eta() << " " << pf.phi() << " " <<  TMath::Abs(pf.charge()) << " " <<  TMath::Abs(pf.pdgId()) << " " <<  pf.pt() << " " <<  pf.isGlobalMuon() << " " <<  pf.hasTrackDetails() << std::endl;
      //	}
	
      if(TMath::Abs(pf.eta()) < 2.4 && 
	 TMath::Abs(pf.charge())==1 &&
	 //	   TMath::Abs(pf.pdgId())==13 && 
	 pf.pt() > 7. &&
	 pf.isGlobalMuon() > 0.5 &&
	 pf.hasTrackDetails() > 0.5
	 ){

	//	  std::cout << "CHECK1: " <<count_dnn<< std::endl;	  

	if(count_dnn < 50){
	  data.tensor<float, 3>()(0, count_dnn, 0) = pf.eta();
	  data.tensor<float, 3>()(0, count_dnn, 1) = pf.phi();
	  data.tensor<float, 3>()(0, count_dnn, 2) = TMath::Log(pf.pt());
	  data.tensor<float, 3>()(0, count_dnn, 3) = TMath::Log(pf.energy());
	  data.tensor<float, 3>()(0, count_dnn, 4) = pf.charge();
	  data.tensor<float, 3>()(0, count_dnn, 5) = TMath::Abs(closestVertex.position().z() - pf.pseudoTrack().vz());
	  data.tensor<float, 3>()(0, count_dnn, 6) = TMath::Sqrt( TMath::Power((closestVertex.position().x() - pf.pseudoTrack().vx()), 2) + TMath::Power((closestVertex.position().y() - pf.pseudoTrack().vy()), 2));
	  data.tensor<float, 3>()(0, count_dnn, 7) = pf.isGlobalMuon();
	    
	  label.matrix<int>()(0, count_dnn) = 0;
	  norm.matrix<float>()(0, count_dnn) = float(1);

	  count_dnn_muon++;
	  count_dnn++;
	}
      }

	
      if(!pf.hasTrackDetails()) continue;
      Float_t precut_dz = pf.vz() - closestVertex.position().z();
      if(TMath::Abs(precut_dz) > c_dz) continue;
	
      npf_qr++;
	

      if(pf.pt() < 0.5) continue;
	
      // use the PF candidates that come from closestVertex
      //      if(pf.vertexRef()->z()!=closestVertex.position().z()) continue;
	
      //      Float_t precut_dz = pf.vertexRef()->z() - closestVertex.position().z();
	
      Bool_t hpflag = pf.trackHighPurity();
      if(!hpflag) continue;
      if(pf.pseudoTrack().hitPattern().numberOfValidPixelHits() < 0) continue;
      if(pf.pseudoTrack().hitPattern().numberOfValidHits() < 3) continue;
      if(pf.pseudoTrack().normalizedChi2() > 100) continue;
	
      if(TMath::Abs(pf.pdgId())!=211) continue; 
      if(TMath::Abs(pf.eta()) > 2.5) continue; 

      //	pfcollection.push_back(pf);
      //	reco::TransientTrack  tt_track = (*builder).build(pf.pseudoTrack());
      //	mytracks.push_back(tt_track);
	
      npf_before_dnn++;	

      pfcand _cand_ = {
	(Int_t)ii,
	(Float_t) abs(precut_dz)
      };
	  
      //	std::cout << cands.size() << std::endl;
      pfcands.push_back(_cand_);
    }

    //sorting by distance to the vertex
    sort(pfcands.begin(), pfcands.end());
      
    // filling information in for evaluation

    //      for(int imu = 0; imu < (int)muoncollection.size() && imu<50; imu++){
    //	const pat::Muon mu = muoncollection[imu];
    //
    //
    //      }

    for(size_t ic = 0; ic < pfcands.size(); ic++){
      Int_t idx = pfcands[ic].cand_idx;

      pat::PackedCandidate pf = (*packedpfcandidates_)[idx];

      //	std::cout << "CHECK2: " <<count_dnn<< " " << pfcands[ic].cand_absdz <<std::endl;

      if(count_dnn < 50){
	data.tensor<float, 3>()(0, count_dnn, 0) = pf.eta();
	data.tensor<float, 3>()(0, count_dnn, 1) = pf.phi();
	data.tensor<float, 3>()(0, count_dnn, 2) = TMath::Log(pf.pt());
	data.tensor<float, 3>()(0, count_dnn, 3) = TMath::Log(pf.energy());
	data.tensor<float, 3>()(0, count_dnn, 4) = pf.charge();
	data.tensor<float, 3>()(0, count_dnn, 5) = TMath::Abs(closestVertex.position().z() - pf.pseudoTrack().vz());
	data.tensor<float, 3>()(0, count_dnn, 6) = TMath::Sqrt( TMath::Power((closestVertex.position().x() - pf.pseudoTrack().vx()), 2) + TMath::Power((closestVertex.position().y() - pf.pseudoTrack().vy()), 2));
	data.tensor<float, 3>()(0, count_dnn, 7) = pf.isGlobalMuon();
	  
	label.matrix<int>()(0, count_dnn) = 0;
	norm.matrix<float>()(0, count_dnn) = float(1);
	count_dnn++;

	pfcollection.push_back(pf);
	reco::TransientTrack  tt_track = (*builder).build(pf.pseudoTrack());
	mytracks.push_back(tt_track);

      }


    }

      
    for(int ic=count_dnn; ic<50; ic++){
	
      //	std::cout << "CHECK3: " <<ic<< std::endl;
      data.tensor<float, 3>()(0, ic, 0) = 0;
      data.tensor<float, 3>()(0, ic, 1) = 0;
      data.tensor<float, 3>()(0, ic, 2) = 0;
      data.tensor<float, 3>()(0, ic, 3) = 0;
      data.tensor<float, 3>()(0, ic, 4) = 0;
      data.tensor<float, 3>()(0, ic, 5) = 0;
      data.tensor<float, 3>()(0, ic, 6) = 0;
      data.tensor<float, 3>()(0, ic, 7) = 0;
	
      label.matrix<int>()(0, ic) = 0;
      norm.matrix<float>()(0, ic) = float(1);

    }
      

    add_global.matrix<float>()(0, 0) = float(count_dnn_muon); // Number of muons around 0.5 cm from PV
    add_global.matrix<float>()(0, 1) = float(count_dnn/100); //Number of charged pf candidates around 0.5 cm from PV
    isTraining.scalar<bool>()() = false; //Number of charged pf candidates around 0.5 cm from PV


    std::vector<tensorflow::Tensor> outputs;
    tensorflow::run(session, {  { "Placeholder:0", data },  { "Placeholder_1:0", label }, { "Placeholder_2:0", add_global } , {"Placeholder_3:0", isTraining}, {"Placeholder_4:0", norm}}, { "Reshape_13:0" }, &outputs);
      
    auto finalOutputTensor = outputs[0].tensor<float, 3>();
      
    //      std::cout << "bg:" << finalOutputTensor(0, 0, 0) << std::endl;
    //      std::cout << "count_muon, count_dnn = " << count_dnn_muon << " " << count_dnn << std::endl;

    for(int ic=count_dnn_muon; ic<count_dnn; ic++){
      //	std::cout << "pf dnn = " << ic << " " << finalOutputTensor(0, ic, 0) << " " << finalOutputTensor(0, ic, 1) << std::endl;
      mydnn.push_back(finalOutputTensor(0, ic, 1));
    }






    /////////////////////////////////
    /////////////////////////////////


    //      std::vector<Int_t> dnnidx;
    //      std::vector<Float_t> dnnval;
    //
    //      if(pfcollection_pre.size()>=1){
    //	dnnidx = DNNidx[std::make_tuple((Int_t)event.id().event(), (Float_t)pfcollection_pre[0].pt(), (Float_t)pfcollection_pre[0].eta(), (Float_t)pfcollection_pre[0].phi())];
    //	dnnval = DNNval[std::make_tuple((Int_t)event.id().event(), (Float_t)pfcollection_pre[0].pt(), (Float_t)pfcollection_pre[0].eta(), (Float_t)pfcollection_pre[0].phi())];
    //	std::cout << "check:" << (Int_t)event.id().event() << " " << (Float_t)pfcollection_pre[0].pt() << " " <<  (Float_t)pfcollection_pre[0].eta() << " " << (Float_t)pfcollection_pre[0].phi() << " " << dnnidx.size() << " " << dnnval.size() << std::endl;
    //      }
    //
    //      for(unsigned int idx=0; idx < dnnidx.size(); idx++){
    //	Int_t dnn_idx = dnnidx[idx]; 
    //	Float_t dnn_val = dnnval[idx];
    //	//	if(iii==-1) continue;
    //	
    //	std::cout << dnn_idx << " " << dnn_val << std::endl;
    //	//	std::cout << iii << " " << std::endl; 
    //	
    //	const pat::PackedCandidate prefilter_pf = (*packedpfcandidates_)[dnn_idx];
    //	//	if(prefilter_pf.pt() < 0.5) continue;	      
    //	//	if(!prefilter_pf.hasTrackDetails()) continue;
    //    
    //	//	std::cout << "\t found list but removed due to no track details ..." << iii << " " << prefilter_pf.hasTrack Details()<< std::endl; 
    //	
    //
    //	Float_t precut_dz = prefilter_pf.vz() - closestVertex.position().z();
    //	if(TMath::Abs(precut_dz) > c_dz) continue;
    //	
    //	Bool_t hpflag = prefilter_pf.trackHighPurity();
    //	if(!hpflag) continue;
    //	
    //	//	if(prefilter_pf.pseudoTrack().hitPattern().numberOfValidHits() < 3) continue;
    //	//	if(prefilter_pf.pseudoTrack().normalizedChi2() > 100) continue;
    //	
    //	//	if(TMath::Abs(pf.pdgId())!=211) continue; 
    //
    //	if(TMath::Abs(prefilter_pf.eta()) > 2.3) continue; 
    //	
    //	//	pfcollection.push_back(prefilter_pf);
    //	reco::TransientTrack  tt_track = (*builder).build(prefilter_pf.pseudoTrack());
    //	mytracks.push_back(tt_track);
    //	mydnn.push_back(dnn_val);
    //
    //      }
      
    /////////////////////////////////
    /////////////////////////////////



  }else{


    for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){   
      
      pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
	
      if(pf.pt() < 0.5) continue;
      if(!pf.hasTrackDetails()) continue;
	
      // use the PF candidates that come from closestVertex
      //      if(pf.vertexRef()->z()!=closestVertex.position().z()) continue;
	
      //      Float_t precut_dz = pf.vertexRef()->z() - closestVertex.position().z();
      Float_t precut_dz = pf.vz() - closestVertex.position().z();
      if(TMath::Abs(precut_dz) > c_dz) continue;
	
      Bool_t hpflag = pf.trackHighPurity();
      if(!hpflag) continue;
      if(pf.pseudoTrack().hitPattern().numberOfValidPixelHits() < 0) continue;
      if(pf.pseudoTrack().hitPattern().numberOfValidHits() < 3) continue;
      if(pf.pseudoTrack().normalizedChi2() > 100) continue;
	
      if(TMath::Abs(pf.pdgId())!=211) continue; 
      if(TMath::Abs(pf.eta()) > 2.5) continue; 

      //      Float_t _dR1 = reco::deltaR(pf.eta(), pf.phi(), 
      //				  mu1_fit->eta(), mu1_fit->phi());
      //      
      //      Float_t _dR2 = reco::deltaR(pf.eta(), pf.phi(), 
      //				  mu2_fit->eta(), mu2_fit->phi());
      //
      //
      //      if(_dR1 < 0.1 || _dR2 < 0.1){
      //	if(TMath::Abs(pf.pdgId()) == 13) continue;
      //      }

      pfcollection.push_back(pf);
      reco::TransientTrack  tt_track = (*builder).build(pf.pseudoTrack());
      mytracks.push_back(tt_track);
	
    }
  }

  // retrieve gen. information 
  Int_t numOfch_pre = (size_t)pfcollection_pre.size();
  Int_t numOfch = (size_t)pfcollection.size();


  std::vector<std::vector<TLorentzVector>> gps;
  std::vector<Int_t> ppdgId;
  std::vector<Int_t> vec_gentaudm;
  std::vector<Int_t> vec_ppdgId;
  std::vector<TLorentzVector> vec_gentaup4;
  std::vector<TLorentzVector> vec_gentau3pp4;
  Int_t isgen3 = 0;
  Int_t isgen3matched = 0;

  bool isMC = runOnMC_;

  if(isMC){
    event.getByToken(genParticlesToken_ , genParticles_); 
    event.getByToken(genTauToken_, genTaus_);

    for( unsigned p=0; p < genParticles_->size(); ++p){
      
      if(TMath::Abs((*genParticles_)[p].pdgId())!=15) continue;
      if(TMath::Abs((*genParticles_)[p].status())!=2) continue;
      
      std::cout << "\t Tau found with # of daughters = " << (*genParticles_)[p].numberOfDaughters() << " with mother = " << (*genParticles_)[p].mother(0)->pdgId() << std::endl;
      
      
      // calculate visible pt ... 

      TLorentzVector genvis;
      std::vector<TLorentzVector> gp;
      Bool_t matched = true;
      Int_t nprong = 0;
      
      for(int idd = 0; idd < (int)(*genParticles_)[p].numberOfDaughters(); idd++){
	
	std::cout << "\t\t -> " << (*genParticles_)[p].daughter(idd)->pdgId() << " (pT, eta, phi) = " 
		  << (*genParticles_)[p].daughter(idd)->pt() << " " 
		  << (*genParticles_)[p].daughter(idd)->eta() << " " 
		  << (*genParticles_)[p].daughter(idd)->phi() << std::endl;


	if(
	   TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==12 ||
	   TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==14 || 
	   TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==16
	   ) continue;


	TLorentzVector _genvis_;
	_genvis_.SetPtEtaPhiM((*genParticles_)[p].daughter(idd)->pt(),
			      (*genParticles_)[p].daughter(idd)->eta(),
			      (*genParticles_)[p].daughter(idd)->phi(),
			      (*genParticles_)[p].daughter(idd)->mass()
			      );
	
	genvis += _genvis_;

	if(TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==211){

	  nprong += 1;
	  
	  // check matching to reco PF objects
	  Float_t min_dr = 999;
	  
	  for(int kkk = 0; kkk < numOfch; kkk ++){
	    
	    pat::PackedCandidate _pf = pfcollection[kkk];
	    
	    if(_pf.pdgId()!=(*genParticles_)[p].daughter(idd)->pdgId()) continue;
	    
	    Float_t _dR = reco::deltaR(
				       _genvis_.Eta(), _genvis_.Phi(),
				       _pf.eta(), _pf.phi()
				       );
	    if(_dR < min_dr && _dR < 0.015 && _pf.pt()/_genvis_.Pt() < 1.15 && _pf.pt()/_genvis_.Pt() > 0.85){
	      //	    if(_dR < min_dr && _dR < 0.1){
	      min_dr = _dR;
	    }
	  }
	  
	  //	  Float_t min_dr2 = 999;
	  //	  for( size_t iii = 0; iii < packedpfcandidates_->size(); ++iii ){   
	  //      
	  //	    pat::PackedCandidate pf = (*packedpfcandidates_)[iii];
	  //	    
	  //	    if(pf.pdgId()!=(*genParticles_)[p].daughter(idd)->pdgId()) continue;
	  //	    
	  //	    Float_t _dR = reco::deltaR(
	  //				       _genvis_.Eta(), _genvis_.Phi(),
	  //				       pf.eta(), pf.phi()
	  //				       );
	  //	    if(_dR < min_dr2 && _dR < 0.1){
	  //	      min_dr2 = _dR;
	  //	    }
	  //	  }
	  //
	  //	  if(min_dr2!=999) std::cout << "pf matched !!!" << std::endl;


	  //////////////////////////////////////
	  if(min_dr == 999) matched = false;
	  //	  else std::cout << "matched!" << std::endl;
	  //	  else{
	  gp.push_back(_genvis_);
	  //	  }
	}
      }


      if(nprong==3) isgen3 += 1;

      //////////////////////////
      // check decay mod. To do this, take matching with tau-genjet. 
      //////////////////////////
      Float_t min_gendr = 999;
      Int_t taugendm = -999;

      for(size_t i = 0; i < genTaus_->size(); ++ i){      
	
	const reco::GenJet & TauCand = (*genTaus_)[i];
	
	reco::Particle::LorentzVector visibleP4 = ((*genTaus_)[i]).p4();

	TLorentzVector visp4;
	visp4.SetPtEtaPhiM(visibleP4.pt(),
			   visibleP4.eta(),
			   visibleP4.phi(),
			   visibleP4.mass());
	
	Float_t dRgen = genvis.DeltaR(visp4);
	
	if(dRgen < min_gendr && dRgen < 0.1){
	  min_gendr = dRgen;
	  taugendm = decaymode_id(JetMCTagUtils::genTauDecayMode(TauCand));
	}
      }
      
      vec_ppdgId.push_back((*genParticles_)[p].mother(0)->pdgId());
      vec_gentaudm.push_back(taugendm);
      vec_gentaup4.push_back(genvis);


  

      
      if(gp.size()==3){
	std::cout << "\t -----> This has been registered with mother = " << (*genParticles_)[p].mother(0)->pdgId() << std::endl;
	gps.push_back(gp);
	ppdgId.push_back((*genParticles_)[p].mother(0)->pdgId());
	vec_gentau3pp4.push_back(genvis);
	
	//	if(TMath::Abs((*genParticles_)[p].mother(0)->pdgId())==541){
	//isgen3matched = matched;
	if(matched) isgen3matched += 1;

	//	}
	//      }else{

	//      }
      }

      //    std::cout << "\t # of gen. taus with 3prong = " << gps.size() << std::endl;

    }
  }
  //////////////////////////////
  // constructing phi meson ... 


  //  std::cout << "Starts to build D* meson out of " << numOfch_pre << " pion candidates" << std::endl;

  //   particle_cand Phi_cand;
  RefCountedKinematicParticle D0_part_rec;
  RefCountedKinematicVertex D0_vertex_rec;
  particle_cand D0cand_rec; 
  Float_t D0_max_pt = -1;
  TLorentzVector D0_tlv_highest;
  //  Float_t D0_vprob_highest = -9;
  Int_t kidx = -1;
  Int_t pidx = -1;
  Int_t qk = -9;
  Int_t qp = -9;
    
  for(int iii = 0; iii < numOfch_pre; iii ++){
     
    pat::PackedCandidate pf1 = pfcollection_pre[iii];
    //    const reco::Track pf1_track = pf1.pseudoTrack();
     
    for(int jjj = iii+1; jjj < numOfch_pre; jjj ++){
       
      pat::PackedCandidate pf2 = pfcollection_pre[jjj];
      //      const reco::Track pf2_track = pf2.pseudoTrack();
       
      Int_t charge = pf1.charge() + pf2.charge();
      
      if(TMath::Abs(charge)!=0) continue; 

      //      bool massflag = false;
      
      for(int kp=0; kp<2; kp++){

	TLorentzVector tlv_kaon;
	TLorentzVector tlv_pion;

	std::vector<RefCountedKinematicParticle> D0Particles;
	
	if(kp==0){
	  D0Particles.push_back(pFactory.particle(mytracks_pre[iii], kaon_mass, chi, ndf, kaon_sigma));
	  D0Particles.push_back(pFactory.particle(mytracks_pre[jjj], pion_mass, chi, ndf, pion_sigma));

	  tlv_kaon.SetPtEtaPhiM(pf1.pt(), pf1.eta(), pf1.phi(), mass_kaon);
	  tlv_pion.SetPtEtaPhiM(pf2.pt(), pf2.eta(), pf2.phi(), mass_pion);

	}else{
	  D0Particles.push_back(pFactory.particle(mytracks_pre[jjj], kaon_mass, chi, ndf, kaon_sigma));
	  D0Particles.push_back(pFactory.particle(mytracks_pre[iii], pion_mass, chi, ndf, pion_sigma));

	  tlv_pion.SetPtEtaPhiM(pf1.pt(), pf1.eta(), pf1.phi(), mass_pion);
	  tlv_kaon.SetPtEtaPhiM(pf2.pt(), pf2.eta(), pf2.phi(), mass_kaon);
	}

	TLorentzVector tlv_D0 = tlv_kaon + tlv_pion;
	//	Float_t D0_mass = D0_part->currentState().mass();
	Float_t D0_mass = tlv_D0.M();
	//	std::cout << "D0 mass = "<< D0_mass << std::endl;
	if(!(D0_mass < (mass_D0 + 0.04) && D0_mass > (mass_D0 - 0.04))) continue;

	
	//reconstructing a J/Psi decay
	RefCountedKinematicTree D0Tree = kpvFitter.fit(D0Particles);
	
	if(D0Tree->isEmpty() || !D0Tree->isValid() || !D0Tree->isConsistent()) continue;
	
	//creating the particle fitter
	KinematicParticleFitter csFitter;
	
	// creating the constraint
	KinematicConstraint* D0_constraint = new MassKinematicConstraint(d0_mass, d0_sigma);
	//the constrained fit
	D0Tree = csFitter.fit(D0_constraint, D0Tree);
	
	//getting the J/Psi KinematicParticle
	D0Tree->movePointerToTheTop();
	RefCountedKinematicParticle D0_part = D0Tree->currentParticle();
	if(!D0_part->currentState().isValid()) continue;
	
	RefCountedKinematicVertex D0_vertex = D0Tree->currentDecayVertex();
	if(!D0_vertex->vertexIsValid()) continue;
	
	if(TMath::Prob(D0_vertex->chiSquared(), D0_vertex->degreesOfFreedom()) <= 0.1) continue;
	
	particle_cand D0cand; 
	D0cand = calculateIPvariables(extrapolator, D0_part, D0_vertex, closestVertex);
	
	//	  if(Taucand.fls3d < 3) continue;
	
	// 6.1.2020 commented out
	if(D0cand.fls3d < 2) continue;
	
	Float_t D0_pt = D0_part->currentState().globalMomentum().perp();



//	std::vector< RefCountedKinematicParticle > children = D0Tree->finalStateParticles();
//	  
//	math::PtEtaPhiMLorentzVector part1 = daughter_p4(children, 0);
//	math::PtEtaPhiMLorentzVector part2 = daughter_p4(children, 1);
//	
//	math::PtEtaPhiMLorentzVector tlv_D0 = part1 + part2;
	


	if(D0_max_pt < D0_pt){
	  D0_max_pt = D0_pt;
	  D0cand_rec = D0cand;
	  D0_vertex_rec = D0_vertex;
	  kidx = (kp==0 ? iii : jjj);
	  pidx = (kp==0 ? jjj : iii);
	  D0_tlv_highest = tlv_D0;
	  //	  D0_vprob_highest = TMath::Prob(D0_vertex->chiSquared(), D0_vertex->degreesOfFreedom());
	  D0_part_rec = D0_part;
	  qk = (kp==0 ? pf1.charge() : pf2.charge());
	  qp = (kp==0 ? pf2.charge() : pf1.charge());
	}
      }
    }
  }

  if(D0_max_pt==-1) return false;
  nBranches_->cutflow_perevt->Fill(3);

  //  std::cout << "D0 candidate found with mass = " << D0_tlv_highest.M() << std::endl;

  /********************************************************************
   *
   * Step4: Construct D* meson
   *
   ********************************************************************/

  RefCountedKinematicParticle Ds_part_rec;
  RefCountedKinematicVertex Ds_vertex_rec;
  particle_cand Dscand_rec; 
  Float_t Ds_max_pt = -1;
  TLorentzVector Ds_tlv_highest;
  //  Float_t Ds_vprob_highest = -9;
  Int_t sidx = -1;
  Int_t qsp = -9;

  for(int kkk = 0; kkk < numOfch_pre; kkk ++){
    
    if(kkk==kidx || kkk == pidx) continue;

    pat::PackedCandidate pf = pfcollection_pre[kkk];

    if(pf.charge()!=qp) continue;
    
    std::vector<RefCountedKinematicParticle> DstarParticles;
    
    //    DstarParticles.push_back(pFactory.particle(mytracks_pre[kidx], kaon_mass, chi, ndf, kaon_sigma));
    //    DstarParticles.push_back(pFactory.particle(mytracks_pre[pidx], pion_mass, chi, ndf, pion_sigma));
    DstarParticles.push_back(pFactory.particle(mytracks_pre[kkk], pion_mass, chi, ndf, pion_sigma));
    DstarParticles.push_back(D0_part_rec);


    //    TLorentzVector tlv_D0;
    TLorentzVector tlv_spion;
//    tlv_D0.SetPtEtaPhiM(D0_part_rec->currentState().globalMomentum().perp(),
//			D0_part_rec->currentState().globalMomentum().eta(),
//			D0_part_rec->currentState().globalMomentum().phi(),
//			D0_part_rec->currentState().mass());
			
    
    tlv_spion.SetPtEtaPhiM(pf.pt(), pf.eta(), pf.phi(), pion_mass);
    //    TLorentzVector tlv_Ds = tlv_D0 + tlv_spion;
    TLorentzVector tlv_Ds = D0_tlv_highest + tlv_spion;
    Float_t Ds_mass = tlv_Ds.M();

    //    Float_t Dsmass = Ds_part->currentState().mass();
    //    std::cout << "Ds mass = " << Ds_mass << std::endl;
    if(!(Ds_mass < (mass_Dstar + 0.03) && Ds_mass > (mass_Dstar - 0.03))) continue;

    //reconstructing a J/Psi decay
    RefCountedKinematicTree DsTree = kpvFitter.fit(DstarParticles);
    
    if(DsTree->isEmpty() || !DsTree->isValid() || !DsTree->isConsistent()) continue;
    
    //creating the particle fitter
    KinematicParticleFitter csFitter;
    
    // creating the constraint
    KinematicConstraint* Ds_constraint = new MassKinematicConstraint(ds_mass, ds_sigma);
    //the constrained fit
    DsTree = csFitter.fit(Ds_constraint, DsTree);
    
    //getting the J/Psi KinematicParticle
    DsTree->movePointerToTheTop();
    RefCountedKinematicParticle Ds_part = DsTree->currentParticle();
    if(!Ds_part->currentState().isValid()) continue;
    
    RefCountedKinematicVertex Ds_vertex = DsTree->currentDecayVertex();
    if(!Ds_vertex->vertexIsValid()) continue;
    
    if(TMath::Prob(Ds_vertex->chiSquared(), Ds_vertex->degreesOfFreedom()) <= 0.1) continue;
    
    particle_cand Dscand; 
    Dscand = calculateIPvariables(extrapolator, Ds_part, Ds_vertex, closestVertex);
    
    //	  if(Taucand.fls3d < 3) continue;
    
    // 6.1.2020 commented out
    if(Dscand.fls3d < 2) continue;
    
    Float_t Ds_pt = Ds_part->currentState().globalMomentum().perp();
    

    //    std::vector< RefCountedKinematicParticle > children = DsTree->finalStateParticles();

    //    std::cout << "# of sizes for Dstar (should be 2) = " << children.size() << std::endl;
	  
    //    math::PtEtaPhiMLorentzVector part1 = daughter_p4(children, 0);
    //math::PtEtaPhiMLorentzVector part2 = daughter_p4(children, 1);
    
    //    math::PtEtaPhiMLorentzVector tlv_Ds = part1 + part2;
    
    if(Ds_max_pt < Ds_pt){
      Ds_max_pt = Ds_pt;
      Dscand_rec = Dscand;
      Ds_part_rec = Ds_part;
      Ds_vertex_rec = Ds_vertex;
      Ds_tlv_highest = tlv_Ds;
      //      Ds_vprob_highest = TMath::Prob(Ds_vertex->chiSquared(), Ds_vertex->degreesOfFreedom());
      sidx = kkk;
      qsp = pf.charge();
    }    
  }

  if(Ds_max_pt==-1) return false;
  nBranches_->cutflow_perevt->Fill(4);

  std::cout << "(D0, Ds) candidates found with mass = (" << D0_tlv_highest.M() << " " << Ds_tlv_highest.M() << ")" << std::endl;
  std::cout << "Starts to build tau candidate out of " << numOfch << " pion candidates" << std::endl;

  std::vector<taucand> cands;
  Int_t npf_after_dnn = 0;
    
  for(int iii = 0; iii < numOfch; iii ++){
      
    if(iii==kidx || iii==pidx || iii==sidx) continue;
    
    pat::PackedCandidate pf1 = pfcollection[iii];

    if(mydnn[iii] < c_dnn) continue;
    npf_after_dnn++;

    //      const reco::Track pf1_track = pf1.pseudoTrack();

    for(int jjj = iii+1; jjj < numOfch; jjj ++){

      if(jjj==kidx || jjj==pidx || jjj==sidx) continue;
	
      pat::PackedCandidate pf2 = pfcollection[jjj];
      if(mydnn[jjj] < c_dnn) continue;
      //	const reco::Track pf2_track = pf2.pseudoTrack();

      for(int kkk = jjj+1; kkk < numOfch; kkk ++){

	if(kkk==kidx || kkk==pidx || kkk==sidx) continue;

	pat::PackedCandidate pf3 = pfcollection[kkk];
	if(mydnn[kkk] < c_dnn) continue;

	//	  const reco::Track pf3_track = pf3.pseudoTrack();


	Int_t tau_charge = pf1.charge() + pf2.charge() + pf3.charge(); 

	if(TMath::Abs(tau_charge)!=1) continue; 


	TLorentzVector tlv_pion1; 
	TLorentzVector tlv_pion2;
	TLorentzVector tlv_pion3;
	  
	tlv_pion1.SetPtEtaPhiM(pf1.pt(), pf1.eta(), pf1.phi(), pf1.mass());
	tlv_pion2.SetPtEtaPhiM(pf2.pt(), pf2.eta(), pf2.phi(), pf2.mass());
	tlv_pion3.SetPtEtaPhiM(pf3.pt(), pf3.eta(), pf3.phi(), pf3.mass());
	  
	TLorentzVector tlv_tau_pre = tlv_pion1 + tlv_pion2 + tlv_pion3;
	if(tlv_tau_pre.Pt() < 3) continue;
	
	//	Float_t taumass = tlv_tau_pre.M();
	//	if(!(taumass > 0.2 && taumass < 1.5)) continue;


	/* reconstruct taus*/

	std::vector<RefCountedKinematicParticle> tauParticles;

	tauParticles.push_back(pFactory.particle(mytracks[iii], pion_mass, chi, ndf, pion_sigma));
	tauParticles.push_back(pFactory.particle(mytracks[jjj], pion_mass, chi, ndf, pion_sigma));
	tauParticles.push_back(pFactory.particle(mytracks[kkk], pion_mass, chi, ndf, pion_sigma));

  
	//reconstructing a tau decay
	RefCountedKinematicTree tauTree = kpvFitter.fit(tauParticles);


	if(tauTree->isEmpty() || !tauTree->isValid() || !tauTree->isConsistent()) continue;

	//creating the particle fitter
	//	  KinematicParticleFitter csFitter;

	// creating the constraint
	//	  KinematicConstraint* tau_constraint = new MassKinematicConstraint(tau_mass, tau_m_sigma);

	//    //the constrained fit
	//    tauTree = csFitter.fit(jpsi_constraint, tauTree);

	//getting the J/Psi KinematicParticle
	tauTree->movePointerToTheTop();


	RefCountedKinematicParticle tau_part = tauTree->currentParticle();
	if(!tau_part->currentState().isValid()) continue;
	RefCountedKinematicVertex tau_vertex = tauTree->currentDecayVertex();
	if(!tau_vertex->vertexIsValid()) continue; 

	// 6.1.2020 commented out
	if(TMath::Prob(tau_vertex->chiSquared(), tau_vertex->degreesOfFreedom()) <= c_vprob) continue;

	  
	std::vector< RefCountedKinematicParticle > tau_children = tauTree->finalStateParticles();
	  
	math::PtEtaPhiMLorentzVector tau1_fit = daughter_p4(tau_children, 0);
	math::PtEtaPhiMLorentzVector tau2_fit = daughter_p4(tau_children, 1);
	math::PtEtaPhiMLorentzVector tau3_fit = daughter_p4(tau_children, 2);
	//
	//	  math::PtEtaPhiMLorentzVector tlv_tau = tau1_fit + tau2_fit + tau3_fit;

	particle_cand Taucand; 
	Taucand = calculateIPvariables(extrapolator, tau_part, tau_vertex, closestVertex);

	//	  if(Taucand.fls3d < 3) continue;

	// 6.1.2020 commented out
	if(Taucand.fls3d < c_fsig) continue;



	//	  
	//	  Float_t vprob_3 = -9;
	//	  TransientVertex vertex_3;
	//	  
	//	  std::vector<reco::TransientTrack> transient_tracks; 
	//	  transient_tracks.push_back(mytracks[iii]);
	//	  transient_tracks.push_back(mytracks[jjj]);
	//	  transient_tracks.push_back(mytracks[kkk]);
	//	  
	//	  std::tie(vprob_3, vertex_3) = vertexProb(transient_tracks);
	//
	//	  if(vprob_3 == -9) continue;
	  


	std::vector<RefCountedKinematicParticle> allParticles;
	allParticles.push_back(pFactory.particle(mytracks[iii], pion_mass, chi, ndf, pion_sigma));
	allParticles.push_back(pFactory.particle(mytracks[jjj], pion_mass, chi, ndf, pion_sigma));
	allParticles.push_back(pFactory.particle(mytracks[kkk], pion_mass, chi, ndf, pion_sigma));

	allParticles.push_back(Ds_part_rec);

	RefCountedKinematicTree bcTree = kpvFitter.fit(allParticles);

	if(bcTree->isEmpty() || !bcTree->isValid() || !bcTree->isConsistent()) continue;

	RefCountedKinematicParticle bc_part = bcTree->currentParticle();
	if(!bc_part->currentState().isValid()) continue;


	RefCountedKinematicVertex bc_vertex = bcTree->currentDecayVertex();
	if(!bc_vertex->vertexIsValid()) continue;
 
	//	  if(TMath::Prob(bc_vertex->chiSquared(), bc_vertex->degreesOfFreedom()) <=0.1) continue; 


	//	  const auto& Taucand = bc_children.at(0); 
	//	  const auto& JPcand = bc_children.at(1); 

	particle_cand Bcand; 
	Bcand = calculateIPvariables(extrapolator, bc_part, bc_vertex, closestVertex);

	  
	std::vector< RefCountedKinematicParticle > bc_children = bcTree->finalStateParticles();

	//	  const auto& state = fitted_children.at(i)->currentState();

	math::PtEtaPhiMLorentzVector tt1_fit = daughter_p4(bc_children, 0);
	math::PtEtaPhiMLorentzVector tt2_fit = daughter_p4(bc_children, 1);
	math::PtEtaPhiMLorentzVector tt3_fit = daughter_p4(bc_children, 2);

	math::PtEtaPhiMLorentzVector tlv_tau_fit = tt1_fit + tt2_fit + tt3_fit;

	// calculation of the isolation 

	if(tlv_tau_fit.Pt() < 3.) continue;


	Float_t iso = 0;
	Int_t ntracks = 0;
	Float_t iso_mindoca = 999; 
  
	  
	for(int itrk = 0; itrk < numOfch;  itrk++){
    
	  if(itrk==iii || itrk==jjj || itrk==kkk) continue;
	  if(itrk==kidx || itrk==pidx || itrk==sidx) continue;

	  iso += pfcollection[itrk].pt();

	  TrajectoryStateOnSurface tsos_pf = extrapolator.extrapolate(mytracks[itrk].impactPointState(), bc_vertex->position());
     
    
	  VertexDistance3D a3d_pf;  

	  std::pair<bool,Measurement1D> cur3DIP_pf = BsDstarTauNuNtuplizer::absoluteImpactParameter(tsos_pf, bc_vertex, a3d_pf);

	  Float_t pvip_pf = cur3DIP_pf.second.value();

	  //    std::cout << itrk << ": Distance of closest apporach to the bc vertex = " << pvip << std::endl;
    
	  if(pvip_pf < 0.03) ntracks+=1;

	  if(iso_mindoca > pvip_pf) iso_mindoca = pvip_pf;
	}



	Float_t max_dr_3prong = -1;

	Float_t dR_12 = reco::deltaR(tau1_fit.Eta(), tau1_fit.Phi(), tau2_fit.Eta(), tau2_fit.Phi());
	Float_t dR_13 = reco::deltaR(tau1_fit.Eta(), tau1_fit.Phi(), tau3_fit.Eta(), tau3_fit.Phi());
	Float_t dR_23 = reco::deltaR(tau2_fit.Eta(), tau2_fit.Phi(), tau3_fit.Eta(), tau3_fit.Phi());

	if(max_dr_3prong < dR_12) max_dr_3prong = dR_12;
	if(max_dr_3prong < dR_13) max_dr_3prong = dR_13;
	if(max_dr_3prong < dR_23) max_dr_3prong = dR_23;


	//	std::cout << "\t tau1 (pt, eta, phi) = "  << tau1_fit.Pt()  << "  " << tau1_fit.Eta() <<  " " << tau1_fit.Phi() << std::endl;
	//	std::cout << "\t tau2 (pt, eta, phi) = "  << tau2_fit.Pt()  << "  " << tau2_fit.Eta() <<  " " << tau2_fit.Phi() << std::endl;
	//	std::cout << "\t tau3 (pt, eta, phi) = "  << tau3_fit.Pt()  << "  " << tau3_fit.Eta() <<  " " << tau3_fit.Phi() << std::endl;
	// check gen. info. 


	Bool_t isRight = false; 
	Bool_t isRight1 = false; 
	Bool_t isRight2 = false; 
	Bool_t isRight3 = false; 
	Float_t dr1 = 999;
	Float_t dr2 = 999;
	Float_t dr3 = 999;
	Float_t ptres1 = 999;
	Float_t ptres2 = 999;
	Float_t ptres3 = 999;

	Int_t pid = -999;
	Float_t matched_gentaupt = -999;

	
	if(isMC){

	  for(unsigned int mmm=0; mmm < gps.size(); mmm++){
	    
	    Bool_t isRight1_ = false;
	    Bool_t isRight2_ = false;
	    Bool_t isRight3_ = false;
	    
	    std::vector<TLorentzVector> tlvs = gps[mmm];
	    
	    for(unsigned int nnn=0; nnn < tlvs.size(); nnn++){

	      if(
		 reco::deltaR(tau1_fit.Eta(), tau1_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.015 &&
		 tau1_fit.Pt()/tlvs[nnn].Pt() > 0.85 && 
		 tau1_fit.Pt()/tlvs[nnn].Pt() < 1.15
		 ){
		isRight1_ = true; 
		dr1 = reco::deltaR(tau1_fit.Eta(), tau1_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi());
		ptres1 = tau1_fit.Pt()/tlvs[nnn].Pt();
		//		std::cout << "\t\t\t gen. Nr" << mmm << " (" << tlvs[nnn].Pt() << ", " << tlvs[nnn].Eta()  << ", " << tlvs[nnn].Phi() << ") is matched with PF index = " << iii << ", with (pt,eta,phi) = " << tau1_fit.Pt() << ", " << tau1_fit.Eta() << ", " << tau1_fit.Phi() <<  std::endl;
		//		std::cout << "\t\t\t fitted pt respoinse = " << (tau1_fit.Pt()/tlvs[nnn].Pt()) << ", delta_R = " << reco::deltaR(tau1_fit.Eta(), tau1_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi())  << std::endl;
		//		std::cout << "\t\t\t original pt respoinse = " << (pf1.pt()/tlvs[nnn].Pt()) << ", delta_R = " << reco::deltaR(pf1.eta(), pf1.phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi())  << std::endl;
	      }
	      
	      if(
		 reco::deltaR(tau2_fit.Eta(), tau2_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.015 &&
		 tau2_fit.Pt()/tlvs[nnn].Pt() > 0.85 && 
		 tau2_fit.Pt()/tlvs[nnn].Pt() < 1.15
		 ){
		isRight2_ = true; 
		dr2 = reco::deltaR(tau2_fit.Eta(), tau2_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi());
		ptres2 = tau2_fit.Pt()/tlvs[nnn].Pt(); 
		//		std::cout << "\t\t\t gen. Nr" << mmm << " (" << tlvs[nnn].Pt() << ", " << tlvs[nnn].Eta()  << ", " << tlvs[nnn].Phi() << ") is matched with PF index = " << jjj << ", with (pt,eta,phi) = " << tau2_fit.Pt() << ", " << tau2_fit.Eta() << ", " << tau2_fit.Phi() <<  std::endl;
		//		std::cout << "\t\t\t fitted pt respoinse = " << (tau2_fit.Pt()/tlvs[nnn].Pt()) << ", delta_R = " << reco::deltaR(tau2_fit.Eta(), tau2_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi())  << std::endl;
		//		std::cout << "\t\t\t original pt respoinse = " << (pf2.pt()/tlvs[nnn].Pt()) << ", delta_R = " << reco::deltaR(pf2.eta(), pf2.phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi())  << std::endl;
	      }

	      
	      if(
		 reco::deltaR(tau3_fit.Eta(), tau3_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.015 &&
		 tau3_fit.Pt()/tlvs[nnn].Pt() > 0.85 && 
		 tau3_fit.Pt()/tlvs[nnn].Pt() < 1.15
		 ){
		isRight3_ = true; 
		dr3 = reco::deltaR(tau3_fit.Eta(), tau3_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi());
		ptres3 = tau3_fit.Pt()/tlvs[nnn].Pt(); 
		//		std::cout << "bgen. Nr" << mmm << " is matched with PF index = " << kkk << std::endl;

		//		std::cout << "\t\t\t gen. Nr" << mmm << " (" << tlvs[nnn].Pt() << ", " << tlvs[nnn].Eta()  << ", " << tlvs[nnn].Phi() << ") is matched with PF index = " << kkk << ", with (pt,eta,phi) = " << tau3_fit.Pt() << ", " << tau3_fit.Eta() << ", " << tau3_fit.Phi() <<  std::endl;
		//		std::cout << "\t\t\t fitted pt respoinse = " << (tau3_fit.Pt()/tlvs[nnn].Pt()) << ", delta_R = " << reco::deltaR(tau3_fit.Eta(), tau3_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi())  << std::endl;
		//		std::cout << "\t\t\t original pt respoinse = " << (pf3.pt()/tlvs[nnn].Pt()) << ", delta_R = " << reco::deltaR(pf3.eta(), pf3.phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi())  << std::endl;

	      }

	      
	    }
	    
	    Bool_t isRight_ = isRight1_ && isRight2_ && isRight3_;
	    if(isRight1_) isRight1 = true;
	    if(isRight2_) isRight2 = true;
	    if(isRight3_) isRight3 = true;

	    if(isRight_){
	      isRight = true;
	      pid = ppdgId[mmm];
	      matched_gentaupt = vec_gentau3pp4[mmm].Pt();
	    }
	  }	
	}
	

	Float_t sumofdnn = -1;
	if(useDNN_){
	  sumofdnn = mydnn[iii] + mydnn[jjj] + mydnn[kkk];
	  //	  std::cout << "dnn output: "<< sumofdnn << " " <<  mydnn[iii] << " " << mydnn[jjj] << " " << mydnn[kkk] << std::endl;
	}



	taucand _cand_ = {
	  iii,
	  jjj,
	  kkk,
	  (Float_t) tlv_tau_fit.Pt(),
	  (Float_t) tlv_tau_fit.Eta(),
	  (Float_t) tlv_tau_fit.Phi(),
	  (Float_t) tlv_tau_fit.M(),
	  (Float_t) Taucand.lip, 
	  (Float_t) Taucand.lips, 
	  (Float_t) Taucand.pvip, 
	  (Float_t) Taucand.pvips, 
	  (Float_t) Taucand.fl3d,
	  (Float_t) Taucand.fls3d, 
	  (Float_t) Taucand.alpha,
	  (Float_t) TMath::Prob(tau_vertex->chiSquared(), tau_vertex->degreesOfFreedom()),
	  (Float_t) tau_vertex->vertexState().position().x(), 
	  (Float_t) tau_vertex->vertexState().position().y(), 
	  (Float_t) tau_vertex->vertexState().position().z(), 
	  (Float_t) max_dr_3prong, 
	  (Int_t) tau_charge,
	  (Bool_t) isRight,
	  (Bool_t) isRight1,
	  (Bool_t) isRight2,
	  (Bool_t) isRight3,
	  (Float_t) dr1,
	  (Float_t) dr2,
	  (Float_t) dr3,
	  (Float_t) ptres1,
	  (Float_t) ptres2,
	  (Float_t) ptres3,
	  (Int_t) pid,
	  (Float_t) matched_gentaupt, 
	  (Float_t) sumofdnn,
	  (Float_t) mydnn[iii],
	  (Float_t) mydnn[jjj],
	  (Float_t) mydnn[kkk],
	  (Float_t) TMath::Prob(bc_vertex->chiSquared(), bc_vertex->degreesOfFreedom()),
	  (Float_t) bc_vertex->vertexState().position().x(),
	  (Float_t) bc_vertex->vertexState().position().y(),
	  (Float_t) bc_vertex->vertexState().position().z(),
	  (Float_t) bc_part->currentState().globalMomentum().perp(),
	  (Float_t) bc_part->currentState().globalMomentum().eta(),
	  (Float_t) bc_part->currentState().globalMomentum().phi(),
	  (Float_t) bc_part->currentState().mass(),
	  (Float_t) Bcand.lip, 
	  (Float_t) Bcand.lips, 
	  (Float_t) Bcand.pvip, 
	  (Float_t) Bcand.pvips, 
	  (Float_t) Bcand.fl3d,
	  (Float_t) Bcand.fls3d, 
	  (Float_t) Bcand.alpha,
	  (Float_t) iso,
	  (Float_t) ntracks,
	  (Float_t) iso_mindoca,
	};
	  
	//	std::cout << cands.size() << std::endl;
	cands.push_back(_cand_);
      }
    }
  }

  sort(cands.begin(), cands.end());

  //  std::vector<Int_t> dict_idx;
  Int_t ncomb = 0;

//  bool isRight_bS = false;
//  bool isRight_aS = false;
//  int isRight_bS_ith = -1;
//  int isRight_aS_ith = -1;
//  int isRight_bS_n = 0;
//  int isRight_aS_n = 0;

  if(cands.size()==0) return false;
  nBranches_->cutflow_perevt->Fill(5);

  //    std::cout << "step1" << std::endl;

  for(int ic=0; ic < (int)cands.size() && ic < 1; ic++){
      
    Int_t _idx1 = cands[ic].cand_tau_id1;
    Int_t _idx2 = cands[ic].cand_tau_id2;
    Int_t _idx3 = cands[ic].cand_tau_id3;
    //    
    //      bool flag_overlap = false;
    //      for(int idc=0; idc<(int) dict_idx.size(); idc++){
    //	
    //	if(_idx1 == dict_idx[idc] || 
    //	   _idx2 == dict_idx[idc] || 
    //	   _idx3 == dict_idx[idc])
    //	  
    //	  flag_overlap = true;
    //      }
    //
    //      if(cands[ic].cand_tau_isRight==true){
    //	isRight_bS = true;
    //	isRight_bS_ith = ic;
    //	isRight_bS_n += 1;
    //      }
    //      
    //      
    //      if(flag_overlap) continue; 
    //
    //      if(cands[ic].cand_tau_isRight==true){
    //	isRight_aS = true;
    //	isRight_aS_ith = ncomb;
    //	isRight_aS_n += 1;
    //      }
      
    //    std::cout << "\t ----> kept ! " << std::endl;
    //      dict_idx.push_back(_idx1);
    //      dict_idx.push_back(_idx2);
    //      dict_idx.push_back(_idx3);
      
    ncomb += 1;


    /********************************************************************
     *
     * Step9: Filling normal branches
     *
     ********************************************************************/

    // if the candidate has less than 2 GeV, remove it.
    //      if(cands[ic].cand_tau_pt < 2.) continue;


    nBranches_->BsDstarTauNu_tau_pt.push_back(cands[ic].cand_tau_pt);
    nBranches_->BsDstarTauNu_tau_eta.push_back(cands[ic].cand_tau_eta);
    nBranches_->BsDstarTauNu_tau_phi.push_back(cands[ic].cand_tau_phi);
    nBranches_->BsDstarTauNu_tau_mass.push_back(cands[ic].cand_tau_mass);

    std::vector<Float_t> rhomass;
    pat::PackedCandidate pf1 = pfcollection[_idx1];
    pat::PackedCandidate pf2 = pfcollection[_idx2];
    pat::PackedCandidate pf3 = pfcollection[_idx3];
      
    TLorentzVector tlv_pion1; 
    TLorentzVector tlv_pion2;
    TLorentzVector tlv_pion3;
      
    tlv_pion1.SetPtEtaPhiM(pf1.pt(), pf1.eta(), pf1.phi(), pf1.mass());
    tlv_pion2.SetPtEtaPhiM(pf2.pt(), pf2.eta(), pf2.phi(), pf2.mass());
    tlv_pion3.SetPtEtaPhiM(pf3.pt(), pf3.eta(), pf3.phi(), pf3.mass());
      
      
    if(pf1.charge()*pf2.charge() == -1){
      //	std::cout << "this is 1" << std::endl;
      TLorentzVector tlv_rho = tlv_pion1 + tlv_pion2;
      rhomass.push_back(tlv_rho.M());
    }
      
    if(pf1.charge()*pf3.charge() == -1){
      //	std::cout << "this is 2" << std::endl;
      TLorentzVector tlv_rho = tlv_pion1 + tlv_pion3;
      rhomass.push_back(tlv_rho.M());
    }
      
    if(pf2.charge()*pf3.charge() == -1){
      //	std::cout << "this is 3" << std::endl;
      TLorentzVector tlv_rho = tlv_pion2 + tlv_pion3;
      rhomass.push_back(tlv_rho.M());
    }
      
    //      std::cout << "rho masses size = " << rhomass.size() << " "  << rhomass.at(0) << " " << rhomass.at(1) << std::endl;
      
    nBranches_->BsDstarTauNu_tau_rhomass1.push_back(rhomass.at(0));
    nBranches_->BsDstarTauNu_tau_rhomass2.push_back(rhomass.at(1));
      

    //      nBranches_->BsDstarTauNu_tau_rhomass1.push_back(-1);
    //      nBranches_->BsDstarTauNu_tau_rhomass2.push_back(-1);

    nBranches_->BsDstarTauNu_tau_q.push_back(cands[ic].cand_tau_charge);
    nBranches_->BsDstarTauNu_tau_vx.push_back(cands[ic].cand_tau_vx);
    nBranches_->BsDstarTauNu_tau_vy.push_back(cands[ic].cand_tau_vy);
    nBranches_->BsDstarTauNu_tau_vz.push_back(cands[ic].cand_tau_vz);

    nBranches_->BsDstarTauNu_tau_max_dr_3prong.push_back(cands[ic].cand_tau_max_dr_3prong);
    nBranches_->BsDstarTauNu_tau_lip.push_back(cands[ic].cand_tau_lip);
    nBranches_->BsDstarTauNu_tau_lips.push_back(cands[ic].cand_tau_lips);
    nBranches_->BsDstarTauNu_tau_pvip.push_back(cands[ic].cand_tau_pvip);
    nBranches_->BsDstarTauNu_tau_pvips.push_back(cands[ic].cand_tau_pvips);
    nBranches_->BsDstarTauNu_tau_fl3d.push_back(cands[ic].cand_tau_fl3d);
    nBranches_->BsDstarTauNu_tau_fls3d.push_back(cands[ic].cand_tau_fls3d);
    nBranches_->BsDstarTauNu_tau_alpha.push_back(cands[ic].cand_tau_alpha);
    nBranches_->BsDstarTauNu_tau_vprob.push_back(cands[ic].cand_tau_vprob);
    nBranches_->BsDstarTauNu_tau_isRight.push_back(cands[ic].cand_tau_isRight);
    nBranches_->BsDstarTauNu_tau_matched_ppdgId.push_back(cands[ic].cand_tau_matched_ppdgId);
    nBranches_->BsDstarTauNu_tau_matched_gentaupt.push_back(cands[ic].cand_tau_matched_gentaupt);
    nBranches_->BsDstarTauNu_tau_sumofdnn.push_back(cands[ic].cand_tau_sumofdnn);
    nBranches_->BsDstarTauNu_tau_pfidx1.push_back(cands[ic].cand_tau_id1);
    nBranches_->BsDstarTauNu_tau_pfidx2.push_back(cands[ic].cand_tau_id2);
    nBranches_->BsDstarTauNu_tau_pfidx3.push_back(cands[ic].cand_tau_id3);


    nBranches_->BsDstarTauNu_tau_pi1_pt.push_back(tlv_pion1.Pt());
    nBranches_->BsDstarTauNu_tau_pi1_eta.push_back(tlv_pion1.Eta());
    nBranches_->BsDstarTauNu_tau_pi1_phi.push_back(tlv_pion1.Phi());
    nBranches_->BsDstarTauNu_tau_pi1_mass.push_back(tlv_pion1.M());
      
    nBranches_->BsDstarTauNu_tau_pi2_pt.push_back(tlv_pion2.Pt());
    nBranches_->BsDstarTauNu_tau_pi2_eta.push_back(tlv_pion2.Eta());
    nBranches_->BsDstarTauNu_tau_pi2_phi.push_back(tlv_pion2.Phi());
    nBranches_->BsDstarTauNu_tau_pi2_mass.push_back(tlv_pion2.M());
      
    nBranches_->BsDstarTauNu_tau_pi3_pt.push_back(tlv_pion3.Pt());
    nBranches_->BsDstarTauNu_tau_pi3_eta.push_back(tlv_pion3.Eta());
    nBranches_->BsDstarTauNu_tau_pi3_phi.push_back(tlv_pion3.Phi());
    nBranches_->BsDstarTauNu_tau_pi3_mass.push_back(tlv_pion3.M());


    nBranches_->BsDstarTauNu_B_pt.push_back(cands[ic].cand_b_pt);
    nBranches_->BsDstarTauNu_B_eta.push_back(cands[ic].cand_b_eta);
    nBranches_->BsDstarTauNu_B_phi.push_back(cands[ic].cand_b_phi);
    nBranches_->BsDstarTauNu_B_mass.push_back(cands[ic].cand_b_mass);
    nBranches_->BsDstarTauNu_B_vprob.push_back(cands[ic].cand_b_vprob);
    nBranches_->BsDstarTauNu_B_lip.push_back(cands[ic].cand_b_lip);
    nBranches_->BsDstarTauNu_B_lips.push_back(cands[ic].cand_b_lips);
    nBranches_->BsDstarTauNu_B_pvip.push_back(cands[ic].cand_b_pvip);
    nBranches_->BsDstarTauNu_B_pvips.push_back(cands[ic].cand_b_pvips);
    nBranches_->BsDstarTauNu_B_fls3d.push_back(cands[ic].cand_b_fls3d);
    nBranches_->BsDstarTauNu_B_fl3d.push_back(cands[ic].cand_b_fl3d);
    nBranches_->BsDstarTauNu_B_alpha.push_back(cands[ic].cand_b_alpha);

    std::vector<RefCountedKinematicParticle> allParticles4doc;
      
    allParticles4doc.push_back(pFactory.particle(tt_muon, muon_mass, chi, ndf, muon_sigma));
    allParticles4doc.push_back(pFactory.particle(mytracks[cands[ic].cand_tau_id1], pion_mass, chi, ndf, pion_sigma));
    allParticles4doc.push_back(pFactory.particle(mytracks[cands[ic].cand_tau_id2], pion_mass, chi, ndf, pion_sigma));
    allParticles4doc.push_back(pFactory.particle(mytracks[cands[ic].cand_tau_id3], pion_mass, chi, ndf, pion_sigma));


    nBranches_->BsDstarTauNu_B_maxdoca.push_back(getMaxDoca(allParticles4doc));
    nBranches_->BsDstarTauNu_B_mindoca.push_back(getMinDoca(allParticles4doc));
    nBranches_->BsDstarTauNu_B_vx.push_back(cands[ic].cand_b_vx);
    nBranches_->BsDstarTauNu_B_vy.push_back(cands[ic].cand_b_vy);
    nBranches_->BsDstarTauNu_B_vz.push_back(cands[ic].cand_b_vz);

    nBranches_->BsDstarTauNu_B_iso.push_back(cands[ic].cand_b_iso);
    nBranches_->BsDstarTauNu_B_iso_ntracks.push_back(cands[ic].cand_b_iso_ntracks);
    nBranches_->BsDstarTauNu_B_iso_mindoca.push_back(cands[ic].cand_b_iso_mindoca);


    TLorentzVector Tlv_B;
    TLorentzVector Tlv_Ds;
    TLorentzVector Tlv_tau;

    nBranches_->BsDstarTauNu_B_pt.push_back(cands[ic].cand_b_pt);
    nBranches_->BsDstarTauNu_B_eta.push_back(cands[ic].cand_b_eta);
    nBranches_->BsDstarTauNu_B_phi.push_back(cands[ic].cand_b_phi);
    nBranches_->BsDstarTauNu_B_mass.push_back(cands[ic].cand_b_mass);

    Tlv_B.SetPtEtaPhiM(cands[ic].cand_b_pt,
		       cands[ic].cand_b_eta,
		       cands[ic].cand_b_phi,
		       cands[ic].cand_b_mass);
    
    Tlv_Ds.SetPtEtaPhiM(Ds_part_rec->currentState().globalMomentum().perp(),
			Ds_part_rec->currentState().globalMomentum().eta(),
			Ds_part_rec->currentState().globalMomentum().phi(),
			Ds_part_rec->currentState().mass());


    Tlv_tau.SetPtEtaPhiM(cands[ic].cand_tau_pt, 
			 cands[ic].cand_tau_eta, 
			 cands[ic].cand_tau_phi, 
			 cands[ic].cand_tau_mass);


    Float_t mm2 = (Tlv_B - Tlv_Ds - Tlv_tau).M2();
    //    std::cout << "check:" << (Tlv_B - Tlv_Ds - Tlv_tau).M() << " " << mm2 << std::endl;

    Float_t q2 = (Tlv_B - Tlv_Ds).M2();
    nBranches_->BsDstarTauNu_B_mm2.push_back(mm2);
    nBranches_->BsDstarTauNu_B_q2.push_back(q2);

    //    std::cout << "before:" << Tlv_tau.Px() << " " << Tlv_tau.Py() << " " << Tlv_tau.Pz() << " " << Tlv_tau.M() << " " << Tlv_tau.E() << std::endl;
    Tlv_tau.Boost( -Tlv_B.BoostVector() );
    //    std::cout << "after:" << Tlv_tau.Px() << " " << Tlv_tau.Py() << " " << Tlv_tau.Pz() << " " << Tlv_tau.M() << " " << Tlv_tau.E()<< std::endl;

    nBranches_->BsDstarTauNu_B_Es.push_back(Tlv_tau.E()); 

    Float_t ptback = cands[ic].cand_b_pt*mass_B0/cands[ic].cand_b_mass;
    nBranches_->BsDstarTauNu_B_ptback.push_back(ptback); 

  
      
  }





  nBranches_->BsDstarTauNu_mu1_pt.push_back(muoncollection[0].pt());
  nBranches_->BsDstarTauNu_mu1_eta.push_back(muoncollection[0].eta());
  nBranches_->BsDstarTauNu_mu1_phi.push_back(muoncollection[0].phi());
  nBranches_->BsDstarTauNu_mu1_mass.push_back(muoncollection[0].mass());
  nBranches_->BsDstarTauNu_mu1_q.push_back(muoncollection[0].charge());
  nBranches_->BsDstarTauNu_mu1_isLoose.push_back(muoncollection[0].isLooseMuon());
  nBranches_->BsDstarTauNu_mu1_isTight.push_back(muoncollection[0].isTightMuon(closestVertex));
  nBranches_->BsDstarTauNu_mu1_isPF.push_back(muoncollection[0].isPFMuon());
  nBranches_->BsDstarTauNu_mu1_isGlobal.push_back(muoncollection[0].isGlobalMuon());
  nBranches_->BsDstarTauNu_mu1_isTracker.push_back(muoncollection[0].isTrackerMuon());
  nBranches_->BsDstarTauNu_mu1_isSoft.push_back(muoncollection[0].isSoftMuon(closestVertex));
  nBranches_->BsDstarTauNu_mu1_vx.push_back(muoncollection[0].vx());
  nBranches_->BsDstarTauNu_mu1_vy.push_back(muoncollection[0].vy());
  nBranches_->BsDstarTauNu_mu1_vz.push_back(muoncollection[0].vz());
  nBranches_->BsDstarTauNu_mu1_iso.push_back(1.);
  nBranches_->BsDstarTauNu_mu1_dbiso.push_back(MuonPFIso(muoncollection[0]));
  
  nBranches_->BsDstarTauNu_PV_vx.push_back(vertices_->begin()->position().x());
  nBranches_->BsDstarTauNu_PV_vy.push_back(vertices_->begin()->position().y());
  nBranches_->BsDstarTauNu_PV_vz.push_back(vertices_->begin()->position().z());

  //      if(myVertex.isValid()){
  //	nBranches_->BsDstarTauNu_bbPV_refit_vx.push_back(myVertex.position().x());
  //	nBranches_->BsDstarTauNu_bbPV_refit_vy.push_back(myVertex.position().y());
  //	nBranches_->BsDstarTauNu_bbPV_refit_vz.push_back(myVertex.position().z());
  //      }else{
  nBranches_->BsDstarTauNu_bbPV_refit_vx.push_back(-1);
  nBranches_->BsDstarTauNu_bbPV_refit_vy.push_back(-1);
  nBranches_->BsDstarTauNu_bbPV_refit_vz.push_back(-1);
  //      }

  nBranches_->BsDstarTauNu_bbPV_vx.push_back(closestVertex.position().x());
  nBranches_->BsDstarTauNu_bbPV_vy.push_back(closestVertex.position().y());
  nBranches_->BsDstarTauNu_bbPV_vz.push_back(closestVertex.position().z());


  //////////////////////////////

  nBranches_->BsDstarTauNu_D0_pt.push_back(D0_part_rec->currentState().globalMomentum().perp());
  nBranches_->BsDstarTauNu_D0_eta.push_back(D0_part_rec->currentState().globalMomentum().eta());
  nBranches_->BsDstarTauNu_D0_phi.push_back(D0_part_rec->currentState().globalMomentum().phi());
  nBranches_->BsDstarTauNu_D0_mass.push_back(D0_part_rec->currentState().mass());
  nBranches_->BsDstarTauNu_D0_vprob.push_back(TMath::Prob(D0_part_rec->chiSquared(), D0_part_rec->degreesOfFreedom()));
  nBranches_->BsDstarTauNu_D0_lip.push_back(D0cand_rec.lip);
  nBranches_->BsDstarTauNu_D0_lips.push_back(D0cand_rec.lips);
  nBranches_->BsDstarTauNu_D0_pvip.push_back(D0cand_rec.pvip);
  nBranches_->BsDstarTauNu_D0_pvips.push_back(D0cand_rec.pvips);
  nBranches_->BsDstarTauNu_D0_fl3d.push_back(D0cand_rec.fl3d);
  nBranches_->BsDstarTauNu_D0_fls3d.push_back(D0cand_rec.fls3d);
  nBranches_->BsDstarTauNu_D0_alpha.push_back(D0cand_rec.alpha);
  nBranches_->BsDstarTauNu_D0_vx.push_back(D0_vertex_rec->vertexState().position().x());
  nBranches_->BsDstarTauNu_D0_vy.push_back(D0_vertex_rec->vertexState().position().y());
  nBranches_->BsDstarTauNu_D0_vz.push_back(D0_vertex_rec->vertexState().position().z());  
  nBranches_->BsDstarTauNu_D0_unfit_pt.push_back(D0_tlv_highest.Pt());
  nBranches_->BsDstarTauNu_D0_unfit_mass.push_back(D0_tlv_highest.M());
  nBranches_->BsDstarTauNu_D0_ptfrac.push_back(D0_part_rec->currentState().globalMomentum().perp()/ptsum);

  nBranches_->BsDstarTauNu_Ds_pt.push_back(Ds_part_rec->currentState().globalMomentum().perp());
  nBranches_->BsDstarTauNu_Ds_eta.push_back(Ds_part_rec->currentState().globalMomentum().eta());
  nBranches_->BsDstarTauNu_Ds_phi.push_back(Ds_part_rec->currentState().globalMomentum().phi());
  nBranches_->BsDstarTauNu_Ds_mass.push_back(Ds_part_rec->currentState().mass());
  nBranches_->BsDstarTauNu_Ds_vprob.push_back(TMath::Prob(Ds_part_rec->chiSquared(), Ds_part_rec->degreesOfFreedom()));
  nBranches_->BsDstarTauNu_Ds_lip.push_back(Dscand_rec.lip);
  nBranches_->BsDstarTauNu_Ds_lips.push_back(Dscand_rec.lips);
  nBranches_->BsDstarTauNu_Ds_pvip.push_back(Dscand_rec.pvip);
  nBranches_->BsDstarTauNu_Ds_pvips.push_back(Dscand_rec.pvips);
  nBranches_->BsDstarTauNu_Ds_fl3d.push_back(Dscand_rec.fl3d);
  nBranches_->BsDstarTauNu_Ds_fls3d.push_back(Dscand_rec.fls3d);
  nBranches_->BsDstarTauNu_Ds_alpha.push_back(Dscand_rec.alpha);
  nBranches_->BsDstarTauNu_Ds_vx.push_back(Ds_vertex_rec->vertexState().position().x());
  nBranches_->BsDstarTauNu_Ds_vy.push_back(Ds_vertex_rec->vertexState().position().y());
  nBranches_->BsDstarTauNu_Ds_vz.push_back(Ds_vertex_rec->vertexState().position().z());  
  nBranches_->BsDstarTauNu_Ds_unfit_pt.push_back(Ds_tlv_highest.Pt());
  nBranches_->BsDstarTauNu_Ds_unfit_mass.push_back(Ds_tlv_highest.M());
  nBranches_->BsDstarTauNu_Ds_ptfrac.push_back(Ds_part_rec->currentState().globalMomentum().perp()/ptsum);

  nBranches_->BsDstarTauNu_k_charge.push_back(qk);
  nBranches_->BsDstarTauNu_pi_charge.push_back(qp);
  nBranches_->BsDstarTauNu_spi_charge.push_back(qsp);

  /********************************************************************
   *
   * Step10: check gen-matching and fill them
   *
   ********************************************************************/

  TVector3 genvertex(-9.,-9.,-9.);
  
  std::vector<const reco::Candidate*> gen_nr_mu;
  

  if(isMC){
    
    for( unsigned p=0; p < genParticles_->size(); ++p){

      //      std::cout << "gen: " << (*genParticles_)[p].pdgId() << " " << (*genParticles_)[p].status() << std::endl;

      // Bc daughters loop
      if(TMath::Abs((*genParticles_)[p].pdgId())==511 && (*genParticles_)[p].status()==2){

	// retrieve production vertex
	genvertex = getVertex((*genParticles_)[p]);
	    
	for(int idd = 0; idd < (int)(*genParticles_)[p].numberOfDaughters(); idd++){
	  Int_t dpid = (*genParticles_)[p].daughter(idd)->pdgId();
	  //	  std::cout << "\t -> " <<  << " " << (*genParticles_)[p].daughter(idd)->status()<< std::endl;
	  if(TMath::Abs(dpid)==13) gen_nr_mu.push_back((*genParticles_)[p].daughter(idd));
	}
      }
	  
	  
    }
  }
      
  

  
  // -9 if there is no Bc found 
  nBranches_->BsDstarTauNu_genPV_vx.push_back(genvertex.x());
  nBranches_->BsDstarTauNu_genPV_vy.push_back(genvertex.y());
  nBranches_->BsDstarTauNu_genPV_vz.push_back(genvertex.z());
  nBranches_->BsDstarTauNu_ngenmuons.push_back(gen_nr_mu.size());





  nBranches_->BsDstarTauNu_isgen3.push_back(isgen3);
  nBranches_->BsDstarTauNu_isgen3matched.push_back(isgen3matched);
  nBranches_->BsDstarTauNu_nch.push_back(numOfch);
  nBranches_->BsDstarTauNu_nch_after_dnn.push_back(npf_after_dnn);
  nBranches_->BsDstarTauNu_nch_before_dnn.push_back(npf_before_dnn);
  nBranches_->BsDstarTauNu_nch_qr.push_back(npf_qr);
  nBranches_->BsDstarTauNu_ngentau3.push_back(gps.size());
  nBranches_->BsDstarTauNu_ngentau.push_back(vec_gentaudm.size());

  if(vec_gentaudm.size() >=1){
    nBranches_->BsDstarTauNu_gentaupt.push_back(vec_gentaup4[0].Pt());
    nBranches_->BsTauTau_gentaudm.push_back(vec_gentaudm[0]);
  }else{
    nBranches_->BsDstarTauNu_gentaupt.push_back(-1);
    nBranches_->BsTauTau_gentaudm.push_back(-1);
  }

  nBranches_->IsBsDstarTauNu.push_back(1.);
  nBranches_->BsDstarTauNu_nCandidates.push_back(ncomb);





  // B Candidate Kinematic fit passed

  return true;
  



}


void BsDstarTauNuNtuplizer::printout(const RefCountedKinematicVertex& myVertex){
  std::cout << "Vertex:" << std::endl;
  if (myVertex->vertexIsValid()) {
    std::cout << "\t Decay vertex: " << myVertex->position() << myVertex->chiSquared() << " " << myVertex->degreesOfFreedom()
	      << std::endl;
  } else
    std::cout << "\t Decay vertex Not valid\n";
}

void BsDstarTauNuNtuplizer::printout(const RefCountedKinematicParticle& myParticle){
  std::cout << "Particle:" << std::endl;
  //accessing the reconstructed Bs meson parameters:
  //SK: uncomment if needed  AlgebraicVector7 bs_par = myParticle->currentState().kinematicParameters().vector();

  //and their joint covariance matrix:
  //SK:uncomment if needed  AlgebraicSymMatrix77 bs_er = myParticle->currentState().kinematicParametersError().matrix();
  std::cout << "\t Momentum at vertex: " << myParticle->currentState().globalMomentum() << std::endl;
  std::cout << "\t Parameters at vertex: " << myParticle->currentState().kinematicParameters().vector() << std::endl;
}

void BsDstarTauNuNtuplizer::printout(const RefCountedKinematicTree& myTree){
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
