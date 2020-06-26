#include "../interface/BsTauTauFHNtuplizer_mr.h"


//===================================================================================================================
BsTauTauFHNtuplizer_mr::BsTauTauFHNtuplizer_mr( edm::EDGetTokenT<pat::MuonCollection>    muonToken   ,
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
  , isTruth_   (runFlags["isTruth"])
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
  }    


}

//===================================================================================================================
BsTauTauFHNtuplizer_mr::~BsTauTauFHNtuplizer_mr( void )
{

}


Int_t BsTauTauFHNtuplizer_mr::decaymode_id(std::string str){
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


TVector3 BsTauTauFHNtuplizer_mr::getVertex(const reco::GenParticle& part){
  return TVector3(part.vx(),part.vy(),part.vz());
}

float BsTauTauFHNtuplizer_mr::MuonPFIso(pat::Muon muon){

  float sumChargedHadronPt = muon.pfIsolationR04().sumChargedHadronPt;
  float sumNeutralHadronEt = muon.pfIsolationR04().sumNeutralHadronEt;
  float sumPhotonEt = muon.pfIsolationR04().sumPhotonEt;
  float sumPUPt = muon.pfIsolationR04().sumPUPt;
  float iso = (sumChargedHadronPt + std::max( 0. ,sumNeutralHadronEt + sumPhotonEt - 0.5 * sumPUPt));// / muon.pt()
 
  return iso;
}




Float_t BsTauTauFHNtuplizer_mr::getMaxDoca(std::vector<RefCountedKinematicParticle> &kinParticles){

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



Float_t BsTauTauFHNtuplizer_mr::getMinDoca(std::vector<RefCountedKinematicParticle> &kinParticles) {

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




std::tuple<Float_t, TransientVertex> BsTauTauFHNtuplizer_mr::vertexProb( const std::vector<reco::TransientTrack>& tracks){

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
std::pair<bool, Measurement1D> BsTauTauFHNtuplizer_mr::absoluteImpactParameter(const TrajectoryStateOnSurface& tsos,
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




particle_cand BsTauTauFHNtuplizer_mr::calculateIPvariables(
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


math::PtEtaPhiMLorentzVector BsTauTauFHNtuplizer_mr::daughter_p4(std::vector< RefCountedKinematicParticle > fitted_children, size_t i){
  const auto& state = fitted_children.at(i)->currentState();

  return math::PtEtaPhiMLorentzVector(
				      state.globalMomentum().perp(), 
				      state.globalMomentum().eta() ,
				      state.globalMomentum().phi() ,
				      state.mass()
				      );
}


bool BsTauTauFHNtuplizer_mr::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  
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


  //    std::cout << "number of matched muon = " << muoncollection.size() << std::endl;
  if(!( muoncollection.size() >= 1)) return false;

    
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
    
  //    std::cout << "number of matched muon = " << muoncollection.size() << std::endl;

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
    
  std::vector<pat::PackedCandidate> pfcollection; 
  //    std::vector<pat::PackedCandidate> pfcollection_pre; 
  std::vector<reco::TransientTrack> mytracks;
  std::vector<Float_t> mydnn;
    


  Int_t npf_before_dnn = 0;
  Int_t npf_qr = 0;


  if(useDNN_){

    //    std::cout << "check1" << std::endl;

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
	  //	  data.tensor<float, 3>()(0, count_dnn, 7) = pf.isGlobalMuon();
	  data.tensor<float, 3>()(0, count_dnn, 7) = (TMath::Abs(pf.pdgId())==13);
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
      npf_before_dnn++;
	  
      //	pfcollection.push_back(pf);
      //	reco::TransientTrack  tt_track = (*builder).build(pf.pseudoTrack());
      //	mytracks.push_back(tt_track);
	

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
	//	data.tensor<float, 3>()(0, count_dnn, 7) = pf.isGlobalMuon();
	data.tensor<float, 3>()(0, count_dnn, 7) = (TMath::Abs(pf.pdgId())==13);
	  
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
  Int_t numOfch = (size_t)pfcollection.size();


  std::vector<std::vector<TLorentzVector>> gps;
  std::vector<Int_t> ppdgId;
  std::vector<Int_t> vec_gentaudm;
  std::vector<Int_t> vec_ppdgId;
  std::vector<TLorentzVector> vec_gentaup4;
  std::vector<TLorentzVector> vec_gentau3pp4;
  std::vector<TLorentzVector> vec_gentau3pp4_full;
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
      TLorentzVector genvis_full;
      std::vector<TLorentzVector> gp;
      Bool_t matched = true;
      Int_t nprong = 0;
      
      for(int idd = 0; idd < (int)(*genParticles_)[p].numberOfDaughters(); idd++){
	
	std::cout << "\t\t -> " << (*genParticles_)[p].daughter(idd)->pdgId() << " (pT, eta, phi) = " 
		  << (*genParticles_)[p].daughter(idd)->pt() << " " 
		  << (*genParticles_)[p].daughter(idd)->eta() << " " 
		  << (*genParticles_)[p].daughter(idd)->phi() << std::endl;

	TLorentzVector _genvis_full;
	_genvis_full.SetPtEtaPhiM((*genParticles_)[p].daughter(idd)->pt(),
			      (*genParticles_)[p].daughter(idd)->eta(),
			      (*genParticles_)[p].daughter(idd)->phi(),
			      (*genParticles_)[p].daughter(idd)->mass()
			      );
	
	genvis_full += _genvis_full;


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
	  Int_t min_index = 999;
	  
	  for(int kkk = 0; kkk < numOfch; kkk ++){
	    
	    pat::PackedCandidate _pf = pfcollection[kkk];
	    
	    if(_pf.pdgId()!=(*genParticles_)[p].daughter(idd)->pdgId()) continue;
	    
	    Float_t _dR = reco::deltaR(
				       _genvis_.Eta(), _genvis_.Phi(),
				       _pf.eta(), _pf.phi()
				       );

	    //	    std::cout << "\t\t\t\t" << kkk  << ", (pt, eta, phi) = " << _pf.pt() << ", " << _pf.eta() << ", " << _pf.phi() << " with dR = " << _dR << " and reco/gen. pT = " << _pf.pt()/_genvis_.Pt() << std::endl;

	    if(_dR < min_dr && _dR < 0.015 && _pf.pt()/_genvis_.Pt() < 1.15 && _pf.pt()/_genvis_.Pt() > 0.85){
	      min_dr = _dR;
	      min_index = kkk; 
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
	  else std::cout << "\t\t\t -----> PF index = " << min_index << " matched!" << std::endl;

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
	vec_gentau3pp4_full.push_back(genvis_full);
	
	//	if(TMath::Abs((*genParticles_)[p].mother(0)->pdgId())==541){
	//	isgen3matched += matched;
	if(matched) isgen3matched += 1;
	//	}
      }//else{
	//	isgen3matched = false;
      //      }
    }

    //    std::cout << "\t # of gen. taus with 3prong = " << gps.size() << std::endl;

  }

  if(vec_gentau3pp4.size()==0) return false;

  //////////////////////////////
  std::cout << "Starts to build tau candidate out of " << numOfch << " pion candidates" << std::endl;

  std::vector<taucandgen> cands;
  Int_t npf_after_dnn = 0;
    
  //  std::cout << "check0" << std::endl;

  for(int iii = 0; iii < numOfch; iii ++){
      
    pat::PackedCandidate pf1 = pfcollection[iii];
    if(useDNN_==true && mydnn[iii] < c_dnn) continue;
    npf_after_dnn++;

    for(int jjj = iii+1; jjj < numOfch; jjj ++){
	
      pat::PackedCandidate pf2 = pfcollection[jjj];
      if(useDNN_==true && mydnn[jjj] < c_dnn) continue;

      for(int kkk = jjj+1; kkk < numOfch; kkk ++){

	pat::PackedCandidate pf3 = pfcollection[kkk];
	if(useDNN_==true && mydnn[kkk] < c_dnn) continue;

	

	//  std::cout << "check1" << std::endl;
	Int_t tau_charge = pf1.charge() + pf2.charge() + pf3.charge(); 

	if(TMath::Abs(tau_charge)!=1) continue; 

	//  std::cout << "check2" << std::endl;

///	std::vector<RefCountedKinematicParticle> tauParticles;
///
///	tauParticles.push_back(pFactory.particle(mytracks[iii], pion_mass, chi, ndf, pion_sigma));
///	tauParticles.push_back(pFactory.particle(mytracks[jjj], pion_mass, chi, ndf, pion_sigma));
///	tauParticles.push_back(pFactory.particle(mytracks[kkk], pion_mass, chi, ndf, pion_sigma));
///
///  
///	//reconstructing a tau decay
///	RefCountedKinematicTree tauTree = kpvFitter.fit(tauParticles);
///
///
///	if(tauTree->isEmpty() || !tauTree->isValid() || !tauTree->isConsistent()) continue;
///
///	//getting the J/Psi KinematicParticle
///	tauTree->movePointerToTheTop();
///
///	RefCountedKinematicParticle tau_part = tauTree->currentParticle();
///	if(!tau_part->currentState().isValid()) continue;
///	RefCountedKinematicVertex tau_vertex = tauTree->currentDecayVertex();
///	if(!tau_vertex->vertexIsValid()) continue; 

	// 6.1.2020 commented out
	//	if(TMath::Prob(tau_vertex->chiSquared(), tau_vertex->degreesOfFreedom()) <= c_vprob) continue;

	  
	//	std::vector< RefCountedKinematicParticle > tau_children = tauTree->finalStateParticles();
	  
	//	math::PtEtaPhiMLorentzVector tau1_fit = daughter_p4(tau_children, 0);
	//	math::PtEtaPhiMLorentzVector tau2_fit = daughter_p4(tau_children, 1);
	//	math::PtEtaPhiMLorentzVector tau3_fit = daughter_p4(tau_children, 2);

	//	math::PtEtaPhiMLorentzVector tlv_tau = tau1_fit + tau2_fit + tau3_fit;

	//	particle_cand Taucand; 
	//	Taucand = calculateIPvariables(extrapolator, tau_part, tau_vertex, closestVertex);


	// 6.1.2020 commented out
	//	if(Taucand.fls3d < c_fsig) continue;

	//	Float_t taumass = tau_part->currentState().mass();
	//	if(!(0.2 < taumass && taumass < 1.5)) continue;

	//	std::cout << "sanity check" << tlv_tau.Pt() << " " << tau_part->currentState().globalMomentum().perp() << std::endl;






	Bool_t isRight = false; 
	//	Bool_t isRight1 = false; 
	//	Bool_t isRight2 = false; 
	//	Bool_t isRight3 = false; 
//	Float_t dr1 = 999;
//	Float_t dr2 = 999;
//	Float_t dr3 = 999;
//	Float_t ptres1 = 999;
//	Float_t ptres2 = 999;
//	Float_t ptres3 = 999;

//	Int_t pid = -999;
	Float_t matched_gentaupt = -999;
	Float_t matched_gentaueta = -999;
	Float_t matched_gentauphi = -999;
	Float_t matched_gentaumass = -999;
	Float_t matched_gentaupt_bd = -999;
	Float_t matched_gentaueta_bd = -999;
	Float_t matched_gentauphi_bd = -999;
	Float_t matched_gentaumass_bd = -999;
	

	//	std::cout << "test1" << std::endl;
	if(isMC){

	  for(unsigned int mmm=0; mmm < gps.size(); mmm++){
	    
	    Bool_t isRight1_ = false;
	    Bool_t isRight2_ = false;
	    Bool_t isRight3_ = false;
	    
	    std::vector<TLorentzVector> tlvs = gps[mmm];
	    
	    for(unsigned int nnn=0; nnn < tlvs.size(); nnn++){

	      if(
		 reco::deltaR(pf1.eta(), pf1.phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.015 &&
		 pf1.pt()/tlvs[nnn].Pt() > 0.85 && 
		 pf1.pt()/tlvs[nnn].Pt() < 1.15
		 ){
		isRight1_ = true; 
		//		dr1 = reco::deltaR(tau1_fit.Eta(), tau1_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi());
		//		ptres1 = tau1_fit.Pt()/tlvs[nnn].Pt();
		//		std::cout << "\t\t\t gen. Nr" << mmm << " (" << tlvs[nnn].Pt() << ", " << tlvs[nnn].Eta()  << ", " << tlvs[nnn].Phi() << ") is matched with PF index = " << iii << ", with (pt,eta,phi) = " << tau1_fit.Pt() << ", " << tau1_fit.Eta() << ", " << tau1_fit.Phi() <<  std::endl;
		//		std::cout << "\t\t\t fitted pt respoinse = " << (tau1_fit.Pt()/tlvs[nnn].Pt()) << ", delta_R = " << reco::deltaR(tau1_fit.Eta(), tau1_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi())  << std::endl;
		//		std::cout << "\t\t\t original pt respoinse = " << (pf1.pt()/tlvs[nnn].Pt()) << ", delta_R = " << reco::deltaR(pf1.eta(), pf1.phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi())  << std::endl;
	      }
	      
	      if(
		 reco::deltaR(pf2.eta(), pf2.phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.015 &&
		 pf2.pt()/tlvs[nnn].Pt() > 0.85 && 
		 pf2.pt()/tlvs[nnn].Pt() < 1.15
		 ){
		isRight2_ = true; 
		//		dr2 = reco::deltaR(tau2_fit.Eta(), tau2_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi());
		//		ptres2 = tau2_fit.Pt()/tlvs[nnn].Pt(); 
		//		std::cout << "\t\t\t gen. Nr" << mmm << " (" << tlvs[nnn].Pt() << ", " << tlvs[nnn].Eta()  << ", " << tlvs[nnn].Phi() << ") is matched with PF index = " << jjj << ", with (pt,eta,phi) = " << tau2_fit.Pt() << ", " << tau2_fit.Eta() << ", " << tau2_fit.Phi() <<  std::endl;
		//		std::cout << "\t\t\t fitted pt respoinse = " << (tau2_fit.Pt()/tlvs[nnn].Pt()) << ", delta_R = " << reco::deltaR(tau2_fit.Eta(), tau2_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi())  << std::endl;
		//		std::cout << "\t\t\t original pt respoinse = " << (pf2.pt()/tlvs[nnn].Pt()) << ", delta_R = " << reco::deltaR(pf2.eta(), pf2.phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi())  << std::endl;
	      }

	      
	      if(
		 reco::deltaR(pf3.eta(), pf3.phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.015 &&
		 pf3.pt()/tlvs[nnn].Pt() > 0.85 && 
		 pf3.pt()/tlvs[nnn].Pt() < 1.15
		 ){
		isRight3_ = true; 
		//		dr3 = reco::deltaR(tau3_fit.Eta(), tau3_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi());
		//		ptres3 = tau3_fit.Pt()/tlvs[nnn].Pt(); 
		//		std::cout << "bgen. Nr" << mmm << " is matched with PF index = " << kkk << std::endl;

		//		std::cout << "\t\t\t gen. Nr" << mmm << " (" << tlvs[nnn].Pt() << ", " << tlvs[nnn].Eta()  << ", " << tlvs[nnn].Phi() << ") is matched with PF index = " << kkk << ", with (pt,eta,phi) = " << tau3_fit.Pt() << ", " << tau3_fit.Eta() << ", " << tau3_fit.Phi() <<  std::endl;
		//		std::cout << "\t\t\t fitted pt respoinse = " << (tau3_fit.Pt()/tlvs[nnn].Pt()) << ", delta_R = " << reco::deltaR(tau3_fit.Eta(), tau3_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi())  << std::endl;
		//		std::cout << "\t\t\t original pt respoinse = " << (pf3.pt()/tlvs[nnn].Pt()) << ", delta_R = " << reco::deltaR(pf3.eta(), pf3.phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi())  << std::endl;

	      }

	      
	    }
	    
	    Bool_t isRight_ = isRight1_ && isRight2_ && isRight3_;
//	    if(isRight1_) isRight1 = true;
//	    if(isRight2_) isRight2 = true;
//	    if(isRight3_) isRight3 = true;
//	    std::cout << "test2" << std::endl;
	    if(isRight_){
	      isRight = true;
	      //	      pid = ppdgId[mmm];
	      matched_gentaupt = vec_gentau3pp4[mmm].Pt();
	      matched_gentaueta = vec_gentau3pp4[mmm].Eta();
	      matched_gentauphi = vec_gentau3pp4[mmm].Phi();
	      matched_gentaumass = vec_gentau3pp4[mmm].M();
	      matched_gentaupt_bd = vec_gentau3pp4_full[mmm].Pt();
	      matched_gentaueta_bd = vec_gentau3pp4_full[mmm].Eta();
	      matched_gentauphi_bd = vec_gentau3pp4_full[mmm].Phi();
	      matched_gentaumass_bd = vec_gentau3pp4_full[mmm].M();
	    }
	  }	
	}

	//	std::cout << "test3" << std::endl;
	//	if(isTruth_ && isMC){
	if(!isRight) continue;
	  //	}
	//	std::cout << "test4" << std::endl;

	taucandgen _cand_ = {
	  iii,
	  jjj,
	  kkk,
	  matched_gentaupt,
	  matched_gentaueta,
	  matched_gentauphi,
	  matched_gentaumass,
	  matched_gentaupt_bd,
	  matched_gentaueta_bd,
	  matched_gentauphi_bd,
	  matched_gentaumass_bd
	};
	//	std::cout << cands.size() << std::endl;
	cands.push_back(_cand_);
      }
    }
  }
  //	std::cout << "test5" << std::endl;
  sort(cands.begin(), cands.end());

  //  std::vector<Int_t> dict_idx;
  //  std::vector<Int_t> dict_charges;
  //  Int_t ncomb = 0;
    
  //    bool isRight_bS = false;
  //    bool isRight_aS = false;
  //    int isRight_bS_ith = -1;
  //    int isRight_aS_ith = -1;
  //    int isRight_bS_n = 0;
  //    int isRight_aS_n = 0;

  if(cands.size()==0) return false;
  //    for(int ic=0; ic < (int)cands.size(); ic++){
  //
  //      
  //      Int_t _idx1 = cands[ic].cand_tau_id1;
  //      Int_t _idx2 = cands[ic].cand_tau_id2;
  //      Int_t _idx3 = cands[ic].cand_tau_id3;
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
  ////      if(cands[ic].cand_tau_isRight==true){
  ////	isRight_bS = true;
  ////	isRight_bS_ith = ic;
  ////	isRight_bS_n += 1;
  ////      }
  //      
  //      
  //      if(flag_overlap) continue; 
  //
  ////      if(cands[ic].cand_tau_isRight==true){
  ////	isRight_aS = true;
  ////	isRight_aS_ith = ncomb;
  ////	isRight_aS_n += 1;
  ////      }
  //
  //      dict_idx.push_back(_idx1);
  //      dict_idx.push_back(_idx2);
  //      dict_idx.push_back(_idx3);
  //
  //      dict_charges.push_back(cands[ic].cand_tau_charge);
  //      
  //      ncomb += 1;
  //      if(ncomb==2) break;
  //    }
  //    
  //    if(ncomb!=2) return false;
    
  std::cout << "Nr of candidates = " << cands.size() << std::endl;

  for(int idx=0; idx < (int)cands.size(); idx++){

    
    Int_t tau_idx1 = cands[idx].cand_tau_id1;
    Int_t tau_idx2 = cands[idx].cand_tau_id2;
    Int_t tau_idx3 = cands[idx].cand_tau_id3;


    pat::PackedCandidate tau_pf1 = pfcollection[tau_idx1];
    pat::PackedCandidate tau_pf2 = pfcollection[tau_idx2];
    pat::PackedCandidate tau_pf3 = pfcollection[tau_idx3];
    
    TLorentzVector tlv_pion1; 
    TLorentzVector tlv_pion2;
    TLorentzVector tlv_pion3;
	
    tlv_pion1.SetPtEtaPhiM(tau_pf1.pt(), tau_pf1.eta(), tau_pf1.phi(), tau_pf1.mass());
    tlv_pion2.SetPtEtaPhiM(tau_pf2.pt(), tau_pf2.eta(), tau_pf2.phi(), tau_pf2.mass());
    tlv_pion3.SetPtEtaPhiM(tau_pf3.pt(), tau_pf3.eta(), tau_pf3.phi(), tau_pf3.mass());
	
    nBranches_->BsTauTauFH_mr_tau_pi1_pt.push_back(tlv_pion1.Pt());
    nBranches_->BsTauTauFH_mr_tau_pi1_eta.push_back(tlv_pion1.Eta());
    nBranches_->BsTauTauFH_mr_tau_pi1_phi.push_back(tlv_pion1.Phi());
    nBranches_->BsTauTauFH_mr_tau_pi1_mass.push_back(tlv_pion1.M());
    nBranches_->BsTauTauFH_mr_tau_pi2_pt.push_back(tlv_pion2.Pt());
    nBranches_->BsTauTauFH_mr_tau_pi2_eta.push_back(tlv_pion2.Eta());
    nBranches_->BsTauTauFH_mr_tau_pi2_phi.push_back(tlv_pion2.Phi());
    nBranches_->BsTauTauFH_mr_tau_pi2_mass.push_back(tlv_pion2.M());
    nBranches_->BsTauTauFH_mr_tau_pi3_pt.push_back(tlv_pion3.Pt());
    nBranches_->BsTauTauFH_mr_tau_pi3_eta.push_back(tlv_pion3.Eta());
    nBranches_->BsTauTauFH_mr_tau_pi3_phi.push_back(tlv_pion3.Phi());
    nBranches_->BsTauTauFH_mr_tau_pi3_mass.push_back(tlv_pion3.M());

    nBranches_->BsTauTauFH_mr_tau_genpt.push_back(cands[idx].cand_gentau_pt);
    nBranches_->BsTauTauFH_mr_tau_geneta.push_back(cands[idx].cand_gentau_eta);
    nBranches_->BsTauTauFH_mr_tau_genphi.push_back(cands[idx].cand_gentau_phi);
    nBranches_->BsTauTauFH_mr_tau_genmass.push_back(cands[idx].cand_gentau_mass);

    nBranches_->BsTauTauFH_mr_tau_genpt_bd.push_back(cands[idx].cand_gentau_pt_bd);
    nBranches_->BsTauTauFH_mr_tau_geneta_bd.push_back(cands[idx].cand_gentau_eta_bd);
    nBranches_->BsTauTauFH_mr_tau_genphi_bd.push_back(cands[idx].cand_gentau_phi_bd);
    nBranches_->BsTauTauFH_mr_tau_genmass_bd.push_back(cands[idx].cand_gentau_mass_bd);

  }

  return true;


}


void BsTauTauFHNtuplizer_mr::printout(const RefCountedKinematicVertex& myVertex){
  std::cout << "Vertex:" << std::endl;
  if (myVertex->vertexIsValid()) {
    std::cout << "\t Decay vertex: " << myVertex->position() << myVertex->chiSquared() << " " << myVertex->degreesOfFreedom()
	      << std::endl;
  } else
    std::cout << "\t Decay vertex Not valid\n";
}

void BsTauTauFHNtuplizer_mr::printout(const RefCountedKinematicParticle& myParticle){
  std::cout << "Particle:" << std::endl;
  //accessing the reconstructed Bs meson parameters:
  //SK: uncomment if needed  AlgebraicVector7 bs_par = myParticle->currentState().kinematicParameters().vector();

  //and their joint covariance matrix:
  //SK:uncomment if needed  AlgebraicSymMatrix77 bs_er = myParticle->currentState().kinematicParametersError().matrix();
  std::cout << "\t Momentum at vertex: " << myParticle->currentState().globalMomentum() << std::endl;
  std::cout << "\t Parameters at vertex: " << myParticle->currentState().kinematicParameters().vector() << std::endl;
}

void BsTauTauFHNtuplizer_mr::printout(const RefCountedKinematicTree& myTree){
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
