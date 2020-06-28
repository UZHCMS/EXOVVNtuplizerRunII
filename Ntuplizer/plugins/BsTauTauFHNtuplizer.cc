#include "../interface/BsTauTauFHNtuplizer.h"


//===================================================================================================================
BsTauTauFHNtuplizer::BsTauTauFHNtuplizer( edm::EDGetTokenT<pat::MuonCollection>    muonToken   ,
					  edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
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
  std::cout << "Confirm (dzcut, fsigcut, vprobcut, dnn cut) = " << c_dz << " " << c_fsig << " " << c_vprob << " " << c_dnn<< std::endl;
  
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
BsTauTauFHNtuplizer::~BsTauTauFHNtuplizer( void )
{

}


//Int_t BsTauTauFHNtuplizer::decaymode_id(std::string str){
//  if(str=="electron") return -2;
//  else if(str=="muon") return -1;
//  else if(str=="oneProng0Pi0") return 0;
//  else if(str=="oneProng1Pi0") return 1;
//  else if(str=="oneProng2Pi0") return 2;
//  else if(str=="oneProng3Pi0") return 3;
//  else if(str=="oneProngOther") return 4;  
//  else if(str=="threeProng0Pi0") return 10;
//  else if(str=="threeProng1Pi0") return 11;
//  else if(str=="threeProngOther") return 14;
//  else if(str=="rare") return 15;
//  else return -9;
//}
//
//
//TVector3 BsTauTauFHNtuplizer::getVertex(const reco::GenParticle& part){
//  return TVector3(part.vx(),part.vy(),part.vz());
//}
//
//float BsTauTauFHNtuplizer::MuonPFIso(pat::Muon muon){
//
//  float sumChargedHadronPt = muon.pfIsolationR04().sumChargedHadronPt;
//  float sumNeutralHadronEt = muon.pfIsolationR04().sumNeutralHadronEt;
//  float sumPhotonEt = muon.pfIsolationR04().sumPhotonEt;
//  float sumPUPt = muon.pfIsolationR04().sumPUPt;
//  float iso = (sumChargedHadronPt + std::max( 0. ,sumNeutralHadronEt + sumPhotonEt - 0.5 * sumPUPt));// / muon.pt()
// 
//  return iso;
//}
//
//
//
//
//Float_t BsTauTauFHNtuplizer::getMaxDoca(std::vector<RefCountedKinematicParticle> &kinParticles){
//
//  double maxDoca = -1.0;
//
//  TwoTrackMinimumDistance md;
//  std::vector<RefCountedKinematicParticle>::iterator in_it, out_it;
//
//  for (out_it = kinParticles.begin(); out_it != kinParticles.end(); ++out_it) {
//    for (in_it = out_it + 1; in_it != kinParticles.end(); ++in_it) {
//      md.calculate((*out_it)->currentState().freeTrajectoryState(),(*in_it)->currentState().freeTrajectoryState());
//      if (md.distance() > maxDoca)
//	maxDoca = md.distance();
//    }
//  }
//
//  return maxDoca;
//}
//
//
//
//Float_t BsTauTauFHNtuplizer::getMinDoca(std::vector<RefCountedKinematicParticle> &kinParticles) {
//
//  double minDoca = 99999.9;
//
//  TwoTrackMinimumDistance md;
//  unsigned j,k,n;
//
//  n = kinParticles.size();
//  for (j = 0; j < n; j++) {
//    for (k = j+1; k < n; k++) {
//      md.calculate(kinParticles[j]->currentState().freeTrajectoryState(),kinParticles[k]->currentState().freeTrajectoryState());
//      if (md.distance() < minDoca)
//	minDoca = md.distance();
//    }
//  }
//
//  return minDoca;
//}
//
//
//
//
//std::tuple<Float_t, TransientVertex> BsTauTauFHNtuplizer::vertexProb( const std::vector<reco::TransientTrack>& tracks){
//
//  Float_t vprob = -1;
//  
//  KalmanVertexFitter kalman_fitter;
//  TransientVertex vertex;
//
//  try{
//    vertex = kalman_fitter.vertex(tracks);
//  }catch(std::exception e){
//    std::cout << "No vertex found ... return" << std::endl;
//    return std::forward_as_tuple(-9, vertex);
//  }
//
//  if(vertex.isValid()){
//
//    vprob =  TMath::Prob(vertex.totalChiSquared(), vertex.degreesOfFreedom());
//
//    //    vx = vertex.position().x();
//    //    vy = vertex.position().y();
//    //    vz = vertex.position().z();
//    
//    return std::forward_as_tuple(vprob, vertex);
//
//  }else{
//
//    return std::forward_as_tuple(-9, vertex);
//
//  }
//}
//
//
////adapt absoluteImpactParameter functionality for RefCountedKinematicVertex
//std::pair<bool, Measurement1D> BsTauTauFHNtuplizer::absoluteImpactParameter(const TrajectoryStateOnSurface& tsos,
//									    RefCountedKinematicVertex vertex,
//									    VertexDistance& distanceComputer){
//  if (!tsos.isValid()) {
//    return std::pair<bool, Measurement1D>(false, Measurement1D(0., 0.));
//  }
//  GlobalPoint refPoint = tsos.globalPosition();
//  GlobalError refPointErr = tsos.cartesianError().position();
//  GlobalPoint vertexPosition = vertex->vertexState().position();
//  GlobalError vertexPositionErr = RecoVertex::convertError(vertex->vertexState().error());
//  return std::pair<bool, Measurement1D>(
//					true,
//					distanceComputer.distance(VertexState(vertexPosition, vertexPositionErr), VertexState(refPoint, refPointErr)));
//}
//
//
//
//
//particle_cand BsTauTauFHNtuplizer::calculateIPvariables(
//							AnalyticalImpactPointExtrapolator extrapolator,
//							RefCountedKinematicParticle particle,
//							RefCountedKinematicVertex vertex,
//							reco::Vertex wrtVertex
//							){
//
//  TrajectoryStateOnSurface tsos = extrapolator.extrapolate(particle->currentState().freeTrajectoryState(),
//							   RecoVertex::convertPos(wrtVertex.position()));
//
//
//  VertexDistance3D a3d;  
//
//  std::pair<bool,Measurement1D> currentIp = IPTools::signedDecayLength3D(tsos, GlobalVector(0,0,1), wrtVertex);
//  std::pair<bool,Measurement1D> cur3DIP = IPTools::absoluteImpactParameter(tsos, wrtVertex, a3d);
//  
//  // flight length
//  Float_t fl3d = a3d.distance(wrtVertex, vertex->vertexState()).value();
//  Float_t fl3de = a3d.distance(wrtVertex, vertex->vertexState()).error();
//  Float_t fls3d = -1;
//
//  if(fl3de!=0) fls3d = fl3d/fl3de;
//
//  // longitudinal impact parameters
//  Float_t lip = currentIp.second.value();
//  Float_t lipe = currentIp.second.error();
//  Float_t lips = -1;
//  
//  if(lipe!=0) lips = lip/lipe;
//
//  // impact parameter to the PV
//  Float_t pvip = cur3DIP.second.value();
//  Float_t pvipe = cur3DIP.second.error();
//  Float_t pvips = -1;
//  
//  if(pvipe!=0) pvips = pvip/pvipe;
//
//  // opening angle
//  TVector3 plab = TVector3(particle->currentState().globalMomentum().x(),
//			   particle->currentState().globalMomentum().y(),
//			   particle->currentState().globalMomentum().z());
//
//  const TVector3 tv3diff = TVector3(vertex->vertexState().position().x() - wrtVertex.position().x(),
//				    vertex->vertexState().position().y() - wrtVertex.position().y(),
//				    vertex->vertexState().position().z() - wrtVertex.position().z()
//				    );
//
//  Float_t alpha = -1;
//
//  if(plab.Mag() != 0. && tv3diff.Mag()!=0){
//    alpha = plab.Dot(tv3diff) / (plab.Mag() * tv3diff.Mag());
//  }
//
//  particle_cand cand = {
//    lip,
//    lips,
//    pvip, 
//    pvips,
//    fl3d,
//    fls3d,
//    alpha
//  };
//
//
//  return cand;
//}
//
//
//math::PtEtaPhiMLorentzVector BsTauTauFHNtuplizer::daughter_p4(std::vector< RefCountedKinematicParticle > fitted_children, size_t i){
//  const auto& state = fitted_children.at(i)->currentState();
//
//  return math::PtEtaPhiMLorentzVector(
//				      state.globalMomentum().perp(), 
//				      state.globalMomentum().eta() ,
//				      state.globalMomentum().phi() ,
//				      state.mass()
//				      );
//}


bool BsTauTauFHNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  
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


  //    std::cout << "number of matched muon = " << muoncollection.size() << std::endl;
  if(!( muoncollection.size() >= 1)) return false;
  nBranches_->cutflow_perevt->Fill(2);
    
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

	//	std::cout << "CHECK1_muon: " <<count_dnn<< std::endl;	  

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
      //      std::cout << "CHECK1_pf: " <<count_dnn<< std::endl;	  

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
      //      std::cout << "pf dnn = " << ic << " " << finalOutputTensor(0, ic, 0) << " " << finalOutputTensor(0, ic, 1) << std::endl;
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
  Int_t isgen3 = 0;
  Int_t isgen3matched = 0;

  bool isMC = runOnMC_;

  if(isMC){
    event.getByToken(genParticlesToken_ , genParticles_); 
    event.getByToken(genTauToken_, genTaus_);

    for( unsigned p=0; p < genParticles_->size(); ++p){
      
      if(TMath::Abs((*genParticles_)[p].pdgId())!=15) continue;
      if(TMath::Abs((*genParticles_)[p].status())!=2) continue;
      
      //      std::cout << "\t Tau found with # of daughters = " << (*genParticles_)[p].numberOfDaughters() << " with mother = " << (*genParticles_)[p].mother(0)->pdgId() << std::endl;
      
      
      // calculate visible pt ... 

      TLorentzVector genvis;
      std::vector<TLorentzVector> gp;
      Bool_t matched = true;
      Int_t nprong = 0;
      
      for(int idd = 0; idd < (int)(*genParticles_)[p].numberOfDaughters(); idd++){
	
//	std::cout << "\t\t -> " << (*genParticles_)[p].daughter(idd)->pdgId() << " (pT, eta, phi) = " 
//		  << (*genParticles_)[p].daughter(idd)->pt() << " " 
//		  << (*genParticles_)[p].daughter(idd)->eta() << " " 
//		  << (*genParticles_)[p].daughter(idd)->phi() << std::endl;


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
	  //	  Int_t min_index = 999;
	  
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
	      //	      min_index = kkk; 
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
	  //	  else std::cout << "\t\t\t -----> PF index = " << min_index << " matched!" << std::endl;

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
	  taugendm = aux.decaymode_id(JetMCTagUtils::genTauDecayMode(TauCand));
	}
      }
      
      vec_ppdgId.push_back((*genParticles_)[p].mother(0)->pdgId());
      vec_gentaudm.push_back(taugendm);
      vec_gentaup4.push_back(genvis);


  

      
      if(gp.size()==3){
	//	std::cout << "\t -----> This has been registered with mother = " << (*genParticles_)[p].mother(0)->pdgId() << std::endl;
	gps.push_back(gp);
	ppdgId.push_back((*genParticles_)[p].mother(0)->pdgId());
	vec_gentau3pp4.push_back(genvis);
	
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


  //////////////////////////////
  std::cout << "Starts to build tau candidate out of " << numOfch << " pion candidates" << std::endl;

  std::vector<taucandsimple> cands;
  Int_t npf_after_dnn = 0;
    
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


	//	std::cout << iii << " " << jjj << " " << kkk << std::endl;

	Int_t tau_charge = pf1.charge() + pf2.charge() + pf3.charge(); 

	if(TMath::Abs(tau_charge)!=1) continue; 



	std::vector<RefCountedKinematicParticle> tauParticles;

	tauParticles.push_back(pFactory.particle(mytracks[iii], aux.pion_mass, chi, ndf, aux.pion_sigma));
	tauParticles.push_back(pFactory.particle(mytracks[jjj], aux.pion_mass, chi, ndf, aux.pion_sigma));
	tauParticles.push_back(pFactory.particle(mytracks[kkk], aux.pion_mass, chi, ndf, aux.pion_sigma));

  
	//reconstructing a tau decay
	RefCountedKinematicTree tauTree = kpvFitter.fit(tauParticles);


	if(tauTree->isEmpty() || !tauTree->isValid() || !tauTree->isConsistent()) continue;

	//getting the J/Psi KinematicParticle
	tauTree->movePointerToTheTop();

	RefCountedKinematicParticle tau_part = tauTree->currentParticle();
	if(!tau_part->currentState().isValid()) continue;
	RefCountedKinematicVertex tau_vertex = tauTree->currentDecayVertex();
	if(!tau_vertex->vertexIsValid()) continue; 

	// 6.1.2020 commented out
	if(TMath::Prob(tau_vertex->chiSquared(), tau_vertex->degreesOfFreedom()) <= c_vprob) continue;

	  
	//	std::vector< RefCountedKinematicParticle > tau_children = tauTree->finalStateParticles();
	  
	//	math::PtEtaPhiMLorentzVector tau1_fit = daughter_p4(tau_children, 0);
	//	math::PtEtaPhiMLorentzVector tau2_fit = daughter_p4(tau_children, 1);
	//	math::PtEtaPhiMLorentzVector tau3_fit = daughter_p4(tau_children, 2);

	//	math::PtEtaPhiMLorentzVector tlv_tau = tau1_fit + tau2_fit + tau3_fit;

	particle_cand Taucand; 
	Taucand = aux.calculateIPvariables(extrapolator, tau_part, tau_vertex, closestVertex);


	// 6.1.2020 commented out
	if(Taucand.fls3d < c_fsig) continue;

	//	Float_t taumass = tau_part->currentState().mass();
	//	if(!(0.2 < taumass && taumass < 1.5)) continue;

	//	std::cout << "sanity check" << tlv_tau.Pt() << " " << tau_part->currentState().globalMomentum().perp() << std::endl;

	taucandsimple _cand_ = {
	  iii,
	  jjj,
	  kkk,
	  (Float_t)tau_part->currentState().globalMomentum().perp(),
	  (Int_t)tau_charge
	};
	//	std::cout << cands.size() << std::endl;
	cands.push_back(_cand_);
      }
    }
  }

  sort(cands.begin(), cands.end());

  //  std::vector<Int_t> dict_idx;
  //  std::vector<Int_t> dict_charges;
  Int_t ncomb = 0;
    
  //    bool isRight_bS = false;
  //    bool isRight_aS = false;
  //    int isRight_bS_ith = -1;
  //    int isRight_aS_ith = -1;
  //    int isRight_bS_n = 0;
  //    int isRight_aS_n = 0;

  if(cands.size()<=1) return false;
  nBranches_->cutflow_perevt->Fill(3);
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

  for(int leg1=0; leg1 < (int)cands.size(); leg1++){
    for(int leg2=leg1+1; leg2 < (int)cands.size(); leg2++){


      //      Int_t _idx1 = cands[ic].cand_tau_id1;
      //      Int_t _idx2 = cands[ic].cand_tau_id2;
      //      Int_t _idx3 = cands[ic].cand_tau_id3;

      //      std::cout << "check1" << std::endl;
    
      Int_t tau1_idx1 = cands[leg1].cand_tau_id1;
      Int_t tau1_idx2 = cands[leg1].cand_tau_id2;
      Int_t tau1_idx3 = cands[leg1].cand_tau_id3;

      //      std::cout << "check2" << std::endl;

      Int_t tau2_idx1 = cands[leg2].cand_tau_id1;
      Int_t tau2_idx2 = cands[leg2].cand_tau_id2;
      Int_t tau2_idx3 = cands[leg2].cand_tau_id3;

      //      std::cout << "check3-1 " << tau1_idx1 << " " << tau1_idx2 << " " << tau1_idx3 << std::endl;
      //      std::cout << "check3-2 " << tau2_idx1 << " " << tau2_idx2 << " " << tau2_idx3 << std::endl;

      if(!isTruth_ || !isMC){
	if(tau1_idx1==tau2_idx1 || tau1_idx1==tau2_idx2 || tau1_idx1==tau2_idx3) continue;
	if(tau1_idx2==tau2_idx1 || tau1_idx2==tau2_idx2 || tau1_idx2==tau2_idx3) continue;
	if(tau1_idx3==tau2_idx1 || tau1_idx3==tau2_idx2 || tau1_idx3==tau2_idx3) continue;
	//	std::cout << "This is to remove overlap" << std::endl;
      }


      std::vector<RefCountedKinematicParticle> tauParticles1;
      tauParticles1.push_back(pFactory.particle(mytracks[tau1_idx1], aux.pion_mass, chi, ndf, aux.pion_sigma));
      tauParticles1.push_back(pFactory.particle(mytracks[tau1_idx2], aux.pion_mass, chi, ndf, aux.pion_sigma));
      tauParticles1.push_back(pFactory.particle(mytracks[tau1_idx3], aux.pion_mass, chi, ndf, aux.pion_sigma));

      std::vector<RefCountedKinematicParticle> tauParticles2;          
      tauParticles2.push_back(pFactory.particle(mytracks[tau2_idx1], aux.pion_mass, chi, ndf, aux.pion_sigma));
      tauParticles2.push_back(pFactory.particle(mytracks[tau2_idx2], aux.pion_mass, chi, ndf, aux.pion_sigma));
      tauParticles2.push_back(pFactory.particle(mytracks[tau2_idx3], aux.pion_mass, chi, ndf, aux.pion_sigma));


      //      std::cout << "check4" << std::endl;    
//      allParticles.push_back(pFactory.particle(mytracks[tau1_idx1], pion_mass, chi, ndf, pion_sigma));
//      allParticles.push_back(pFactory.particle(mytracks[tau1_idx2], pion_mass, chi, ndf, pion_sigma));
//      allParticles.push_back(pFactory.particle(mytracks[tau1_idx3], pion_mass, chi, ndf, pion_sigma));
//      allParticles.push_back(pFactory.particle(mytracks[tau2_idx1], pion_mass, chi, ndf, pion_sigma));
//      allParticles.push_back(pFactory.particle(mytracks[tau2_idx2], pion_mass, chi, ndf, pion_sigma));
//      allParticles.push_back(pFactory.particle(mytracks[tau2_idx3], pion_mass, chi, ndf, pion_sigma));
    
	
      //reconstructing a tau decay
      RefCountedKinematicTree tauTree1 = kpvFitter.fit(tauParticles1);
      
      if(tauTree1->isEmpty() || !tauTree1->isValid() || !tauTree1->isConsistent()) continue;
      
      tauTree1->movePointerToTheTop();
      
      RefCountedKinematicParticle tau_part1 = tauTree1->currentParticle();
      if(!tau_part1->currentState().isValid()) continue;
      RefCountedKinematicVertex tau_vertex1 = tauTree1->currentDecayVertex();
      if(!tau_vertex1->vertexIsValid()) continue; 


      //      std::cout << "check5" << std::endl;    
      //	// 6.1.2020 commented out
      //	//	if(TMath::Prob(tau_vertex1->chiSquared(), tau_vertex1->degreesOfFreedom()) <= c_vprob) continue;
      //
      //	  
      //	std::vector< RefCountedKinematicParticle > tau_children1 = tauTree1->finalStateParticles();
      //	
      //	math::PtEtaPhiMLorentzVector tau1_pi1 = daughter_p4(tau_children1, 0);
      //	math::PtEtaPhiMLorentzVector tau1_pi2 = daughter_p4(tau_children1, 1);
      //	math::PtEtaPhiMLorentzVector tau1_pi3 = daughter_p4(tau_children1, 2);
      //	
      //	math::PtEtaPhiMLorentzVector tlv_tau1_fit = tau1_pi1 + tau1_pi2 + tau1_pi3;
      //	
      //	particle_cand Taucand1; 
      //	Taucand1 = calculateIPvariables(extrapolator, tau_part1, tau_vertex1, closestVertex);
      //	
      //	if(Taucand1.fls3d < c_fsig) continue;
      //	
      //	
      //	if(tlv_tau1_fit.Pt() < 3.) continue;
	

      ///////////////////////////////


      //    std::vector<RefCountedKinematicParticle> tauParticles2;
      //    
      //    tauParticles2.push_back(pFactory.particle(mytracks[tau2_idx1], pion_mass, chi, ndf, pion_sigma));
      //    tauParticles2.push_back(pFactory.particle(mytracks[tau2_idx2], pion_mass, chi, ndf, pion_sigma));
      //    tauParticles2.push_back(pFactory.particle(mytracks[tau2_idx3], pion_mass, chi, ndf, pion_sigma));
      //    
      //    
      //    //reconstructing a tau decay
          RefCountedKinematicTree tauTree2 = kpvFitter.fit(tauParticles2);
              
          if(tauTree2->isEmpty() || !tauTree2->isValid() || !tauTree2->isConsistent()) continue;
          
          tauTree2->movePointerToTheTop();
                    
          RefCountedKinematicParticle tau_part2 = tauTree2->currentParticle();
          if(!tau_part2->currentState().isValid()) continue;
          RefCountedKinematicVertex tau_vertex2 = tauTree2->currentDecayVertex();
          if(!tau_vertex2->vertexIsValid()) continue; 

	  //      std::cout << "check6" << std::endl;          
      //    // 6.1.2020 commented out
      //    if(TMath::Prob(tau_vertex2->chiSquared(), tau_vertex2->degreesOfFreedom()) <= c_vprob) continue;
      //
      //	  
      //    std::vector< RefCountedKinematicParticle > tau_children2 = tauTree2->finalStateParticles();
      //    
      //    math::PtEtaPhiMLorentzVector tau2_pi1 = daughter_p4(tau_children2, 0);
      //    math::PtEtaPhiMLorentzVector tau2_pi2 = daughter_p4(tau_children2, 1);
      //    math::PtEtaPhiMLorentzVector tau2_pi3 = daughter_p4(tau_children2, 2);
      //    
      //    math::PtEtaPhiMLorentzVector tlv_tau2_fit = tau2_pi1 + tau2_pi2 + tau2_pi3;
      //    
      //    particle_cand Taucand2; 
      //    Taucand2 = calculateIPvariables(extrapolator, tau_part2, tau_vertex2, closestVertex);
      //
      //    if(Taucand2.fls3d < c_fsig) continue;
      //	  
      //
      //    if(tlv_tau2_fit.Pt() < 3.) continue;

      ///////////////////////////////////////////////////////////

      std::vector<RefCountedKinematicParticle> allParticles;  
      allParticles.push_back(tau_part1);
      allParticles.push_back(tau_part2);

      //      std::cout << "check7" << std::endl;                

      RefCountedKinematicTree bcTree = kpvFitter.fit(allParticles);
	
      if(bcTree->isEmpty() || !bcTree->isValid() || !bcTree->isConsistent()) continue;	
	
      RefCountedKinematicParticle bc_part = bcTree->currentParticle();
      if(!bc_part->currentState().isValid()) continue;
	
      RefCountedKinematicVertex bc_vertex = bcTree->currentDecayVertex();
      if(!bc_vertex->vertexIsValid()) continue;
	
      //      std::cout << "check8" << std::endl;                

      particle_cand Bcand; 
      Bcand = aux.calculateIPvariables(extrapolator, bc_part, bc_vertex, closestVertex);
	
      std::vector< RefCountedKinematicParticle > tau_children = bcTree->finalStateParticles();
	
      math::PtEtaPhiMLorentzVector tau1_pi1 = aux.daughter_p4(tau_children, 0);
      math::PtEtaPhiMLorentzVector tau1_pi2 = aux.daughter_p4(tau_children, 1);
      math::PtEtaPhiMLorentzVector tau1_pi3 = aux.daughter_p4(tau_children, 2);
      math::PtEtaPhiMLorentzVector tau2_pi1 = aux.daughter_p4(tau_children, 3);
      math::PtEtaPhiMLorentzVector tau2_pi2 = aux.daughter_p4(tau_children, 4);
      math::PtEtaPhiMLorentzVector tau2_pi3 = aux.daughter_p4(tau_children, 5);

      math::PtEtaPhiMLorentzVector tlv_tau1_fit = tau1_pi1 + tau1_pi2 + tau1_pi3;
      math::PtEtaPhiMLorentzVector tlv_tau2_fit = tau2_pi1 + tau2_pi2 + tau2_pi3;


      //      std::cout << "check9" << std::endl;

      particle_cand Taucand1; 
      particle_cand Taucand2; 

      Taucand1 = aux.calculateIPvariables(extrapolator, tau_part1, tau_vertex1, closestVertex);
      Taucand2 = aux.calculateIPvariables(extrapolator, tau_part2, tau_vertex2, closestVertex);

      //      std::cout << "check10" << std::endl;     
      Float_t iso = 0;
      Int_t ntracks = 0;
      Float_t iso_mindoca = 999; 
	
	
      for(int itrk = 0; itrk < numOfch;  itrk++){
	  
	if(itrk==tau1_idx1 || itrk==tau1_idx2 || itrk==tau1_idx3) continue;
	if(itrk==tau2_idx1 || itrk==tau2_idx2 || itrk==tau2_idx3) continue;
      
	iso += pfcollection[itrk].pt();
      
	TrajectoryStateOnSurface tsos_pf = extrapolator.extrapolate(mytracks[itrk].impactPointState(), bc_vertex->position());
	
	  
	VertexDistance3D a3d_pf;  
	  
	std::pair<bool,Measurement1D> cur3DIP_pf = aux.absoluteImpactParameter(tsos_pf, bc_vertex, a3d_pf);
	  
	Float_t pvip_pf = cur3DIP_pf.second.value();
	  
	//    std::cout << itrk << ": Distance of closest apporach to the bc vertex = " << pvip << std::endl;
	  
	if(pvip_pf < 0.03) ntracks+=1;
	  
	if(iso_mindoca > pvip_pf) iso_mindoca = pvip_pf;
      }
	
    
      //      std::cout << "check11" << std::endl;     	
      Float_t max_dr1 = -1;
	
      Float_t tau1_dR_12 = reco::deltaR(tau1_pi1.Eta(), tau1_pi1.Phi(), tau1_pi2.Eta(), tau1_pi2.Phi());
      Float_t tau1_dR_13 = reco::deltaR(tau1_pi1.Eta(), tau1_pi1.Phi(), tau1_pi3.Eta(), tau1_pi3.Phi());
      Float_t tau1_dR_23 = reco::deltaR(tau1_pi2.Eta(), tau1_pi2.Phi(), tau1_pi3.Eta(), tau1_pi3.Phi());
    
      if(max_dr1 < tau1_dR_12) max_dr1 = tau1_dR_12;
      if(max_dr1 < tau1_dR_13) max_dr1 = tau1_dR_13;
      if(max_dr1 < tau1_dR_23) max_dr1 = tau1_dR_23;
	

      Float_t max_dr2 = -1;
	
      Float_t tau2_dR_12 = reco::deltaR(tau2_pi1.Eta(), tau2_pi1.Phi(), tau2_pi2.Eta(), tau2_pi2.Phi());
      Float_t tau2_dR_13 = reco::deltaR(tau2_pi1.Eta(), tau2_pi1.Phi(), tau2_pi3.Eta(), tau2_pi3.Phi());
      Float_t tau2_dR_23 = reco::deltaR(tau2_pi2.Eta(), tau2_pi2.Phi(), tau2_pi3.Eta(), tau2_pi3.Phi());
	
      if(max_dr2 < tau2_dR_12) max_dr2 = tau2_dR_12;
      if(max_dr2 < tau2_dR_13) max_dr2 = tau2_dR_13;
      if(max_dr2 < tau2_dR_23) max_dr2 = tau2_dR_23;
	

      //      std::cout << "check12" << std::endl;     	
    
      Bool_t isRight1 = false; 
      Int_t pid1 = -999;
      Float_t matched_gentaupt1 = -999;
	
      Bool_t isRight2 = false; 
      Int_t pid2 = -999;
      Float_t matched_gentaupt2 = -999;
	
      if(isMC){

	//	std::cout << "# of taus =" << gps.size() << std::endl;
	  
	for(unsigned int mmm=0; mmm < gps.size(); mmm++){
	    
	  Bool_t isRight1_ = false;
	  Bool_t isRight2_ = false;
	  Bool_t isRight3_ = false;
	    
	  std::vector<TLorentzVector> tlvs = gps[mmm];


	  //	  std::cout << "\t # of CHs = " << tlvs.size() << std::endl;
	    
	  for(unsigned int nnn=0; nnn < tlvs.size(); nnn++){
	      

//	    std::cout << "\t\t gen. Nr" << mmm << " (pT, eta, phi) = " << tlvs[nnn].Pt() << ", " << tlvs[nnn].Eta()  << ", " << tlvs[nnn].Phi() << std::endl;
//	    std::cout << "\t\t pi1 index = " << tau1_idx1 << ", with (pt,eta,phi) = " << tau1_pi1.Pt() << ", " << tau1_pi1.Eta() << ", " << tau1_pi1.Phi() <<  std::endl;
//	    std::cout << "\t\t pi2 index = " << tau1_idx2 << ", with (pt,eta,phi) = " << tau1_pi2.Pt() << ", " << tau1_pi2.Eta() << ", " << tau1_pi2.Phi() <<  std::endl;
//	    std::cout << "\t\t pi3 index = " << tau1_idx3 << ", with (pt,eta,phi) = " << tau1_pi3.Pt() << ", " << tau1_pi3.Eta() << ", " << tau1_pi3.Phi() <<  std::endl;
//
//	    std::cout << "\t\t pt response = " << (tau1_pi1.Pt()/tlvs[nnn].Pt()) << ", delta_R = " << reco::deltaR(tau1_pi1.Eta(), tau1_pi1.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi())  << std::endl;
//	    std::cout << "\t\t pt response = " << (tau1_pi2.Pt()/tlvs[nnn].Pt()) << ", delta_R = " << reco::deltaR(tau1_pi2.Eta(), tau1_pi2.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi())  << std::endl;
//	    std::cout << "\t\t pt response = " << (tau1_pi3.Pt()/tlvs[nnn].Pt()) << ", delta_R = " << reco::deltaR(tau1_pi3.Eta(), tau1_pi3.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi())  << std::endl;

	    if(
	       //	       reco::deltaR(tau1_pi1.Eta(), tau1_pi1.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.1
	       reco::deltaR(tau1_pi1.Eta(), tau1_pi1.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.015 &&
	       tau1_pi1.Pt()/tlvs[nnn].Pt() > 0.85 && 
	       tau1_pi1.Pt()/tlvs[nnn].Pt() < 1.15

	       ){
	      isRight1_ = true; 

//	      std::cout << "\t\t\t gen. Nr" << mmm << " (" << tlvs[nnn].Pt() << ", " << tlvs[nnn].Eta()  << ", " << tlvs[nnn].Phi() << ") is matched with PF index = " << iii << ", with (pt,eta,phi) = " << tau1_fit.Pt() << ", " << tau1_fit.Eta() << ", " << tau1_fit.Phi() <<  std::endl;
//	      std::cout << "\t\t\t fitted pt respoinse = " << (tau1_fit.Pt()/tlvs[nnn].Pt()) << ", delta_R = " << reco::deltaR(tau1_fit.Eta(), tau1_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi())  << std::endl;
//	      std::cout << "\t\t\t original pt respoinse = " << (pf1.pt()/tlvs[nnn].Pt()) << ", delta_R = " << reco::deltaR(pf1.eta(), pf1.phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi())  << std::endl;

	    }
	      
	    if(
	       //	       reco::deltaR(tau1_pi2.Eta(), tau1_pi2.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.1
	       reco::deltaR(tau1_pi2.Eta(), tau1_pi2.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.015 &&
	       tau1_pi2.Pt()/tlvs[nnn].Pt() > 0.85 && 
	       tau1_pi2.Pt()/tlvs[nnn].Pt() < 1.15

	       ) isRight2_ = true; 
	      
	    if(
	       //	       reco::deltaR(tau1_pi3.Eta(), tau1_pi3.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.1
	       reco::deltaR(tau1_pi3.Eta(), tau1_pi3.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.015 &&
	       tau1_pi3.Pt()/tlvs[nnn].Pt() > 0.85 && 
	       tau1_pi3.Pt()/tlvs[nnn].Pt() < 1.15

	       ) isRight3_ = true; 
	  }
	    
	  Bool_t isRight_ = isRight1_ && isRight2_ && isRight3_;

	  //	  std::cout << "\t\t check1:" << isRight1_ << " " << isRight2_ << " " << isRight3_ << std::endl;
	    
	  if(isRight_){
	    isRight1 = true;
	    pid1 = ppdgId[mmm];
	    matched_gentaupt1 = vec_gentau3pp4[mmm].Pt();
	  }
	}
	  
	/////////////////////////////////////////


	for(unsigned int mmm=0; mmm < gps.size(); mmm++){
	    
	  Bool_t isRight1_ = false;
	  Bool_t isRight2_ = false;
	  Bool_t isRight3_ = false;
	    
	  std::vector<TLorentzVector> tlvs = gps[mmm];
	    
	  for(unsigned int nnn=0; nnn < tlvs.size(); nnn++){

//	    std::cout << "\t\t gen. Nr" << mmm << " (pT, eta, phi) = " << tlvs[nnn].Pt() << ", " << tlvs[nnn].Eta()  << ", " << tlvs[nnn].Phi() << std::endl;
//	    std::cout << "\t\t pi1 index = " << tau2_idx1 << ", with (pt,eta,phi) = " << tau2_pi1.Pt() << ", " << tau2_pi1.Eta() << ", " << tau2_pi1.Phi() <<  std::endl;
//	    std::cout << "\t\t pi2 index = " << tau2_idx2 << ", with (pt,eta,phi) = " << tau2_pi2.Pt() << ", " << tau2_pi2.Eta() << ", " << tau2_pi2.Phi() <<  std::endl;
//	    std::cout << "\t\t pi3 index = " << tau2_idx3 << ", with (pt,eta,phi) = " << tau2_pi3.Pt() << ", " << tau2_pi3.Eta() << ", " << tau2_pi3.Phi() <<  std::endl;
//	    std::cout << "\t\t pt response = " << (tau2_pi1.Pt()/tlvs[nnn].Pt()) << ", delta_R = " << reco::deltaR(tau2_pi1.Eta(), tau2_pi1.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi())  << std::endl;
//	    std::cout << "\t\t pt response = " << (tau2_pi2.Pt()/tlvs[nnn].Pt()) << ", delta_R = " << reco::deltaR(tau2_pi2.Eta(), tau2_pi2.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi())  << std::endl;
//	    std::cout << "\t\t pt response = " << (tau2_pi3.Pt()/tlvs[nnn].Pt()) << ", delta_R = " << reco::deltaR(tau2_pi3.Eta(), tau2_pi3.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi())  << std::endl;

	    if(
	       //	       reco::deltaR(tau2_pi1.Eta(), tau2_pi1.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.1
	       reco::deltaR(tau2_pi1.Eta(), tau2_pi1.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.015 &&
	       tau2_pi1.Pt()/tlvs[nnn].Pt() > 0.85 && 
	       tau2_pi1.Pt()/tlvs[nnn].Pt() < 1.15

	       ) isRight1_ = true; 
	      
	    if(
	       //	       reco::deltaR(tau2_pi2.Eta(), tau2_pi2.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.1
	       reco::deltaR(tau2_pi2.Eta(), tau2_pi2.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.015 &&
	       tau2_pi2.Pt()/tlvs[nnn].Pt() > 0.85 && 
	       tau2_pi2.Pt()/tlvs[nnn].Pt() < 1.15

	       ) isRight2_ = true; 
	      
	    if(
	       //	       reco::deltaR(tau2_pi3.Eta(), tau2_pi3.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.1
	       reco::deltaR(tau2_pi3.Eta(), tau2_pi3.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.015 &&
	       tau2_pi3.Pt()/tlvs[nnn].Pt() > 0.85 && 
	       tau2_pi3.Pt()/tlvs[nnn].Pt() < 1.15

	       ) isRight3_ = true; 
	  }
	    
	  Bool_t isRight_ = isRight1_ && isRight2_ && isRight3_;

	  //	  std::cout << "\t\t check2:" << isRight1_ << " " << isRight2_ << " " << isRight3_ << std::endl;
	    
	  if(isRight_){
	    isRight2 = true;
	    pid2 = ppdgId[mmm];
	    matched_gentaupt2 = vec_gentau3pp4[mmm].Pt();
	  }
	}

      }


      if(isTruth_ && isMC){
	//	std::cout << "\t\t isRight1, isRigh2 = " << isRight1 << " " << isRight2 << std::endl;
	if(!(isRight1 && isRight2)) continue;
      }

      //      std::cout << "check14" << std::endl;     	

      /********************************************************************
       *
       * Step9: Filling normal branches
       *
       ********************************************************************/

      // if the candidate has less than 2 GeV, remove it.
      //      if(cands[ic].cand_tau_fullfit_pt < 2.) continue;

    
      nBranches_->BsTauTauFH_tau1_pt.push_back(tlv_tau1_fit.Pt());
      nBranches_->BsTauTauFH_tau1_eta.push_back(tlv_tau1_fit.Eta());
      nBranches_->BsTauTauFH_tau1_phi.push_back(tlv_tau1_fit.Phi());
      nBranches_->BsTauTauFH_tau1_mass.push_back(tlv_tau1_fit.M());
    
      std::vector<Float_t> rhomass1;
      pat::PackedCandidate tau1_pf1 = pfcollection[tau1_idx1];
      pat::PackedCandidate tau1_pf2 = pfcollection[tau1_idx2];
      pat::PackedCandidate tau1_pf3 = pfcollection[tau1_idx3];
    
      TLorentzVector tlv1_pion1; 
      TLorentzVector tlv1_pion2;
      TLorentzVector tlv1_pion3;
	
      tlv1_pion1.SetPtEtaPhiM(tau1_pf1.pt(), tau1_pf1.eta(), tau1_pf1.phi(), tau1_pf1.mass());
      tlv1_pion2.SetPtEtaPhiM(tau1_pf2.pt(), tau1_pf2.eta(), tau1_pf2.phi(), tau1_pf2.mass());
      tlv1_pion3.SetPtEtaPhiM(tau1_pf3.pt(), tau1_pf3.eta(), tau1_pf3.phi(), tau1_pf3.mass());
	
	
      if(tau1_pf1.charge()*tau1_pf2.charge() == -1){
	TLorentzVector tlv_rho = tlv1_pion1 + tlv1_pion2;
	rhomass1.push_back(tlv_rho.M());
      }
	
      if(tau1_pf1.charge()*tau1_pf3.charge() == -1){
	TLorentzVector tlv_rho = tlv1_pion1 + tlv1_pion3;
	rhomass1.push_back(tlv_rho.M());
      }
	
      if(tau1_pf2.charge()*tau1_pf3.charge() == -1){
	TLorentzVector tlv_rho = tlv1_pion2 + tlv1_pion3;
	rhomass1.push_back(tlv_rho.M());
      }

      //      std::cout << "check15" << std::endl;     		
      //	std::cout << "rho masses size = " << rhomass.size() << std::endl;
      if(rhomass1.size()==2){    
	nBranches_->BsTauTauFH_tau1_rhomass1.push_back(rhomass1.at(0));
	nBranches_->BsTauTauFH_tau1_rhomass2.push_back(rhomass1.at(1));
      }else{
	nBranches_->BsTauTauFH_tau1_rhomass1.push_back(rhomass1.at(-1));
	nBranches_->BsTauTauFH_tau1_rhomass2.push_back(rhomass1.at(-1));
      }

      nBranches_->BsTauTauFH_tau1_q.push_back(cands[leg1].cand_tau_charge);
      nBranches_->BsTauTauFH_tau1_vx.push_back(tau_vertex1->vertexState().position().x());
      nBranches_->BsTauTauFH_tau1_vy.push_back(tau_vertex1->vertexState().position().y());
      nBranches_->BsTauTauFH_tau1_vz.push_back(tau_vertex1->vertexState().position().z());
    
      nBranches_->BsTauTauFH_tau1_max_dr_3prong.push_back(max_dr1);
      nBranches_->BsTauTauFH_tau1_lip.push_back(Taucand1.lip);
      nBranches_->BsTauTauFH_tau1_lips.push_back(Taucand1.lips);
      nBranches_->BsTauTauFH_tau1_pvip.push_back(Taucand1.pvip);
      nBranches_->BsTauTauFH_tau1_pvips.push_back(Taucand1.pvips);
      nBranches_->BsTauTauFH_tau1_fl3d.push_back(Taucand1.fl3d);
      nBranches_->BsTauTauFH_tau1_fls3d.push_back(Taucand1.fls3d);
      nBranches_->BsTauTauFH_tau1_alpha.push_back(Taucand1.alpha);
      nBranches_->BsTauTauFH_tau1_vprob.push_back(TMath::Prob(tau_vertex1->chiSquared(), tau_vertex1->degreesOfFreedom()));
      nBranches_->BsTauTauFH_tau1_isRight.push_back(isRight1);
      nBranches_->BsTauTauFH_tau1_matched_ppdgId.push_back(pid1);
      nBranches_->BsTauTauFH_tau1_matched_gentaupt.push_back(matched_gentaupt1);
      nBranches_->BsTauTauFH_tau1_sumofdnn.push_back(mydnn[tau1_idx1] + mydnn[tau1_idx2] + mydnn[tau1_idx3]);
      nBranches_->BsTauTauFH_tau1_pfidx1.push_back(tau1_idx1);
      nBranches_->BsTauTauFH_tau1_pfidx2.push_back(tau1_idx2);
      nBranches_->BsTauTauFH_tau1_pfidx3.push_back(tau1_idx3);
      nBranches_->BsTauTauFH_tau1_pi1_dnn.push_back(mydnn[tau1_idx1]);
      nBranches_->BsTauTauFH_tau1_pi2_dnn.push_back(mydnn[tau1_idx2]);
      nBranches_->BsTauTauFH_tau1_pi3_dnn.push_back(mydnn[tau1_idx3]);

      nBranches_->BsTauTauFH_tau1_pi1_pt.push_back(tlv1_pion1.Pt());
      nBranches_->BsTauTauFH_tau1_pi1_eta.push_back(tlv1_pion1.Eta());
      nBranches_->BsTauTauFH_tau1_pi1_phi.push_back(tlv1_pion1.Phi());
      nBranches_->BsTauTauFH_tau1_pi1_mass.push_back(tlv1_pion1.M());

      nBranches_->BsTauTauFH_tau1_pi2_pt.push_back(tlv1_pion2.Pt());
      nBranches_->BsTauTauFH_tau1_pi2_eta.push_back(tlv1_pion2.Eta());
      nBranches_->BsTauTauFH_tau1_pi2_phi.push_back(tlv1_pion2.Phi());
      nBranches_->BsTauTauFH_tau1_pi2_mass.push_back(tlv1_pion2.M());

      nBranches_->BsTauTauFH_tau1_pi3_pt.push_back(tlv1_pion3.Pt());
      nBranches_->BsTauTauFH_tau1_pi3_eta.push_back(tlv1_pion3.Eta());
      nBranches_->BsTauTauFH_tau1_pi3_phi.push_back(tlv1_pion3.Phi());
      nBranches_->BsTauTauFH_tau1_pi3_mass.push_back(tlv1_pion3.M());
    

      //      std::cout << "check16" << std::endl;     		    


      nBranches_->BsTauTauFH_tau2_pt.push_back(tlv_tau2_fit.Pt());
      nBranches_->BsTauTauFH_tau2_eta.push_back(tlv_tau2_fit.Eta());
      nBranches_->BsTauTauFH_tau2_phi.push_back(tlv_tau2_fit.Phi());
      nBranches_->BsTauTauFH_tau2_mass.push_back(tlv_tau2_fit.M());
    
      std::vector<Float_t> rhomass2;
      pat::PackedCandidate tau2_pf1 = pfcollection[tau2_idx1];
      pat::PackedCandidate tau2_pf2 = pfcollection[tau2_idx2];
      pat::PackedCandidate tau2_pf3 = pfcollection[tau2_idx3];
    
      TLorentzVector tlv2_pion1; 
      TLorentzVector tlv2_pion2;
      TLorentzVector tlv2_pion3;
    
      tlv2_pion1.SetPtEtaPhiM(tau2_pf1.pt(), tau2_pf1.eta(), tau2_pf1.phi(), tau2_pf1.mass());
      tlv2_pion2.SetPtEtaPhiM(tau2_pf2.pt(), tau2_pf2.eta(), tau2_pf2.phi(), tau2_pf2.mass());
      tlv2_pion3.SetPtEtaPhiM(tau2_pf3.pt(), tau2_pf3.eta(), tau2_pf3.phi(), tau2_pf3.mass());
    
    
      if(tau2_pf1.charge()*tau2_pf2.charge() == -1){
	TLorentzVector tlv_rho = tlv2_pion1 + tlv2_pion2;
	rhomass2.push_back(tlv_rho.M());
      }
    
      if(tau2_pf1.charge()*tau2_pf3.charge() == -1){
	TLorentzVector tlv_rho = tlv2_pion1 + tlv2_pion3;
	rhomass2.push_back(tlv_rho.M());
      }
    
      if(tau2_pf2.charge()*tau2_pf3.charge() == -1){
	TLorentzVector tlv_rho = tlv2_pion2 + tlv2_pion3;
	rhomass2.push_back(tlv_rho.M());
      }
    
      //	std::cout << "rho masses size = " << rhomass.size() << std::endl;


      if(rhomass1.size()==2){    
	nBranches_->BsTauTauFH_tau2_rhomass1.push_back(rhomass2.at(0));
	nBranches_->BsTauTauFH_tau2_rhomass2.push_back(rhomass2.at(1));
      }else{
	nBranches_->BsTauTauFH_tau2_rhomass1.push_back(rhomass2.at(-1));
	nBranches_->BsTauTauFH_tau2_rhomass2.push_back(rhomass2.at(-1));
      }

      nBranches_->BsTauTauFH_tau2_q.push_back(cands[leg2].cand_tau_charge);
      nBranches_->BsTauTauFH_tau2_vx.push_back(tau_vertex2->vertexState().position().x());
      nBranches_->BsTauTauFH_tau2_vy.push_back(tau_vertex2->vertexState().position().y());
      nBranches_->BsTauTauFH_tau2_vz.push_back(tau_vertex2->vertexState().position().z());
    
      nBranches_->BsTauTauFH_tau2_max_dr_3prong.push_back(max_dr2);
      nBranches_->BsTauTauFH_tau2_lip.push_back(Taucand2.lip);
      nBranches_->BsTauTauFH_tau2_lips.push_back(Taucand2.lips);
      nBranches_->BsTauTauFH_tau2_pvip.push_back(Taucand2.pvip);
      nBranches_->BsTauTauFH_tau2_pvips.push_back(Taucand2.pvips);
      nBranches_->BsTauTauFH_tau2_fl3d.push_back(Taucand2.fl3d);
      nBranches_->BsTauTauFH_tau2_fls3d.push_back(Taucand2.fls3d);
      nBranches_->BsTauTauFH_tau2_alpha.push_back(Taucand2.alpha);
      nBranches_->BsTauTauFH_tau2_vprob.push_back(TMath::Prob(tau_vertex2->chiSquared(), tau_vertex2->degreesOfFreedom()));
      nBranches_->BsTauTauFH_tau2_isRight.push_back(isRight2);
      nBranches_->BsTauTauFH_tau2_matched_ppdgId.push_back(pid2);
      nBranches_->BsTauTauFH_tau2_matched_gentaupt.push_back(matched_gentaupt2);
      nBranches_->BsTauTauFH_tau2_sumofdnn.push_back(mydnn[tau2_idx1] + mydnn[tau2_idx2] + mydnn[tau2_idx3]);
      nBranches_->BsTauTauFH_tau2_pfidx1.push_back(tau2_idx1);
      nBranches_->BsTauTauFH_tau2_pfidx2.push_back(tau2_idx2);
      nBranches_->BsTauTauFH_tau2_pfidx3.push_back(tau2_idx3);
      nBranches_->BsTauTauFH_tau2_pi1_dnn.push_back(mydnn[tau2_idx1]);
      nBranches_->BsTauTauFH_tau2_pi2_dnn.push_back(mydnn[tau2_idx2]);
      nBranches_->BsTauTauFH_tau2_pi3_dnn.push_back(mydnn[tau2_idx3]);
    
      nBranches_->BsTauTauFH_tau2_pi1_pt.push_back(tlv2_pion1.Pt());
      nBranches_->BsTauTauFH_tau2_pi1_eta.push_back(tlv2_pion1.Eta());
      nBranches_->BsTauTauFH_tau2_pi1_phi.push_back(tlv2_pion1.Phi());
      nBranches_->BsTauTauFH_tau2_pi1_mass.push_back(tlv2_pion1.M());

      nBranches_->BsTauTauFH_tau2_pi2_pt.push_back(tlv2_pion2.Pt());
      nBranches_->BsTauTauFH_tau2_pi2_eta.push_back(tlv2_pion2.Eta());
      nBranches_->BsTauTauFH_tau2_pi2_phi.push_back(tlv2_pion2.Phi());
      nBranches_->BsTauTauFH_tau2_pi2_mass.push_back(tlv2_pion2.M());

      nBranches_->BsTauTauFH_tau2_pi3_pt.push_back(tlv2_pion3.Pt());
      nBranches_->BsTauTauFH_tau2_pi3_eta.push_back(tlv2_pion3.Eta());
      nBranches_->BsTauTauFH_tau2_pi3_phi.push_back(tlv2_pion3.Phi());
      nBranches_->BsTauTauFH_tau2_pi3_mass.push_back(tlv2_pion3.M());

    
      nBranches_->BsTauTauFH_B_pt.push_back(bc_part->currentState().globalMomentum().perp());
      nBranches_->BsTauTauFH_B_eta.push_back(bc_part->currentState().globalMomentum().eta());
      nBranches_->BsTauTauFH_B_phi.push_back(bc_part->currentState().globalMomentum().phi());
      nBranches_->BsTauTauFH_B_mass.push_back(bc_part->currentState().mass());
      nBranches_->BsTauTauFH_B_vprob.push_back(TMath::Prob(bc_vertex->chiSquared(), bc_vertex->degreesOfFreedom()));
      nBranches_->BsTauTauFH_B_lip.push_back(Bcand.lip);
      nBranches_->BsTauTauFH_B_lips.push_back(Bcand.lips);
      nBranches_->BsTauTauFH_B_pvip.push_back(Bcand.pvip);
      nBranches_->BsTauTauFH_B_pvips.push_back(Bcand.pvips);
      nBranches_->BsTauTauFH_B_fls3d.push_back(Bcand.fls3d);
      nBranches_->BsTauTauFH_B_fl3d.push_back(Bcand.fl3d);
      nBranches_->BsTauTauFH_B_alpha.push_back(Bcand.alpha);
    
      //	allParticles4doc.push_back(pFactory.particle(tt_muon, muon_mass, chi, ndf, muon_sigma));
      //	allParticles4doc.push_back(pFactory.particle(mytracks[_idx1], pion_mass, chi, ndf, pion_sigma));
      //	allParticles4doc.push_back(pFactory.particle(mytracks[_idx2], pion_mass, chi, ndf, pion_sigma));
      //	allParticles4doc.push_back(pFactory.particle(mytracks[_idx3], pion_mass, chi, ndf, pion_sigma));


      nBranches_->BsTauTauFH_B_maxdoca.push_back(-1);
      nBranches_->BsTauTauFH_B_mindoca.push_back(-1);
      nBranches_->BsTauTauFH_B_vx.push_back(bc_vertex->vertexState().position().x());
      nBranches_->BsTauTauFH_B_vy.push_back(bc_vertex->vertexState().position().y());
      nBranches_->BsTauTauFH_B_vz.push_back(bc_vertex->vertexState().position().z());
    
      nBranches_->BsTauTauFH_B_iso.push_back(iso);
      nBranches_->BsTauTauFH_B_iso_ntracks.push_back(ntracks);
      nBranches_->BsTauTauFH_B_iso_mindoca.push_back(iso_mindoca);

      ncomb++;
    }
  }

  //  if(ncomb==0) return false;
  nBranches_->cutflow_perevt->Fill(4);
    
    
  nBranches_->BsTauTauFH_mu1_pt.push_back(muoncollection[0].pt());
  nBranches_->BsTauTauFH_mu1_eta.push_back(muoncollection[0].eta());
  nBranches_->BsTauTauFH_mu1_phi.push_back(muoncollection[0].phi());
  nBranches_->BsTauTauFH_mu1_mass.push_back(muoncollection[0].mass());
  nBranches_->BsTauTauFH_mu1_q.push_back(muoncollection[0].charge());
  nBranches_->BsTauTauFH_mu1_isLoose.push_back(muoncollection[0].isLooseMuon());
  nBranches_->BsTauTauFH_mu1_isTight.push_back(muoncollection[0].isTightMuon(closestVertex));
  nBranches_->BsTauTauFH_mu1_isPF.push_back(muoncollection[0].isPFMuon());
  nBranches_->BsTauTauFH_mu1_isGlobal.push_back(muoncollection[0].isGlobalMuon());
  nBranches_->BsTauTauFH_mu1_isTracker.push_back(muoncollection[0].isTrackerMuon());
  nBranches_->BsTauTauFH_mu1_isSoft.push_back(muoncollection[0].isSoftMuon(closestVertex));
  nBranches_->BsTauTauFH_mu1_vx.push_back(muoncollection[0].vx());
  nBranches_->BsTauTauFH_mu1_vy.push_back(muoncollection[0].vy());
  nBranches_->BsTauTauFH_mu1_vz.push_back(muoncollection[0].vz());
  nBranches_->BsTauTauFH_mu1_iso.push_back(1.);
  nBranches_->BsTauTauFH_mu1_dbiso.push_back(aux.MuonPFIso(muoncollection[0]));
    
  nBranches_->BsTauTauFH_PV_vx.push_back(vertices_->begin()->position().x());
  nBranches_->BsTauTauFH_PV_vy.push_back(vertices_->begin()->position().y());
  nBranches_->BsTauTauFH_PV_vz.push_back(vertices_->begin()->position().z());
    
  //      if(myVertex.isValid()){
  //	nBranches_->BsTauTauFH_bbPV_refit_vx.push_back(myVertex.position().x());
  //	nBranches_->BsTauTauFH_bbPV_refit_vy.push_back(myVertex.position().y());
  //	nBranches_->BsTauTauFH_bbPV_refit_vz.push_back(myVertex.position().z());
  //      }else{
  nBranches_->BsTauTauFH_bbPV_refit_vx.push_back(-1);
  nBranches_->BsTauTauFH_bbPV_refit_vy.push_back(-1);
  nBranches_->BsTauTauFH_bbPV_refit_vz.push_back(-1);
  //      }
    
  nBranches_->BsTauTauFH_bbPV_vx.push_back(closestVertex.position().x());
  nBranches_->BsTauTauFH_bbPV_vy.push_back(closestVertex.position().y());
  nBranches_->BsTauTauFH_bbPV_vz.push_back(closestVertex.position().z());
    

  //////////////////////////////




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
      if(TMath::Abs((*genParticles_)[p].pdgId())==531 && (*genParticles_)[p].status()==2){
	  
	// retrieve production vertex
	genvertex = aux.getVertex((*genParticles_)[p]);
	  
	for(int idd = 0; idd < (int)(*genParticles_)[p].numberOfDaughters(); idd++){
	  Int_t dpid = (*genParticles_)[p].daughter(idd)->pdgId();
	  //	  std::cout << "\t -> " <<  << " " << (*genParticles_)[p].daughter(idd)->status()<< std::endl;
	  if(TMath::Abs(dpid)==15) gen_nr_mu.push_back((*genParticles_)[p].daughter(idd));
	}
      }
	
	
    }
  }
    
  

  
  // -9 if there is no Bc found 
  nBranches_->BsTauTauFH_genPV_vx.push_back(genvertex.x());
  nBranches_->BsTauTauFH_genPV_vy.push_back(genvertex.y());
  nBranches_->BsTauTauFH_genPV_vz.push_back(genvertex.z());
  nBranches_->BsTauTauFH_ngenmuons.push_back(gen_nr_mu.size());
    
    
    
  nBranches_->BsTauTauFH_isgen3.push_back(isgen3);
  nBranches_->BsTauTauFH_isgen3matched.push_back(isgen3matched);
  nBranches_->BsTauTauFH_nch.push_back(numOfch);
  nBranches_->BsTauTauFH_nch_after_dnn.push_back(npf_after_dnn);
  nBranches_->BsTauTauFH_nch_before_dnn.push_back(npf_before_dnn);
  nBranches_->BsTauTauFH_nch_qr.push_back(npf_qr);
  nBranches_->BsTauTauFH_ngentau3.push_back(gps.size());
  nBranches_->BsTauTauFH_ngentau.push_back(vec_gentaudm.size());

  if(vec_gentaudm.size() >=1){

    for(size_t ii = 0; ii < vec_gentaup4.size(); ii++){
      nBranches_->BsTauTauFH_gentaupt.push_back(vec_gentaup4[ii].Pt());
      nBranches_->BsTauTauFH_gentaudm.push_back(vec_gentaudm[ii]);
    }

  }else{
    nBranches_->BsTauTauFH_gentaupt.push_back(-1);
    nBranches_->BsTauTauFH_gentaudm.push_back(-1);
  }

  nBranches_->IsBsTauTauFH.push_back(1.);
  nBranches_->BsTauTauFH_nCandidates.push_back(ncomb);
  nBranches_->BsTauTauFH_ntaus.push_back(cands.size());

  //  nBranches_->cutflow_perevt->Fill(5);
  return true;


}


//void BsTauTauFHNtuplizer::printout(const RefCountedKinematicVertex& myVertex){
//  std::cout << "Vertex:" << std::endl;
//  if (myVertex->vertexIsValid()) {
//    std::cout << "\t Decay vertex: " << myVertex->position() << myVertex->chiSquared() << " " << myVertex->degreesOfFreedom()
//	      << std::endl;
//  } else
//    std::cout << "\t Decay vertex Not valid\n";
//}
//
//void BsTauTauFHNtuplizer::printout(const RefCountedKinematicParticle& myParticle){
//  std::cout << "Particle:" << std::endl;
//  //accessing the reconstructed Bs meson parameters:
//  //SK: uncomment if needed  AlgebraicVector7 bs_par = myParticle->currentState().kinematicParameters().vector();
//
//  //and their joint covariance matrix:
//  //SK:uncomment if needed  AlgebraicSymMatrix77 bs_er = myParticle->currentState().kinematicParametersError().matrix();
//  std::cout << "\t Momentum at vertex: " << myParticle->currentState().globalMomentum() << std::endl;
//  std::cout << "\t Parameters at vertex: " << myParticle->currentState().kinematicParameters().vector() << std::endl;
//}
//
//void BsTauTauFHNtuplizer::printout(const RefCountedKinematicTree& myTree){
//  if (!myTree->isValid()) {
//    std::cout << "Tree is invalid. Fit failed.\n";
//    return;
//  }
//
//  //accessing the tree components, move pointer to top
//  myTree->movePointerToTheTop();
//
//  //We are now at the top of the decay tree getting the B_s reconstructed KinematicPartlcle
//  RefCountedKinematicParticle b_s = myTree->currentParticle();
//  printout(b_s);
//
//  // The B_s decay vertex
//  RefCountedKinematicVertex b_dec_vertex = myTree->currentDecayVertex();
//  printout(b_dec_vertex);
//
//  // Get all the children of Bs:
//  //In this way, the pointer is not moved
//  std::vector<RefCountedKinematicParticle> bs_children = myTree->finalStateParticles();
//
//  for (unsigned int i = 0; i < bs_children.size(); ++i) {
//    printout(bs_children[i]);
//  }
//
//  std::cout << "\t ------------------------------------------" << std::endl;
//
//  //Now navigating down the tree , pointer is moved:
//  bool child = myTree->movePointerToTheFirstChild();
//
//  if (child)
//    while (myTree->movePointerToTheNextChild()) {
//      RefCountedKinematicParticle aChild = myTree->currentParticle();
//      printout(aChild);
//    }
//}
