#include "../interface/JpsiMuNtuplizer.h"


//===================================================================================================================
JpsiMuNtuplizer::JpsiMuNtuplizer( edm::EDGetTokenT<pat::MuonCollection>    muonToken   ,
                                  edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
                                  edm::EDGetTokenT<reco::BeamSpot> bsToken,
                                  edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken,
                                  edm::EDGetTokenT<pat::PackedCandidateCollection> losttrackToken,
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
    , losttrackToken_(losttrackToken) 
    , HLTtriggersToken_	( triggertoken )
    , triggerObjects_	( triggerobject )
    , genParticlesToken_( genptoken )
    , runOnMC_   (runFlags["runOnMC"])
    , doCutFlow_ (runFlags["doCutFlow"])
   
{
}

//===================================================================================================================
JpsiMuNtuplizer::~JpsiMuNtuplizer( void )
{

}


TVector3 JpsiMuNtuplizer::getVertex(const reco::GenParticle& part){
    return TVector3(part.vx(),part.vy(),part.vz());
}

float JpsiMuNtuplizer::MuonPFIso(pat::Muon muon){

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


//adapt absoluteImpactParameter functionality for RefCountedKinematicVertex
std::pair<bool, Measurement1D> JpsiMuNtuplizer::absoluteImpactParameter(const TrajectoryStateOnSurface& tsos,
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


math::PtEtaPhiMLorentzVector JpsiMuNtuplizer::daughter_p4(std::vector< RefCountedKinematicParticle > fitted_children, size_t i){
  const auto& state = fitted_children.at(i)->currentState();

  return math::PtEtaPhiMLorentzVector(
				      state.globalMomentum().perp(), 
				      state.globalMomentum().eta() ,
				      state.globalMomentum().phi() ,
				      state.mass()
				      );
}


void JpsiMuNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  
    // std::cout << "---------------- event, run, lumi = " << event.id().event() << " " << event.id().run() << " " << event.id().luminosityBlock() << "----------------" << std::endl;
  
    /********************************************************************
     *
     * Step1: check if the J/psi trigger is fired.
     * Namely, HLT_DoubleMu4_JpsiTrk_Displaced_v
     * and  HLT_Dimuon0_Jpsi3p5_Muon2_v
     ********************************************************************/

    event.getByToken(HLTtriggersToken_, HLTtriggers_);

    bool isTriggered = false;
    const edm::TriggerNames& trigNames = event.triggerNames(*HLTtriggers_);
    std::string finalTriggerName="";
    std::string finalTriggerFilterObjName="";

    for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {

        // if(trigNames.triggerName(i).find("HLT_DoubleMu4_JpsiTrk_Displaced_v")!= std::string::npos || trigNames.triggerName(i).find("HLT_Dimuon0_Jpsi3p5_Muon2_v")!= std::string::npos ){
           
        //     nBranches_->HLT_BPH_isFired[trigNames.triggerName(i)] = HLTtriggers_->accept(i);
        if(trigNames.triggerName(i).find("HLT_DoubleMu4_JpsiTrk_Displaced_v")!= std::string::npos){
            nBranches_->HLT_BPH_isFired[trigNames.triggerName(i)] = HLTtriggers_->accept(i);
            if(HLTtriggers_->accept(i)){
                isTriggered = true;
                finalTriggerName=trigNames.triggerName(i);  
                finalTriggerFilterObjName="hltJpsiTkVertexFilter";
                // std::cout << "finalTriggerName = "  << finalTriggerName << std::endl;
               
            }
        }

    }
    // Second trigger if fist one didn't fire
    if (!isTriggered){
        for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {
            if(trigNames.triggerName(i).find("HLT_Dimuon0_Jpsi3p5_Muon2_v")!= std::string::npos){
                nBranches_->HLT_BPH_isFired[trigNames.triggerName(i)] = HLTtriggers_->accept(i);
                if(HLTtriggers_->accept(i)){
                    isTriggered = true;
                    finalTriggerName=trigNames.triggerName(i);  
                    finalTriggerFilterObjName="hltVertexmumuFilterJpsiMuon3p5";
                    // std::cout << "finalTriggerName = "  << finalTriggerName << std::endl;
                    
                }
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
    std::vector<int> muoncollection_id;
    muoncollection.clear();
    muoncollection_id.clear();

    // bools keeping track of muon 1 or 2 for cutflow purposes
    bool mu1passpteta = false;
    bool mu2passpteta = false;
    bool mu3passpteta = false;
    bool mu1passtrigmatch = false;
    bool mu2passtrigmatch = false;
 // Precut
    if( doCutFlow_) {
       nBranches_->cutflow_perevt->Fill(0);
      
       for(size_t imuon = 0; imuon < muons_->size(); ++ imuon){

           const pat::Muon & muon = (*muons_)[imuon];
           // right now we're just checking that there are at least two muons passing our kinematic cuts
           // so we don't care if there are more than two muons that pass right now
           if(mu2passpteta) continue;
           if(muon.pt() < 4) continue;
           if(TMath::Abs(muon.eta()) > 2.4) continue;
           if(!(muon.track().isNonnull())) continue;
           // If we're on the Second muon, then the second muon is what passed 
           if (mu1passpteta && !mu2passpteta) { mu2passpteta = true;}
           // If we're on the First muon, then the first muon is what passed
           if (!mu1passpteta) {mu1passpteta = true;}
           }

       if (mu1passpteta) {nBranches_->cutflow_perevt->Fill(1);}
       if (mu2passpteta) {nBranches_->cutflow_perevt->Fill(2);}
    }
    if(!isTriggered) return;
 // evt Triggered
    if( doCutFlow_) {
       nBranches_->cutflow_perevt->Fill(3);
      }
//Now time to begin the cutflow in earnest
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
	
                if(obj.filterLabels()[hh].find(finalTriggerFilterObjName) != std::string::npos){
                    isFilterExist = true;
                }
            }
      
            if(!isFilterExist) continue;
      
            Float_t trigger_dR = reco::deltaR(obj.eta(), obj.phi(),
                                              muon.eta(), muon.phi());
      
            if(trigger_dR < 0.1) trigMatch = true;
        }

        if(!trigMatch) continue;
        // If we're on the Second muon, then the second muon is what passed 
        if (mu1passtrigmatch && !mu2passtrigmatch) { mu2passtrigmatch = true;}
        // If we're on the First muon, then the first muon is what passed
        if (!mu1passtrigmatch) {mu1passtrigmatch = true;}
        muoncollection.push_back(muon);
        muoncollection_id.push_back(imuon);
    }
  if( doCutFlow_) {
    if (mu1passtrigmatch) {nBranches_->cutflow_perevt->Fill(4);}
    if (mu2passtrigmatch) {nBranches_->cutflow_perevt->Fill(5);}
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
            // Jpsi mass cut passed
      
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

    if ( doCutFlow_) {
       nBranches_->cutflow_perevt->Fill(6);
       }


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

    std::vector< RefCountedKinematicParticle > jpsi_children = jpTree->finalStateParticles();

    math::PtEtaPhiMLorentzVector mu1_fit = daughter_p4(jpsi_children, 0);
    math::PtEtaPhiMLorentzVector mu2_fit = daughter_p4(jpsi_children, 1);


//    std::cout << "before fit1: " << muoncollection[mcidx_mu1].pt() << " " << muoncollection[mcidx_mu1].eta() << " " << muoncollection[mcidx_mu1].phi() << " " << muoncollection[mcidx_mu1].mass() << " " << muoncollection[mcidx_mu1].charge()  << std::endl;
//    std::cout << "before fit2: " << muoncollection[mcidx_mu2].pt() << " " << muoncollection[mcidx_mu2].eta() << " " << muoncollection[mcidx_mu2].phi() << " " << muoncollection[mcidx_mu2].mass() << " " << muoncollection[mcidx_mu2].charge() << std::endl;
//
//    std::cout << "after fit1: " << mu1_fit.pt() << " " << mu1_fit.eta() << " " << mu1_fit.phi() << " " << mu1_fit.mass() << std::endl;
//    std::cout << "after fit2: " << mu2_fit.pt() << " " << mu2_fit.eta() << " " << mu2_fit.phi() << " " << mu2_fit.mass() << std::endl;

    // Jpsi kinematic fit passed
    if( doCutFlow_) {
      nBranches_->cutflow_perevt->Fill(7);
      }

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
     * and found that 1) or 2) is the best.
     * 
     * Note: we can refit the vertex each time, by excluding the muons from J/psi
     * but since the efficiency of selecting the right vertex is already quite high (w/o vertex refit)
     * and I don't think it is necessary to complicate things at this stage 
     * (we do the refit once we select B candidate)
     *
     ********************************************************************/

    // define extrapolator
    edm::ESHandle<MagneticField> fieldHandle;
    iSetup.get<IdealMagneticFieldRecord>().get(fieldHandle);
    fMagneticField = fieldHandle.product();

    AnalyticalImpactPointExtrapolator extrapolator(fMagneticField);

    Float_t max_criteria = 999;
    reco::Vertex closestVertex; 
    int counter = 0;

    for( reco::VertexCollection::const_iterator vtx = vertices_->begin(); vtx != vertices_->end(); ++vtx){

        //    std::cout << "vtx:" << vtx->position().z() << std::endl;
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

        counter += 1;
    }


  
    /********************************************************************
     *
     * Step6: 3rd muon selection
     *        Just select highest in pT but there might be better selection ... 
     *
     ********************************************************************/
    bool doMultipleMuon3 = true;  // false: to slect jus 1 J/psi, 1 mu3 and 1 B Meson candidate; true: to select multiple mu3 and B meson candidates

    std::vector<pat::Muon> mu3_vec;
    Float_t max_pt3 = -1;

    for(size_t imuon = 0; imuon < muons_->size(); ++ imuon){

        const pat::Muon & muon = (*muons_)[imuon];

        if(muon.pt() < 2) continue; // there was a 4 but maybe it's best to optimize it.
        if(TMath::Abs(muon.eta()) > 2.4) continue;
        if(!(muon.track().isNonnull())) continue;
        if(imuon==idx_mu1 || imuon==idx_mu2) continue;
        mu3passpteta = true;
        // shall we put delta z cut here maybe? 

	
	//	std::cout << "muon3 pt = " << muon.pt() << std::endl;
	
        max_pt3=2;
        if(muon.pt() > max_pt3){
          if(!doMultipleMuon3) { max_pt3 = muon.pt();}
            mu3_vec.push_back(muon);
            if (!doMultipleMuon3) break;
      
        }
    }

    if(max_pt3==-1) return;
    //mu3 passed pT and eta cut
    if( doCutFlow_ && mu3passpteta ) {
      nBranches_->cutflow_perevt->Fill(8);
      }

    //    int nmuo =0;
    for( auto mu3: mu3_vec){
      //        nmuo++;
     
        const reco::TrackRef track3_muon = mu3.muonBestTrack();
        reco::TransientTrack tt3_muon = (*builder).build(track3_muon);

        /********************************************************************
         *
         * Step7: bbPV vertex refit by excluding 3 muon candidates
         *
         ********************************************************************/

        event.getByToken( packedpfcandidatesToken_               , packedpfcandidates_      ); 
        event.getByToken( losttrackToken_               , losttrack_      ); 

        std::vector<pat::PackedCandidate> mypfcollection; 
        std::vector<reco::TransientTrack> mytracks;

        for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){   
    
            pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
    
            //    if(pf.pt() < 0.5) continue;
    
            //    std::cout << pf.pdgId() << std::endl;
            if(!pf.hasTrackDetails()) continue;
            //    std::cout << "\t -> passed" << std::endl;

            //    bool ismatch = false;
            //    for( reco::VertexCollection::const_iterator vtx = vertices_->begin(); vtx != vertices_->end(); ++vtx){
            //      //      std::cout << counter << " " << vtx->position().z() << std::endl;
            //      if(pf.vertexRef()->z()==vtx->position().z()) ismatch = true;
            //    }
            //
            //    if(!ismatch) std::cout << "This does not belong any PV !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;

            //    std::cout << pf.vertexRef()->z() << " " << closestVertex.position().z() << std::endl;
            if(pf.vertexRef()->z()!=closestVertex.position().z()) continue;

            Float_t _dR1 = reco::deltaR(pf.eta(), pf.phi(), 
                                        track1_muon->eta(), track1_muon->phi());

            Float_t _dR2 = reco::deltaR(pf.eta(), pf.phi(), 
                                        track2_muon->eta(), track2_muon->phi());

            Float_t _dR3 = reco::deltaR(pf.eta(), pf.phi(), 
                                        track3_muon->eta(), track3_muon->phi());


            if(_dR1 < 0.1 || _dR2 < 0.1 || _dR3 < 0.1){
                //      std::cout << "This track is matched to the muon: " << pf.pdgId() << " " << pf.pt() << " " << track1_muon->pt() << " " << track2_muon->pt() << " " << track3_muon->pt() << std::endl;

                if(TMath::Abs(pf.pdgId()) == 13) continue;
            }

            reco::TransientTrack  tt_track = (*builder).build(pf.pseudoTrack());
            //    transientTrack.setBeamSpot(beamspot_);
            mytracks.push_back(tt_track);
    
            //   Float_t precut_dz = closestVertex.position().z() - pf.vz();
            //    if(TMath::Abs(precut_dz) > 0.5) continue;
            //      Float_t precut_dz = vz_jpsi - pf.vz();
            //      if(TMath::Abs(precut_dz) > 0.5) continue;

            //    if(!pf.trackHighPurity()) continue;    
            //    if(pf.pseudoTrack().hitPattern().numberOfValidHits() < 3) continue;
            //    if(pf.pseudoTrack().normalizedChi2() > 100) continue;
    
            //    if(TMath::Abs(pf.pdgId())==211 && TMath::Abs(pf.eta()) < 2.3){
            if(pf.pt() > 0.5) mypfcollection.push_back(pf);
            //    }
        }



        TransientVertex myVertex; 
  
        if( mytracks.size()>1 ){
            AdaptiveVertexFitter  theFitter;
    
            myVertex = theFitter.vertex(mytracks);  // if you don't want the beam constraint
            //    TransientVertex myVertex = theFitter.vertex(mytracks, beamspot_);  // if you don't want the beam constraint
            //    std::cout << "TEST   " << closestVertex.position().z() << " " << myVertex.position().z() << std::endl;

        }
  







        //  for( size_t ii = 0; ii < losttrack_->size(); ++ii ){   
        //    
        //    pat::PackedCandidate pf = (*losttrack_)[ii];
        //    
        //    //    if(pf.pt() < 0.5) continue;
        //    
        //    if(!pf.hasTrackDetails()) continue;
        //
        //    std::cout << "lost track" << pf.vertexRef()->z() << " " << closestVertex.position().z() << std::endl;
        //    if(pf.vertexRef()->z()!=closestVertex.position().z()) continue;
        //
        //  }





        /********************************************************************
         *
         * Step7: Kinematic fit for B candidate
         *
         ********************************************************************/

        std::vector<RefCountedKinematicParticle> allParticles;

        allParticles.push_back(pFactory.particle(tt3_muon, muon_mass, chi, ndf, muon_sigma));
        allParticles.push_back(jpsi_part);

        RefCountedKinematicTree bcTree = kpvFitter.fit(allParticles);

        if(bcTree->isEmpty() || !bcTree->isValid() || !bcTree->isConsistent()) return;

        RefCountedKinematicParticle bc_part = bcTree->currentParticle();
        if(!bc_part->currentState().isValid()) return;

        RefCountedKinematicVertex bc_vertex = bcTree->currentDecayVertex();
        if(!bc_vertex->vertexIsValid()) return; 

        particle_cand JPcand;
        particle_cand Bcand;

        if(myVertex.isValid()){
            //    std::cout << "refittex vertex available" << std::endl;
            JPcand = calculateIPvariables(extrapolator, jpsi_part, jpsi_vertex, myVertex);
            Bcand = calculateIPvariables(extrapolator, bc_part, bc_vertex, myVertex);
        }else{
            //    std::cout << "no refittex vertex available" << std::endl;
            JPcand = calculateIPvariables(extrapolator, jpsi_part, jpsi_vertex, closestVertex);
            Bcand = calculateIPvariables(extrapolator, bc_part, bc_vertex, closestVertex);
        }



	std::vector< RefCountedKinematicParticle > bc_children = bcTree->finalStateParticles();
	math::PtEtaPhiMLorentzVector mu3_fit = daughter_p4(bc_children, 0);


	//	std::cout << "before fit3: " << mu3.pt() << " " << mu3.eta() << " " << mu3.phi() << " " << mu3.mass() << std::endl;
	//	std::cout << "after fit3: " << mu3_fit.pt() << " " << mu3_fit.eta() << " " << mu3_fit.phi() << " " << mu3_fit.mass() << std::endl;



        // std::cout << "J/psi candidate (lip, lips, pvip, pvips, fl3d, fls3d, alpha) = " << JPcand.lip << " " << JPcand.lips << " " << JPcand.pvip << " " << JPcand.pvips << " " << JPcand.fl3d << " " << JPcand.fls3d << " " << JPcand.alpha << std::endl;

        //std::cout << "B candidate (lip, lips, pvip, pvips, fl3d, fls3d, alpha) = " << Bcand.lip << " " << Bcand.lips << " " << Bcand.pvip << " " << Bcand.pvips << " " << Bcand.fl3d << " " << Bcand.fls3d << " " << Bcand.alpha << std::endl;



        /********************************************************************
         *
         * Step8: calculation of the isolation quantities
         *
         ********************************************************************/

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


        Float_t iso = 0;
        Float_t iso_mu1 = 0;
        Float_t iso_mu2 = 0;
        Float_t iso_mu3 = 0;
        Int_t ntracks = 0;
        Float_t iso_mindoca = 999; 
  
  
        for(size_t itrk = 0; itrk < mypfcollection.size();  itrk++){
    
            Float_t _dR = reco::deltaR(mypfcollection[itrk].eta(), mypfcollection[itrk].phi(),
                                       bc_part->currentState().globalMomentum().eta(), bc_part->currentState().globalMomentum().phi());

            if(_dR < 0.7) iso += mypfcollection[itrk].pt();

            Float_t _dR1 = reco::deltaR(mypfcollection[itrk].eta(), mypfcollection[itrk].phi(),
                                        mu1_fit.eta(), mu1_fit.phi());
					//                                        muoncollection[mcidx_mu1].eta(), muoncollection[mcidx_mu1].phi());

            if(_dR1 < 0.7) iso_mu1 += mypfcollection[itrk].pt();

            Float_t _dR2 = reco::deltaR(mypfcollection[itrk].eta(), mypfcollection[itrk].phi(),
                                        mu2_fit.eta(), mu2_fit.phi());
	    //                                        muoncollection[mcidx_mu2].eta(), muoncollection[mcidx_mu2].phi());

            if(_dR2 < 0.7) iso_mu2 += mypfcollection[itrk].pt();

            Float_t _dR3 = reco::deltaR(mypfcollection[itrk].eta(), mypfcollection[itrk].phi(),
					mu3_fit.eta(), mu3_fit.phi());
					//                                        mu3.eta(), mu3.phi());

            if(_dR3 < 0.7) iso_mu3 += mypfcollection[itrk].pt();


            //    TrajectoryStateOnSurface tsos = extrapolator.extrapolate(mytracks[itrk].impactPointState());
            TrajectoryStateOnSurface tsos_pf = extrapolator.extrapolate(mytracks[itrk].impactPointState(), bc_vertex->position());

    
            //    TrajectoryStateOnSurface tsos = extrapolator.extrapolate(particle->currentState().freeTrajectoryState(),
            //							     RecoVertex::convertPos(wrtVertex.position()));
    

            VertexDistance3D a3d_pf;  
            // use own function here ... 
            std::pair<bool,Measurement1D> cur3DIP_pf = JpsiMuNtuplizer::absoluteImpactParameter(tsos_pf, bc_vertex, a3d_pf);
    
            Float_t pvip_pf = cur3DIP_pf.second.value();

            //    std::cout << itrk << ": Distance of closest apporach to the bc vertex = " << pvip << std::endl;
    
            if(pvip_pf < 0.03) ntracks+=1;

            if(iso_mindoca > pvip_pf) iso_mindoca = pvip_pf;
        }

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
         * Step9: Filling normal branches
         *
         ********************************************************************/

        nBranches_->Jpsi_mu1_pt.push_back(mu1_fit.pt());
        nBranches_->Jpsi_mu1_eta.push_back(mu1_fit.eta());
        nBranches_->Jpsi_mu1_phi.push_back(mu1_fit.phi());
        nBranches_->Jpsi_mu1_mass.push_back(mu1_fit.mass());
        nBranches_->Jpsi_mu1_unfit_pt.push_back(muoncollection[mcidx_mu1].pt());
        nBranches_->Jpsi_mu1_unfit_eta.push_back(muoncollection[mcidx_mu1].eta());
        nBranches_->Jpsi_mu1_unfit_phi.push_back(muoncollection[mcidx_mu1].phi());
        nBranches_->Jpsi_mu1_unfit_mass.push_back(muoncollection[mcidx_mu1].mass());
        nBranches_->Jpsi_mu1_q.push_back(muoncollection[mcidx_mu1].charge());
        nBranches_->Jpsi_mu1_isLoose.push_back(muoncollection[mcidx_mu1].isLooseMuon());
        nBranches_->Jpsi_mu1_isTight.push_back(muoncollection[mcidx_mu1].isTightMuon(closestVertex));
        nBranches_->Jpsi_mu1_isPF.push_back(muoncollection[mcidx_mu1].isPFMuon());
        nBranches_->Jpsi_mu1_isGlobal.push_back(muoncollection[mcidx_mu1].isGlobalMuon());
        nBranches_->Jpsi_mu1_isTracker.push_back(muoncollection[mcidx_mu1].isTrackerMuon());
        nBranches_->Jpsi_mu1_isSoft.push_back(muoncollection[mcidx_mu1].isSoftMuon(closestVertex));
        nBranches_->Jpsi_mu1_vx.push_back(muoncollection[mcidx_mu1].vx());
        nBranches_->Jpsi_mu1_vy.push_back(muoncollection[mcidx_mu1].vy());
        nBranches_->Jpsi_mu1_vz.push_back(muoncollection[mcidx_mu1].vz());
        nBranches_->Jpsi_mu1_iso.push_back(iso_mu1);
        nBranches_->Jpsi_mu1_dbiso.push_back(MuonPFIso(muoncollection[mcidx_mu1]));
  
        nBranches_->Jpsi_mu2_pt.push_back(mu2_fit.pt());
        nBranches_->Jpsi_mu2_eta.push_back(mu2_fit.eta());
        nBranches_->Jpsi_mu2_phi.push_back(mu2_fit.phi());
        nBranches_->Jpsi_mu2_mass.push_back(mu2_fit.mass());
        nBranches_->Jpsi_mu2_unfit_pt.push_back(muoncollection[mcidx_mu2].pt());
        nBranches_->Jpsi_mu2_unfit_eta.push_back(muoncollection[mcidx_mu2].eta());
        nBranches_->Jpsi_mu2_unfit_phi.push_back(muoncollection[mcidx_mu2].phi());
        nBranches_->Jpsi_mu2_unfit_mass.push_back(muoncollection[mcidx_mu2].mass());
        nBranches_->Jpsi_mu2_q.push_back(muoncollection[mcidx_mu2].charge());
        nBranches_->Jpsi_mu2_isLoose.push_back(muoncollection[mcidx_mu2].isLooseMuon());
        nBranches_->Jpsi_mu2_isTight.push_back(muoncollection[mcidx_mu2].isTightMuon(closestVertex));
        nBranches_->Jpsi_mu2_isPF.push_back(muoncollection[mcidx_mu2].isPFMuon());
        nBranches_->Jpsi_mu2_isGlobal.push_back(muoncollection[mcidx_mu2].isGlobalMuon());
        nBranches_->Jpsi_mu2_isTracker.push_back(muoncollection[mcidx_mu2].isTrackerMuon());
        nBranches_->Jpsi_mu2_isSoft.push_back(muoncollection[mcidx_mu2].isSoftMuon(closestVertex));
        nBranches_->Jpsi_mu2_vx.push_back(muoncollection[mcidx_mu2].vx());
        nBranches_->Jpsi_mu2_vy.push_back(muoncollection[mcidx_mu2].vy());
        nBranches_->Jpsi_mu2_vz.push_back(muoncollection[mcidx_mu2].vz());
        nBranches_->Jpsi_mu2_iso.push_back(iso_mu2);
        nBranches_->Jpsi_mu2_dbiso.push_back(MuonPFIso(muoncollection[mcidx_mu2]));

        nBranches_->Jpsi_mu3_unfit_pt.push_back(mu3.pt());
        nBranches_->Jpsi_mu3_unfit_eta.push_back(mu3.eta());
        nBranches_->Jpsi_mu3_unfit_phi.push_back(mu3.phi());
        nBranches_->Jpsi_mu3_unfit_mass.push_back(mu3.mass());
        nBranches_->Jpsi_mu3_pt.push_back(mu3_fit.pt());
        nBranches_->Jpsi_mu3_eta.push_back(mu3_fit.eta());
        nBranches_->Jpsi_mu3_phi.push_back(mu3_fit.phi());
        nBranches_->Jpsi_mu3_mass.push_back(mu3_fit.mass());
        nBranches_->Jpsi_mu3_q.push_back(mu3.charge());
        nBranches_->Jpsi_mu3_isLoose.push_back(mu3.isLooseMuon());
        nBranches_->Jpsi_mu3_isTight.push_back(mu3.isTightMuon(closestVertex));
        nBranches_->Jpsi_mu3_isPF.push_back(mu3.isPFMuon());
        nBranches_->Jpsi_mu3_isGlobal.push_back(mu3.isGlobalMuon());
        nBranches_->Jpsi_mu3_isTracker.push_back(mu3.isTrackerMuon());
        nBranches_->Jpsi_mu3_isSoft.push_back(mu3.isSoftMuon(closestVertex));
        nBranches_->Jpsi_mu3_vx.push_back(mu3.vx());
        nBranches_->Jpsi_mu3_vy.push_back(mu3.vy());
        nBranches_->Jpsi_mu3_vz.push_back(mu3.vz());
        nBranches_->Jpsi_mu3_iso.push_back(iso_mu3);
        nBranches_->Jpsi_mu3_dbiso.push_back(MuonPFIso(mu3));

        nBranches_->Jpsi_PV_vx.push_back(vertices_->begin()->position().x());
        nBranches_->Jpsi_PV_vy.push_back(vertices_->begin()->position().y());
        nBranches_->Jpsi_PV_vz.push_back(vertices_->begin()->position().z());

        if(myVertex.isValid()){
            nBranches_->Jpsi_bbPV_refit_vx.push_back(myVertex.position().x());
            nBranches_->Jpsi_bbPV_refit_vy.push_back(myVertex.position().y());
            nBranches_->Jpsi_bbPV_refit_vz.push_back(myVertex.position().z());
        }else{
            nBranches_->Jpsi_bbPV_refit_vx.push_back(-1);
            nBranches_->Jpsi_bbPV_refit_vy.push_back(-1);
            nBranches_->Jpsi_bbPV_refit_vz.push_back(-1);
        }

        nBranches_->Jpsi_bbPV_vx.push_back(closestVertex.position().x());
        nBranches_->Jpsi_bbPV_vy.push_back(closestVertex.position().y());
        nBranches_->Jpsi_bbPV_vz.push_back(closestVertex.position().z());

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
        nBranches_->Jpsi_unfit_pt.push_back(jpsi_tlv_highest.Pt());
        nBranches_->Jpsi_unfit_mass.push_back(jpsi_tlv_highest.M());
        nBranches_->Jpsi_unfit_vprob.push_back(jpsi_vprob_highest);

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

        nBranches_->Jpsi_trimu_iso.push_back(iso);
        nBranches_->Jpsi_trimu_iso_ntracks.push_back(ntracks);
        nBranches_->Jpsi_trimu_iso_mindoca.push_back(iso_mindoca);


        nBranches_->Jpsi_trimu_unfit_pt.push_back(tlv_B.Pt());
        nBranches_->Jpsi_trimu_unfit_mass.push_back(tlv_B.M());
        nBranches_->Jpsi_trimu_unfit_vprob.push_back(vprob_bc);
  
        if(vprob_bc!=-9){
            nBranches_->Jpsi_trimu_unfit_vx.push_back(vertex_bc.position().x());
            nBranches_->Jpsi_trimu_unfit_vy.push_back(vertex_bc.position().y());
            nBranches_->Jpsi_trimu_unfit_vz.push_back(vertex_bc.position().z());
        }




        /********************************************************************
         *
         * Step10: check gen-matching and fill them
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

  
        bool flag_nr_match = false;
        if(gen_nr_mu.size()==1){
            Float_t _dR = reco::deltaR(gen_nr_mu[0]->eta(), gen_nr_mu[0]->phi(), 
                                       mu3_fit.eta(), mu3_fit.phi());

            if(_dR < 0.1) flag_nr_match = true;
        }

        bool flag_jpsi_match = false;
        if(gen_jpsi_mu.size()==2){

            Float_t _dR_11 = reco::deltaR(gen_jpsi_mu[0]->eta(), gen_jpsi_mu[0]->phi(), 
					  mu1_fit.eta(), mu1_fit.phi());
					  //                                          muoncollection[mcidx_mu1].eta(), muoncollection[mcidx_mu1].phi());
    
            Float_t _dR_22 = reco::deltaR(gen_jpsi_mu[1]->eta(), gen_jpsi_mu[1]->phi(), 
					  mu2_fit.eta(), mu2_fit.phi());
					  //                                          muoncollection[mcidx_mu2].eta(), muoncollection[mcidx_mu2].phi());

            Float_t _dR_21 = reco::deltaR(gen_jpsi_mu[1]->eta(), gen_jpsi_mu[1]->phi(), 
					  mu1_fit.eta(), mu1_fit.phi());
					  //                                          muoncollection[mcidx_mu1].eta(), muoncollection[mcidx_mu1].phi());
    
            Float_t _dR_12 = reco::deltaR(gen_jpsi_mu[0]->eta(), gen_jpsi_mu[0]->phi(), 
					  mu2_fit.eta(), mu2_fit.phi());
					  //                                          muoncollection[mcidx_mu2].eta(), muoncollection[mcidx_mu2].phi());
      

            if(_dR_11 < 0.1 && _dR_22 < 0.1) flag_jpsi_match = true;
            if(_dR_21 < 0.1 && _dR_12 < 0.1) flag_jpsi_match = true;
        }

  
        // -9 if there is no Bc found 
        nBranches_->Jpsi_genPV_vx.push_back(genvertex.x());
        nBranches_->Jpsi_genPV_vy.push_back(genvertex.y());
        nBranches_->Jpsi_genPV_vz.push_back(genvertex.z());
        nBranches_->Jpsi_ngenmuons.push_back(gen_nr_mu.size() + gen_jpsi_mu.size());
        nBranches_->Jpsi_isgenmatched.push_back((int)flag_jpsi_match);
        nBranches_->Jpsi_mu3_isgenmatched.push_back((int)flag_nr_match);

        nBranches_->IsJpsiMu.push_back(1.);
    }
    // B Candidate Kinematic fit passed
    if( doCutFlow_) {
      nBranches_->cutflow_perevt->Fill(9);
      }
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
