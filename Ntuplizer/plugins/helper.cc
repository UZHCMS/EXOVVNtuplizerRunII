#include "../interface/helper.h"

TVector3 helper::getVertex(const reco::GenParticle& part){
  return TVector3(part.vx(),part.vy(),part.vz());
}

Int_t helper::decaymode_id(std::string str){
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


float helper::MuonPFIso(pat::Muon muon){
  
  float sumChargedHadronPt = muon.pfIsolationR04().sumChargedHadronPt;
  float sumNeutralHadronEt = muon.pfIsolationR04().sumNeutralHadronEt;
  float sumPhotonEt = muon.pfIsolationR04().sumPhotonEt;
  float sumPUPt = muon.pfIsolationR04().sumPUPt;
  float iso = (sumChargedHadronPt + std::max( 0. ,sumNeutralHadronEt + sumPhotonEt - 0.5 * sumPUPt));// / muon.pt()
  
  return iso;
}



Float_t helper::getMaxDoca(std::vector<RefCountedKinematicParticle> &kinParticles){

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



Float_t helper::getMinDoca(std::vector<RefCountedKinematicParticle> &kinParticles) {

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


void helper::printout(const RefCountedKinematicVertex& myVertex){
    std::cout << "Vertex:" << std::endl;
    if (myVertex->vertexIsValid()) {
        std::cout << "\t Decay vertex: " << myVertex->position() << myVertex->chiSquared() << " " << myVertex->degreesOfFreedom()
                  << std::endl;
    } else
        std::cout << "\t Decay vertex Not valid\n";
}

void helper::printout(const RefCountedKinematicParticle& myParticle){
    std::cout << "Particle:" << std::endl;
    //accessing the reconstructed Bs meson parameters:
    //SK: uncomment if needed  AlgebraicVector7 bs_par = myParticle->currentState().kinematicParameters().vector();

    //and their joint covariance matrix:
    //SK:uncomment if needed  AlgebraicSymMatrix77 bs_er = myParticle->currentState().kinematicParametersError().matrix();
    std::cout << "\t Momentum at vertex: " << myParticle->currentState().globalMomentum() << std::endl;
    std::cout << "\t Parameters at vertex: " << myParticle->currentState().kinematicParameters().vector() << std::endl;
}

void helper::printout(const RefCountedKinematicTree& myTree){
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


particle_cand helper::calculateIPvariables(
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


std::pair<bool, Measurement1D> helper::absoluteImpactParameter(const TrajectoryStateOnSurface& tsos,
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






math::PtEtaPhiMLorentzVector helper::daughter_p4(std::vector< RefCountedKinematicParticle > fitted_children, size_t i){
  const auto& state = fitted_children.at(i)->currentState();

  return math::PtEtaPhiMLorentzVector(
				      state.globalMomentum().perp(), 
				      state.globalMomentum().eta() ,
				      state.globalMomentum().phi() ,
				      state.mass()
				      );
}


std::tuple<Float_t, TransientVertex> helper::vertexProb( const std::vector<reco::TransientTrack>& tracks){

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