#include "../interface/JpsiMuNtuplizer.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include <CondFormats/JetMETObjects/interface/JetResolutionObject.h>
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include <cmath>
#include <vector>

//===================================================================================================================
JpsiMuNtuplizer::JpsiMuNtuplizer( edm::EDGetTokenT<pat::MuonCollection>    muonToken   ,
				  edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
				  edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken,
				  edm::EDGetTokenT<edm::TriggerResults> triggertoken,
				  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobject,
				  NtupleBranches* nBranches )
  : CandidateNtuplizer ( nBranches )
  , muonToken_	        ( muonToken )
  , verticeToken_          ( verticeToken )
  , packedpfcandidatesToken_(packedpfcandidatesToken) 
  , HLTtriggersToken_	( triggertoken )
  , triggerObjects_	( triggerobject )

   
{
}

//===================================================================================================================
JpsiMuNtuplizer::~JpsiMuNtuplizer( void )
{

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

  const double pCut(pmin), coneSize(r);
  const bool verbose(false);

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









void JpsiMuNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  
  /*
   *
   * Step1: check if the J/psi trigger is fired.
   * Namely, HLT_DoubleMu4_JpsiTrk_Displaced_v
   *
   */

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

  std::cout << "finalTriggerName = "  << finalTriggerName << std::endl;



  /*
   *
   * Step2: pre-select muons for building J/psi candidates ... 
   *
   */

  event.getByToken(verticeToken_   , vertices_     );
  event.getByToken(muonToken_	, muons_    );
  event.getByToken(triggerObjects_  , triggerObjects);

  reco::VertexCollection::const_iterator firstGoodVertex = vertices_->begin();


  std::vector<pat::Muon> muoncollection;
  std::vector<int> muoncollection_id;
  muoncollection.clear();
  muoncollection_id.clear();

  for(size_t imuon = 0; imuon < muons_->size(); ++ imuon){

    const pat::Muon & muon = (*muons_)[imuon];

    if(muon.pt() < 4) continue;
    if(TMath::Abs(muon.eta()) > 2.4) continue;
    if(!(muon.track().isNonnull())) continue;

    bool isSoft = muon.isSoftMuon(*firstGoodVertex);
    bool isGlobal = muon.isGlobalMuon();
    //    bool isTracker = muon.isTrackerMuon();
    //    bool isLoose = muon.isLooseMuon();
    //    bool isTight =  muon.isTightMuon(*firstGoodVertex);
    //    bool isPF = muon.isPFMuon();

    if(!(isSoft && isGlobal)) continue;

    // remove muons, that is far away from PV ?
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


  /*
   *
   * Step3: building J/psi candidates 
   *
   */
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

  Float_t Jpsi_max_pt = -1;
  unsigned int idx_mu1 = -1;
  unsigned int idx_mu2 = -1;
  unsigned int mcidx_mu1 = -1;
  unsigned int mcidx_mu2 = -1;
  TLorentzVector Jpsi_tlv_highest;
  Float_t Jpsi_vprob_highest;
  TransientVertex Jpsi_vertex_highest;

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
      const reco::TrackRef track1_muon_ = muoncollection[imu].muonBestTrack();
      const reco::TrackRef track2_muon_ = muoncollection[jmu].muonBestTrack();
      
      transient_tracks_dimuon.push_back((*builder).build(track1_muon_));
      transient_tracks_dimuon.push_back((*builder).build(track2_muon_));

      Float_t vprob_jpsi = -9;
      TransientVertex vertex_jpsi;
      std::tie(vprob_jpsi, vertex_jpsi) = vertexProb(transient_tracks_dimuon);

      if(!(vprob_jpsi > 0)) continue;

      if(Jpsi_max_pt < jpsi_pt){
	Jpsi_max_pt = jpsi_pt;
	idx_mu1 = muoncollection_id[imu];
	idx_mu2 = muoncollection_id[jmu];
	mcidx_mu1 = imu;
	mcidx_mu2 = jmu;
	Jpsi_tlv_highest = tlv_jpsi;
	Jpsi_vprob_highest = vprob_jpsi;
	Jpsi_vertex_highest = vertex_jpsi;
      }
    }
  }

  if(Jpsi_max_pt < 8) return;

  std::cout << "J/psi candidate found with (vprob, pt) = (" << Jpsi_vprob_highest << ", " << Jpsi_max_pt << ")" << std::endl;



  /*
   *
   * Step4: 3rd muon selection
   *
   */

  //  event.getByToken( packedpfcandidatesToken_               , packedpfcandidates_      ); 

  pat::Muon mu3;
  Float_t max_pt3 = -1;

  for(size_t imuon = 0; imuon < muons_->size(); ++ imuon){

    const pat::Muon & muon = (*muons_)[imuon];

    if(muon.pt() < 4) continue;
    if(TMath::Abs(muon.eta()) > 2.4) continue;
    if(!(muon.track().isNonnull())) continue;
    if(imuon==idx_mu1 || imuon==idx_mu2) continue;

    bool isSoft = muon.isSoftMuon(*firstGoodVertex);
    bool isTight =  muon.isTightMuon(*firstGoodVertex);

    if(!(isSoft && isTight)) continue;

    /// whatever selection you want to have for the 3rd muon ... 
    
    if(muon.pt() > max_pt3){
      max_pt3 = muon.pt();
      mu3 = muon;
      std::cout << "This is the one: " << imuon  << std::endl;
    }
  }

  if(max_pt3==-1) return;

  std::cout << "3rd muon found with pt= " << max_pt3 << std::endl;

  /*
   *
   * Step5: filling up branches ...
   *
   */


  /*
   * J/psi kinematics
   */
  GlobalVector Jpsi_direction_highest(Jpsi_tlv_highest.Px(), Jpsi_tlv_highest.Py(), Jpsi_tlv_highest.Pz()); //To compute sign of IP

  double Jpsi_flightSig3D = reco::SecondaryVertex::computeDist3d(*firstGoodVertex,Jpsi_vertex_highest,Jpsi_direction_highest,true).significance();
  double Jpsi_flightLength3D = reco::SecondaryVertex::computeDist3d(*firstGoodVertex,Jpsi_vertex_highest,Jpsi_direction_highest,true).value();
  double Jpsi_flightLengthErr3D = reco::SecondaryVertex::computeDist3d(*firstGoodVertex,Jpsi_vertex_highest,Jpsi_direction_highest,true).error();

  double Jpsi_flightSig2D = reco::SecondaryVertex::computeDist2d(*firstGoodVertex,Jpsi_vertex_highest,Jpsi_direction_highest,true).significance();
  double Jpsi_flightLength2D = reco::SecondaryVertex::computeDist2d(*firstGoodVertex,Jpsi_vertex_highest,Jpsi_direction_highest,true).value();
  double Jpsi_flightLengthErr2D = reco::SecondaryVertex::computeDist2d(*firstGoodVertex,Jpsi_vertex_highest,Jpsi_direction_highest,true).error();


  TLorentzVector tlv_muon3;
  tlv_muon3.SetPtEtaPhiM(mu3.pt(), mu3.eta(), mu3.phi(), mu3.mass());

  TLorentzVector tlv_B = Jpsi_tlv_highest + tlv_muon3;

  std::vector<reco::TransientTrack> transient_tracks_trimuon;

  const reco::TrackRef track1_muon = muoncollection[mcidx_mu1].muonBestTrack();
  const reco::TrackRef track2_muon = muoncollection[mcidx_mu2].muonBestTrack();
  const reco::TrackRef track3_muon = mu3.muonBestTrack();
  transient_tracks_trimuon.push_back((*builder).build(track1_muon));
  transient_tracks_trimuon.push_back((*builder).build(track2_muon));
  transient_tracks_trimuon.push_back((*builder).build(track3_muon));
  
  Float_t vprob_B = -9;
  TransientVertex B_vertex;
  std::tie(vprob_B, B_vertex) = vertexProb(transient_tracks_trimuon);

  
  if(vprob_B==-9) return;

  GlobalVector direction(tlv_B.Px(), tlv_B.Py(), tlv_B.Pz()); //To compute sign of IP

  double B_flightSig3D = -999;
  double B_flightLength3D =  -999; 
  double B_flightLengthErr3D = -999; 
	
  double B_flightSig2D =  -999; 
  double B_flightLength2D = -999; 
  double B_flightLengthErr2D = -999; 

  B_flightSig3D = reco::SecondaryVertex::computeDist3d(*firstGoodVertex,B_vertex,direction,true).significance();
  B_flightLength3D = reco::SecondaryVertex::computeDist3d(*firstGoodVertex,B_vertex,direction,true).value();
  B_flightLengthErr3D = reco::SecondaryVertex::computeDist3d(*firstGoodVertex,B_vertex,direction,true).error();
	    
  B_flightSig2D = reco::SecondaryVertex::computeDist2d(*firstGoodVertex,B_vertex,direction,true).significance();
  B_flightLength2D = reco::SecondaryVertex::computeDist2d(*firstGoodVertex,B_vertex,direction,true).value();
  B_flightLengthErr2D = reco::SecondaryVertex::computeDist2d(*firstGoodVertex,B_vertex,direction,true).error();

  // distance to closest apporach
  double doca = -1;
  TrajectoryStateClosestToPoint mu1TS = (*builder).build(track1_muon).impactPointTSCP();
  TrajectoryStateClosestToPoint mu2TS = (*builder).build(track2_muon).impactPointTSCP();
  if (mu1TS.isValid() && mu2TS.isValid()) {
    ClosestApproachInRPhi cApp;
    cApp.calculate(mu1TS.theState(), mu2TS.theState());
    if (cApp.status()){
      std::cout << "distance of closest apporach: " << cApp.distance() << std::endl;
      doca = cApp.distance();
    }
  }
  


  // impact parameter calculation
  double d03D = IPTools::absoluteImpactParameter3D((*builder).build(track1_muon), *firstGoodVertex).second.value();


  // opening angle between (PV, SV) vector and the B's flight direction
  math::XYZVector vtx_vec(B_vertex.position().x() - firstGoodVertex->position().x(),
			  B_vertex.position().y() - firstGoodVertex->position().y(),
			  B_vertex.position().z() - firstGoodVertex->position().z());
  

  math::XYZVector pt_vec(tlv_B.Px(),
			 tlv_B.Py(),
			 tlv_B.Pz());
  
  double den = (vtx_vec.R() * pt_vec.R());
  double alpha = -1; 

  if(den != 0.){
    alpha = vtx_vec.Dot(pt_vec)/den; 
  }

  //  std::cout << "opening angle:" << alpha << std::endl;




  // Isolation calculation
  

//  GlobalPoint point = fitter.fitted_vtx();
//  auto bs_pos = bs.position(point.z());
//  math::XYZVector delta(point.x() - bs_pos.x(), point.y() - bs_pos.y(), 0.);
//  math::XYZVector pt(p4.px(), p4.py(), 0.);
//  double den = (delta.R() * pt.R());
//  return (den != 0.) ? delta.Dot(pt)/den : -2;

//
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


  nBranches_->Jpsi_flightSig3D.push_back(Jpsi_flightSig3D);
  nBranches_->Jpsi_flightLength3D.push_back(Jpsi_flightLength3D);
  nBranches_->Jpsi_flightLengthErr3D.push_back(Jpsi_flightLengthErr3D);
  nBranches_->Jpsi_flightSig2D.push_back(Jpsi_flightSig2D);
  nBranches_->Jpsi_flightLength2D.push_back(Jpsi_flightLength2D);
  nBranches_->Jpsi_flightLengthErr2D.push_back(Jpsi_flightLengthErr2D);

  nBranches_->Jpsi_trimu_flightSig3D.push_back(B_flightSig3D);
  nBranches_->Jpsi_trimu_flightLength3D.push_back(B_flightLength3D);
  nBranches_->Jpsi_trimu_flightLengthErr3D.push_back(B_flightLengthErr3D);
  nBranches_->Jpsi_trimu_flightSig2D.push_back(B_flightSig2D);
  nBranches_->Jpsi_trimu_flightLength2D.push_back(B_flightLength2D);
  nBranches_->Jpsi_trimu_flightLengthErr2D.push_back(B_flightLengthErr2D);

  //  isJpsiMu_=1;
  //save muon branches
  nBranches_->Jpsi_mu1_isLoose.push_back(muoncollection[mcidx_mu1].isLooseMuon());
  nBranches_->Jpsi_mu1_isTight.push_back(muoncollection[mcidx_mu1].isTightMuon(*firstGoodVertex));
  nBranches_->Jpsi_mu1_isPF.push_back(muoncollection[mcidx_mu1].isPFMuon());
  nBranches_->Jpsi_mu1_isGlobal.push_back(muoncollection[mcidx_mu1].isGlobalMuon());
  nBranches_->Jpsi_mu1_isTracker.push_back(muoncollection[mcidx_mu1].isTrackerMuon());
  nBranches_->Jpsi_mu1_isSoft.push_back(muoncollection[mcidx_mu1].isSoftMuon(*firstGoodVertex));
  nBranches_->Jpsi_mu1_pt.push_back(muoncollection[mcidx_mu1].pt());
  nBranches_->Jpsi_mu1_eta.push_back(muoncollection[mcidx_mu1].eta());
  nBranches_->Jpsi_mu1_phi.push_back(muoncollection[mcidx_mu1].phi());
  nBranches_->Jpsi_mu1_ch.push_back(muoncollection[mcidx_mu1].charge());
  
  nBranches_->Jpsi_mu2_isLoose.push_back(muoncollection[mcidx_mu2].isLooseMuon());
  nBranches_->Jpsi_mu2_isTight.push_back(muoncollection[mcidx_mu2].isTightMuon(*firstGoodVertex));
  nBranches_->Jpsi_mu2_isPF.push_back(muoncollection[mcidx_mu2].isPFMuon());
  nBranches_->Jpsi_mu2_isGlobal.push_back(muoncollection[mcidx_mu2].isGlobalMuon());
  nBranches_->Jpsi_mu2_isTracker.push_back(muoncollection[mcidx_mu2].isTrackerMuon());
  nBranches_->Jpsi_mu2_isSoft.push_back(muoncollection[mcidx_mu2].isSoftMuon(*firstGoodVertex));
  nBranches_->Jpsi_mu2_pt.push_back(muoncollection[mcidx_mu2].pt());
  nBranches_->Jpsi_mu2_eta.push_back(muoncollection[mcidx_mu2].eta());
  nBranches_->Jpsi_mu2_phi.push_back(muoncollection[mcidx_mu2].phi());
  nBranches_->Jpsi_mu2_ch.push_back(muoncollection[mcidx_mu2].charge());

  nBranches_->Jpsi_mu3_isLoose.push_back(mu3.isLooseMuon());
  nBranches_->Jpsi_mu3_isTight.push_back(mu3.isTightMuon(*firstGoodVertex));
  nBranches_->Jpsi_mu3_isPF.push_back(mu3.isPFMuon());
  nBranches_->Jpsi_mu3_isGlobal.push_back(mu3.isGlobalMuon());
  nBranches_->Jpsi_mu3_isTracker.push_back(mu3.isTrackerMuon());
  nBranches_->Jpsi_mu3_isSoft.push_back(mu3.isSoftMuon(*firstGoodVertex));
  nBranches_->Jpsi_mu3_pt.push_back(mu3.pt());
  nBranches_->Jpsi_mu3_eta.push_back(mu3.eta());
  nBranches_->Jpsi_mu3_phi.push_back(mu3.phi());
  nBranches_->Jpsi_mu3_ch.push_back(mu3.charge());
  nBranches_->Jpsi_mu3_x.push_back(mu3.vx());
  nBranches_->Jpsi_mu3_y.push_back(mu3.vy());
  nBranches_->Jpsi_mu3_z.push_back(mu3.vz());

  nBranches_->Jpsi_PV_x.push_back(firstGoodVertex->position().x());
  nBranches_->Jpsi_PV_y.push_back(firstGoodVertex->position().y());
  nBranches_->Jpsi_PV_z.push_back(firstGoodVertex->position().z());

  nBranches_->Jpsi_dx.push_back(Jpsi_vertex_highest.position().x()-firstGoodVertex->position().x());
  nBranches_->Jpsi_dy.push_back(Jpsi_vertex_highest.position().y()-firstGoodVertex->position().y());
  nBranches_->Jpsi_dz.push_back(Jpsi_vertex_highest.position().z()-firstGoodVertex->position().z());	  

  nBranches_->Jpsi_pt.push_back(Jpsi_tlv_highest.Pt());
  nBranches_->Jpsi_eta.push_back(Jpsi_tlv_highest.Eta());
  nBranches_->Jpsi_phi.push_back(Jpsi_tlv_highest.Phi());
  nBranches_->Jpsi_mass.push_back(Jpsi_tlv_highest.M());
  nBranches_->Jpsi_vtxprob.push_back(Jpsi_vprob_highest);
  nBranches_->Jpsi_vtxz.push_back(Jpsi_vertex_highest.position().z());

  nBranches_->Jpsi_trimu_pt.push_back(tlv_B.Pt());
  nBranches_->Jpsi_trimu_eta.push_back(tlv_B.Eta());
  nBranches_->Jpsi_trimu_phi.push_back(tlv_B.Phi());
  nBranches_->Jpsi_trimu_mass.push_back(tlv_B.M());

  nBranches_->Jpsi_trimu_vtxprob.push_back(vprob_B);
  nBranches_->Jpsi_trimu_vtxz.push_back(B_vertex.position().z());

  nBranches_->Jpsi_trimu_dx.push_back(B_vertex.position().x()-firstGoodVertex->position().x());
  nBranches_->Jpsi_trimu_dy.push_back(B_vertex.position().y()-firstGoodVertex->position().y());
  nBranches_->Jpsi_trimu_dz.push_back(B_vertex.position().z()-firstGoodVertex->position().z());
  nBranches_->IsJpsiMu.push_back(1.);

  nBranches_->Jpsi_doca.push_back(doca);
  nBranches_->Jpsi_d03D.push_back(d03D);
  nBranches_->Jpsi_alpha.push_back(alpha);



  //if(isJpsiMu_==0){
  //std::cout<<"this should be working and i am returning"<<std::endl;
  //    return;
  //  }

}


