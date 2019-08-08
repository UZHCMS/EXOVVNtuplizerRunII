#include "../interface/JpsiEleNtuplizer.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include <CondFormats/JetMETObjects/interface/JetResolutionObject.h>
#include <cmath>

//===================================================================================================================
JpsiEleNtuplizer::JpsiEleNtuplizer( edm::EDGetTokenT<edm::View<pat::Electron>>    electronToken   ,
			      //			      std::vector<edm::EDGetTokenT<reco::VertexCollection> > tokens,
			      edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
			      NtupleBranches* nBranches )
  : CandidateNtuplizer ( nBranches )
  , electronToken_         ( electronToken )
  , verticeToken_          ( verticeToken )
    
   
{
}

//===================================================================================================================
JpsiEleNtuplizer::~JpsiEleNtuplizer( void )
{

}
std::tuple<Float_t, Float_t> JpsiEleNtuplizer::vertexProb( const std::vector<reco::TransientTrack>& tracks) {
  
  Float_t vprob = -1;
  Float_t vz = -1;
  //  Float_t xy = -1;

  KalmanVertexFitter kalman_fitter;
  TransientVertex vertex;
  vertex = kalman_fitter.vertex(tracks);

  if(vertex.isValid()){

    vprob =  TMath::Prob(vertex.totalChiSquared(), vertex.degreesOfFreedom());

    //    if(vprob==0) vprob = -1;
    vz = vertex.position().z();
  }  

  
  return std::forward_as_tuple(vprob, vz);
  //  return vprob;
}
//===================================================================================================================
void JpsiEleNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
  //std::cout<<"!!!--->Jpsi Ntuplizer<---!!!"<<std::endl;
  nBranches_->IsJpsiEle.push_back(1);
  event.getByToken(verticeToken_   , vertices_     );
  event.getByToken(electronToken_	, electrons_    );
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);  
  std::cout<<"im getting to the JpsiEle part"<<std::endl;
  
  Float_t vprobJpsi=-1; 
  Float_t vzJpsi=-999; 
  Float_t vprobB=-1;
  Float_t vzB=-999; 
  std::vector<reco::TransientTrack> transient_tracks_muon_pair;
  std::vector<reco::TransientTrack> transient_tracks_trimuon;
  unsigned int nmus = 0;
  unsigned int isJpsimunu_=0;
  for (const pat::Electron &ele : *electrons_){
    nmus++;
    std::cout<<nmus<<std::endl;
  }
  /*
  for (const pat::Muon &mu : *muons_) {
    nmus++;
    if(mu.pt()<2 ||fabs(mu.eta())>2.4) continue;
    //if( mu.isMuon()==0 || mu.isGlobalMuon()==0 || mu.isTrackerMuon()==0 ) continue;
    //Muon softness
    reco::VertexCollection::const_iterator firstGoodVertex = vertices_->begin();
    Int_t isSoft1=mu.isSoftMuon(*firstGoodVertex);
    if(isSoft1<0.5) continue;

    unsigned int nmus2=0;
    for (const pat::Muon &mu2 : *muons_){
      nmus2++;
      if (nmus2>=nmus) continue;
      if(mu2.pt()<2 || fabs(mu2.eta())>2.4) continue;
      //if(mu2.isMuon()==0 || mu2.isGlobalMuon()==0 || mu2.isTrackerMuon()==0) continue;
      Int_t isSoft2=mu2.isSoftMuon(*firstGoodVertex);
      if(isSoft2<0.5) continue;
      if (((mu.charge() + mu2.charge()) !=0)) continue;
      TLorentzVector tlv_muon1;
      TLorentzVector tlv_muon2;
      tlv_muon1.SetPtEtaPhiM(mu.pt(),mu.eta(),mu.phi(), mu.mass());
      tlv_muon2.SetPtEtaPhiM(mu2.pt(), mu2.eta(), mu2.phi(), mu2.mass());
      TLorentzVector tlv_jpsi = tlv_muon1 + tlv_muon2;
      double mass=tlv_jpsi.M();
      if(mass<2.6 || mass>3.6) continue;      
      
      //Vertex requirement
      transient_tracks_muon_pair.clear();
      const reco::TrackRef    track1 = mu.muonBestTrack(); 
      const reco::TrackRef    track2 = mu2.muonBestTrack();
      transient_tracks_muon_pair.push_back((*builder).build(track1.get()));
      transient_tracks_muon_pair.push_back((*builder).build(track2.get()));
      KalmanVertexFitter kalman_fitter;
      TransientVertex Jpsi_vertex;
      Jpsi_vertex = kalman_fitter.vertex(transient_tracks_muon_pair);
      std::tie(vprobJpsi, vzJpsi) = vertexProb(transient_tracks_muon_pair);

      //std::cout<<"-------->VPROB:"<<vprobJpsi<<std::endl;
      unsigned int nmus3=0;
      for (const pat::Muon &mu3 : *muons_){
	nmus3++;
	if(mu3.pt()<2) continue;
	if(nmus3==nmus || nmus3==nmus2) continue;
	//std::cout<<"1:"<<nmus<<"\t, 2:"<<nmus2<<"\t ,3:"<<nmus3<<std::endl;
	//std::cout<<"1:"<<mu.pt()<<"\t, 2:"<<mu2.pt()<<"\t ,3:"<<mu3.pt()<<std::endl;
	TLorentzVector tlv_muon3;
	TransientVertex B_vertex;
	tlv_muon3.SetPtEtaPhiM(mu3.pt(), mu3.eta(), mu3.phi(), mu3.mass());
	TLorentzVector tlv_B = tlv_jpsi + tlv_muon3;
	double b_mass =tlv_B.M();
	//std::cout<<"trimuon scenario, bmass: "<< b_mass<<std::endl;
	//vtx prob part
	transient_tracks_trimuon.clear();
	const reco::TrackRef    track3 = mu3.muonBestTrack(); 
	transient_tracks_trimuon.push_back((*builder).build(track1.get()));
	transient_tracks_trimuon.push_back((*builder).build(track2.get()));
	transient_tracks_trimuon.push_back((*builder).build(track3.get()));
	KalmanVertexFitter B_kalman_fitter;
	B_vertex = B_kalman_fitter.vertex(transient_tracks_trimuon);
	//std::cout<<"before getting to the probability thing: "<<std::endl;

	std::tie(vprobB, vzB) = vertexProb(transient_tracks_trimuon);
	//std::cout<<"after getting to the probability thing: "<<std::endl;

	isJpsimunu_=1;
	//save muon branches
	nBranches_->Jpsi_mu1_isLoose.push_back(mu.isLooseMuon());
	nBranches_->Jpsi_mu1_isTight.push_back(mu.isTightMuon(*firstGoodVertex));
	nBranches_->Jpsi_mu1_isPF.push_back(mu.isPFMuon());
	nBranches_->Jpsi_mu1_isGlobal.push_back(mu.isGlobalMuon());
	nBranches_->Jpsi_mu1_isTracker.push_back(mu.isTrackerMuon());
	nBranches_->Jpsi_mu1_isSoft.push_back(mu.isSoftMuon(*firstGoodVertex));

	nBranches_->Jpsi_mu2_isLoose.push_back(mu2.isLooseMuon());
	nBranches_->Jpsi_mu2_isTight.push_back(mu2.isTightMuon(*firstGoodVertex));
	nBranches_->Jpsi_mu2_isPF.push_back(mu.isPFMuon());
	nBranches_->Jpsi_mu2_isGlobal.push_back(mu2.isGlobalMuon());
	nBranches_->Jpsi_mu2_isTracker.push_back(mu2.isTrackerMuon());
	nBranches_->Jpsi_mu2_isSoft.push_back(mu2.isSoftMuon(*firstGoodVertex));


	nBranches_->Jpsi_mu1_pt.push_back(mu.pt());
	nBranches_->Jpsi_mu1_eta.push_back(mu.eta());
	nBranches_->Jpsi_mu1_phi.push_back(mu.phi());
	nBranches_->Jpsi_mu1_ch.push_back(mu.charge());

	nBranches_->Jpsi_mu2_pt.push_back(mu2.pt());
	nBranches_->Jpsi_mu2_eta.push_back(mu2.eta());
	nBranches_->Jpsi_mu2_phi.push_back(mu2.phi());
	nBranches_->Jpsi_mu2_ch.push_back(mu.charge());

	nBranches_->Jpsi_mu3_pt.push_back(mu3.pt());
	nBranches_->Jpsi_mu3_eta.push_back(mu3.eta());
	nBranches_->Jpsi_mu3_phi.push_back(mu3.phi());
	nBranches_->Jpsi_mu3_ch.push_back(mu3.charge());

	nBranches_->Jpsi_mu3_isLoose.push_back(mu3.isLooseMuon());
	nBranches_->Jpsi_mu3_isTight.push_back(mu3.isTightMuon(*firstGoodVertex));
	nBranches_->Jpsi_mu3_isPF.push_back(mu3.isPFMuon());
	nBranches_->Jpsi_mu3_isGlobal.push_back(mu3.isGlobalMuon());
	nBranches_->Jpsi_mu3_isTracker.push_back(mu3.isTrackerMuon());
	nBranches_->Jpsi_mu3_isSoft.push_back(mu3.isSoftMuon(*firstGoodVertex));
	//std::cout<<firstGoodVertex->position().z()<<"<-PV, Jpsi->"<< Jpsi_vertex.position().z()<< std::endl;
	//jpsi variables
	nBranches_->Jpsi_PV_x.push_back(firstGoodVertex->position().x());
	nBranches_->Jpsi_PV_y.push_back(firstGoodVertex->position().y());
	nBranches_->Jpsi_PV_z.push_back(firstGoodVertex->position().z());
	if(Jpsi_vertex.isValid()){
	  nBranches_->Jpsi_dx.push_back(Jpsi_vertex.position().x()-firstGoodVertex->position().x());
	  nBranches_->Jpsi_dy.push_back(Jpsi_vertex.position().y()-firstGoodVertex->position().y());
	  nBranches_->Jpsi_dz.push_back(Jpsi_vertex.position().z()-firstGoodVertex->position().z());	  
	}else{
	  nBranches_->Jpsi_dx.push_back(-999);
	  nBranches_->Jpsi_dy.push_back(-999);
	  nBranches_->Jpsi_dz.push_back(-999);
	}
	nBranches_->Jpsi_pt.push_back(tlv_jpsi.Pt());
	nBranches_->Jpsi_eta.push_back(tlv_jpsi.Eta());
	nBranches_->Jpsi_phi.push_back(tlv_jpsi.Phi());
	nBranches_->Jpsi_mass.push_back(mass);
	nBranches_->Jpsi_vtxprob.push_back(vprobJpsi);
	nBranches_->Jpsi_vtxz.push_back(vzJpsi);
	//std::cout<<"this is before the trimuon part"<< std::endl;

        //Trimuon (B) variables
	nBranches_->Jpsi_trimu_pt.push_back(tlv_B.Pt());
	nBranches_->Jpsi_trimu_eta.push_back(tlv_B.Eta());
	nBranches_->Jpsi_trimu_phi.push_back(tlv_B.Phi());
	nBranches_->Jpsi_trimu_mass.push_back(b_mass);

	//std::cout<<"this is before the trimuon vtx"<< std::endl;

	nBranches_->Jpsi_trimu_vtxprob.push_back(vprobB);
	//std::cout<<"this is before the trimuon vtx z direction"<< std::endl;
	//std::cout<<vprobB<<"<-yutas, normal->"<<vprobB2<< std::endl;
	//std::cout<<vzB<<"<-yutas, normal->"<<B_vertex.position().z()<< std::endl;

	//std::cout<<sizeof(B_vertex.position())<<std::endl;
	//std::cout<<B_vertex.position()<<std::endl;
	//std::cout<<B_vertex.position().z()<<std::endl;
	nBranches_->Jpsi_trimu_vtxz.push_back(vzB);
	//std::cout<<"this is before the trimuon position"<< std::endl;
	if(B_vertex.isValid()){
	  nBranches_->Jpsi_trimu_dx.push_back(B_vertex.position().x()-firstGoodVertex->position().x());
	  nBranches_->Jpsi_trimu_dy.push_back(B_vertex.position().y()-firstGoodVertex->position().y());
	  nBranches_->Jpsi_trimu_dz.push_back(B_vertex.position().z()-firstGoodVertex->position().z());
	}else{
	   nBranches_->Jpsi_trimu_dx.push_back(-999);
	   nBranches_->Jpsi_trimu_dy.push_back(-999);
	   nBranches_->Jpsi_trimu_dz.push_back(-999);
	}
	//std::cout<<"this should be the end of the Jpsi part"<< std::endl;
    
      }
    }
  }
  //
  nBranches_->IsJpsimunu.push_back(isJpsimunu_);
  // within the analyze()
  
  


  if(isJpsimunu_==0){
    //std::cout<<"this should be working and i am returning"<<std::endl;
    return;

    }*/
}
