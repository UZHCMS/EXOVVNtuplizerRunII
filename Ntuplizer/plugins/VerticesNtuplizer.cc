#include "../interface/VerticesNtuplizer.h"

#include <cmath>

//===================================================================================================================
VerticesNtuplizer::VerticesNtuplizer( std::vector<edm::EDGetTokenT<reco::VertexCollection>> tokens,
				      edm::EDGetTokenT<reco::BeamSpot>             beamToken,
				      NtupleBranches* nBranches,
				      std::map< std::string, bool >& runFlags) 
   : CandidateNtuplizer( nBranches )
   , vtxToken_ ( tokens[0])
   , beamToken_( beamToken)
   , isJpsiMu_( runFlags["doJpsiMu"]  )
   , isJpsiEle_( runFlags["doJpsiEle"]  )
   , isJpsiTau_( runFlags["doJpsiTau"]  )
{

}

//===================================================================================================================
VerticesNtuplizer::~VerticesNtuplizer( void )
{

}

//===================================================================================================================
bool VerticesNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){
  
    //Skip events with no jspi if that analysis is chosen
   
//    bool rflag = false;
//    
//    if(isJpsiEle_ || isJpsiMu_){
//      if ( nBranches_->JpsiMu_B_pt.size() >=1) rflag = true;
//    }
//    
//    if(isJpsiTau_){
//      if ( nBranches_->JpsiTau_B_pt.size() >=1)  rflag = true;
//    }
//
//    if(rflag==false) return;

  event.getByToken(vtxToken_, vertices_);
  event.getByToken(beamToken_, beamSpot_);
  //std::cout<<beamSpot_->z0()<<std::endl;
  nBranches_->BeamSpot_x0.push_back(beamSpot_->x0());
  nBranches_->BeamSpot_y0.push_back(beamSpot_->y0());
  nBranches_->BeamSpot_z0.push_back(beamSpot_->z0());
  //reco::BeamSpot beamSpot;
  //double z0 = beamSpot.z0();
  //std::cout<<"this is z0 "<< z0 << std::endl;
  //edm::Handle<reco::BeamSpot> beamSpotHandle;
  //event.getByLabel("offlineBeamSpot", beamSpotHandle);
  nBranches_->PV_N = vertices_->size();
  
  nBranches_->PV_filter = true;
  
  reco::VertexCollection::const_iterator firstGoodVertex = vertices_->end();
  int firstGoodVertexIdx = 0;
  for( reco::VertexCollection::const_iterator vtx = vertices_->begin(); vtx != vertices_->end(); ++vtx, ++firstGoodVertexIdx){
  
    nBranches_->PV_chi2.push_back(vtx->chi2());
    nBranches_->PV_ndof.push_back(vtx->ndof());
    nBranches_->PV_rho.push_back(vtx->position().Rho());
    nBranches_->PV_z.push_back(vtx->position().Z());
    
    bool isFake = (vtx->chi2()==0 && vtx->ndof()==0);
    if( !isFake && vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0) {
      firstGoodVertex = vtx;
      break;
    }
    
  }
  
  if ( firstGoodVertex==vertices_->end() ) nBranches_->PV_filter = false;
  

  /*
  if ( beamSpotHandle.isValid() ){
    beamSpot = *beamSpotHandle;

  } else{
    edm::LogInfo("MyAnalyzer")
      << "No beam spot available from EventSetup \n";
  }

  double x0 = beamSpot.x0();
  double y0 = beamSpot.y0();
  double z0 = beamSpot.z0();
  double sigmaz = beamSpot.sigmaZ();
  double dxdz = beamSpot.dxdz();
  double BeamWidthX = beamSpot.BeamWidthX();
  double BeamWidthY = beamSpot.BeamWidthY();

  // print the beam spot object
  // cout << beamSpot << endl;
  std::cout<<"this is z0 "<< z0 << std::endl;*/

  return true;
}
