#include "../interface/GenEventNtuplizer.h"

//===================================================================================================================
GenEventNtuplizer::GenEventNtuplizer( std::vector< edm::EDGetTokenT< GenEventInfoProduct > > tokens, NtupleBranches* nBranches,  std::vector< edm::EDGetTokenT< LHEEventProduct > > tokens_lhe )
   : CandidateNtuplizer( nBranches )
   , geneventToken_( tokens[0] )
   , lheEventProductToken_( tokens_lhe[0])
{

}

//===================================================================================================================
GenEventNtuplizer::~GenEventNtuplizer( void )
{

}

//===================================================================================================================
void GenEventNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  event.getByToken(geneventToken_, geneventInfo_);  
  
  nBranches_->genWeight=geneventInfo_->weight();
  nBranches_->qScale=geneventInfo_->qScale();
  nBranches_->PDF_x.push_back((geneventInfo_->pdf()->x).first);
  nBranches_->PDF_x.push_back((geneventInfo_->pdf()->x).second);
  nBranches_->PDF_xPDF.push_back((geneventInfo_->pdf()->xPDF).first);
  nBranches_->PDF_xPDF.push_back((geneventInfo_->pdf()->xPDF).second);
  nBranches_->PDF_id.push_back((geneventInfo_->pdf()->id).first);
  nBranches_->PDF_id.push_back((geneventInfo_->pdf()->id).second);
  


  //gen Parton HT
  //taken from https://github.com/IHEP-CMS/BSMFramework/blob/CMSSW_805p1/BSM3G_TNT_Maker/src/EventInfoSelector.cc#L63-L80
  //Definition taken from https://github.com/cmkuo/ggAnalysis/blob/a24edc65be23b402d761c75545192ce79cddf316/ggNtuplizer/plugins/ggNtuplizer_genParticles.cc#L201 
  //Zaixing has a somehow different, but likely equivalent implementation
  //https://github.com/zaixingmao/FSA/blob/miniAOD_dev_7_4_14/DataFormats/src/PATFinalStateEvent.cc#L153


  event.getByToken(lheEventProductToken_, lheEventProduct_);
  float lheHt_ = 0.;
  int nParton = 0;

  std::vector<TLorentzVector> tlv;

  if(lheEventProduct_.isValid()){
  
    const lhef::HEPEUP& lheEvent = lheEventProduct_->hepeup();
    std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
    size_t numParticles = lheParticles.size();

    for(size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle){
      int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
      int status = lheEvent.ISTUP[idxParticle];



      if(status == 1 && ((absPdgId >= 1 && absPdgId <= 6) || absPdgId == 21)){ // quarks and gluons
	lheHt_ += TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.) + TMath::Power(lheParticles[idxParticle][1], 2.)); // first entry is px, second py
      } 

      // Yuta added 7 Sep : For gen weighting
      if(status == 1 && ((absPdgId >= 1 && absPdgId <= 5) || absPdgId == 21)){ // quarks and gluons, excluding top
	nParton++;
      }

      if(status == 1 && (absPdgId == 11 || absPdgId == 13 || absPdgId == 15)){ // leptons (needed for DY, W stitching)
	TLorentzVector tl;
	tl.SetPxPyPzE(lheParticles[idxParticle][0],
		      lheParticles[idxParticle][1],
		      lheParticles[idxParticle][2],
		      lheParticles[idxParticle][3]);						

	tlv.push_back(tl);
      }


    }
  }


  // Yuta added 7 Sep : For gen weighting
  if(tlv.size()==2){
    Float_t genboson_mass = (tlv[0] + tlv[1]).M();
    Float_t genboson_pt = (tlv[0] + tlv[1]).Pt();

    nBranches_->lheV_mass = genboson_mass;
    nBranches_->lheV_pt = genboson_pt;
  }else{
    nBranches_->lheV_mass = -1;
    nBranches_->lheV_pt = -1;
  }

  //  nBranches_->lheV_pt = 0.;
  //  nBranches_->lheNj = 0; // Does anybody use this ?
  nBranches_->lheNj = nParton;
  nBranches_->lheHT = lheHt_;


 
}
