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
  int nLeptons = 0;
  int nParton = 0;
  int nBParton = 0;
  
  double weightFacUp(0.), weightFacDown(0.), weightRenUp(0.), weightRenDown(0.), weightFacRenUp(0.), weightFacRenDown(0.);
  double pdfRMS(0.);

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
	      if(absPdgId==5) nBParton++;
      }

      if(status == 1 && (absPdgId == 11 || absPdgId == 13 || absPdgId == 15)){ // leptons (needed for DY, W stitching)
	        nLeptons++;
	        TLorentzVector tl;
	        tl.SetPxPyPzE(lheParticles[idxParticle][0],
		      lheParticles[idxParticle][1],
		      lheParticles[idxParticle][2],
		      lheParticles[idxParticle][3]);						
        	tlv.push_back(tl);
      }

    }
    
    
    // Gen weights for QCD scales and PDF
    //  https://indico.cern.ch/event/459797/contributions/1961581/attachments/1181555/1800214/mcaod-Feb15-2016.pdf
    const LHEEventProduct* Product = lheEventProduct_.product();
    
    if(Product->weights().size() >= 9) {
        weightFacUp = Product->weights()[1].wgt / Product->originalXWGTUP();
        weightFacDown = Product->weights()[2].wgt / Product->originalXWGTUP();
        weightRenUp = Product->weights()[3].wgt / Product->originalXWGTUP();
        weightRenDown = Product->weights()[6].wgt / Product->originalXWGTUP();
        weightFacRenUp = Product->weights()[4].wgt / Product->originalXWGTUP();
        weightFacRenDown = Product->weights()[8].wgt / Product->originalXWGTUP();
    }    
    std::vector<double> pdfWeights;
    for(unsigned int i = 10; i <= 110 && i < Product->weights().size(); i++) {
        pdfWeights.push_back(Product->weights()[i].wgt / Product->originalXWGTUP());
    }
    
    // Calculate RMS
    if(pdfWeights.size() > 0) {
        double sum = std::accumulate(pdfWeights.begin(), pdfWeights.end(), 0.0);
        double mean = sum / pdfWeights.size();

        double sq_sum = std::inner_product(pdfWeights.begin(), pdfWeights.end(), pdfWeights.begin(), 0.0);
        pdfRMS = std::sqrt(sq_sum / pdfWeights.size() - mean * mean);
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

  nBranches_->lheNl = nLeptons;
  nBranches_->lheNj = nParton;
  nBranches_->lheNb = nBParton;
  nBranches_->lheHT = lheHt_;


  nBranches_->genFacWeightUp = weightFacUp;
  nBranches_->genFacWeightDown = weightFacDown;
  nBranches_->genRenWeightUp = weightRenUp;
  nBranches_->genRenWeightDown = weightRenDown;
  nBranches_->genFacRenWeightUp = weightFacRenUp;
  nBranches_->genFacRenWeightDown = weightFacRenDown;
  nBranches_->PDF_rms = 1. + pdfRMS;

}
