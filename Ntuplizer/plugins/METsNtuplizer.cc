#include "../interface/METsNtuplizer.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include <TFormula.h>
#include "FWCore/Framework/interface/ConsumesCollector.h"

//===================================================================================================================        
METsNtuplizer::METsNtuplizer( 	edm::EDGetTokenT<pat::METCollection>     mettoken    , 
				edm::EDGetTokenT<pat::METCollection>     metpuppitoken    , 
				edm::EDGetTokenT<pat::METCollection>     metmvatoken    , 
 	 			edm::EDGetTokenT<pat::JetCollection>	 jettoken    ,
				edm::EDGetTokenT<pat::MuonCollection> 	 muontoken   ,
				edm::EDGetTokenT<double> 		 rhotoken    ,
				edm::EDGetTokenT<reco::VertexCollection> vtxtoken    ,
				edm::EDGetTokenT<double>	      metSigtoken    ,
				edm::EDGetTokenT<math::Error<2>::type> metCovtoken ,

				std::vector<std::string> 		 jecAK4labels,
				std::vector<std::string>	  	 corrformulas,
				NtupleBranches* 			 nBranches   , 
				std::map< std::string, bool >&                        runFlags 	)
									
: CandidateNtuplizer ( nBranches    )
, metInputToken_     ( mettoken     )
, metpuppiInputToken_( metpuppitoken)
, metmvaInputToken_  ( metmvatoken  )
, jetInputToken_     ( jettoken     )
, muonInputToken_    ( muontoken    )	    
, rhoToken_	     ( rhotoken     )	    
, verticeToken_	     ( vtxtoken     )	
, metSigToken_       ( metSigtoken  )
, metCovToken_       ( metCovtoken  )
, jetCorrLabel_	     ( jecAK4labels )
, corrFormulas_	     ( corrformulas )	
, doMETSVFIT_        ( runFlags["doMETSVFIT"]  )					    
, doMVAMET_        ( runFlags["doMVAMET"]  )					    

{
        if( jetCorrLabel_.size() != 0 ){
	   offsetCorrLabel_.push_back(jetCorrLabel_[0]);
	   initJetCorrFactors();
	   doCorrOnTheFly_ = true;
	}
	else doCorrOnTheFly_ = false;	 
}

//===================================================================================================================
METsNtuplizer::~METsNtuplizer( void )
{
	TypeICorrMap_.clear();
}

//===================================================================================================================
void METsNtuplizer::initJetCorrFactors( void ){

  std::vector<JetCorrectorParameters> vPar;
  for ( std::vector<std::string>::const_iterator payloadBegin = jetCorrLabel_.begin(), payloadEnd = jetCorrLabel_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vPar.push_back(pars);
  }
  
  // Make the FactorizedJetCorrector
  jecAK4_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
  
  vPar.clear();
  for ( std::vector<std::string>::const_iterator payloadBegin = offsetCorrLabel_.begin(), payloadEnd = offsetCorrLabel_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vPar.push_back(pars);
  }
  
  jecOffset_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );  
}

//===================================================================================================================
double METsNtuplizer::getJEC( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ ){
   
   double jetCorrFactor = 1.;
   if ( fabs(rawJetP4.eta()) < jetCorrEtaMax ){     
      jecAK4_->setJetEta( rawJetP4.eta() );
      jecAK4_->setJetPt ( rawJetP4.pt()	);
      jecAK4_->setJetE  ( rawJetP4.energy() );
      jecAK4_->setJetPhi( rawJetP4.phi()    );
      jecAK4_->setJetA  ( jet.jetArea()	);
      jecAK4_->setRho   ( *(rho_.product()) );
      jecAK4_->setNPV   ( vertices_->size() );
      jetCorrFactor = jecAK4_->getCorrection();	  
   }

   reco::Candidate::LorentzVector corrJetP4 = rawJetP4;
   corrJetP4 *= jetCorrFactor;
   
   return jetCorrFactor;
 
}

//===================================================================================================================
double METsNtuplizer::getJECOffset( reco::Candidate::LorentzVector& rawJetP4, const pat::Jet& jet, double& jetCorrEtaMax, std::vector<std::string> jecPayloadNames_ ){
   
   double jetCorrFactor = 1.;
   if ( fabs(rawJetP4.eta()) < jetCorrEtaMax ){     
      jecOffset_->setJetEta( rawJetP4.eta()	);
      jecOffset_->setJetPt ( rawJetP4.pt()	);
      jecOffset_->setJetE  ( rawJetP4.energy()  );
      jecOffset_->setJetPhi( rawJetP4.phi()     );
      jecOffset_->setJetA  ( jet.jetArea()	);
      jecOffset_->setRho   ( *(rho_.product())  );
      jecOffset_->setNPV   ( vertices_->size()  );
      jetCorrFactor = jecOffset_->getCorrection();	  
   }

   reco::Candidate::LorentzVector corrJetP4 = rawJetP4;
   corrJetP4 *= jetCorrFactor;
   
   return jetCorrFactor;
 
}

//===================================================================================================================
void METsNtuplizer::addTypeICorr( edm::Event const & event ){

   TypeICorrMap_.clear();
      
   event.getByToken(jetInputToken_    	, jets_    );
   event.getByToken(rhoToken_          	, rho_     );
   event.getByToken(verticeToken_      	, vertices_);   
   event.getByToken(muonInputToken_   	, muons_   );
        
   bool skipEM_                    = true; 
   double skipEMfractionThreshold_ = 0.9;
   bool skipMuons_                 = true;
   std::string skipMuonSelection_string = "isGlobalMuon | isStandAloneMuon";
   StringCutObjectSelector<reco::Candidate>* skipMuonSelection_ = new StringCutObjectSelector<reco::Candidate>(skipMuonSelection_string,true);
  
   double jetCorrEtaMax_           = 9.9;
   double type1JetPtThreshold_     = 15.0;
   
   double corrEx    = 0;
   double corrEy    = 0;
   double corrSumEt = 0;

   for (const pat::Jet &jet : *jets_) {
	     	   
     double emEnergyFraction = jet.chargedEmEnergyFraction() + jet.neutralEmEnergyFraction();
     if ( skipEM_ && emEnergyFraction > skipEMfractionThreshold_ ) continue;
     
     reco::Candidate::LorentzVector rawJetP4 = jet.correctedP4(0); 
     double corr = getJEC(rawJetP4, jet, jetCorrEtaMax_, jetCorrLabel_);    
              
     if ( skipMuons_ ) {
       const std::vector<reco::CandidatePtr> & cands = jet.daughterPtrVector();
       for ( std::vector<reco::CandidatePtr>::const_iterator cand = cands.begin();
             cand != cands.end(); ++cand ) {
     	 const reco::PFCandidate *pfcand = dynamic_cast<const reco::PFCandidate *>(cand->get());
     	 const reco::Candidate *mu = (pfcand != 0 ? ( pfcand->muonRef().isNonnull() ? pfcand->muonRef().get() : 0) : cand->get());
         if ( mu != 0 && (*skipMuonSelection_)(*mu) ) {
           reco::Candidate::LorentzVector muonP4 = (*cand)->p4();
           rawJetP4 -= muonP4;
         }
       }
     }
           
     reco::Candidate::LorentzVector corrJetP4 = corr*rawJetP4;     
       
     if ( corrJetP4.pt() > type1JetPtThreshold_ ) {

	reco::Candidate::LorentzVector tmpP4 = jet.correctedP4(0);
	corr = getJECOffset(tmpP4, jet, jetCorrEtaMax_, offsetCorrLabel_);
	reco::Candidate::LorentzVector rawJetP4offsetCorr = corr*rawJetP4;
	
	corrEx    -= (corrJetP4.px() - rawJetP4offsetCorr.px());
	corrEy    -= (corrJetP4.py() - rawJetP4offsetCorr.py());
	corrSumEt += (corrJetP4.Et() - rawJetP4offsetCorr.Et());
	
     }
     
 }
 
 TypeICorrMap_["corrEx"]    = corrEx;
 TypeICorrMap_["corrEy"]    = corrEy;
 TypeICorrMap_["corrSumEt"] = corrSumEt;
 
}

//===================================================================================================================
void METsNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){



  // PFMET
	
  event.getByToken(metInputToken_, METs_ );

  if( doCorrOnTheFly_ ) addTypeICorr(event);

  for (const pat::MET &met : *METs_) {
    //const float rawPt	= met.shiftedPt(pat::MET::NoShift, pat::MET::Raw);
    //const float rawPhi  = met.shiftedPhi(pat::MET::NoShift, pat::MET::Raw);
    //const float rawSumEt= met.shiftedSumEt(pat::MET::NoShift, pat::MET::Raw);
    const float rawPt	 = met.uncorPt();//met.shiftedPt(pat::MET::METUncertainty::NoShift, pat::MET::METUncertaintyLevel::Raw);
    const float rawPhi   = met.uncorPhi();//met.shiftedPhi(pat::MET::METUncertainty::NoShift, pat::MET::METUncertaintyLevel::Raw);
    const float rawSumEt = met.uncorSumEt();//met.shiftedSumEt(pat::MET::METUncertainty::NoShift, pat::MET::METUncertaintyLevel::Raw);
        
    TVector2 rawMET_;
    rawMET_.SetMagPhi (rawPt, rawPhi );

    Double_t rawPx = rawMET_.Px();
    Double_t rawPy = rawMET_.Py();
    Double_t rawEt = std::hypot(rawPx,rawPy);  
    
    nBranches_->METraw_et .push_back(rawEt); 
    nBranches_->METraw_phi.push_back(rawPhi);
    nBranches_->METraw_sumEt.push_back(rawSumEt);

    double pxcorr = met.px();
    double pycorr = met.py();
    double et	  = met.et();
    double sumEtcorr = met.sumEt();

    if( doCorrOnTheFly_ ){
       pxcorr = rawPx+TypeICorrMap_["corrEx"];
       pycorr = rawPy+TypeICorrMap_["corrEy"];
       et     = std::hypot(pxcorr,pycorr);
       sumEtcorr = rawSumEt+TypeICorrMap_["corrSumEt"];
    }
    
    TLorentzVector corrmet; corrmet.SetPxPyPzE(pxcorr,pycorr,0.,et);
    nBranches_->MET_et .push_back(et);
    nBranches_->MET_phi.push_back(corrmet.Phi());
    nBranches_->MET_sumEt.push_back(sumEtcorr);
    nBranches_->MET_corrPx.push_back(TypeICorrMap_["corrEx"]);
    nBranches_->MET_corrPy.push_back(TypeICorrMap_["corrEy"]); 	 
 
   
  } 


  // Y.T added 6 Sep. For puppi MET
  event.getByToken(metpuppiInputToken_, METspuppi_ );

  for (const pat::MET &met : *METspuppi_) {
    nBranches_->MET_puppi_et.push_back(met.pt());
    nBranches_->MET_puppi_phi.push_back(met.phi());

  }

  
  // Y.T added 10 Sep. For MVA MET -> add options here 

  std::cout << "doMVAMET_ = " << doMVAMET_ << std::endl;

  if (doMVAMET_) {

    event.getByToken(metmvaInputToken_, METsmva_ );
    
    int nmva = 0;
    
    for (const pat::MET &met : *METsmva_) {
      
      nmva++;
      
      std::vector<float> recoil_pt;
      std::vector<float> recoil_eta;
      std::vector<float> recoil_phi;
      std::vector<int> recoil_pdgId;
      
      recoil_pt.clear();
      recoil_eta.clear();
      recoil_phi.clear();
      recoil_pdgId.clear();

      //    std::cout << "met = " << met.pt() << std::endl;
      
      for(auto name: met.userCandNames()){
	reco::CandidatePtr aRecoCand = met.userCand(name);
	
	recoil_pt.push_back(aRecoCand->p4().Pt());
	recoil_eta.push_back(aRecoCand->p4().Eta());
	recoil_phi.push_back(aRecoCand->p4().Phi());
	recoil_pdgId.push_back(aRecoCand->pdgId());
      }

      nBranches_->MET_mva_et.push_back(met.pt());
      nBranches_->MET_mva_phi.push_back(met.phi());
      nBranches_->MET_mva_cov00.push_back(met.getSignificanceMatrix()(0,0));
      nBranches_->MET_mva_cov10.push_back(met.getSignificanceMatrix()(1,0));
      nBranches_->MET_mva_cov11.push_back(met.getSignificanceMatrix()(1,1));
      
      nBranches_->MET_mva_recoil_pt.push_back(recoil_pt);
      nBranches_->MET_mva_recoil_eta.push_back(recoil_eta);
      nBranches_->MET_mva_recoil_phi.push_back(recoil_phi);
      nBranches_->MET_mva_recoil_pdgId.push_back(recoil_pdgId);
      
    }

    nBranches_->MET_Nmva.push_back(nmva);
  }



  if (doMETSVFIT_) {
 

    event.getByToken (metSigToken_, significanceHandle);
    event.getByToken (metCovToken_, covHandle);
  //event.getByLabel ("METSignificance", "METSignificance", significanceHandle);
  //event.getByLabel ("METSignificance", "METCovariance", covHandle);
 
    nBranches_->MET_significance.push_back( (*significanceHandle));
    nBranches_->MET_cov00.push_back(    (*covHandle)(0,0));
    nBranches_->MET_cov10.push_back(    (*covHandle)(1,0));
    nBranches_->MET_cov11.push_back(    (*covHandle)(1,1));
  //FROM LOW MASS ANALYSIS
  // covMET[0][0] = (*covHandle)(0,0);
  // covMET[1][0] = (*covHandle)(1,0);
  // covMET[0][1] = covMET[1][0]; // (1,0) is the only one saved
  // covMET[1][1] = (*covHandle)(1,1);
  
  }


}

