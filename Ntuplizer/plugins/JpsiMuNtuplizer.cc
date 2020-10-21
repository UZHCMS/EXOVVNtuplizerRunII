#include "../interface/JpsiMuNtuplizer.h"


//===================================================================================================================
JpsiMuNtuplizer::JpsiMuNtuplizer( edm::EDGetTokenT<pat::MuonCollection>    muonToken   ,
                                  edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
                                  edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken,
                                  edm::EDGetTokenT<edm::TriggerResults> triggertoken,
                                  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobject,
                                  edm::EDGetTokenT<reco::GenParticleCollection> genptoken,
                                  std::map< std::string, bool >& runFlags,
                                  NtupleBranches* nBranches )
: CandidateNtuplizer ( nBranches )
  , muonToken_	        ( muonToken )
  , verticeToken_          ( verticeToken )
  , packedpfcandidatesToken_(packedpfcandidatesToken) 
  , HLTtriggersToken_	( triggertoken )
  , triggerObjects_	( triggerobject )
  , genParticlesToken_( genptoken )
  , runOnMC_   (runFlags["runOnMC"])
  , useHammer_  (runFlags["useHammer"])
  , verbose_   (runFlags["verbose"])
   
{

  if(verbose_){
    std::cout << "[JpsiMuNtuplizer] runOnMC    = " << runOnMC_ << std::endl;
    std::cout << "[JpsiMuNtuplizer] UseHammer  = " << useHammer_ << std::endl;
  }



  if(runOnMC_ && useHammer_){

    ran = new TRandom3();
    ran->SetSeed(1);

    if(verbose_) std::cout << "[JpsiMuNtuplizer] Setting up Hammer" << std::endl;

    hammer.setUnits("GeV");

    std::vector<std::string> processes = {"BcJpsiMuNu", "BcJpsiTauNu"};
  
    for(auto proc : processes) {
      if(verbose_) std::cout << "[JpsiMuNtuplizer] \t Hammer added: " << proc << std::endl;
      hammer.includeDecay(proc);
    }


    hammer.setFFInputScheme({{"BcJpsi", "Kiselev"}, {"TauPiPiPi","RCT"}});
    hammer.addFFScheme("Scheme1", {
	{"BcJpsi", "BGLVar"}, 
	  {"TauPiPiPi", "RCT"}
      });

    hammer.addPurePSVertices({"TauPiPiPiNu"}, Hammer::WTerm::DENOMINATOR);
    hammer.addPurePSVertices({"TauPiPiPiNu"}, Hammer::WTerm::NUMERATOR);
    //    hammer.addPurePSVertices({"TauPiNu"}, Hammer::WTerm::DENOMINATOR);
    //    hammer.addPurePSVertices({"TauPiNu"}, Hammer::WTerm::NUMERATOR);
    hammer.initRun();

    if(verbose_) std::cout << "[JpsiMuNtuplizer] Finish setting up Hammer" << std::endl;

//    std::string centralValuesOpt = "BctoJpsiBGLVar: {abcdmatrix: [";
//    centralValuesOpt += "[-0.000435761, -0.00128207, 0.00367311, 0.00449467, -0.0063562, 0.0021275, 0.00271668, 0.00210176, 0.299812, -0.95395, 0.],"; 
//    centralValuesOpt += "[-0.00875884, 0.0361598, -0.11538, -0.152852, 0.817625, 0.530242, -0.0934758, 0.058003, -0.00971328, -0.0086655, 0.],"; 
//    centralValuesOpt += "[0.701822, 0.703935, -0.0695183, -0.0742765, -0.0323615, -0.0215506, 0.00727573, -0.00202158, 0.000970731, -0.00139538, 0.],"; 
//    centralValuesOpt += "[-0.000281186, 0.000638864, 0.000796885, 0.00243559, -0.0146715, 0.00715778, -0.00127426, 0.0820568, -0.302581, -0.094792, -0.944696],"; 
//    centralValuesOpt += "[-0.0178757, -0.0131949, -0.117564, -0.146799, 0.436528, -0.745198, 0.225176, 0.408729, 0.0127131, -0.000151653, 0.018235],"; 
//    centralValuesOpt += "[0.684286, -0.705427, -0.172203, -0.0603275, -0.0136353, 0.024984, -0.00418088, -0.00158171, -0.000972362, -0.000486225, -0.000351983],"; 
//    centralValuesOpt += "[-0.00288391, 0.00160327, -0.0201438, 0.00736481, -0.0154784, 0.0332032, -0.0122592, 0.0665072, 0.90324, 0.28412, -0.311523],"; 
//    centralValuesOpt += "[0.066627, 0.022617, -0.187896, 0.954851, 0.200779, -0.0664202, 0.0196037, -0.0534551, -0.00372263, 0.000996768, -0.00489996],"; 
//    centralValuesOpt += "[0.00379569, 0.00957017, -0.0101066, 0.123275, -0.248312, 0.299884, -0.0916655, 0.901238, -0.0462642, -0.00996493, 0.100667],"; 
//    centralValuesOpt += "[-0.185266, 0.0689267, -0.949214, -0.136689, -0.187331, 0.0447062, 0.0210228, -0.0576221, -0.0172231, -0.00843905, 0.00352652],"; 
//    centralValuesOpt += "[0.00400356, -0.0028094, 0.0393078, 0.0151223, -0.0462619, 0.254824, 0.964935, -0.000850206, 0.00236844, 0.00459156, 0.00012354]]}";
//
//
//    std::cout << "[Hammer]: Central values\n\t" << centralValuesOpt << std::endl;
//    hammer.setOptions(centralValuesOpt);
    hammer.saveOptionCard("Opts.yml", false);

    //    std::cout << "... finishes " << std::endl;

    // add central setting, if needed. 

    //// from https://arxiv.org/pdf/1909.10691.pdf
    //    std::map<std::string, std::vector<double>> paramsBGL {
    //      {"avec", {0.004698, -0.02218, 0.1503}},
    //      {"bvec", {0.003424, -0.02595, 0.3897}},
    //      {"cvec", {-0.003164, 0.08731}},
    //      {"dvec", {0.04011, -0.2132, 0.008157}}
    //    };
    //    std::string centralValuesOpt = "\"BctoJpsiBGL: {";
    //    
    //    for(auto elem: paramsBGL) {
    //      centralValuesOpt += Form("%s: {", elem.first.c_str());
    //      
    //      for(size_t ii=0; ii < elem.second.size(); ii++){
    //	if(ii==elem.second.size()-1){
    //	  centralValuesOpt += Form("%f ", elem.second[ii]);
    //	}else{
    //	  centralValuesOpt += Form("%f, ", elem.second[ii]);
    //	}
    //
    //      }
    //      if(elem.first.c_str()==std::string("dvec")) centralValuesOpt += "}";
    //      else centralValuesOpt += "}, ";
    //    }
    //    centralValuesOpt += "}\"";
    //    std::cout << "[Hammer]: Central values\n\t" << centralValuesOpt << std::endl;
    //    hammer.setOptions(centralValuesOpt);
    //    std::cout << "... finishes " << std::endl;
   

    // Generate FF toys ... 

    for(int imc=0; imc < numberofToys;imc++){

      vector<double> deltas; 
      int idx1 = 0;
      Float_t chi2 = 0;
      
      for(auto pars1: _FFErrNames) {

	if(idx1==10) break;
	Float_t mean = ran->Gaus(_FFmean[idx1], _FFErr[idx1]);
	
	Float_t _chi2 = (mean - _FFmean[idx1])*Inv[idx1]*(mean - _FFmean[idx1]);
	chi2 += _chi2;

	deltas.push_back(mean - _FFmean[idx1]);
	
	//	if(imc<=1) std::cout << imc << " "  << pars1 << " " << mean << std::endl;

	idx1+=1; 
      }

      //      if(chi2 > 11.536) continue;
      
      map<string, double> settings;
      int idx_err = 0;
      
      for(auto pars1: _FFErrNames) {
	
	Float_t newerr = 0;
	
	for(int j=0; j<10; j++) {
	  newerr += deltas[j]*eigVec[idx_err][j];
	}
	

	if(idx_err==10){
	  settings[pars1] = 0;
	}else{
	  settings[pars1] = newerr;
	}

	idx_err += 1;
      }


      //      if(imc <= 1) std::cout << imc << " settings of " << "delta_a0" << ": " << settings["delta_a0"] << std::endl;

      FFdict.push_back(settings);



    }
    
    
    if(verbose_) std::cout << "Saved " << FFdict.size() << " FF variations" << std::endl;



 
  }
}

//===================================================================================================================
JpsiMuNtuplizer::~JpsiMuNtuplizer( void )
{

  if(runOnMC_ && useHammer_){
  
    if(verbose_) std::cout <<"[JpsiMuNtuplizer] Evaluate partial width" << std::endl;
  
    std::vector<std::string> processes = {"BcJpsiTau+Nu", "BcJpsiMu+Nu"};

    for(auto proc : processes) {
      std::map<std::string, double> outRate;
      if(verbose_) std::cout << "\t Process: " << proc << std::endl;
      
      outRate["den"] = hammer.getDenominatorRate(proc);
      
      if(outRate["den"] == 0) {
	std::cout <<"[JpsiMuNtuplizer]: ERROR Failed to get default partial width" << std::endl;
	continue;
      }else{
	if(verbose_) std::cout << Form("[JpsiMuNtuplizer] Default rate: %1.3e", outRate["den"]) << std::endl;
      }
      
      std::map<std::string, double> settings;
      for(auto n: parName) settings["delta_" + n] = 0;      
      //      for(auto pars: _FFErrNames) {
      //	settings[pars] = 0;
      //      }
      
      hammer.setFFEigenvectors("BctoJpsi", "BGLVar", settings);
      
      outRate["Central"] = hammer.getRate(proc, "Scheme1");
      std::cout << Form("[JpsiMuNtuplizer] New rate: %1.3e (ratio = %.3f)", outRate["Central"], outRate["Central"]/outRate["den"]) << std::endl;
      
      nBranches_->hammer_width->SetBinContent(1, outRate["den"]);
      nBranches_->hammer_width->SetBinContent(2, outRate["Central"]);
    
    
      for(int i=0; i<11; i++) { //Loop over eigenVar
	for (int j=0; j<2; j++) { //Loop over pos/neg direction
	  map<string, double> settings;
	  for (int k=0; k<11; k++) { //Loop over parameters
	    settings["delta_" + parName[k]] = eigVar[i][k][j];
	  }

	  hammer.setFFEigenvectors("BctoJpsi", "BGLVar", settings);
	  auto rate = hammer.getRate(proc, "Scheme1");
	  
	  std::string var_name = "eig";
	  var_name += std::to_string(i);
	  var_name += j==0? "_Up" : "_Down";
	  
	  if(verbose_) std::cout << "[JpsiMuNtuplizer] " << var_name << " --> " << Form(": %1.3e", rate) << std::endl;
	  
	  if(var_name==std::string("eig0_Up")) nBranches_->hammer_width->SetBinContent(3, rate);
	  else if(var_name==std::string("eig0_Down")) nBranches_->hammer_width->SetBinContent(4, rate);
	  else if(var_name==std::string("eig1_Up")) nBranches_->hammer_width->SetBinContent(5, rate);
	  else if(var_name==std::string("eig1_Down")) nBranches_->hammer_width->SetBinContent(6, rate);
	  else if(var_name==std::string("eig2_Up")) nBranches_->hammer_width->SetBinContent(7, rate);
	  else if(var_name==std::string("eig2_Down")) nBranches_->hammer_width->SetBinContent(8, rate);
	  
	  else if(var_name==std::string("eig3_Up")) nBranches_->hammer_width->SetBinContent(9, rate);
	  else if(var_name==std::string("eig3_Down")) nBranches_->hammer_width->SetBinContent(10, rate);
	  else if(var_name==std::string("eig4_Up")) nBranches_->hammer_width->SetBinContent(11, rate);
	  else if(var_name==std::string("eig4_Down")) nBranches_->hammer_width->SetBinContent(12, rate);
	  else if(var_name==std::string("eig5_Up")) nBranches_->hammer_width->SetBinContent(13, rate);
	  else if(var_name==std::string("eig5_Down")) nBranches_->hammer_width->SetBinContent(14, rate);
	  
	  else if(var_name==std::string("eig6_Up")) nBranches_->hammer_width->SetBinContent(15, rate);
	  else if(var_name==std::string("eig6_Down")) nBranches_->hammer_width->SetBinContent(16, rate);
	  else if(var_name==std::string("eig7_Up")) nBranches_->hammer_width->SetBinContent(17, rate);
	  else if(var_name==std::string("eig7_Down")) nBranches_->hammer_width->SetBinContent(18, rate);
	  
	  else if(var_name==std::string("eig8_Up")) nBranches_->hammer_width->SetBinContent(19, rate);
	  else if(var_name==std::string("eig8_Down")) nBranches_->hammer_width->SetBinContent(20, rate);
	  else if(var_name==std::string("eig9_Up")) nBranches_->hammer_width->SetBinContent(21, rate);
	  else if(var_name==std::string("eig9_Down")) nBranches_->hammer_width->SetBinContent(22, rate);
	  else if(var_name==std::string("eig10_Up")) nBranches_->hammer_width->SetBinContent(23, rate);
	  else if(var_name==std::string("eig10_Down")) nBranches_->hammer_width->SetBinContent(24, rate);
	  



//	  string var_name = "BGL" + varName[i];
//	  var_name += j==0? "Up" : "Down";
//	  outRate[var_name] = rate;
//	  if(verbose) {cout << var_name << Form(": %1.3e", rate) << endl;}
	}
      }


    }
  }
}


bool JpsiMuNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  if(verbose_) std::cout << "[JpsiMuNtuplizer] ---------------- event, run, lumi = " << event.id().event() << " " << event.id().run() << " " << event.id().luminosityBlock() << "----------------" << std::endl;


  // for gen. matching
  TVector3 genvertex(-9.,-9.,-9.);

  std::vector<TLorentzVector> gen_nr_mu;
  std::vector<TLorentzVector> gen_jpsi_mu;

  TLorentzVector pB_gen;
  TLorentzVector pJpsi_gen;

  if(runOnMC_){

    event.getByToken(genParticlesToken_ , genParticles_); 
    
    for( unsigned p=0; p < genParticles_->size(); ++p){
	  
      if(!(TMath::Abs((*genParticles_)[p].pdgId())==541 && (*genParticles_)[p].status()==2)) continue;
    
      auto _part_ = (*genParticles_)[p];
    
      bool isJpsi = false;
      bool isTau = false;
    
      for(auto d : _part_.daughterRefVector()) {
	if(TMath::Abs(d->pdgId()) == 443){
	  isJpsi = true;
	}
	if(TMath::Abs(d->pdgId()) == 13 || TMath::Abs(d->pdgId()) == 15){
	  isTau = true;
	}
      }
      
      if(!(isJpsi && isTau)) continue;
    
      pB_gen.SetPtEtaPhiM(_part_.pt(), _part_.eta(), _part_.phi(), _part_.mass());
    
      genvertex = aux.getVertex((*genParticles_)[p]);
    
      for(auto d : _part_.daughterRefVector()) {
	if(TMath::Abs(d->pdgId()) == 443){
	
	  pJpsi_gen.SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
	
	  for (auto dd : d->daughterRefVector()) {
	    if(TMath::Abs(dd->pdgId())==13){
	      TLorentzVector jpsi_muon_tlv;
	      jpsi_muon_tlv.SetPtEtaPhiM(dd->pt(), dd->eta(), dd->phi(), dd->mass());
	      gen_jpsi_mu.push_back(jpsi_muon_tlv);
	    }
	  }
	
	}
	if(TMath::Abs(d->pdgId())==13){
	  TLorentzVector nojpsi_muon_tlv;
	  nojpsi_muon_tlv.SetPtEtaPhiM(d->pt(), d->eta(), d->phi(), d->mass());
	
	  gen_nr_mu.push_back(nojpsi_muon_tlv);
	}


	if(TMath::Abs(d->pdgId())==15){

	  //	  std::cout << "gen:" << d->pdgId() << std::endl;
	  for (auto dd : d->daughterRefVector()) {
	    //	    std::cout << "\t --> " << dd->pdgId() << std::endl;
	    if(TMath::Abs(dd->pdgId())==13){

	      TLorentzVector nojpsi_muon_tlv;
	      nojpsi_muon_tlv.SetPtEtaPhiM(dd->pt(), dd->eta(), dd->phi(), dd->mass());
	      
	      gen_nr_mu.push_back(nojpsi_muon_tlv);
	      
	    }
	  }
	}


      }
    }
  }


  TLorentzVector q2_gen; 
  q2_gen = pB_gen - pJpsi_gen;

  if(runOnMC_){ 
    nBranches_->q2_nocut->Fill(q2_gen.M2());
  }
  /********************************************************************
   *
   * Step1: check if the J/psi trigger is fired.
   * Namely, HLT_DoubleMu4_JpsiTrk_Displaced_v
   *
   ********************************************************************/


  event.getByToken(HLTtriggersToken_, HLTtriggers_);
  nBranches_->cutflow_perevt->Fill(0);
  
  bool isTriggered = false;
  const edm::TriggerNames& trigNames = event.triggerNames(*HLTtriggers_);
  std::string finalTriggerName="";
  std::string finalTriggerFilterObjName="";

  for (unsigned int i = 0, n = HLTtriggers_->size(); i < n; ++i) {

    if(trigNames.triggerName(i).find("HLT_DoubleMu4_JpsiTrk_Displaced_v")!= std::string::npos){
      nBranches_->HLT_BPH_isFired[trigNames.triggerName(i)] = HLTtriggers_->accept(i);
      if(HLTtriggers_->accept(i)){
	isTriggered = true;
	finalTriggerName=trigNames.triggerName(i);  
	finalTriggerFilterObjName="hltJpsiTkVertexFilter";
      }
    }
  }

  if(!isTriggered) return false;
  nBranches_->cutflow_perevt->Fill(1);

  if(verbose_) std::cout << "[JpsiMuNtuplizer] Trigger fired" << std::endl;

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
  std::vector<int> muoncollection_id;
  muoncollection.clear();
  muoncollection_id.clear();

  for(size_t imuon = 0; imuon < muons_->size(); ++ imuon){
      
    const pat::Muon & muon = (*muons_)[imuon];
      
    if(muon.pt() < 4) continue;
    if(TMath::Abs(muon.eta()) > 2.4) continue;
    if(!(muon.track().isNonnull())) continue;
      
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

    muoncollection.push_back(muon);
    muoncollection_id.push_back(imuon);
  }


  if(!( muoncollection.size() >= 2)) return false;
  nBranches_->cutflow_perevt->Fill(2);
    
  if(verbose_) std::cout << "[JpsiMuNtuplizer] At least 2 muons are found" << std::endl;

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
      std::tie(vprob_jpsi, vertex_jpsi) = aux.vertexProb(transient_tracks_dimuon);
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

  if(jpsi_max_pt == -1) return false;
  nBranches_->cutflow_perevt->Fill(3);

  if(verbose_) std::cout << "[JpsiMuNtuplizer] J/psi found" << std::endl;

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

  muonParticles.push_back(pFactory.particle(tt1_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));
  muonParticles.push_back(pFactory.particle(tt2_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));
  
  //creating the vertex fitter
  KinematicParticleVertexFitter kpvFitter;

  //reconstructing a J/Psi decay
  RefCountedKinematicTree jpTree = kpvFitter.fit(muonParticles);

  if(jpTree->isEmpty() || !jpTree->isValid() || !jpTree->isConsistent()) return false;

  nBranches_->cutflow_perevt->Fill(4);

  //creating the particle fitter
  KinematicParticleFitter csFitter;

  // creating the constraint
  KinematicConstraint* jpsi_constraint = new MassKinematicConstraint(aux.jpsi_mass, aux.jp_m_sigma);

  //the constrained fit
  jpTree = csFitter.fit(jpsi_constraint, jpTree);

  //getting the J/Psi KinematicParticle
  jpTree->movePointerToTheTop();
  RefCountedKinematicParticle jpsi_part = jpTree->currentParticle();
  if(!jpsi_part->currentState().isValid()) return false;
  nBranches_->cutflow_perevt->Fill(5);

  RefCountedKinematicVertex jpsi_vertex = jpTree->currentDecayVertex();
  if(!jpsi_vertex->vertexIsValid()) return false; 
  nBranches_->cutflow_perevt->Fill(6);

  if(TMath::Prob(jpsi_vertex->chiSquared(), jpsi_vertex->degreesOfFreedom()) <=0) return false;
  nBranches_->cutflow_perevt->Fill(7);

  std::vector< RefCountedKinematicParticle > jpsi_children = jpTree->finalStateParticles();

  math::PtEtaPhiMLorentzVector mu1_fit = aux.daughter_p4(jpsi_children, 0);
  math::PtEtaPhiMLorentzVector mu2_fit = aux.daughter_p4(jpsi_children, 1);



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
   * and found that 1) or 2) is the best. Here, we use method (1).
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
    
    particle_cand cand = aux.calculateIPvariables(extrapolator, jpsi_part, jpsi_vertex, *vtx);

    //    if(TMath::Abs(cand.pvip) < max_criteria){    
    //      max_criteria = TMath::Abs(cand.pvip);

    if(TMath::Abs(cand.lip) < max_criteria){
      max_criteria = TMath::Abs(cand.lip);
      closestVertex = *vtx;
    }

    counter += 1;
  }


  if(!(muoncollection[mcidx_mu1].isSoftMuon(closestVertex) > 0.5 && muoncollection[mcidx_mu2].isSoftMuon(closestVertex) > 0.5)) return false;
  nBranches_->cutflow_perevt->Fill(8);

  if(verbose_) std::cout << "[JpsiMuNtuplizer] J/psi soft muon ID passed" << std::endl;

  particle_cand JPcand;
  JPcand = aux.calculateIPvariables(extrapolator, jpsi_part, jpsi_vertex, closestVertex);
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
    //        mu3passpteta = true;
    // shall we put delta z cut here maybe? 

	
    //	std::cout << "muon3 pt = " << muon.pt() << std::endl;
	
    max_pt3=2;
    if(muon.pt() > max_pt3){
      if(!doMultipleMuon3) { max_pt3 = muon.pt();}
      mu3_vec.push_back(muon);
      if (!doMultipleMuon3) break;
    }
  }

  if(max_pt3==-1) return false;
  nBranches_->cutflow_perevt->Fill(9);


  for( auto mu3: mu3_vec){

    const reco::TrackRef track3_muon = mu3.muonBestTrack();
    reco::TransientTrack tt3_muon = (*builder).build(track3_muon);
    /********************************************************************
     *
     * Step7: bbPV vertex refit by excluding 3 muon candidates
     *
     ********************************************************************/

    event.getByToken( packedpfcandidatesToken_               , packedpfcandidates_      ); 

    std::vector<pat::PackedCandidate> mypfcollection; 
    std::vector<reco::TransientTrack> mytracks;

    for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){   
    
      pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
      if(!pf.hasTrackDetails()) continue;

      if(pf.vertexRef()->z()!=closestVertex.position().z()) continue;

      Float_t _dR1 = reco::deltaR(pf.eta(), pf.phi(), 
				  track1_muon->eta(), track1_muon->phi());

      Float_t _dR2 = reco::deltaR(pf.eta(), pf.phi(), 
				  track2_muon->eta(), track2_muon->phi());

      Float_t _dR3 = reco::deltaR(pf.eta(), pf.phi(), 
				  track3_muon->eta(), track3_muon->phi());


      if(_dR1 < 0.1 || _dR2 < 0.1 || _dR3 < 0.1){

	if(TMath::Abs(pf.pdgId()) == 13) continue;
      }

      reco::TransientTrack  tt_track = (*builder).build(pf.pseudoTrack());

      mytracks.push_back(tt_track);

      if(pf.pt() > 0.5) mypfcollection.push_back(pf);
    }
  


    /********************************************************************
     *
     * Step7: Kinematic fit for B candidate
     *
     ********************************************************************/
    std::vector<RefCountedKinematicParticle> allParticles;

    allParticles.push_back(pFactory.particle(tt3_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));
    allParticles.push_back(jpsi_part);

    RefCountedKinematicTree bcTree = kpvFitter.fit(allParticles);

    if(bcTree->isEmpty() || !bcTree->isValid() || !bcTree->isConsistent()) continue;

    RefCountedKinematicParticle bc_part = bcTree->currentParticle();
    if(!bc_part->currentState().isValid()) continue;

    RefCountedKinematicVertex bc_vertex = bcTree->currentDecayVertex();
    if(!bc_vertex->vertexIsValid()) continue;

    particle_cand Bcand;

    Bcand = aux.calculateIPvariables(extrapolator, bc_part, bc_vertex, closestVertex);

    std::vector< RefCountedKinematicParticle > bc_children = bcTree->finalStateParticles();
    math::PtEtaPhiMLorentzVector mu3_fit = aux.daughter_p4(bc_children, 0);

    /********************************************************************
     *
     * Step8: calculation of the isolation quantities
     *
     ********************************************************************/

    TLorentzVector tlv_mu3;
    tlv_mu3.SetPtEtaPhiM(mu3.pt(), mu3.eta(), mu3.phi(), mu3.mass());

    TLorentzVector tlv_B = jpsi_tlv_highest + tlv_mu3;

    std::vector<reco::TransientTrack> transient_tracks_trimuon;
    transient_tracks_trimuon.push_back(tt1_muon);
    transient_tracks_trimuon.push_back(tt2_muon);
    transient_tracks_trimuon.push_back(tt3_muon);

    Float_t vprob_bc = -9;
    TransientVertex vertex_bc;
    std::tie(vprob_bc, vertex_bc) = aux.vertexProb(transient_tracks_trimuon);

    Float_t iso = 0;
    Int_t ntracks = 0;
    Float_t iso_mindoca = 999; 
  

    for(size_t itrk = 0; itrk < mypfcollection.size();  itrk++){
    
      Float_t _dR = reco::deltaR(mypfcollection[itrk].eta(), mypfcollection[itrk].phi(),
				 bc_part->currentState().globalMomentum().eta(), bc_part->currentState().globalMomentum().phi());

      if(_dR < 0.7) iso += mypfcollection[itrk].pt();


      TrajectoryStateOnSurface tsos_pf = extrapolator.extrapolate(mytracks[itrk].impactPointState(), bc_vertex->position());

      VertexDistance3D a3d_pf;  

      std::pair<bool,Measurement1D> cur3DIP_pf = aux.absoluteImpactParameter(tsos_pf, bc_vertex, a3d_pf);
    
      Float_t pvip_pf = cur3DIP_pf.second.value();

   
      if(pvip_pf < 0.03) ntracks+=1;

      if(iso_mindoca > pvip_pf) iso_mindoca = pvip_pf;
    }

    /********************************************************************
     *
     * Step9: Filling normal branches
     *
     ********************************************************************/

    nBranches_->JpsiMu_mu1_pt.push_back(mu1_fit.pt());
    nBranches_->JpsiMu_mu1_eta.push_back(mu1_fit.eta());
    nBranches_->JpsiMu_mu1_phi.push_back(mu1_fit.phi());
    nBranches_->JpsiMu_mu1_mass.push_back(mu1_fit.mass());
    nBranches_->JpsiMu_mu1_unfit_pt.push_back(muoncollection[mcidx_mu1].pt());
    nBranches_->JpsiMu_mu1_unfit_eta.push_back(muoncollection[mcidx_mu1].eta());
    nBranches_->JpsiMu_mu1_unfit_phi.push_back(muoncollection[mcidx_mu1].phi());
    nBranches_->JpsiMu_mu1_unfit_mass.push_back(muoncollection[mcidx_mu1].mass());
    nBranches_->JpsiMu_mu1_q.push_back(muoncollection[mcidx_mu1].charge());
    nBranches_->JpsiMu_mu1_isLoose.push_back(muoncollection[mcidx_mu1].isLooseMuon());
    nBranches_->JpsiMu_mu1_isTight.push_back(muoncollection[mcidx_mu1].isTightMuon(closestVertex));
    nBranches_->JpsiMu_mu1_isPF.push_back(muoncollection[mcidx_mu1].isPFMuon());
    nBranches_->JpsiMu_mu1_isGlobal.push_back(muoncollection[mcidx_mu1].isGlobalMuon());
    nBranches_->JpsiMu_mu1_isTracker.push_back(muoncollection[mcidx_mu1].isTrackerMuon());
    nBranches_->JpsiMu_mu1_isSoft.push_back(muoncollection[mcidx_mu1].isSoftMuon(closestVertex));
    nBranches_->JpsiMu_mu1_vx.push_back(muoncollection[mcidx_mu1].vx());
    nBranches_->JpsiMu_mu1_vy.push_back(muoncollection[mcidx_mu1].vy());
    nBranches_->JpsiMu_mu1_vz.push_back(muoncollection[mcidx_mu1].vz());
    nBranches_->JpsiMu_mu1_dbiso.push_back(aux.MuonPFIso(muoncollection[mcidx_mu1]));

    nBranches_->JpsiMu_mu2_pt.push_back(mu2_fit.pt());
    nBranches_->JpsiMu_mu2_eta.push_back(mu2_fit.eta());
    nBranches_->JpsiMu_mu2_phi.push_back(mu2_fit.phi());
    nBranches_->JpsiMu_mu2_mass.push_back(mu2_fit.mass());
    nBranches_->JpsiMu_mu2_unfit_pt.push_back(muoncollection[mcidx_mu2].pt());
    nBranches_->JpsiMu_mu2_unfit_eta.push_back(muoncollection[mcidx_mu2].eta());
    nBranches_->JpsiMu_mu2_unfit_phi.push_back(muoncollection[mcidx_mu2].phi());
    nBranches_->JpsiMu_mu2_unfit_mass.push_back(muoncollection[mcidx_mu2].mass());
    nBranches_->JpsiMu_mu2_q.push_back(muoncollection[mcidx_mu2].charge());
    nBranches_->JpsiMu_mu2_isLoose.push_back(muoncollection[mcidx_mu2].isLooseMuon());
    nBranches_->JpsiMu_mu2_isTight.push_back(muoncollection[mcidx_mu2].isTightMuon(closestVertex));
    nBranches_->JpsiMu_mu2_isPF.push_back(muoncollection[mcidx_mu2].isPFMuon());
    nBranches_->JpsiMu_mu2_isGlobal.push_back(muoncollection[mcidx_mu2].isGlobalMuon());
    nBranches_->JpsiMu_mu2_isTracker.push_back(muoncollection[mcidx_mu2].isTrackerMuon());
    nBranches_->JpsiMu_mu2_isSoft.push_back(muoncollection[mcidx_mu2].isSoftMuon(closestVertex));
    nBranches_->JpsiMu_mu2_vx.push_back(muoncollection[mcidx_mu2].vx());
    nBranches_->JpsiMu_mu2_vy.push_back(muoncollection[mcidx_mu2].vy());
    nBranches_->JpsiMu_mu2_vz.push_back(muoncollection[mcidx_mu2].vz());
    nBranches_->JpsiMu_mu2_dbiso.push_back(aux.MuonPFIso(muoncollection[mcidx_mu2]));

    nBranches_->JpsiMu_mu3_unfit_pt.push_back(mu3.pt());
    nBranches_->JpsiMu_mu3_unfit_eta.push_back(mu3.eta());
    nBranches_->JpsiMu_mu3_unfit_phi.push_back(mu3.phi());
    nBranches_->JpsiMu_mu3_unfit_mass.push_back(mu3.mass());
    nBranches_->JpsiMu_mu3_pt.push_back(mu3_fit.pt());
    nBranches_->JpsiMu_mu3_eta.push_back(mu3_fit.eta());
    nBranches_->JpsiMu_mu3_phi.push_back(mu3_fit.phi());
    nBranches_->JpsiMu_mu3_mass.push_back(mu3_fit.mass());
    nBranches_->JpsiMu_mu3_q.push_back(mu3.charge());
    nBranches_->JpsiMu_mu3_isLoose.push_back(mu3.isLooseMuon());
    nBranches_->JpsiMu_mu3_isTight.push_back(mu3.isTightMuon(closestVertex));
    nBranches_->JpsiMu_mu3_isPF.push_back(mu3.isPFMuon());
    nBranches_->JpsiMu_mu3_isGlobal.push_back(mu3.isGlobalMuon());
    nBranches_->JpsiMu_mu3_isTracker.push_back(mu3.isTrackerMuon());
    nBranches_->JpsiMu_mu3_isSoft.push_back(mu3.isSoftMuon(closestVertex));
    nBranches_->JpsiMu_mu3_vx.push_back(mu3.vx());
    nBranches_->JpsiMu_mu3_vy.push_back(mu3.vy());
    nBranches_->JpsiMu_mu3_vz.push_back(mu3.vz());
    //        nBranches_->JpsiMu_mu3_iso.push_back(iso_mu3);
    nBranches_->JpsiMu_mu3_dbiso.push_back(aux.MuonPFIso(mu3));

    std::vector<RefCountedKinematicParticle> mu13;
    mu13.push_back(pFactory.particle(tt1_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));
    mu13.push_back(pFactory.particle(tt3_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));

    std::vector<RefCountedKinematicParticle> mu23;
    mu23.push_back(pFactory.particle(tt2_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));
    mu23.push_back(pFactory.particle(tt3_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));

    nBranches_->JpsiMu_mu3_doca2mu1.push_back(aux.getMaxDoca(mu13));
    nBranches_->JpsiMu_mu3_doca2mu2.push_back(aux.getMaxDoca(mu23));

    nBranches_->JpsiMu_B_pt.push_back(bc_part->currentState().globalMomentum().perp());
    nBranches_->JpsiMu_B_eta.push_back(bc_part->currentState().globalMomentum().eta());
    nBranches_->JpsiMu_B_phi.push_back(bc_part->currentState().globalMomentum().phi());
    nBranches_->JpsiMu_B_mass.push_back(bc_part->currentState().mass());
    nBranches_->JpsiMu_B_vprob.push_back(TMath::Prob(bc_part->chiSquared(), bc_part->degreesOfFreedom()));
    nBranches_->JpsiMu_B_lip.push_back(Bcand.lip);
    nBranches_->JpsiMu_B_lips.push_back(Bcand.lips);
    nBranches_->JpsiMu_B_pvip.push_back(Bcand.pvip);
    nBranches_->JpsiMu_B_pvips.push_back(Bcand.pvips);
    nBranches_->JpsiMu_B_fl3d.push_back(Bcand.fl3d);
    nBranches_->JpsiMu_B_fls3d.push_back(Bcand.fls3d);
    nBranches_->JpsiMu_B_alpha.push_back(Bcand.alpha);


    std::vector<RefCountedKinematicParticle> allParticles4doc;

    allParticles4doc.push_back(pFactory.particle(tt1_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));
    allParticles4doc.push_back(pFactory.particle(tt2_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));
    allParticles4doc.push_back(pFactory.particle(tt3_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));

    nBranches_->JpsiMu_B_maxdoca.push_back(aux.getMaxDoca(allParticles4doc));
    nBranches_->JpsiMu_B_mindoca.push_back(aux.getMinDoca(allParticles4doc));
    nBranches_->JpsiMu_B_vx.push_back(bc_vertex->vertexState().position().x());
    nBranches_->JpsiMu_B_vy.push_back(bc_vertex->vertexState().position().y());
    nBranches_->JpsiMu_B_vz.push_back(bc_vertex->vertexState().position().z());  

    nBranches_->JpsiMu_B_iso.push_back(iso);
    nBranches_->JpsiMu_B_iso_ntracks.push_back(ntracks);
    nBranches_->JpsiMu_B_iso_mindoca.push_back(iso_mindoca);


    nBranches_->JpsiMu_B_unfit_pt.push_back(tlv_B.Pt());
    nBranches_->JpsiMu_B_unfit_mass.push_back(tlv_B.M());
    nBranches_->JpsiMu_B_unfit_vprob.push_back(vprob_bc);
  
    if(vprob_bc!=-9){
      nBranches_->JpsiMu_B_unfit_vx.push_back(vertex_bc.position().x());
      nBranches_->JpsiMu_B_unfit_vy.push_back(vertex_bc.position().y());
      nBranches_->JpsiMu_B_unfit_vz.push_back(vertex_bc.position().z());
    }

    TLorentzVector Tlv_B;
    TLorentzVector Tlv_Jpsi;
    TLorentzVector Tlv_mu3;

    Tlv_B.SetPtEtaPhiM(bc_part->currentState().globalMomentum().perp(),
		       bc_part->currentState().globalMomentum().eta(),
		       bc_part->currentState().globalMomentum().phi(),
		       bc_part->currentState().mass()
		       );

    // calculation of the corrected mass
    TVector3 *bvector = new TVector3(Tlv_B.Px(), Tlv_B.Py(), Tlv_B.Pz());
      
    Float_t pperp = bvector->Mag()*TMath::Sin(TMath::ACos(Bcand.alpha));
    
    Float_t mcorr = pperp + TMath::Sqrt(pperp*pperp + Tlv_B.M()*Tlv_B.M());
    
    nBranches_->JpsiMu_B_mcorr.push_back(mcorr);
    //


    Tlv_B *= aux.mass_Bc/bc_part->currentState().mass();
	
    Tlv_Jpsi.SetPtEtaPhiM(jpsi_part->currentState().globalMomentum().perp(),
			  jpsi_part->currentState().globalMomentum().eta(),
			  jpsi_part->currentState().globalMomentum().phi(),
			  jpsi_part->currentState().mass());
      
    Tlv_mu3.SetPtEtaPhiM(mu3_fit.pt(),
			 mu3_fit.eta(),
			 mu3_fit.phi(),
			 mu3_fit.mass());

    Float_t q2 = (Tlv_B - Tlv_Jpsi).M2();


    nBranches_->JpsiMu_B_q2.push_back(q2);


    Float_t mm2 = (Tlv_B - Tlv_Jpsi - Tlv_mu3).M2();
    Float_t ptmiss = (Tlv_B - Tlv_Jpsi - Tlv_mu3).Pt();

    nBranches_->JpsiMu_B_mm2.push_back(mm2);
    nBranches_->JpsiMu_B_ptmiss.push_back(ptmiss);

    Tlv_mu3.Boost( -Tlv_B.BoostVector() );

    nBranches_->JpsiMu_B_Es.push_back(Tlv_mu3.E()); 

    nBranches_->JpsiMu_B_ptback.push_back(Tlv_B.Pt()); 



  
    bool flag_nr_match = false;
    if(gen_nr_mu.size()==1){
      Float_t _dR = reco::deltaR(gen_nr_mu[0].Eta(), gen_nr_mu[0].Phi(), 
				 mu3_fit.eta(), mu3_fit.phi());

      if(_dR < 0.1) flag_nr_match = true;
    }

    bool flag_jpsi_match = false;
    if(gen_jpsi_mu.size()==2){

      Float_t _dR_11 = reco::deltaR(gen_jpsi_mu[0].Eta(), gen_jpsi_mu[0].Phi(), 
				    mu1_fit.eta(), mu1_fit.phi());
      //                                          muoncollection[mcidx_mu1].eta(), muoncollection[mcidx_mu1].phi());
    
      Float_t _dR_22 = reco::deltaR(gen_jpsi_mu[1].Eta(), gen_jpsi_mu[1].Phi(), 
				    mu2_fit.eta(), mu2_fit.phi());
      //                                          muoncollection[mcidx_mu2].eta(), muoncollection[mcidx_mu2].phi());

      Float_t _dR_21 = reco::deltaR(gen_jpsi_mu[1].Eta(), gen_jpsi_mu[1].Phi(), 
				    mu1_fit.eta(), mu1_fit.phi());
      //                                          muoncollection[mcidx_mu1].eta(), muoncollection[mcidx_mu1].phi());
    
      Float_t _dR_12 = reco::deltaR(gen_jpsi_mu[0].Eta(), gen_jpsi_mu[0].Phi(), 
				    mu2_fit.eta(), mu2_fit.phi());
      //                                          muoncollection[mcidx_mu2].eta(), muoncollection[mcidx_mu2].phi());
      

      if(_dR_11 < 0.1 && _dR_22 < 0.1) flag_jpsi_match = true;
      if(_dR_21 < 0.1 && _dR_12 < 0.1) flag_jpsi_match = true;
    }
    nBranches_->JpsiMu_isgenmatched.push_back((int)flag_jpsi_match);
    nBranches_->JpsiMu_mu3_isgenmatched.push_back((int)flag_nr_match);


  }





  nBranches_->JpsiMu_PV_vx.push_back(vertices_->begin()->position().x());
  nBranches_->JpsiMu_PV_vy.push_back(vertices_->begin()->position().y());
  nBranches_->JpsiMu_PV_vz.push_back(vertices_->begin()->position().z());
    
  nBranches_->JpsiMu_bbPV_vx.push_back(closestVertex.position().x());
  nBranches_->JpsiMu_bbPV_vy.push_back(closestVertex.position().y());
  nBranches_->JpsiMu_bbPV_vz.push_back(closestVertex.position().z());

  nBranches_->JpsiMu_Jpsi_pt.push_back(jpsi_part->currentState().globalMomentum().perp());
  nBranches_->JpsiMu_Jpsi_eta.push_back(jpsi_part->currentState().globalMomentum().eta());
  nBranches_->JpsiMu_Jpsi_phi.push_back(jpsi_part->currentState().globalMomentum().phi());
  nBranches_->JpsiMu_Jpsi_mass.push_back(jpsi_part->currentState().mass());
  nBranches_->JpsiMu_Jpsi_vprob.push_back(TMath::Prob(jpsi_part->chiSquared(), jpsi_part->degreesOfFreedom()));
  nBranches_->JpsiMu_Jpsi_lip.push_back(JPcand.lip);
  nBranches_->JpsiMu_Jpsi_lips.push_back(JPcand.lips);
  nBranches_->JpsiMu_Jpsi_pvip.push_back(JPcand.pvip);
  nBranches_->JpsiMu_Jpsi_pvips.push_back(JPcand.pvips);
  nBranches_->JpsiMu_Jpsi_fl3d.push_back(JPcand.fl3d);
  nBranches_->JpsiMu_Jpsi_fls3d.push_back(JPcand.fls3d);
  nBranches_->JpsiMu_Jpsi_alpha.push_back(JPcand.alpha);
  nBranches_->JpsiMu_Jpsi_maxdoca.push_back(aux.getMaxDoca(muonParticles));
  nBranches_->JpsiMu_Jpsi_mindoca.push_back(aux.getMinDoca(muonParticles));
  nBranches_->JpsiMu_Jpsi_vx.push_back(jpsi_vertex->vertexState().position().x());
  nBranches_->JpsiMu_Jpsi_vy.push_back(jpsi_vertex->vertexState().position().y());
  nBranches_->JpsiMu_Jpsi_vz.push_back(jpsi_vertex->vertexState().position().z());  
  nBranches_->JpsiMu_Jpsi_unfit_pt.push_back(jpsi_tlv_highest.Pt());
  nBranches_->JpsiMu_Jpsi_unfit_mass.push_back(jpsi_tlv_highest.M());
  nBranches_->JpsiMu_Jpsi_unfit_vprob.push_back(jpsi_vprob_highest);
    
  if(jpsi_vprob_highest!=-9){
    nBranches_->JpsiMu_Jpsi_unfit_vx.push_back(jpsi_vertex_highest.position().x());
    nBranches_->JpsiMu_Jpsi_unfit_vy.push_back(jpsi_vertex_highest.position().y());
    nBranches_->JpsiMu_Jpsi_unfit_vz.push_back(jpsi_vertex_highest.position().z());
  }
    


  nBranches_->IsJpsiMu.push_back(1.);
  nBranches_->JpsiMu_nCandidates.push_back(nBranches_->JpsiMu_mu1_pt.size());


  if(!runOnMC_) return true;

  if(useHammer_){
    hammer.initEvent();
    
    Hammer::Process Bc2JpsiLNu;
    
    int idxB = -1;
    std::vector<size_t> Bvtx_idxs;
    int idxTau = -1;
    std::vector<size_t> Tauvtx_idxs;
    int idxJpsi = -1;
    std::vector<size_t> Jpsivtx_idxs;

//    for( unsigned p=0; p < genParticles_->size(); ++p){
//	  
//
//      if(!(TMath::Abs((*genParticles_)[p].pdgId())==541 && (*genParticles_)[p].status()==2)) continue;
//
//      auto _part_ = (*genParticles_)[p];
//      std::cout << "[JpsiMuNtuplizer] Hammer: " << _part_.pdgId() << " (" << _part_.status() << ")" << std::endl;
//
//      bool isJpsi = false;
//      bool isTau = false;
//      
//
//      for(auto d : _part_.daughterRefVector()) {	
//	if(TMath::Abs(d->pdgId()) == 443) isJpsi = true;
//	if(TMath::Abs(d->pdgId()) == 15 || TMath::Abs(d->pdgId()) == 13) isTau = true;
//	
//	std::cout << "[JpsiMuNtuplizer] Hammer: ---> " << d->pdgId() << " (" << d->status() << ")" << std::endl;
//
//	for (auto dd : d->daughterRefVector()) {
//	  std::cout << "[JpsiMuNtuplizer] Hammer: -----> " << dd->pdgId() << " (" << dd->status()  << ")" << std::endl;
//	}
//      }
//      std::cout << isJpsi << " " << isTau  << " " << (isJpsi==true && isTau==true) << std::endl;
//      if(!(isJpsi==true && isTau==true)){
//	std::cout << "[JpsiMuNtuplizer] Hammer: This B is rejected !!!!" << std::endl;
//      }      
//    }


    for( unsigned p=0; p < genParticles_->size(); ++p){
	  
      // Bc daughters loop
      // only allow Bc+ as the MC is produced as such
      // ie if we take Bc-, we take probe side by mistake ... 
      //      if(!((*genParticles_)[p].pdgId()==541 && (*genParticles_)[p].status()==2)) continue;

      
      if(!(TMath::Abs((*genParticles_)[p].pdgId())==541 && (*genParticles_)[p].status()==2)) continue;
	
      auto _part_ = (*genParticles_)[p];
      //      if(verbose_) std::cout << "[JpsiMuNtuplizer] Hammer: " << _part_.pdgId() << " (" << _part_.status() << ")" << std::endl;
      
      bool isJpsi = false;
      bool isTau = false;

      
      // check if the Bc candidate has both J/psi and tau
      for(auto d : _part_.daughterRefVector()) {
	if(TMath::Abs(d->pdgId()) == 443) isJpsi = true;
	if(TMath::Abs(d->pdgId()) == 15 || TMath::Abs(d->pdgId())==13) isTau = true;
      }

      //      std::cout << isJpsi << " " << isTau  << " " << (isJpsi==true && isTau==true) << std::endl;
      if(!(isJpsi==true && isTau==true)){
	//	if(verbose_) std::cout << "[JpsiMuNtuplizer] Hammer: This B is rejected !!!!" << std::endl;
	continue;
      }
      
      
      Hammer::Particle pB({_part_.energy(), _part_.px(), _part_.py(), _part_.pz()}, _part_.pdgId());
      
      idxB = Bc2JpsiLNu.addParticle(pB);
      
      if(verbose_) std::cout << "[JpsiMuNtuplizer] Hammer: " << (*genParticles_)[p].pdgId() << " " << (*genParticles_)[p].status() << std::endl;
      
      
      for(auto d : _part_.daughterRefVector()) {
	
	Hammer::Particle B_dau({d->energy(), d->px(), d->py(), d->pz()}, d->pdgId());
	
	auto idx_d = Bc2JpsiLNu.addParticle(B_dau);
	Bvtx_idxs.push_back(idx_d);
	
	if(verbose_) std::cout << "[JpsiMuNtuplizer] Hammer: \t gen: " << d->pdgId() << " " << d->status() << std::endl;	  
	
	if(TMath::Abs(d->pdgId()) == 15) {
	  
	  idxTau = idx_d;
	 
	  for (auto dd : d->daughterRefVector()) {
	    Hammer::Particle Tau_dau({dd->energy(), dd->px(), dd->py(), dd->pz()}, dd->pdgId());
	    auto idx_dd = Bc2JpsiLNu.addParticle(Tau_dau);
	    Tauvtx_idxs.push_back(idx_dd);
	    if(verbose_) std::cout << "[JpsiMuNtuplizer] Hammer: \t\t gen: " << dd->pdgId() << " " << dd->status() << std::endl;
	  }
	}
	
	else if(TMath::Abs(d->pdgId()) == 443) {
	  
	  idxJpsi = idx_d;
	  
	  for (auto dd : d->daughterRefVector()) {
	    Hammer::Particle Jpsi_dau({dd->energy(), dd->px(), dd->py(), dd->pz()}, dd->pdgId());
	    auto idx_dd = Bc2JpsiLNu.addParticle(Jpsi_dau);
	    Jpsivtx_idxs.push_back(idx_dd);
	    if(verbose_) std::cout << "[JpsiMuNtuplizer] Hammer: \t\t gen: " << dd->pdgId() << " " << dd->status() << std::endl;
	  }

	}
      }
    }

    if(verbose_) std::cout << "[JpsiMuNtuplizer] Hammer idx (B, tau, Jpsi) = " << idxB << " " << idxTau << " " << idxJpsi << std::endl;
      
    Bc2JpsiLNu.addVertex(idxB, Bvtx_idxs);
    
    if(idxTau != -1) {
      Bc2JpsiLNu.addVertex(idxTau, Tauvtx_idxs);
    }
    if(idxJpsi != -1) {
      Bc2JpsiLNu.addVertex(idxJpsi, Jpsivtx_idxs);
    }
	
    hammer.addProcess(Bc2JpsiLNu);
    hammer.processEvent();

    std::map<std::string, double> settings;
    for(auto n: parName) settings["delta_" + n] = 0;
    //    for(auto pars: _FFErrNames) {
    //      settings[pars] = 0;
    //    }

    hammer.setFFEigenvectors("BctoJpsi", "BGLVar", settings);
	
    auto weights = hammer.getWeights("Scheme1");

    Float_t weight = -1;

    if(!weights.empty()){
      for(auto elem: weights) {
	if(isnan(elem.second)) {
	  std::cout << "[JpsiMuNtuplizer] ERROR: BGL Central nan weight: " << elem.second << std::endl;
	}else{
	  weight = elem.second;
	}
      }
    }

    nBranches_->JpsiMu_hammer_ebe.push_back(weight);

    //    std::cout << "-----------------------" << std::endl;
    //    std::cout << "base weight = " << weight << std::endl;


    ///////////////////////////////////////////////////////////////////////
    // MC 
    ///////////////////////////////////////////////////////////////////////

    std::vector<float> hweights; 

    for(int imc=0; imc < numberofToys; imc++){
      hammer.setFFEigenvectors("BctoJpsi", "BGLVar", FFdict[imc]);
      auto weights = hammer.getWeights("Scheme1");
      Float_t weight_sys = -1;

      if(!weights.empty()){
	for(auto elem: weights) {
	  if(isnan(elem.second)) {
	    std::cout << "[ERROR]: BGL nan weight: " << elem.second << std::endl;
	  }else{
	    weight_sys = elem.second;
	  }
	}
      }

      if(flag_fill==false){
	std::vector<float> settings_ff;
	for(auto pars: _FFErrNames) {
	  settings_ff.push_back(FFdict[imc][pars]);
	}
      
	nBranches_->JpsiMu_hammer_ff.push_back(settings_ff);      
      }

      hweights.push_back(weight_sys);
    }

    flag_fill = true;


    nBranches_->JpsiMu_hammer_ebe_toy.push_back(hweights);


    //////////////////////////
    // ordinary method 
    //////////////////////////

    for(int i=0; i<11; i++) { //Loop over eigenVar
      for (int j=0; j<2; j++) { //Loop over pos/neg direction

        map<string, double> settings;

        for (int k=0; k<11; k++) { //Loop over parameters
          settings["delta_" + parName[k]] = eigVar[i][k][j];
        }

        hammer.setFFEigenvectors("BctoJpsi", "BGLVar", settings);
        auto weights = hammer.getWeights("Scheme1");
        string var_name = "eig";
	var_name += std::to_string(i);
        var_name += j==0? "_Up" : "_Down";



	Float_t weight_sys = -1;

	if(!weights.empty()){
	  for(auto elem: weights) {
	    if(isnan(elem.second)) {
	      std::cout << "[ERROR]: BGL nan weight: " << elem.second << std::endl;
	    }else{
	      weight_sys = elem.second;
	    }
	  }
	}

	if(var_name==std::string("eig0_Up")) nBranches_->JpsiMu_hammer_ebe_a0_up.push_back(weight_sys);
	else if(var_name==std::string("eig0_Down")) nBranches_->JpsiMu_hammer_ebe_a0_down.push_back(weight_sys);
	else if(var_name==std::string("eig1_Up")) nBranches_->JpsiMu_hammer_ebe_a1_up.push_back(weight_sys);
	else if(var_name==std::string("eig1_Down")) nBranches_->JpsiMu_hammer_ebe_a1_down.push_back(weight_sys);
	else if(var_name==std::string("eig2_Up")) nBranches_->JpsiMu_hammer_ebe_a2_up.push_back(weight_sys);
	else if(var_name==std::string("eig2_Down")) nBranches_->JpsiMu_hammer_ebe_a2_down.push_back(weight_sys);

	else if(var_name==std::string("eig3_Up")) nBranches_->JpsiMu_hammer_ebe_b0_up.push_back(weight_sys);
	else if(var_name==std::string("eig3_Down")) nBranches_->JpsiMu_hammer_ebe_b0_down.push_back(weight_sys);
	else if(var_name==std::string("eig4_Up")) nBranches_->JpsiMu_hammer_ebe_b1_up.push_back(weight_sys);
	else if(var_name==std::string("eig4_Down")) nBranches_->JpsiMu_hammer_ebe_b1_down.push_back(weight_sys);
	else if(var_name==std::string("eig5_Up")) nBranches_->JpsiMu_hammer_ebe_b2_up.push_back(weight_sys);
	else if(var_name==std::string("eig5_Down")) nBranches_->JpsiMu_hammer_ebe_b2_down.push_back(weight_sys);

	else if(var_name==std::string("eig6_Up")) nBranches_->JpsiMu_hammer_ebe_c1_up.push_back(weight_sys);
	else if(var_name==std::string("eig6_Down")) nBranches_->JpsiMu_hammer_ebe_c1_down.push_back(weight_sys);
	else if(var_name==std::string("eig7_Up")) nBranches_->JpsiMu_hammer_ebe_c2_up.push_back(weight_sys);
	else if(var_name==std::string("eig7_Down")) nBranches_->JpsiMu_hammer_ebe_c2_down.push_back(weight_sys);

	else if(var_name==std::string("eig8_Up")) nBranches_->JpsiMu_hammer_ebe_d0_up.push_back(weight_sys);
	else if(var_name==std::string("eig8_Down")) nBranches_->JpsiMu_hammer_ebe_d0_down.push_back(weight_sys);
	else if(var_name==std::string("eig9_Up")) nBranches_->JpsiMu_hammer_ebe_d1_up.push_back(weight_sys);
	else if(var_name==std::string("eig9_Down")) nBranches_->JpsiMu_hammer_ebe_d1_down.push_back(weight_sys);
	else if(var_name==std::string("eig10_Up")) nBranches_->JpsiMu_hammer_ebe_d2_up.push_back(weight_sys);
	else if(var_name==std::string("eig10_Down")) nBranches_->JpsiMu_hammer_ebe_d2_down.push_back(weight_sys);
	
	//	std::cout << "test1 \t " << var_name << " " << weight_sys << std::endl;
      }
    }

    //    nBranches_->JpsiTau_hammer_ebe_up.push_back(TMath::Sqrt(hammer_up));
    //    nBranches_->JpsiTau_hammer_ebe_down.push_back(TMath::Sqrt(hammer_down));

  }


	
  nBranches_->JpsiMu_q2_gen.push_back(q2_gen.M2());
  nBranches_->JpsiMu_B_pt_gen.push_back(pB_gen.Pt());
  nBranches_->JpsiMu_B_eta_gen.push_back(pB_gen.Eta());
  nBranches_->JpsiMu_B_phi_gen.push_back(pB_gen.Phi());
  nBranches_->JpsiMu_B_mass_gen.push_back(pB_gen.M());



  // -9 if there is no Bc found 
  nBranches_->JpsiMu_genPV_vx.push_back(genvertex.x());
  nBranches_->JpsiMu_genPV_vy.push_back(genvertex.y());
  nBranches_->JpsiMu_genPV_vz.push_back(genvertex.z());
  nBranches_->JpsiMu_ngenmuons.push_back(gen_nr_mu.size() + gen_jpsi_mu.size());

  return true;
}


