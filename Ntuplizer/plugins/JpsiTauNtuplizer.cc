#include "../interface/JpsiTauNtuplizer.h"


//===================================================================================================================
JpsiTauNtuplizer::JpsiTauNtuplizer( edm::EDGetTokenT<pat::MuonCollection>    muonToken   ,
				    edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
				    edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken,
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
  , HLTtriggersToken_	( triggertoken )
  , triggerObjects_	( triggerobject )
  , genParticlesToken_( genptoken )
  , genTauToken_( genttoken )
  , runOnMC_   (runFlags["runOnMC"])
  , useDNN_   (runFlags["useDNN"])
  , useHammer_   (runFlags["useHammer"])
  , isTruth_   (runFlags["isTruth"])
  , verbose_   (runFlags["verbose"])
  , c_dz (runValues["dzcut"])
  , c_fsig (runValues["fsigcut"])
  , c_vprob (runValues["vprobcut"])
  , c_dnn (runValues["dnncut"])
  , c_charge (runValues["tau_charge"])
  , dnnfile_ (runStrings["dnnfile"])      
   
{

  if(verbose_){
    std::cout << "[JpsiTauNtuplizer] runOnMC    = " << runOnMC_ << std::endl;
    std::cout << "[JpsiTauNtuplizer] UseDNN     = " << useDNN_ << std::endl;
    std::cout << "[JpsiTauNtuplizer] UseHammer  = " << useHammer_ << std::endl;
    std::cout << "[JpsiTauNtuplizer] isTruth    = " << isTruth_ << std::endl;
    std::cout << "[JpsiTauNtuplizer] dzcut      = " << c_dz << std::endl;
    std::cout << "[JpsiTauNtuplizer] fsigcut    = " << c_fsig << std::endl;
    std::cout << "[JpsiTauNtuplizer] vprob      = " << c_vprob << std::endl;
    std::cout << "[JpsiTauNtuplizer] dnn cut    = " << c_dnn << std::endl;
    std::cout << "[JpsiTauNtuplizer] tau charge = " << c_charge << std::endl;
  }


  if(useDNN_){
    
    std::string dnnfilepath = edm::FileInPath("EXOVVNtuplizerRunII/Ntuplizer/" +  dnnfile_).fullPath();

    if(verbose_) std::cout << "[JpsiTauNtuplizer] DNN file   = " << dnnfilepath << std::endl;

    std::string tbr = "DUMMY";  // to be replaced
    auto pos = dnnfilepath.find(tbr);
    auto len = tbr.length();
    if (pos != std::string::npos) {
      dnnfilepath.replace(pos, len, "");
    }
    
    graphDef = tensorflow::loadMetaGraph(dnnfilepath);
    session = tensorflow::createSession(graphDef, dnnfilepath);
    
    data = tensorflow::Tensor(tensorflow::DT_FLOAT, { 1, 50, 8 }); // single batch of dimension 10
    label = tensorflow::Tensor(tensorflow::DT_INT32, { 1,50}); 
    add_global = tensorflow::Tensor(tensorflow::DT_FLOAT, { 1, 2 }); 
    isTraining = tensorflow::Tensor(tensorflow::DT_BOOL, tensorflow::TensorShape()); 
    norm = tensorflow::Tensor(tensorflow::DT_FLOAT, { 1, 50 }); 
    
  }
  

  if(runOnMC_ && useHammer_){

    if(verbose_) std::cout << "[JpsiTauNtuplizer] Setting up Hammer" << std::endl;

    hammer.setUnits("GeV");

    std::vector<std::string> processes = {"BcJpsiMuNu", "BcJpsiTauNu"};
  
    for(auto proc : processes) {
      if(verbose_) std::cout << "[JpsiTauNtuplizer] \t Hammer added: " << proc << std::endl;
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

    if(verbose_) std::cout << "[JpsiTauNtuplizer] Finish setting up Hammer" << std::endl;

    // add abcdmatrix to be the principal components

    // from https://arxiv.org/pdf/1909.10691.pdf
    
    vector<vector<double>> abcdmat{ 
      {-0.0168285, -0.346957, -0.561914, -0.200744, -0.0357712, 0.488216, -0.490298, -0.112359, 0.172, -0.032536, 0},
      {-0.273517, 0.260029, 0.424574, 0.223675, 0.150618, 0.445983, -0.455354, -0.0707909, -0.438665, -0.0301344, 0},
      {0.32971, 0.294397, 0.343693, 0.137184, -0.0550957, 0.12977, -0.237882, -0.142826, 0.755413, -0.0128438, 0.},
      {-0.117461, 0.394041, -0.424535, 0.3952, -0.222928, -0.418512, -0.37727, 0.126589, -0.0227784, -0.154801, -0.294672},
      {-0.381043, -0.199939, 0.184817, -0.198149, 0.210655, -0.413679, -0.424446, 0.30485, 0.213187, 0.391218, 0.233061},
      {0.370329, -0.173951, 0.200939, -0.341875, 0.038249, -0.416914, -0.364073, -0.451958, -0.253803, -0.290139, -0.121033},
      {-0.37948, 0.252966, -0.0779388, -0.294611, 0.294091, 0.00507043, 0.161783, -0.42904, 0.150762, 0.255879, -0.560504},
      {0.295611, 0.345131, -0.0289752, -0.497994, 0.248978, 0.129983, -0.0623405, 0.619602, -0.0690973, -0.149173, -0.225734},
      {0.173793, 0.527999, -0.299466, -0.180232, 0.0603336, -0.0444077, -0.0241702, -0.274265, -0.154552, 0.325246, 0.598131},
      {-0.433734, 0.131009, -0.0426141, -0.14, 0.213242, -0.0533621, 0.0997566, -0.0717846, 0.213096, -0.736981, 0.348311},
      {0.265219, -0.144152, -0.180805, 0.437778, 0.824543, -0.0637481, 0.00642872, -0.0196135, 0.0125086, -0.00453039, 0.000634603}};


//    std::map<std::string, std::vector<double>> paramsBGL {
//      {"avec", {0.004698, -0.02218, 0.1503}},
//	{"bvec", {0.003424, -0.02595, 0.3897}},
//          {"cvec", {-0.003164, 0.08731}},
//	    {"dvec", {0.04011, -0.2132, 0.008157}}
//    };
    std::string centralValuesOpt = "\"BctoJpsiBGL: {abcdmatrix: [";
    centralValuesOpt += "[-0.0168285, -0.346957, -0.561914, -0.200744, -0.0357712, 0.488216, -0.490298, -0.112359, 0.172, -0.032536, 0],";
    centralValuesOpt += "[-0.273517, 0.260029, 0.424574, 0.223675, 0.150618, 0.445983, -0.455354, -0.0707909, -0.438665, -0.0301344, 0],";
    centralValuesOpt += "[0.32971, 0.294397, 0.343693, 0.137184, -0.0550957, 0.12977, -0.237882, -0.142826, 0.755413, -0.0128438, 0.],";
    centralValuesOpt += "[-0.117461, 0.394041, -0.424535, 0.3952, -0.222928, -0.418512, -0.37727, 0.126589, -0.0227784, -0.154801, -0.294672],";
    centralValuesOpt += "[-0.381043, -0.199939, 0.184817, -0.198149, 0.210655, -0.413679, -0.424446, 0.30485, 0.213187, 0.391218, 0.233061],";
    centralValuesOpt += "[0.370329, -0.173951, 0.200939, -0.341875, 0.038249, -0.416914, -0.364073, -0.451958, -0.253803, -0.290139, -0.121033],";
    centralValuesOpt += "[-0.37948, 0.252966, -0.0779388, -0.294611, 0.294091, 0.00507043, 0.161783, -0.42904, 0.150762, 0.255879, -0.560504],";
    centralValuesOpt += "[0.295611, 0.345131, -0.0289752, -0.497994, 0.248978, 0.129983, -0.0623405, 0.619602, -0.0690973, -0.149173, -0.225734],";
    centralValuesOpt += "[0.173793, 0.527999, -0.299466, -0.180232, 0.0603336, -0.0444077, -0.0241702, -0.274265, -0.154552, 0.325246, 0.598131],";
    centralValuesOpt += "[-0.433734, 0.131009, -0.0426141, -0.14, 0.213242, -0.0533621, 0.0997566, -0.0717846, 0.213096, -0.736981, 0.348311],";
    centralValuesOpt += "[0.265219, -0.144152, -0.180805, 0.437778, 0.824543, -0.0637481, 0.00642872, -0.0196135, 0.0125086, -0.00453039, 0.000634603]]}\"";
    
	  
	
//        for(auto elem: paramsBGL) {
//          centralValuesOpt += Form("%s: {", elem.first.c_str());
//          
//          for(size_t ii=0; ii < elem.second.size(); ii++){
//    	if(ii==elem.second.size()-1){
//    	  centralValuesOpt += Form("%f ", elem.second[ii]);
//    	}else{
//    	  centralValuesOpt += Form("%f, ", elem.second[ii]);
//    	}
//    
//          }
//          if(elem.first.c_str()==std::string("dvec")) centralValuesOpt += "}";
//          else centralValuesOpt += "}, ";
//        }
//        centralValuesOpt += "}\"";
        std::cout << "[Hammer]: Central values\n\t" << centralValuesOpt << std::endl;
        hammer.setOptions(centralValuesOpt);
        std::cout << "... finishes " << std::endl;
    



  }
}




//===================================================================================================================
JpsiTauNtuplizer::~JpsiTauNtuplizer( void )
{

  
  tensorflow::closeSession(session);
  delete graphDef;


  if(runOnMC_ && useHammer_){
  
    if(verbose_) std::cout <<"[JpsiTauNtuplizer] Evaluate partial width" << std::endl;
  
    std::vector<std::string> processes = {"BcJpsiTau+Nu"};
    
    for(auto proc : processes) {
      std::map<std::string, double> outRate;
      if(verbose_) std::cout << "\t Process: " << proc << std::endl;
      
      outRate["den"] = hammer.getDenominatorRate(proc);
      
      if(outRate["den"] == 0) {
	std::cout <<"[JpsiTauNtuplizer]: ERROR Failed to get default partial width" << std::endl;
	continue;
      }else{
	if(verbose_) std::cout << Form("[JpsiTauNtuplizer] Default rate: %1.3e", outRate["den"]) << std::endl;
      }
      
      std::map<std::string, double> settings;
      
      for(auto pars: _FFErrNames) {
	settings[pars] = 0;
      }
      
      hammer.setFFEigenvectors("BctoJpsi", "BGLVar", settings);
      
      outRate["Central"] = hammer.getRate(proc, "Scheme1");
      std::cout << Form("[JpsiTauNtuplizer] New rate: %1.3e (ratio = %.3f)", outRate["Central"], outRate["Central"]/outRate["den"]) << std::endl;
      
      nBranches_->hammer_width->SetBinContent(1, outRate["den"]);
      nBranches_->hammer_width->SetBinContent(2, outRate["Central"]);
    
    
      int idx1 = 0;
      
      for(auto pars1: _FFErrNames) {
	
	for(int isUp=0; isUp < 2; isUp++) { // up, down variation    
	  
	  std::map<std::string, double> settings;
	  
	  int idx2 = 0;
	  
	  for(auto pars2: _FFErrNames) {
	    
	    if(idx1 == idx2){
	      if(isUp==1) settings[pars2] = _FFErr[idx2];
	      else settings[pars2] = -_FFErr[idx2];
	    }else{
	      settings[pars2] = 0;
	    }
	    idx2 += 1;
	  }
	  
	  
	  hammer.setFFEigenvectors("BctoJpsi", "BGLVar", settings);
	  
	  auto rate = hammer.getRate(proc, "Scheme1");
	  
	  std::string var_name = pars1;
	  var_name += isUp==1? "_Up" : "_Down";
	  
	  if(verbose_) std::cout << "[JpsiTauNtuplizer] " << var_name << " --> " << Form(": %1.3e", rate) << std::endl;
	  
	  if(var_name==std::string("delta_a0_Up")) nBranches_->hammer_width->SetBinContent(3, rate);
	  else if(var_name==std::string("delta_a0_Down")) nBranches_->hammer_width->SetBinContent(4, rate);
	  else if(var_name==std::string("delta_a1_Up")) nBranches_->hammer_width->SetBinContent(5, rate);
	  else if(var_name==std::string("delta_a1_Down")) nBranches_->hammer_width->SetBinContent(6, rate);
	  else if(var_name==std::string("delta_a2_Up")) nBranches_->hammer_width->SetBinContent(7, rate);
	  else if(var_name==std::string("delta_a2_Down")) nBranches_->hammer_width->SetBinContent(8, rate);
	  
	  else if(var_name==std::string("delta_b0_Up")) nBranches_->hammer_width->SetBinContent(9, rate);
	  else if(var_name==std::string("delta_b0_Down")) nBranches_->hammer_width->SetBinContent(10, rate);
	  else if(var_name==std::string("delta_b1_Up")) nBranches_->hammer_width->SetBinContent(11, rate);
	  else if(var_name==std::string("delta_b1_Down")) nBranches_->hammer_width->SetBinContent(12, rate);
	  else if(var_name==std::string("delta_b2_Up")) nBranches_->hammer_width->SetBinContent(13, rate);
	  else if(var_name==std::string("delta_b2_Down")) nBranches_->hammer_width->SetBinContent(14, rate);
	  
	  else if(var_name==std::string("delta_c1_Up")) nBranches_->hammer_width->SetBinContent(15, rate);
	  else if(var_name==std::string("delta_c1_Down")) nBranches_->hammer_width->SetBinContent(16, rate);
	  else if(var_name==std::string("delta_c2_Up")) nBranches_->hammer_width->SetBinContent(17, rate);
	  else if(var_name==std::string("delta_c2_Down")) nBranches_->hammer_width->SetBinContent(18, rate);
	  
	  else if(var_name==std::string("delta_d0_Up")) nBranches_->hammer_width->SetBinContent(19, rate);
	  else if(var_name==std::string("delta_d0_Down")) nBranches_->hammer_width->SetBinContent(20, rate);
	  else if(var_name==std::string("delta_d1_Up")) nBranches_->hammer_width->SetBinContent(21, rate);
	  else if(var_name==std::string("delta_d1_Down")) nBranches_->hammer_width->SetBinContent(22, rate);
	  else if(var_name==std::string("delta_d2_Up")) nBranches_->hammer_width->SetBinContent(23, rate);
	  else if(var_name==std::string("delta_d2_Down")) nBranches_->hammer_width->SetBinContent(24, rate);
	  
	  
	}
	idx1 += 1;
	
      }
    }
  }
}



bool JpsiTauNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  if(verbose_) std::cout << "[JpsiTauNtuplizer] ---------------- event, run, lumi = " << event.id().event() << " " << event.id().run() << " " << event.id().luminosityBlock() << "----------------" << std::endl;
  

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
	  //	  std::cout << "J/psi status = " << d->status() << std::endl;
	  isJpsi = true;
	}
	if(TMath::Abs(d->pdgId()) == 15){
	  //	  std::cout << "Tau status = " << d->status() << std::endl;
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
      }
    }
  }

  TLorentzVector q2_gen; 
  q2_gen = pB_gen - pJpsi_gen;

  //std::cout << "test:" << q2_gen.Pt() << " " << q2_gen.M2() << " " << pB_gen.Pt() << " " << pJpsi_gen.Pt() << std::endl;
  nBranches_->q2_nocut->Fill(q2_gen.M2());
  

  
  /********************************************************************
   *
   * Step1: Check if the J/psi trigger is fired.
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

  if(verbose_) std::cout << "[JpsiTauNtuplizer] Trigger fired" << std::endl;

  /********************************************************************
   *
   * Step2: pre-select muons for building J/psi candidates ... 
   * For muons, no requirement applied (soft muon ID is required once ther vertex is chosen)
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

  if(verbose_) std::cout << "[JpsiTauNtuplizer] At least 2 muons are found" << std::endl;

  /********************************************************************
   *
   * Step3: building J/psi candidates 
   *
   ********************************************************************/

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

  Float_t jpsi_max_pt = -1;
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
      
      std::vector<reco::TransientTrack> transient_tracks_dimuon;
      
      transient_tracks_dimuon.push_back((*builder).build(muoncollection[imu].muonBestTrack()));
      transient_tracks_dimuon.push_back((*builder).build(muoncollection[jmu].muonBestTrack()));
      
      Float_t vprob_jpsi = -9;
      TransientVertex vertex_jpsi;
      std::tie(vprob_jpsi, vertex_jpsi) = aux.vertexProb(transient_tracks_dimuon);

      if(jpsi_max_pt < jpsi_pt){
	jpsi_max_pt = jpsi_pt;
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

  if(verbose_) std::cout << "[JpsiTauNtuplizer] J/psi found" << std::endl;

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
   * I don't think it is necessary to complicate things at this stage 
   * (we can do the refit once we select B candidate)
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

  if(verbose_) std::cout << "[JpsiTauNtuplizer] J/psi soft muon ID passed" << std::endl;

  particle_cand JPcand;
  JPcand = aux.calculateIPvariables(extrapolator, jpsi_part, jpsi_vertex, closestVertex);

  /********************************************************************
   *
   * Step6: Tau selection
   *        Just select highest in pT but there might be better selection ... 
   *
   ********************************************************************/
    
  event.getByToken( packedpfcandidatesToken_               , packedpfcandidates_      ); 
    
  std::vector<pat::PackedCandidate> pfcollection; 
  std::vector<reco::TransientTrack> mytracks;
  std::vector<Float_t> mydnn;
    
  Int_t npf_before_dnn = 0;
  Int_t npf_qr = 0;

  if(useDNN_){

    std::vector<pfcand> pfcands;

    Int_t count_dnn = 0;
    Int_t count_dnn_muon = 0;

    for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){   
	
      pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
	
      if(TMath::Abs(pf.eta()) < 2.4 && 
	 TMath::Abs(pf.charge())==1 &&
	 pf.pt() > 4. &&
	 pf.isGlobalMuon() > 0.5 &&
	 pf.hasTrackDetails() > 0.5
	 ){


	if(count_dnn < 50){
	  data.tensor<float, 3>()(0, count_dnn, 0) = pf.eta();
	  data.tensor<float, 3>()(0, count_dnn, 1) = pf.phi();
	  data.tensor<float, 3>()(0, count_dnn, 2) = TMath::Log(pf.pt());
	  data.tensor<float, 3>()(0, count_dnn, 3) = TMath::Log(pf.energy());
	  data.tensor<float, 3>()(0, count_dnn, 4) = pf.charge();
	  data.tensor<float, 3>()(0, count_dnn, 5) = TMath::Abs(closestVertex.position().z() - pf.pseudoTrack().vz());
	  data.tensor<float, 3>()(0, count_dnn, 6) = TMath::Sqrt( TMath::Power((closestVertex.position().x() - pf.pseudoTrack().vx()), 2) + TMath::Power((closestVertex.position().y() - pf.pseudoTrack().vy()), 2));
	  data.tensor<float, 3>()(0, count_dnn, 7) = pf.isGlobalMuon();
	    
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
	
      Bool_t hpflag = pf.trackHighPurity();
      if(!hpflag) continue;
      if(pf.pseudoTrack().hitPattern().numberOfValidPixelHits() < 0) continue;
      if(pf.pseudoTrack().hitPattern().numberOfValidHits() < 3) continue;
      if(pf.pseudoTrack().normalizedChi2() > 100) continue;
	
      if(TMath::Abs(pf.pdgId())!=211) continue; 
      if(TMath::Abs(pf.eta()) > 2.5) continue; 

      npf_before_dnn++;	

      pfcand _cand_ = {
	(Int_t)ii,
	(Float_t) abs(precut_dz)
      };
	  
      pfcands.push_back(_cand_);
    }

    //sorting by distance to the vertex
    sort(pfcands.begin(), pfcands.end());


    for(size_t ic = 0; ic < pfcands.size(); ic++){
      Int_t idx = pfcands[ic].cand_idx;

      pat::PackedCandidate pf = (*packedpfcandidates_)[idx];

      if(count_dnn < 50){
	data.tensor<float, 3>()(0, count_dnn, 0) = pf.eta();
	data.tensor<float, 3>()(0, count_dnn, 1) = pf.phi();
	data.tensor<float, 3>()(0, count_dnn, 2) = TMath::Log(pf.pt());
	data.tensor<float, 3>()(0, count_dnn, 3) = TMath::Log(pf.energy());
	data.tensor<float, 3>()(0, count_dnn, 4) = pf.charge();
	data.tensor<float, 3>()(0, count_dnn, 5) = TMath::Abs(closestVertex.position().z() - pf.pseudoTrack().vz());
	data.tensor<float, 3>()(0, count_dnn, 6) = TMath::Sqrt( TMath::Power((closestVertex.position().x() - pf.pseudoTrack().vx()), 2) + TMath::Power((closestVertex.position().y() - pf.pseudoTrack().vy()), 2));
	data.tensor<float, 3>()(0, count_dnn, 7) = pf.isGlobalMuon();
	  
	label.matrix<int>()(0, count_dnn) = 0;
	norm.matrix<float>()(0, count_dnn) = float(1);
	count_dnn++;

	pfcollection.push_back(pf);
	reco::TransientTrack  tt_track = (*builder).build(pf.pseudoTrack());
	mytracks.push_back(tt_track);

      }
    }

      
    for(int ic=count_dnn; ic<50; ic++){

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

    for(int ic=count_dnn_muon; ic<count_dnn; ic++){

      mydnn.push_back(finalOutputTensor(0, ic, 1));
    }


  }else{


    for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){   
      
      pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
	
      if(pf.pt() < 0.5) continue;
      if(!pf.hasTrackDetails()) continue;
	
      // use the PF candidates that come from closestVertex
      //      if(pf.vertexRef()->z()!=closestVertex.position().z()) continue;
	
      Float_t precut_dz = pf.vz() - closestVertex.position().z();
      if(TMath::Abs(precut_dz) > c_dz) continue;
	
      Bool_t hpflag = pf.trackHighPurity();
      if(!hpflag) continue;
      if(pf.pseudoTrack().hitPattern().numberOfValidPixelHits() < 0) continue;
      if(pf.pseudoTrack().hitPattern().numberOfValidHits() < 3) continue;
      if(pf.pseudoTrack().normalizedChi2() > 100) continue;
	
      if(TMath::Abs(pf.pdgId())!=211) continue; 
      if(TMath::Abs(pf.eta()) > 2.5) continue; 

      pfcollection.push_back(pf);
      reco::TransientTrack  tt_track = (*builder).build(pf.pseudoTrack());
      mytracks.push_back(tt_track);
	
    }
  }



  Int_t numOfch = (size_t)pfcollection.size();


  std::vector<std::vector<TLorentzVector>> gps;
  std::vector<Int_t> ppdgId;
  std::vector<Int_t> vec_gentaudm;
  std::vector<Int_t> vec_ppdgId;
  std::vector<TLorentzVector> vec_gentaup4;
  std::vector<TLorentzVector> vec_gentaup4_vis;
  std::vector<TLorentzVector> vec_gentau3pp4;
  Int_t isgen3 = 0;
  Int_t isgen3matched = 0;

  if(runOnMC_){
    event.getByToken(genTauToken_, genTaus_);
    
    for( unsigned p=0; p < genParticles_->size(); ++p){
      
      if(TMath::Abs((*genParticles_)[p].pdgId())!=15) continue;
      if(TMath::Abs((*genParticles_)[p].status())!=2) continue;
      
      if(verbose_) std::cout << "[JpsiTauNtuplizer] Tau found with # of daughters = " << (*genParticles_)[p].numberOfDaughters() << " with mother = " << (*genParticles_)[p].mother(0)->pdgId() << std::endl;
      
      
      TLorentzVector genvis;
      TLorentzVector genvis_full;

      std::vector<TLorentzVector> gp;
      Bool_t matched = true;
      Int_t nprong = 0;
      
      for(int idd = 0; idd < (int)(*genParticles_)[p].numberOfDaughters(); idd++){
	
	if(verbose_){ std::cout << "[JpsiTauNtuplizer] \t -> " << (*genParticles_)[p].daughter(idd)->pdgId() << " (pT, eta, phi) = " 
			       << (*genParticles_)[p].daughter(idd)->pt() << " " 
			       << (*genParticles_)[p].daughter(idd)->eta() << " " 
			       << (*genParticles_)[p].daughter(idd)->phi() << std::endl;
	}


	TLorentzVector _genvis_full;
	_genvis_full.SetPtEtaPhiM((*genParticles_)[p].daughter(idd)->pt(),
				  (*genParticles_)[p].daughter(idd)->eta(),
				  (*genParticles_)[p].daughter(idd)->phi(),
				  (*genParticles_)[p].daughter(idd)->mass()
				  );
	
	genvis_full += _genvis_full;

	if(
	   TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==12 ||
	   TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==14 || 
	   TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==16
	   ){

	  continue;
	}


	TLorentzVector _genvis_;
	_genvis_.SetPtEtaPhiM((*genParticles_)[p].daughter(idd)->pt(),
			      (*genParticles_)[p].daughter(idd)->eta(),
			      (*genParticles_)[p].daughter(idd)->phi(),
			      (*genParticles_)[p].daughter(idd)->mass()
			      );
	
	genvis += _genvis_;

	if(TMath::Abs((*genParticles_)[p].daughter(idd)->pdgId())==211){

	  nprong += 1;
	  
	  Float_t min_dr = 999;

	  for(int kkk = 0; kkk < numOfch; kkk ++){
	    
	    pat::PackedCandidate _pf = pfcollection[kkk];
	    
	    if(_pf.pdgId()!=(*genParticles_)[p].daughter(idd)->pdgId()) continue;
	    
	    Float_t _dR = reco::deltaR(
				       _genvis_.Eta(), _genvis_.Phi(),
				       _pf.eta(), _pf.phi()
				       );

	    if(_dR < min_dr && _dR < 0.015 && _pf.pt()/_genvis_.Pt() < 1.15 && _pf.pt()/_genvis_.Pt() > 0.85){
	      min_dr = _dR;
	    }
	  }

	  if(min_dr == 999) matched = false;

	  gp.push_back(_genvis_);

	}
      }


      if(nprong==3) isgen3 += 1;


      // check decay mod. To do this, take matching with tau-genjet. 

      Float_t min_gendr = 999;
      Int_t taugendm = -999;
      //      Float_t taugenvis = -999;

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
	  //	  taugenvis = visibleP4.pt();
	  
	  if(verbose_) std::cout << "[JpsiTauNtuplizer] \t matched gen decay mode = " << JetMCTagUtils::genTauDecayMode(TauCand) << " (" <<  taugendm << ")" << std::endl;

	}
      }
      
      vec_ppdgId.push_back((*genParticles_)[p].mother(0)->pdgId());
      vec_gentaudm.push_back(taugendm);
      vec_gentaup4_vis.push_back(genvis);
      vec_gentaup4.push_back(genvis_full);
      
      //      std::cout << "check: " << genvis.Pt() << " " << taugenvis << std::endl;
      
      if(gp.size()==3){
	if(verbose_) std::cout << "[JpsiTauNtuplizer] \t 3prong found with mother = " << (*genParticles_)[p].mother(0)->pdgId() << std::endl;
	gps.push_back(gp);
	ppdgId.push_back((*genParticles_)[p].mother(0)->pdgId());
	vec_gentau3pp4.push_back(genvis);
	
	if(matched) isgen3matched += 1;

      }
    }
  }


  ///////////////////////////////
  // Event display 
  ///////////////////////////////
  
  //  std::cout <<"----------------------------------------------------------" << std::endl;

//  if(runOnMC_ > 0.5 && vec_gentaudm.size()==1 && isgen3matched==1){
//
//    Int_t counter_npf = 0;
//    for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){   
//      
//      pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
//      
//      if(pf.pt() < 0.5) continue;
//      if(!pf.hasTrackDetails()) continue;
//	
//      // use the PF candidates that come from closestVertex
//      //      if(pf.vertexRef()->z()!=closestVertex.position().z()) continue;
//	
//      //      Float_t precut_dz = pf.vz() - closestVertex.position().z();
//      //      if(TMath::Abs(precut_dz) > c_dz) continue;
//	
//      Bool_t hpflag = pf.trackHighPurity();
//      if(!hpflag) continue;
//      if(pf.pseudoTrack().hitPattern().numberOfValidPixelHits() < 0) continue;
//      if(pf.pseudoTrack().hitPattern().numberOfValidHits() < 3) continue;
//      if(pf.pseudoTrack().normalizedChi2() > 100) continue;
//	
//      if(TMath::Abs(pf.pdgId())!=211) continue; 
//      if(TMath::Abs(pf.eta()) > 2.5) continue; 
//
//      counter_npf ++;
//      nBranches_->JpsiTau_ed_pfeta.push_back(pf.eta());
//      nBranches_->JpsiTau_ed_pfphi.push_back(pf.phi());
//      nBranches_->JpsiTau_ed_isRight.push_back(-9);
//      nBranches_->JpsiTau_ed_pfdnn.push_back(-9);
//      nBranches_->JpsiTau_ed_genpt.push_back(vec_gentaup4_vis[0].Pt());
//      nBranches_->JpsiTau_ed_id.push_back(event.id().event());
//	
//    }
//
//    std::cout << "total_pf = " << counter_npf << std::endl;
//
//
//    for(int iii = 0; iii < numOfch; iii ++){
//    
//      pat::PackedCandidate pf = pfcollection[iii];
//      //    std::cout << pf.pt() << " " << pf.eta() << " " << pf.phi() << std::endl;
//      
//      Bool_t isRight_ = false;
//      
//      
//      for(unsigned int mmm=0; mmm < gps.size(); mmm++){
//	
//	std::vector<TLorentzVector> tlvs = gps[mmm];
//	
//	for(unsigned int nnn=0; nnn < tlvs.size(); nnn++){
//	  
//	  if(
//	     reco::deltaR(pf.eta(), pf.phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.015 &&
//	     pf.pt()/tlvs[nnn].Pt() > 0.85 && 
//	     pf.pt()/tlvs[nnn].Pt() < 1.15
//	     ){
//	    isRight_ = true; 
//	    
//	  }
//	}
//      }
//      
//      //      std::cout << "PF id = " << iii << ", (pt,eta,phi) = " << pf.pt() << " " << pf.eta() << " " << pf.phi() << ", isRight = " << isRight_ << ", dnn score = " << mydnn[iii] << ",  gen tau pt = " << vec_gentaup4_vis[0].Pt() << std::endl;
//
//      //      if(isRight_) std::cout << "matched! " << iii << std::endl;
//
//      nBranches_->JpsiTau_ed_pfeta.push_back(pf.eta());
//      nBranches_->JpsiTau_ed_pfphi.push_back(pf.phi());
//      nBranches_->JpsiTau_ed_isRight.push_back(isRight_);
//      nBranches_->JpsiTau_ed_pfdnn.push_back(mydnn[iii]);
//      nBranches_->JpsiTau_ed_genpt.push_back(vec_gentaup4_vis[0].Pt());
//      nBranches_->JpsiTau_ed_id.push_back(event.id().event());
//
//    }
//  }




  if(verbose_) std::cout << "[JpsiTauNtuplizer] Starts to build tau candidate out of " << numOfch << " pion candidates" << std::endl;

  std::vector<taucand> cands;
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

	Int_t tau_charge = pf1.charge() + pf2.charge() + pf3.charge(); 

	if(TMath::Abs(tau_charge)!=(int)c_charge) continue; 


	/* reconstruct taus*/

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


	if(TMath::Prob(tau_vertex->chiSquared(), tau_vertex->degreesOfFreedom()) <= c_vprob) continue;

	  
	std::vector< RefCountedKinematicParticle > tau_children = tauTree->finalStateParticles();
	  
	math::PtEtaPhiMLorentzVector tau1_fit = aux.daughter_p4(tau_children, 0);
	math::PtEtaPhiMLorentzVector tau2_fit = aux.daughter_p4(tau_children, 1);
	math::PtEtaPhiMLorentzVector tau3_fit = aux.daughter_p4(tau_children, 2);


	particle_cand Taucand; 
	Taucand = aux.calculateIPvariables(extrapolator, tau_part, tau_vertex, closestVertex);

	//	std::cout << iii << " " << jjj << " " << kkk << std::endl;

	if(Taucand.fls3d < c_fsig){
	  //	  std::cout <<"remove" << std::endl;
	  continue;
	}

	std::vector<RefCountedKinematicParticle> allParticles;

	allParticles.push_back(pFactory.particle(mytracks[iii], aux.pion_mass, chi, ndf, aux.pion_sigma));
	allParticles.push_back(pFactory.particle(mytracks[jjj], aux.pion_mass, chi, ndf, aux.pion_sigma));
	allParticles.push_back(pFactory.particle(mytracks[kkk], aux.pion_mass, chi, ndf, aux.pion_sigma));
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

	math::PtEtaPhiMLorentzVector tt1_fit = aux.daughter_p4(bc_children, 0);
	math::PtEtaPhiMLorentzVector tt2_fit = aux.daughter_p4(bc_children, 1);
	math::PtEtaPhiMLorentzVector tt3_fit = aux.daughter_p4(bc_children, 2);

	math::PtEtaPhiMLorentzVector tlv_tau_fit = tt1_fit + tt2_fit + tt3_fit;

	if(tlv_tau_fit.Pt() < 3.){
	  //	  std::cout <<"remove pt" << std::endl;
	  continue;
	}

	
	// isolation calculation w.r.t SV
	Float_t iso = 0;
	Int_t ntracks = 0;
	Float_t iso_mindoca = 999; 
	
	for(int itrk = 0; itrk < numOfch;  itrk++){
    
	  if(itrk==iii || itrk==jjj || itrk==kkk) continue;

	  iso += pfcollection[itrk].pt();

	  TrajectoryStateOnSurface tsos_pf = extrapolator.extrapolate(mytracks[itrk].impactPointState(), bc_vertex->position());
    
    
	  VertexDistance3D a3d_pf;  

	  std::pair<bool,Measurement1D> cur3DIP_pf = aux.absoluteImpactParameter(tsos_pf, bc_vertex, a3d_pf);

	  Float_t pvip_pf = cur3DIP_pf.second.value();
    
	  if(pvip_pf < 0.03) ntracks+=1;

	  if(iso_mindoca > pvip_pf) iso_mindoca = pvip_pf;
        }



	Float_t max_dr_3prong = -1;

	Float_t dR_12 = reco::deltaR(tau1_fit.Eta(), tau1_fit.Phi(), tau2_fit.Eta(), tau2_fit.Phi());
	Float_t dR_13 = reco::deltaR(tau1_fit.Eta(), tau1_fit.Phi(), tau3_fit.Eta(), tau3_fit.Phi());
	Float_t dR_23 = reco::deltaR(tau2_fit.Eta(), tau2_fit.Phi(), tau3_fit.Eta(), tau3_fit.Phi());

	if(max_dr_3prong < dR_12) max_dr_3prong = dR_12;
	if(max_dr_3prong < dR_13) max_dr_3prong = dR_13;
	if(max_dr_3prong < dR_23) max_dr_3prong = dR_23;



	Bool_t isRight = false; 
	Bool_t isRight1 = false; 
	Bool_t isRight2 = false; 
	Bool_t isRight3 = false; 

	Int_t pid = -999;
	Float_t matched_gentaupt = -999;
	
	if(runOnMC_){

	  for(unsigned int mmm=0; mmm < gps.size(); mmm++){
	    
	    Bool_t isRight1_ = false;
	    Bool_t isRight2_ = false;
	    Bool_t isRight3_ = false;
	    
	    std::vector<TLorentzVector> tlvs = gps[mmm];
	    
	    for(unsigned int nnn=0; nnn < tlvs.size(); nnn++){

	      if(
		 reco::deltaR(tau1_fit.Eta(), tau1_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.015 &&
		 tau1_fit.Pt()/tlvs[nnn].Pt() > 0.85 && 
		 tau1_fit.Pt()/tlvs[nnn].Pt() < 1.15
		 ){
		isRight1_ = true; 
	      }
	      
	      if(
		 reco::deltaR(tau2_fit.Eta(), tau2_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.015 &&
		 tau2_fit.Pt()/tlvs[nnn].Pt() > 0.85 && 
		 tau2_fit.Pt()/tlvs[nnn].Pt() < 1.15
		 ){
		isRight2_ = true; 
	      }

	      
	      if(
		 reco::deltaR(tau3_fit.Eta(), tau3_fit.Phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.015 &&
		 tau3_fit.Pt()/tlvs[nnn].Pt() > 0.85 && 
		 tau3_fit.Pt()/tlvs[nnn].Pt() < 1.15
		 ){
		isRight3_ = true; 
	      }

	      
	    }
	    
	    Bool_t isRight_ = isRight1_ && isRight2_ && isRight3_;
	    if(isRight1_) isRight1 = true;
	    if(isRight2_) isRight2 = true;
	    if(isRight3_) isRight3 = true;

	    if(isRight_){
	      isRight = true;
	      pid = ppdgId[mmm];
	      matched_gentaupt = vec_gentau3pp4[mmm].Pt();
	    }
	  }	
	}

	if(isTruth_ && runOnMC_){
	  if(!isRight) continue;
	}



	Float_t sumofdnn = -1;
	Float_t dnn1 = -1;
	Float_t dnn2 = -1;
	Float_t dnn3 = -1;

	if(useDNN_){
	  sumofdnn = mydnn[iii] + mydnn[jjj] + mydnn[kkk];
	  dnn1 = mydnn[iii];
	  dnn2 = mydnn[jjj];
	  dnn3 = mydnn[kkk];
	}


	taucand _cand_ = {
	  iii,
	  jjj,
	  kkk,
	  (Float_t) tlv_tau_fit.Pt(),
	  (Float_t) tlv_tau_fit.Eta(),
	  (Float_t) tlv_tau_fit.Phi(),
	  (Float_t) tlv_tau_fit.M(),
	  Taucand,
	  (Float_t) TMath::Prob(tau_vertex->chiSquared(), tau_vertex->degreesOfFreedom()),
	  (Float_t) tau_vertex->vertexState().position().x(), 
	  (Float_t) tau_vertex->vertexState().position().y(), 
	  (Float_t) tau_vertex->vertexState().position().z(), 
	  (Float_t) max_dr_3prong, 
	  (Int_t) tau_charge,
	  (Bool_t) isRight,
	  (Bool_t) isRight1,
	  (Bool_t) isRight2,
	  (Bool_t) isRight3,
	  (Int_t) pid,
	  (Float_t) matched_gentaupt, 
	  (Float_t) sumofdnn,
	  (Float_t) dnn1,
	  (Float_t) dnn2,
	  (Float_t) dnn3,
	  (Float_t) TMath::Prob(bc_vertex->chiSquared(), bc_vertex->degreesOfFreedom()),
	  (Float_t) bc_vertex->vertexState().position().x(),
	  (Float_t) bc_vertex->vertexState().position().y(),
	  (Float_t) bc_vertex->vertexState().position().z(),
	  (Float_t) bc_part->currentState().globalMomentum().perp(),
	  (Float_t) bc_part->currentState().globalMomentum().eta(),
	  (Float_t) bc_part->currentState().globalMomentum().phi(),
	  (Float_t) bc_part->currentState().mass(),
	  Bcand,
	  (Float_t) iso,
	  (Float_t) ntracks,
	  (Float_t) iso_mindoca,
	};
	  
	cands.push_back(_cand_);
      }
    }
  }
    
  sort(cands.begin(), cands.end());


  if(cands.size()==0) return false;

  if(verbose_) std::cout << "[JpsiTauNtuplizer] " << cands.size() << " tau candidates were found" << std::endl;

  nBranches_->cutflow_perevt->Fill(9);

  Int_t ncomb = 0;

  for(int ic=0; ic < (int)cands.size(); ic++){
      
    ncomb += 1;


    /********************************************************************
     *
     * Step9: Filling normal branches
     *
     ********************************************************************/


    nBranches_->JpsiTau_tau_pt.push_back(cands[ic].cand_tau_pt);
    nBranches_->JpsiTau_tau_eta.push_back(cands[ic].cand_tau_eta);
    nBranches_->JpsiTau_tau_phi.push_back(cands[ic].cand_tau_phi);
    nBranches_->JpsiTau_tau_mass.push_back(cands[ic].cand_tau_mass);
    nBranches_->JpsiTau_tau_q.push_back(cands[ic].cand_tau_charge);
    nBranches_->JpsiTau_tau_vx.push_back(cands[ic].cand_tau_vx);
    nBranches_->JpsiTau_tau_vy.push_back(cands[ic].cand_tau_vy);
    nBranches_->JpsiTau_tau_vz.push_back(cands[ic].cand_tau_vz);
    nBranches_->JpsiTau_tau_max_dr_3prong.push_back(cands[ic].cand_tau_max_dr_3prong);
    nBranches_->JpsiTau_tau_lip.push_back(cands[ic].cand_tau.lip);
    nBranches_->JpsiTau_tau_lips.push_back(cands[ic].cand_tau.lips);
    nBranches_->JpsiTau_tau_pvip.push_back(cands[ic].cand_tau.pvip);
    nBranches_->JpsiTau_tau_pvips.push_back(cands[ic].cand_tau.pvips);
    nBranches_->JpsiTau_tau_fl3d.push_back(cands[ic].cand_tau.fl3d);
    nBranches_->JpsiTau_tau_fls3d.push_back(cands[ic].cand_tau.fls3d);
    nBranches_->JpsiTau_tau_alpha.push_back(cands[ic].cand_tau.alpha);
    nBranches_->JpsiTau_tau_vprob.push_back(cands[ic].cand_tau_vprob);
    nBranches_->JpsiTau_tau_isRight.push_back(cands[ic].cand_tau_isRight);
    nBranches_->JpsiTau_tau_isRight1.push_back(cands[ic].cand_tau_isRight1);
    nBranches_->JpsiTau_tau_isRight2.push_back(cands[ic].cand_tau_isRight2);
    nBranches_->JpsiTau_tau_isRight3.push_back(cands[ic].cand_tau_isRight3);
    nBranches_->JpsiTau_tau_matched_ppdgId.push_back(cands[ic].cand_tau_matched_ppdgId);
    nBranches_->JpsiTau_tau_matched_gentaupt.push_back(cands[ic].cand_tau_matched_gentaupt);
    nBranches_->JpsiTau_tau_sumofdnn.push_back(cands[ic].cand_tau_sumofdnn);
    nBranches_->JpsiTau_tau_pfidx1.push_back(cands[ic].cand_tau_id1);
    nBranches_->JpsiTau_tau_pfidx2.push_back(cands[ic].cand_tau_id2);
    nBranches_->JpsiTau_tau_pfidx3.push_back(cands[ic].cand_tau_id3);
    nBranches_->JpsiTau_tau_pi1_dnn.push_back(cands[ic].cand_tau_pi1_dnn);
    nBranches_->JpsiTau_tau_pi2_dnn.push_back(cands[ic].cand_tau_pi2_dnn);
    nBranches_->JpsiTau_tau_pi3_dnn.push_back(cands[ic].cand_tau_pi3_dnn);

    std::vector<Float_t> rhomass;
    pat::PackedCandidate pf1 = pfcollection[cands[ic].cand_tau_id1];
    pat::PackedCandidate pf2 = pfcollection[cands[ic].cand_tau_id2];
    pat::PackedCandidate pf3 = pfcollection[cands[ic].cand_tau_id3];
      
    TLorentzVector tlv_pion1; 
    TLorentzVector tlv_pion2;
    TLorentzVector tlv_pion3;
      
    tlv_pion1.SetPtEtaPhiM(pf1.pt(), pf1.eta(), pf1.phi(), pf1.mass());
    tlv_pion2.SetPtEtaPhiM(pf2.pt(), pf2.eta(), pf2.phi(), pf2.mass());
    tlv_pion3.SetPtEtaPhiM(pf3.pt(), pf3.eta(), pf3.phi(), pf3.mass());

	
    if(pf1.charge()*pf2.charge() == -1){
      TLorentzVector tlv_rho = tlv_pion1 + tlv_pion2;
      rhomass.push_back(tlv_rho.M());
    }
      
    if(pf1.charge()*pf3.charge() == -1){
      TLorentzVector tlv_rho = tlv_pion1 + tlv_pion3;
      rhomass.push_back(tlv_rho.M());
    }
      
    if(pf2.charge()*pf3.charge() == -1){
      TLorentzVector tlv_rho = tlv_pion2 + tlv_pion3;
      rhomass.push_back(tlv_rho.M());
    }
      
    if(rhomass.size()==2){
      nBranches_->JpsiTau_tau_rhomass1.push_back(rhomass.at(0));
      nBranches_->JpsiTau_tau_rhomass2.push_back(rhomass.at(1));
    }else{
      nBranches_->JpsiTau_tau_rhomass1.push_back(-1);
      nBranches_->JpsiTau_tau_rhomass2.push_back(-1);
    }


    nBranches_->JpsiTau_tau_pi1_pt.push_back(tlv_pion1.Pt());
    nBranches_->JpsiTau_tau_pi1_eta.push_back(tlv_pion1.Eta());
    nBranches_->JpsiTau_tau_pi1_phi.push_back(tlv_pion1.Phi());
    nBranches_->JpsiTau_tau_pi1_mass.push_back(tlv_pion1.M());
    nBranches_->JpsiTau_tau_pi1_q.push_back(pf1.charge());
      
    nBranches_->JpsiTau_tau_pi2_pt.push_back(tlv_pion2.Pt());
    nBranches_->JpsiTau_tau_pi2_eta.push_back(tlv_pion2.Eta());
    nBranches_->JpsiTau_tau_pi2_phi.push_back(tlv_pion2.Phi());
    nBranches_->JpsiTau_tau_pi2_mass.push_back(tlv_pion2.M());
    nBranches_->JpsiTau_tau_pi2_q.push_back(pf2.charge());

    nBranches_->JpsiTau_tau_pi3_pt.push_back(tlv_pion3.Pt());
    nBranches_->JpsiTau_tau_pi3_eta.push_back(tlv_pion3.Eta());
    nBranches_->JpsiTau_tau_pi3_phi.push_back(tlv_pion3.Phi());
    nBranches_->JpsiTau_tau_pi3_mass.push_back(tlv_pion3.M());
    nBranches_->JpsiTau_tau_pi3_q.push_back(pf3.charge());

    nBranches_->JpsiTau_B_pt.push_back(cands[ic].cand_b_pt);
    nBranches_->JpsiTau_B_eta.push_back(cands[ic].cand_b_eta);
    nBranches_->JpsiTau_B_phi.push_back(cands[ic].cand_b_phi);
    nBranches_->JpsiTau_B_mass.push_back(cands[ic].cand_b_mass);
    nBranches_->JpsiTau_B_vprob.push_back(cands[ic].cand_b_vprob);
    nBranches_->JpsiTau_B_lip.push_back(cands[ic].cand_b.lip);
    nBranches_->JpsiTau_B_lips.push_back(cands[ic].cand_b.lips);
    nBranches_->JpsiTau_B_pvip.push_back(cands[ic].cand_b.pvip);
    nBranches_->JpsiTau_B_pvips.push_back(cands[ic].cand_b.pvips);
    nBranches_->JpsiTau_B_fls3d.push_back(cands[ic].cand_b.fls3d);
    nBranches_->JpsiTau_B_fl3d.push_back(cands[ic].cand_b.fl3d);
    nBranches_->JpsiTau_B_alpha.push_back(cands[ic].cand_b.alpha);

    std::vector<RefCountedKinematicParticle> allParticles4doc;
      
    allParticles4doc.push_back(pFactory.particle(tt1_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));
    allParticles4doc.push_back(pFactory.particle(tt2_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));
    allParticles4doc.push_back(pFactory.particle(mytracks[cands[ic].cand_tau_id1], aux.pion_mass, chi, ndf, aux.pion_sigma));
    allParticles4doc.push_back(pFactory.particle(mytracks[cands[ic].cand_tau_id2], aux.pion_mass, chi, ndf, aux.pion_sigma));
    allParticles4doc.push_back(pFactory.particle(mytracks[cands[ic].cand_tau_id3], aux.pion_mass, chi, ndf, aux.pion_sigma));


    nBranches_->JpsiTau_B_maxdoca.push_back(aux.getMaxDoca(allParticles4doc));
    nBranches_->JpsiTau_B_mindoca.push_back(aux.getMinDoca(allParticles4doc));
    nBranches_->JpsiTau_B_vx.push_back(cands[ic].cand_b_vx);
    nBranches_->JpsiTau_B_vy.push_back(cands[ic].cand_b_vy);
    nBranches_->JpsiTau_B_vz.push_back(cands[ic].cand_b_vz);

    nBranches_->JpsiTau_B_iso.push_back(cands[ic].cand_b_iso);
    nBranches_->JpsiTau_B_iso_ntracks.push_back(cands[ic].cand_b_iso_ntracks);
    nBranches_->JpsiTau_B_iso_mindoca.push_back(cands[ic].cand_b_iso_mindoca);


      
    TLorentzVector Tlv_B;
    TLorentzVector Tlv_Jpsi;
    TLorentzVector Tlv_tau;
      
    Tlv_B.SetPtEtaPhiM(cands[ic].cand_b_pt,
		       cands[ic].cand_b_eta,
		       cands[ic].cand_b_phi,
		       cands[ic].cand_b_mass);

    // calculation of the corrected mass
    TVector3 *bvector = new TVector3(Tlv_B.Px(), Tlv_B.Py(), Tlv_B.Pz());
      
    Float_t pperp = bvector->Mag()*TMath::Sin(TMath::ACos(cands[ic].cand_b.alpha));
    
    Float_t mcorr = pperp + TMath::Sqrt(pperp*pperp + Tlv_B.M()*Tlv_B.M());
    
    nBranches_->JpsiTau_B_mcorr.push_back(mcorr);


    Tlv_B *= aux.mass_Bc/cands[ic].cand_b_mass;
      
    Tlv_Jpsi.SetPtEtaPhiM(jpsi_part->currentState().globalMomentum().perp(),
			  jpsi_part->currentState().globalMomentum().eta(),
			  jpsi_part->currentState().globalMomentum().phi(),
			  jpsi_part->currentState().mass());
      
    Tlv_tau.SetPtEtaPhiM(cands[ic].cand_tau_pt, 
			 cands[ic].cand_tau_eta, 
			 cands[ic].cand_tau_phi, 
			 cands[ic].cand_tau_mass);

      

    Float_t q2 = (Tlv_B - Tlv_Jpsi).M2();


    nBranches_->JpsiTau_B_q2.push_back(q2);


    Float_t mm2 = (Tlv_B - Tlv_Jpsi - Tlv_tau).M2();
    Float_t ptmiss = (Tlv_B - Tlv_Jpsi - Tlv_tau).Pt();

    nBranches_->JpsiTau_B_mm2.push_back(mm2);
    nBranches_->JpsiTau_B_ptmiss.push_back(ptmiss);

    Tlv_tau.Boost( -Tlv_B.BoostVector() );

    nBranches_->JpsiTau_B_Es.push_back(Tlv_tau.E()); 

    nBranches_->JpsiTau_B_ptback.push_back(Tlv_B.Pt()); 

  }

  nBranches_->JpsiTau_mu1_pt.push_back(mu1_fit.pt());
  nBranches_->JpsiTau_mu1_eta.push_back(mu1_fit.eta());
  nBranches_->JpsiTau_mu1_phi.push_back(mu1_fit.phi());
  nBranches_->JpsiTau_mu1_mass.push_back(mu1_fit.mass());
  nBranches_->JpsiTau_mu1_q.push_back(muoncollection[mcidx_mu1].charge());
  nBranches_->JpsiTau_mu1_isLoose.push_back(muoncollection[mcidx_mu1].isLooseMuon());
  nBranches_->JpsiTau_mu1_isTight.push_back(muoncollection[mcidx_mu1].isTightMuon(closestVertex));
  nBranches_->JpsiTau_mu1_isPF.push_back(muoncollection[mcidx_mu1].isPFMuon());
  nBranches_->JpsiTau_mu1_isGlobal.push_back(muoncollection[mcidx_mu1].isGlobalMuon());
  nBranches_->JpsiTau_mu1_isTracker.push_back(muoncollection[mcidx_mu1].isTrackerMuon());
  nBranches_->JpsiTau_mu1_isSoft.push_back(muoncollection[mcidx_mu1].isSoftMuon(closestVertex));
  nBranches_->JpsiTau_mu1_vx.push_back(muoncollection[mcidx_mu1].vx());
  nBranches_->JpsiTau_mu1_vy.push_back(muoncollection[mcidx_mu1].vy());
  nBranches_->JpsiTau_mu1_vz.push_back(muoncollection[mcidx_mu1].vz());
  nBranches_->JpsiTau_mu1_iso.push_back(1.);
  nBranches_->JpsiTau_mu1_dbiso.push_back(aux.MuonPFIso(muoncollection[mcidx_mu1]));
  
  nBranches_->JpsiTau_mu2_pt.push_back(mu2_fit.pt());
  nBranches_->JpsiTau_mu2_eta.push_back(mu2_fit.eta());
  nBranches_->JpsiTau_mu2_phi.push_back(mu2_fit.phi());
  nBranches_->JpsiTau_mu2_mass.push_back(mu2_fit.mass());
  nBranches_->JpsiTau_mu2_q.push_back(muoncollection[mcidx_mu2].charge());
  nBranches_->JpsiTau_mu2_isLoose.push_back(muoncollection[mcidx_mu2].isLooseMuon());
  nBranches_->JpsiTau_mu2_isTight.push_back(muoncollection[mcidx_mu2].isTightMuon(closestVertex));
  nBranches_->JpsiTau_mu2_isPF.push_back(muoncollection[mcidx_mu2].isPFMuon());
  nBranches_->JpsiTau_mu2_isGlobal.push_back(muoncollection[mcidx_mu2].isGlobalMuon());
  nBranches_->JpsiTau_mu2_isTracker.push_back(muoncollection[mcidx_mu2].isTrackerMuon());
  nBranches_->JpsiTau_mu2_isSoft.push_back(muoncollection[mcidx_mu2].isSoftMuon(closestVertex));
  nBranches_->JpsiTau_mu2_vx.push_back(muoncollection[mcidx_mu2].vx());
  nBranches_->JpsiTau_mu2_vy.push_back(muoncollection[mcidx_mu2].vy());
  nBranches_->JpsiTau_mu2_vz.push_back(muoncollection[mcidx_mu2].vz());
  nBranches_->JpsiTau_mu2_iso.push_back(2.);
  nBranches_->JpsiTau_mu2_dbiso.push_back(aux.MuonPFIso(muoncollection[mcidx_mu2]));

  nBranches_->JpsiTau_PV_vx.push_back(vertices_->begin()->position().x());
  nBranches_->JpsiTau_PV_vy.push_back(vertices_->begin()->position().y());
  nBranches_->JpsiTau_PV_vz.push_back(vertices_->begin()->position().z());

  nBranches_->JpsiTau_bbPV_vx.push_back(closestVertex.position().x());
  nBranches_->JpsiTau_bbPV_vy.push_back(closestVertex.position().y());
  nBranches_->JpsiTau_bbPV_vz.push_back(closestVertex.position().z());

  nBranches_->JpsiTau_Jpsi_pt.push_back(jpsi_part->currentState().globalMomentum().perp());
  nBranches_->JpsiTau_Jpsi_eta.push_back(jpsi_part->currentState().globalMomentum().eta());
  nBranches_->JpsiTau_Jpsi_phi.push_back(jpsi_part->currentState().globalMomentum().phi());
  nBranches_->JpsiTau_Jpsi_mass.push_back(jpsi_part->currentState().mass());
  nBranches_->JpsiTau_Jpsi_vprob.push_back(TMath::Prob(jpsi_part->chiSquared(), jpsi_part->degreesOfFreedom()));
  nBranches_->JpsiTau_Jpsi_lip.push_back(JPcand.lip);
  nBranches_->JpsiTau_Jpsi_lips.push_back(JPcand.lips);
  nBranches_->JpsiTau_Jpsi_pvip.push_back(JPcand.pvip);
  nBranches_->JpsiTau_Jpsi_pvips.push_back(JPcand.pvips);
  nBranches_->JpsiTau_Jpsi_fl3d.push_back(JPcand.fl3d);
  nBranches_->JpsiTau_Jpsi_fls3d.push_back(JPcand.fls3d);
  nBranches_->JpsiTau_Jpsi_alpha.push_back(JPcand.alpha);
  nBranches_->JpsiTau_Jpsi_maxdoca.push_back(aux.getMaxDoca(muonParticles));
  nBranches_->JpsiTau_Jpsi_mindoca.push_back(aux.getMinDoca(muonParticles));
  nBranches_->JpsiTau_Jpsi_vx.push_back(jpsi_vertex->vertexState().position().x());
  nBranches_->JpsiTau_Jpsi_vy.push_back(jpsi_vertex->vertexState().position().y());
  nBranches_->JpsiTau_Jpsi_vz.push_back(jpsi_vertex->vertexState().position().z());  
  nBranches_->JpsiTau_Jpsi_unfit_pt.push_back(jpsi_tlv_highest.Pt());
  nBranches_->JpsiTau_Jpsi_unfit_mass.push_back(jpsi_tlv_highest.M());
  nBranches_->JpsiTau_Jpsi_unfit_vprob.push_back(jpsi_vprob_highest);

  if(jpsi_vprob_highest!=-9){
    nBranches_->JpsiTau_Jpsi_unfit_vx.push_back(jpsi_vertex_highest.position().x());
    nBranches_->JpsiTau_Jpsi_unfit_vy.push_back(jpsi_vertex_highest.position().y());
    nBranches_->JpsiTau_Jpsi_unfit_vz.push_back(jpsi_vertex_highest.position().z());
  }



  /********************************************************************
   *
   * Step10: check gen-matching and fill them
   *
   ********************************************************************/

  nBranches_->JpsiTau_nch.push_back(numOfch);
  nBranches_->JpsiTau_nch_after_dnn.push_back(npf_after_dnn);
  nBranches_->JpsiTau_nch_before_dnn.push_back(npf_before_dnn);
  nBranches_->JpsiTau_nch_qr.push_back(npf_qr);
  nBranches_->IsJpsiTau.push_back(1.);
  nBranches_->JpsiTau_nCandidates.push_back(ncomb);

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
//      std::cout << "[JpsiTauNtuplizer] Hammer: " << _part_.pdgId() << " (" << _part_.status() << ")" << std::endl;
//
//      bool isJpsi = false;
//      bool isTau = false;
//      
//
//      for(auto d : _part_.daughterRefVector()) {	
//	if(TMath::Abs(d->pdgId()) == 443) isJpsi = true;
//	if(TMath::Abs(d->pdgId()) == 15) isTau = true;
//	
//	std::cout << "[JpsiTauNtuplizer] Hammer: ---> " << d->pdgId() << " (" << d->status() << ")" << std::endl;
//
//	for (auto dd : d->daughterRefVector()) {
//	  std::cout << "[JpsiTauNtuplizer] Hammer: -----> " << dd->pdgId() << " (" << dd->status()  << ")" << std::endl;
//	}
//      }
//      std::cout << isJpsi << " " << isTau  << " " << (isJpsi==true && isTau==true) << std::endl;
//      if(!(isJpsi==true && isTau==true)){
//	std::cout << "[JpsiTauNtuplizer] Hammer: This B is rejected !!!!" << std::endl;
//      }      
//    }


    for( unsigned p=0; p < genParticles_->size(); ++p){
	  
      // Bc daughters loop
      // only allow Bc+ as the MC is produced as such
      // ie if we take Bc-, we take probe side by mistake ... 
      //      if(!((*genParticles_)[p].pdgId()==541 && (*genParticles_)[p].status()==2)) continue;

      
      if(!(TMath::Abs((*genParticles_)[p].pdgId())==541 && (*genParticles_)[p].status()==2)) continue;
	
      auto _part_ = (*genParticles_)[p];
      //      if(verbose_) std::cout << "[JpsiTauNtuplizer] Hammer2: " << _part_.pdgId() << " (" << _part_.status() << ")" << std::endl;
      
      bool isJpsi = false;
      bool isTau = false;

      
      // check if the Bc candidate has both J/psi and tau
      for(auto d : _part_.daughterRefVector()) {
	if(TMath::Abs(d->pdgId()) == 443) isJpsi = true;
	if(TMath::Abs(d->pdgId()) == 15) isTau = true;
      }

      //      std::cout << isJpsi << " " << isTau  << " " << (isJpsi==true && isTau==true) << std::endl;
      if(!(isJpsi==true && isTau==true)){
	//	if(verbose_) std::cout << "[JpsiTauNtuplizer] Hammer: This B is rejected !!!!" << std::endl;
	continue;
      }
      
      
      Hammer::Particle pB({_part_.energy(), _part_.px(), _part_.py(), _part_.pz()}, _part_.pdgId());
      
      idxB = Bc2JpsiLNu.addParticle(pB);
      
      if(verbose_) std::cout << "[JpsiTauNtuplizer] Hammer: " << (*genParticles_)[p].pdgId() << " " << (*genParticles_)[p].status() << std::endl;
      
      
      for(auto d : _part_.daughterRefVector()) {
	
	Hammer::Particle B_dau({d->energy(), d->px(), d->py(), d->pz()}, d->pdgId());
	
	auto idx_d = Bc2JpsiLNu.addParticle(B_dau);
	Bvtx_idxs.push_back(idx_d);

	if(verbose_) std::cout << "[JpsiTauNtuplizer] Hammer: \t gen: " << d->pdgId() << " " << d->status() << std::endl;	  
	
	
	if(TMath::Abs(d->pdgId()) == 15) {
	  
	  idxTau = idx_d;
	  
	  for (auto dd : d->daughterRefVector()) {
	    Hammer::Particle Tau_dau({dd->energy(), dd->px(), dd->py(), dd->pz()}, dd->pdgId());
	    auto idx_dd = Bc2JpsiLNu.addParticle(Tau_dau);
	    Tauvtx_idxs.push_back(idx_dd);
	    if(verbose_) std::cout << "[JpsiTauNtuplizer] Hammer: \t\t gen: " << dd->pdgId() << " " << dd->status() << std::endl;
	  }
	}
	
	else if(TMath::Abs(d->pdgId()) == 443) {
	  
	  idxJpsi = idx_d;
	  
	  for (auto dd : d->daughterRefVector()) {
	    Hammer::Particle Jpsi_dau({dd->energy(), dd->px(), dd->py(), dd->pz()}, dd->pdgId());
	    auto idx_dd = Bc2JpsiLNu.addParticle(Jpsi_dau);
	    Jpsivtx_idxs.push_back(idx_dd);
	    if(verbose_) std::cout << "[JpsiTauNtuplizer] Hammer: \t\t gen: " << dd->pdgId() << " " << dd->status() << std::endl;
	  }

	}
      }
    }

    if(verbose_) std::cout << "[JpsiTauNtuplizer] Hammer idx (B, tau, Jpsi) = " << idxB << " " << idxTau << " " << idxJpsi << std::endl;
      
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

    for(auto pars: _FFErrNames) {
      settings[pars] = 0;
    }

    hammer.setFFEigenvectors("BctoJpsi", "BGLVar", settings);
	
    auto weights = hammer.getWeights("Scheme1");

    Float_t weight = -1;

    if(!weights.empty()){
      for(auto elem: weights) {
	if(isnan(elem.second)) {
	  std::cout << "[JpsiTauNtuplizer] ERROR: BGL Central nan weight: " << elem.second << std::endl;
	}else{
	  weight = elem.second;
	}
      }
    }

    nBranches_->JpsiTau_hammer_ebe.push_back(weight);

	
    int idx1 = 0;
	
    for(auto pars1: _FFErrNames) {
	  
      for(int isUp=0; isUp < 2; isUp++) { // up, down variation    
	    
	std::map<std::string, double> settings;
	    
	int idx2 = 0;
	    
	for(auto pars2: _FFErrNames) {
	      
	  if(idx1 == idx2){
	    if(isUp==1) settings[pars2] = _FFErr[idx2];
	    else settings[pars2] = -_FFErr[idx2];
	  }else{
	    settings[pars2] = 0;
	  }
	  idx2 += 1;
	}

 
	hammer.setFFEigenvectors("BctoJpsi", "BGLVar", settings);
	auto weights = hammer.getWeights("Scheme1");
	std::string var_name = pars1;
	var_name += isUp==1? "_Up" : "_Down";

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

	if(var_name==std::string("delta_a0_Up")) nBranches_->JpsiTau_hammer_ebe_a0_up.push_back(weight_sys);
	else if(var_name==std::string("delta_a0_Down")) nBranches_->JpsiTau_hammer_ebe_a0_down.push_back(weight_sys);
	else if(var_name==std::string("delta_a1_Up")) nBranches_->JpsiTau_hammer_ebe_a1_up.push_back(weight_sys);
	else if(var_name==std::string("delta_a1_Down")) nBranches_->JpsiTau_hammer_ebe_a1_down.push_back(weight_sys);
	else if(var_name==std::string("delta_a2_Up")) nBranches_->JpsiTau_hammer_ebe_a2_up.push_back(weight_sys);
	else if(var_name==std::string("delta_a2_Down")) nBranches_->JpsiTau_hammer_ebe_a2_down.push_back(weight_sys);

	else if(var_name==std::string("delta_b0_Up")) nBranches_->JpsiTau_hammer_ebe_b0_up.push_back(weight_sys);
	else if(var_name==std::string("delta_b0_Down")) nBranches_->JpsiTau_hammer_ebe_b0_down.push_back(weight_sys);
	else if(var_name==std::string("delta_b1_Up")) nBranches_->JpsiTau_hammer_ebe_b1_up.push_back(weight_sys);
	else if(var_name==std::string("delta_b1_Down")) nBranches_->JpsiTau_hammer_ebe_b1_down.push_back(weight_sys);
	else if(var_name==std::string("delta_b2_Up")) nBranches_->JpsiTau_hammer_ebe_b2_up.push_back(weight_sys);
	else if(var_name==std::string("delta_b2_Down")) nBranches_->JpsiTau_hammer_ebe_b2_down.push_back(weight_sys);

	else if(var_name==std::string("delta_c1_Up")) nBranches_->JpsiTau_hammer_ebe_c1_up.push_back(weight_sys);
	else if(var_name==std::string("delta_c1_Down")) nBranches_->JpsiTau_hammer_ebe_c1_down.push_back(weight_sys);
	else if(var_name==std::string("delta_c2_Up")) nBranches_->JpsiTau_hammer_ebe_c2_up.push_back(weight_sys);
	else if(var_name==std::string("delta_c2_Down")) nBranches_->JpsiTau_hammer_ebe_c2_down.push_back(weight_sys);

	else if(var_name==std::string("delta_d0_Up")) nBranches_->JpsiTau_hammer_ebe_d0_up.push_back(weight_sys);
	else if(var_name==std::string("delta_d0_Down")) nBranches_->JpsiTau_hammer_ebe_d0_down.push_back(weight_sys);
	else if(var_name==std::string("delta_d1_Up")) nBranches_->JpsiTau_hammer_ebe_d1_up.push_back(weight_sys);
	else if(var_name==std::string("delta_d1_Down")) nBranches_->JpsiTau_hammer_ebe_d1_down.push_back(weight_sys);
	else if(var_name==std::string("delta_d2_Up")) nBranches_->JpsiTau_hammer_ebe_d2_up.push_back(weight_sys);
	else if(var_name==std::string("delta_d2_Down")) nBranches_->JpsiTau_hammer_ebe_d2_down.push_back(weight_sys);
	   
      }
      
      idx1 += 1;
      
    }
  }







  nBranches_->JpsiTau_q2_gen.push_back(q2_gen.M2());
  nBranches_->JpsiTau_B_pt_gen.push_back(pB_gen.Pt());
  nBranches_->JpsiTau_B_eta_gen.push_back(pB_gen.Eta());
  nBranches_->JpsiTau_B_phi_gen.push_back(pB_gen.Phi());
  nBranches_->JpsiTau_B_mass_gen.push_back(pB_gen.M());
  

  
  bool flag_jpsi_match = false;
  if(gen_jpsi_mu.size()==2){

    Float_t _dR_11 = reco::deltaR(gen_jpsi_mu[0].Eta(), gen_jpsi_mu[0].Phi(), 
				  mu1_fit.eta(), mu1_fit.phi());
    
    Float_t _dR_22 = reco::deltaR(gen_jpsi_mu[1].Eta(), gen_jpsi_mu[1].Phi(), 
				  mu2_fit.eta(), mu2_fit.phi());

    Float_t _dR_21 = reco::deltaR(gen_jpsi_mu[1].Eta(), gen_jpsi_mu[1].Phi(), 
				  mu1_fit.eta(), mu1_fit.phi());
    
    Float_t _dR_12 = reco::deltaR(gen_jpsi_mu[0].Eta(), gen_jpsi_mu[0].Phi(), 
				  mu2_fit.eta(), mu2_fit.phi());

    if(_dR_11 < 0.1 && _dR_22 < 0.1) flag_jpsi_match = true;
    else if(_dR_21 < 0.1 && _dR_12 < 0.1) flag_jpsi_match = true;
  }
  
  nBranches_->JpsiTau_genPV_vx.push_back(genvertex.x());
  nBranches_->JpsiTau_genPV_vy.push_back(genvertex.y());
  nBranches_->JpsiTau_genPV_vz.push_back(genvertex.z());
  nBranches_->JpsiTau_ngenmuons.push_back(gen_nr_mu.size() + gen_jpsi_mu.size());
  nBranches_->JpsiTau_isgenmatched.push_back((int)flag_jpsi_match);


  nBranches_->JpsiTau_isgen3.push_back(isgen3);
  nBranches_->JpsiTau_isgen3matched.push_back(isgen3matched);
  nBranches_->JpsiTau_ngentau3.push_back(gps.size());
  nBranches_->JpsiTau_ngentau.push_back(vec_gentaudm.size());


  if(vec_gentaudm.size() >=1){
    nBranches_->JpsiTau_gentaupt.push_back(vec_gentaup4_vis[0].Pt());
    nBranches_->JpsiTau_gentaueta.push_back(vec_gentaup4_vis[0].Eta());
    nBranches_->JpsiTau_gentauphi.push_back(vec_gentaup4_vis[0].Phi());
    nBranches_->JpsiTau_gentaumass.push_back(vec_gentaup4_vis[0].M());
    nBranches_->JpsiTau_gentaupt_bd.push_back(vec_gentaup4[0].Pt());
    nBranches_->JpsiTau_gentaueta_bd.push_back(vec_gentaup4[0].Eta());
    nBranches_->JpsiTau_gentauphi_bd.push_back(vec_gentaup4[0].Phi());
    nBranches_->JpsiTau_gentaumass_bd.push_back(vec_gentaup4[0].M());
    nBranches_->JpsiTau_gentaudm.push_back(vec_gentaudm[0]);
  }else{
    nBranches_->JpsiTau_gentaupt.push_back(-9);
    nBranches_->JpsiTau_gentaueta.push_back(-9);
    nBranches_->JpsiTau_gentauphi.push_back(-9);
    nBranches_->JpsiTau_gentaumass.push_back(-9);
    nBranches_->JpsiTau_gentaupt_bd.push_back(-9);
    nBranches_->JpsiTau_gentaueta_bd.push_back(-9);
    nBranches_->JpsiTau_gentauphi_bd.push_back(-9);
    nBranches_->JpsiTau_gentaumass_bd.push_back(-9);
    nBranches_->JpsiTau_gentaudm.push_back(-9);
  }

  return true;

}


