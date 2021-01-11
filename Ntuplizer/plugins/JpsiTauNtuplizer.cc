#include "../interface/JpsiTauNtuplizer.h"


//===================================================================================================================
JpsiTauNtuplizer::JpsiTauNtuplizer( edm::EDGetTokenT<pat::MuonCollection>    muonToken   ,
				    edm::EDGetTokenT<reco::VertexCollection> verticeToken, 
				    edm::EDGetTokenT<pat::PackedCandidateCollection> packedpfcandidatesToken,
				    edm::EDGetTokenT<edm::TriggerResults> triggertoken,
				    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerobject,
				    edm::EDGetTokenT<reco::GenParticleCollection> genptoken,
				    edm::EDGetTokenT<pat::PackedGenParticleCollection> packedgenptoken,
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
  , packedgenParticlesToken_( packedgenptoken )
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
  , dnnfile_old_ (runStrings["dnnfile_old"])      
  , dnnfile_perPF_ (runStrings["dnnfile_perPF"])      
  , dnnfile_perEVT_ (runStrings["dnnfile_perEVT"])      
   
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
    
    std::string dnnfilepath_old = edm::FileInPath("EXOVVNtuplizerRunII/Ntuplizer/" +  dnnfile_old_).fullPath();
    std::string dnnfilepath_perPF = edm::FileInPath("EXOVVNtuplizerRunII/Ntuplizer/" +  dnnfile_perPF_).fullPath();
    std::string dnnfilepath_perEVT = edm::FileInPath("EXOVVNtuplizerRunII/Ntuplizer/" +  dnnfile_perEVT_).fullPath();


    if(verbose_) std::cout << "[JpsiTauNtuplizer] DNN old  = " << dnnfilepath_old << std::endl;
    if(verbose_) std::cout << "[JpsiTauNtuplizer] DNN per pion  = " << dnnfilepath_perPF << std::endl;
    if(verbose_) std::cout << "[JpsiTauNtuplizer] DNN per event = " << dnnfilepath_perEVT << std::endl;


    std::string tbr_old = "DUMMY";  // to be replaced
    auto pos = dnnfilepath_old.find(tbr_old);
    auto len = tbr_old.length();
    if (pos != std::string::npos) {
      dnnfilepath_old.replace(pos, len, "");
    }


    std::string tbr_perPF = "DUMMY";  // to be replaced
    pos = dnnfilepath_perPF.find(tbr_perPF);
    len = tbr_perPF.length();
    if (pos != std::string::npos) {
      dnnfilepath_perPF.replace(pos, len, "");
    }

    std::string tbr_perEVT = "DUMMY";  // to be replaced
    pos = dnnfilepath_perEVT.find(tbr_perEVT);
    len = tbr_perEVT.length();
    if (pos != std::string::npos) {
      dnnfilepath_perEVT.replace(pos, len, "");
    }
    
    graphDef_old = tensorflow::loadMetaGraph(dnnfilepath_old);
    session_old = tensorflow::createSession(graphDef_old, dnnfilepath_old);

    graphDef_perPF = tensorflow::loadMetaGraph(dnnfilepath_perPF);
    session_perPF = tensorflow::createSession(graphDef_perPF, dnnfilepath_perPF);

    graphDef_perEVT = tensorflow::loadMetaGraph(dnnfilepath_perEVT);
    session_perEVT = tensorflow::createSession(graphDef_perEVT, dnnfilepath_perEVT);

    
    data = tensorflow::Tensor(tensorflow::DT_FLOAT, { 1, numberofDNN, 13 }); // single batch of dimension 10
    label_perPF = tensorflow::Tensor(tensorflow::DT_INT32, { 1,numberofDNN}); 
    label_perEVT = tensorflow::Tensor(tensorflow::DT_INT32, { 1, 1});
    //    add_global = tensorflow::Tensor(tensorflow::DT_FLOAT, { 1, 2 }); 
    isTraining = tensorflow::Tensor(tensorflow::DT_BOOL, tensorflow::TensorShape()); 
    //    norm = tensorflow::Tensor(tensorflow::DT_FLOAT, { 1, numberofDNN }); 


    data_old = tensorflow::Tensor(tensorflow::DT_FLOAT, { 1, 50, 8 }); // single batch of dimension 10
    label_old = tensorflow::Tensor(tensorflow::DT_INT32, { 1,50}); 
    add_global_old = tensorflow::Tensor(tensorflow::DT_FLOAT, { 1, 2 }); 
    isTraining_old = tensorflow::Tensor(tensorflow::DT_BOOL, tensorflow::TensorShape()); 
    norm_old = tensorflow::Tensor(tensorflow::DT_FLOAT, { 1, 50 }); 


    
  }
  

  if(runOnMC_ && useHammer_){

    ran = new TRandom3();
    ran->SetSeed(1);

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
//    std::cout << "[Hammer]: ABCD matrix:\n\t" << centralValuesOpt << std::endl;
//    hammer.setOptions(centralValuesOpt);

    
    // central values ... 
    //    std::string centralValues = "BctoJpsiBGLVar: {avec: [-0.4182681183, -0.1819441298, 0.1140362192], bvec: [-0.0901217995, 0.0074593144, -0.0127818946], cvec: [0.0028349558, 0.0310274035], dvec: [0.0013219585, -0.0043841356, -0.0000000258]}";

    //    std::cout << "[Hammer]: Central values:\n\t" << centralValues << std::endl;
    //    hammer.setOptions(centralValues);



    hammer.saveOptionCard("Opts.yml", false);
    
    //    std::cout << "... finishes " << std::endl;


    //    string centralValuesOpt = "BctoJpsiBGLVar: {";

    //    for(auto i=0; i<4; i++) {
    //    centralValuesOpt += Form("avec: %f, ", avec_cent);
      //}


     
	
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
      //      std::vector<float> settings_ff;

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

	//	settings_ff.push_back(settings[pars1]);

	idx_err += 1;
      }


      //      if(imc <= 1) std::cout << imc << " settings of " << "delta_a0" << ": " << settings["delta_a0"] << std::endl;

      FFdict.push_back(settings);
      
      //      nBranches_->JpsiTau_hammer_ff.push_back(settings_ff);

    }

    
    if(verbose_) std::cout << "Saved " << FFdict.size() << " FF variations" << std::endl;



  }
}




//===================================================================================================================
JpsiTauNtuplizer::~JpsiTauNtuplizer( void )
{

  
  tensorflow::closeSession(session_old);
  tensorflow::closeSession(session_perPF);
  tensorflow::closeSession(session_perEVT);
  delete graphDef_old;
  delete graphDef_perPF;
  delete graphDef_perEVT;


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
      for(auto n: parName) settings["delta_" + n] = 0;
      
      hammer.setFFEigenvectors("BctoJpsi", "BGLVar", settings);
      
      outRate["Central"] = hammer.getRate(proc, "Scheme1");
      std::cout << Form("[JpsiTauNtuplizer] New rate: %1.3e (ratio = %.3f)", outRate["Central"], outRate["Central"]/outRate["den"]) << std::endl;
      
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
  
	  if(verbose_) std::cout << "[JpsiTauNtuplizer] " << var_name << " --> " << Form(": %1.3e", rate) << std::endl;
	  
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


    
//      int idx1 = 0;
//      
//      for(auto pars1: _FFErrNames) {
//	
//	for(int isUp=0; isUp < 2; isUp++) { // up, down variation    
//	  
//	  std::map<std::string, double> settings;
//	  
//	  int idx2 = 0;
//	  
//	  for(auto pars2: _FFErrNames) {
//	    
//	    if(idx1 == idx2){
//	      if(isUp==1) settings[pars2] = _FFErr[idx2];
//	      else settings[pars2] = -_FFErr[idx2];
//	    }else{
//	      settings[pars2] = 0;
//	    }
//	    idx2 += 1;
//	  }
//	  
//	  
//	  hammer.setFFEigenvectors("BctoJpsi", "BGLVar", settings);
//	  
//	  auto rate = hammer.getRate(proc, "Scheme1");
//	  
//	  std::string var_name = pars1;
//	  var_name += isUp==1? "_Up" : "_Down";
//	  
//	  if(verbose_) std::cout << "[JpsiTauNtuplizer] " << var_name << " --> " << Form(": %1.3e", rate) << std::endl;
//	  
//	  if(var_name==std::string("delta_a0_Up")) nBranches_->hammer_width->SetBinContent(3, rate);
//	  else if(var_name==std::string("delta_a0_Down")) nBranches_->hammer_width->SetBinContent(4, rate);
//	  else if(var_name==std::string("delta_a1_Up")) nBranches_->hammer_width->SetBinContent(5, rate);
//	  else if(var_name==std::string("delta_a1_Down")) nBranches_->hammer_width->SetBinContent(6, rate);
//	  else if(var_name==std::string("delta_a2_Up")) nBranches_->hammer_width->SetBinContent(7, rate);
//	  else if(var_name==std::string("delta_a2_Down")) nBranches_->hammer_width->SetBinContent(8, rate);
//	  
//	  else if(var_name==std::string("delta_b0_Up")) nBranches_->hammer_width->SetBinContent(9, rate);
//	  else if(var_name==std::string("delta_b0_Down")) nBranches_->hammer_width->SetBinContent(10, rate);
//	  else if(var_name==std::string("delta_b1_Up")) nBranches_->hammer_width->SetBinContent(11, rate);
//	  else if(var_name==std::string("delta_b1_Down")) nBranches_->hammer_width->SetBinContent(12, rate);
//	  else if(var_name==std::string("delta_b2_Up")) nBranches_->hammer_width->SetBinContent(13, rate);
//	  else if(var_name==std::string("delta_b2_Down")) nBranches_->hammer_width->SetBinContent(14, rate);
//	  
//	  else if(var_name==std::string("delta_c1_Up")) nBranches_->hammer_width->SetBinContent(15, rate);
//	  else if(var_name==std::string("delta_c1_Down")) nBranches_->hammer_width->SetBinContent(16, rate);
//	  else if(var_name==std::string("delta_c2_Up")) nBranches_->hammer_width->SetBinContent(17, rate);
//	  else if(var_name==std::string("delta_c2_Down")) nBranches_->hammer_width->SetBinContent(18, rate);
//	  
//	  else if(var_name==std::string("delta_d0_Up")) nBranches_->hammer_width->SetBinContent(19, rate);
//	  else if(var_name==std::string("delta_d0_Down")) nBranches_->hammer_width->SetBinContent(20, rate);
//	  else if(var_name==std::string("delta_d1_Up")) nBranches_->hammer_width->SetBinContent(21, rate);
//	  else if(var_name==std::string("delta_d1_Down")) nBranches_->hammer_width->SetBinContent(22, rate);
//	  else if(var_name==std::string("delta_d2_Up")) nBranches_->hammer_width->SetBinContent(23, rate);
//	  else if(var_name==std::string("delta_d2_Down")) nBranches_->hammer_width->SetBinContent(24, rate);
//	  
//	  
//	}
//	idx1 += 1;
//	
//      }
    }
  }
}



bool JpsiTauNtuplizer::fillBranches( edm::Event const & event, const edm::EventSetup& iSetup ){

  if(verbose_) std::cout << "[JpsiTauNtuplizer] ---------------- event, run, lumi = " << event.id().event() << " " << event.id().run() << " " << event.id().luminosityBlock() << "----------------" << std::endl;
  

  TVector3 genvertex(-9.,-9.,-9.);

  std::vector<TLorentzVector> gen_nr_mu;
  std::vector<TLorentzVector> gen_jpsi_mu;

  std::vector<float> match_pt;
  std::vector<float> match_eta;
  std::vector<int> match_pdg;
  std::vector<int> match_ppdg;
  std::vector<float> match_phi;
  std::vector<bool> match_isSignal;
  std::vector<int> match_nprong;

  TLorentzVector pB_gen;
  TLorentzVector pJpsi_gen;
   
  if(runOnMC_){ 
    
    event.getByToken(genParticlesToken_ , genParticles_);   
    event.getByToken(packedgenParticlesToken_ , packedgenParticles_);   

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




    // check stable particles !!!


    for(size_t igp=0; igp < genParticles_->size(); ++igp){
      
      bool isB( (abs((*genParticles_)[igp].pdgId())>=511 && abs((*genParticles_)[igp].pdgId())<=545));
      if(!isB) continue;
      
      Bool_t ismother_B = false;
      
      for( unsigned int m=0; m<(*genParticles_)[igp].numberOfMothers(); ++m ){
	
	if( abs((*genParticles_)[igp].mother(m)->pdgId())>=511 && abs((*genParticles_)[igp].mother(m)->pdgId())<=545 ){
	  ismother_B = true;
	}
      }
    
      if( (*genParticles_)[igp].numberOfDaughters()!=1 && 
	  (*genParticles_)[igp].numberOfMothers()==1 && 
	  (*genParticles_)[igp].mother(0)->pdgId()==(*genParticles_)[igp].pdgId() 
	  ){
	ismother_B = false; 
      }
      
      if(! ( (*genParticles_)[igp].numberOfDaughters()==1 && (*genParticles_)[igp].daughter(0)->pdgId()==(*genParticles_)[igp].pdgId() ) && ismother_B==false ){

	
	const reco::Candidate * bMeson = &(*genParticles_)[igp];
	//      std::cout << "igp = " << igp << std::endl;
	//      std::cout << "PdgID: " << bMeson->pdgId() << " pt " << bMeson->pt() << " eta: " << bMeson->eta() << " phi: " << bMeson->phi() << std::endl;
	
	for(size_t jpgp=0; jpgp < packedgenParticles_->size(); ++jpgp){
	  //	std::cout << "jpgp = " << jpgp  << " / " << packedgenParticles_->size()<< std::endl;
	  //	if(!(*packedgenParticles_)[jpgp]) std::cout << "This is bad!!!" << std::endl;
	  //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
	  //	  std::cout << "try-1" << std::endl;
	  const reco::Candidate * motherInPrunedCollection = (*packedgenParticles_)[jpgp].mother(0);
	  //	std::cout << "try0" << std::endl;
	  
	  if(motherInPrunedCollection != nullptr && aux.isAncestor( bMeson , motherInPrunedCollection)){
	    //	  std::cout << "     PdgID: " << (*packedgenParticles_)[jpgp].pdgId() << " pt " << (*packedgenParticles_)[jpgp].pt() << " eta: " << (*packedgenParticles_)[jpgp].eta() << " phi: " << (*packedgenParticles_)[jpgp].phi() << " " << (*packedgenParticles_)[jpgp].charge()  << " mother = " <<  (*packedgenParticles_)[jpgp].mother(0)->pdgId()  << std::endl;
	    
	    if(TMath::Abs((*packedgenParticles_)[jpgp].charge() )==1){
	      
	      match_pt.push_back( (*packedgenParticles_)[jpgp].pt() );
	      match_eta.push_back( (*packedgenParticles_)[jpgp].eta() );
	      match_pdg.push_back( (*packedgenParticles_)[jpgp].pdgId());
	      match_ppdg.push_back((*genParticles_)[igp].pdgId() );
	      match_phi.push_back( (*packedgenParticles_)[jpgp].phi() );
	      
	      bool isSignal = false;
	      if( TMath::Abs((*genParticles_)[igp].pdgId())==541 && TMath::Abs((*packedgenParticles_)[jpgp].pdgId())==211 && TMath::Abs( (*packedgenParticles_)[jpgp].mother(0)->pdgId() )==15){
		if((*packedgenParticles_)[jpgp].mother(0)->numberOfMothers()!=0){
		  if(TMath::Abs( (*packedgenParticles_)[jpgp].mother(0)->mother(0)->pdgId()) == 541) isSignal = true;
		  //		std::cout << "     --> from " <<  (*packedgenParticles_)[jpgp].mother(0)->pdgId() << " from " << (*packedgenParticles_)[jpgp].mother(0)->mother(0)->pdgId()  << std::endl;
		}
	      }
	      
	      Int_t nprong = 0;
	      
	      if( TMath::Abs((*packedgenParticles_)[jpgp].pdgId())==211 && TMath::Abs( (*packedgenParticles_)[jpgp].mother(0)->pdgId() )==15){
		//	      std::cout << " XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX tau " <<  std::endl;
		for(unsigned int jdau=0; jdau < (*packedgenParticles_)[jpgp].mother(0)->numberOfDaughters(); jdau++){
		  //		std::cout << " ---------- " << (*packedgenParticles_)[jpgp].mother(0)->daughter(jdau)->pdgId() <<  std::endl;
		  if(TMath::Abs((*packedgenParticles_)[jpgp].mother(0)->daughter(jdau)->pdgId())==211){
		    nprong ++; 
		  }
		  
		  //	      std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX " <<  std::endl;
		}
	      }
	      
	      //	    std::cout << "prong ... " << nprong << std::endl;
	      
	      match_nprong.push_back(nprong);
	      match_isSignal.push_back(isSignal);
	      
	    }
	  }
	}
      }
    }



  }


  //std::cout << "test:" << q2_gen.Pt() << " " << q2_gen.M2() << " " << pB_gen.Pt() << " " << pJpsi_gen.Pt() << std::endl;
  TLorentzVector q2_gen; 
  q2_gen = pB_gen - pJpsi_gen;
  
  if(runOnMC_){ 
    nBranches_->q2_nocut->Fill(q2_gen.M2());
  }














  
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
  //  std::vector<int> muoncollection_id;
  muoncollection.clear();
  //  muoncollection_id.clear();


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
    //    muoncollection_id.push_back(imuon);
  }

  nBranches_->nmuon->Fill( muoncollection.size() );

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
  
//  //creating the vertex fitter
//  KinematicParticleVertexFitter kpvFitter;
//   
//  //reconstructing a J/Psi decay
//  RefCountedKinematicTree jpTree = kpvFitter.fit(muonParticles);
//   
//  if(jpTree->isEmpty() || !jpTree->isValid() || !jpTree->isConsistent()) return false;
//
//  nBranches_->cutflow_perevt->Fill(4);
//
//  //creating the particle fitter
//  KinematicParticleFitter csFitter;
//
//  // creating the constraint
//  KinematicConstraint* jpsi_constraint = new MassKinematicConstraint(aux.jpsi_mass, aux.jp_m_sigma);
//  //the constrained fit
//  jpTree = csFitter.fit(jpsi_constraint, jpTree);
//
//  //getting the J/Psi KinematicParticle
//  jpTree->movePointerToTheTop();
//  RefCountedKinematicParticle jpsi_part = jpTree->currentParticle();
//  if(!jpsi_part->currentState().isValid()) return false; 
//  nBranches_->cutflow_perevt->Fill(5);
//
//  RefCountedKinematicVertex jpsi_vertex = jpTree->currentDecayVertex();
//  if(!jpsi_vertex->vertexIsValid()) return false; 
//  nBranches_->cutflow_perevt->Fill(6);
//
//  if(TMath::Prob(jpsi_vertex->chiSquared(), jpsi_vertex->degreesOfFreedom()) <=0) return false;
//  nBranches_->cutflow_perevt->Fill(7);



  RefCountedKinematicParticle jpsi_part;
  RefCountedKinematicVertex jpsi_vertex;
  RefCountedKinematicTree jpTree;
  Bool_t jpsifit_flag;
  //  std::tie(jpsifit_flag, jpsi_part, jpsi_vertex, jpTree) = aux.KinematicFit(muonParticles, aux.jpsi_mass, aux.jpsi_m_sigma);
  std::tie(jpsifit_flag, jpsi_part, jpsi_vertex, jpTree) = aux.KinematicFit(muonParticles, -1, -1);

  if(!jpsifit_flag) return false;
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
  TransverseImpactPointExtrapolator extrapolatort(fMagneticField);

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

	  for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){   
	    
	    pat::PackedCandidate _pf = (*packedpfcandidates_)[ii];

	    //	  for(int kkk = 0; kkk < numOfch; kkk ++){
	    
	    //	    pat::PackedCandidate _pf = pfcollection[kkk];
	    if(!aux.basicPFcut(_pf)) continue;
	    Float_t precut_dz = _pf.vz() - closestVertex.position().z();
	    if(TMath::Abs(precut_dz) > c_dz) return false;

	    
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







  /********************************************************************
   *
   * Adding more attributes for each PF candidate
   *
   ********************************************************************/

  event.getByToken( packedpfcandidatesToken_               , packedpfcandidates_      ); 

  std::vector<pfcand_struct> reducedpfcands;
  //  std::vector<attribute> attributes;

  
  for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){   
      
    pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
    
//    if(pf.pt() < 0.5) continue;
//    if(!pf.hasTrackDetails()) continue;
//    
//    // use the PF candidates that come from closestVertex    
//    Float_t precut_dz = pf.vz() - closestVertex.position().z();
//    if(TMath::Abs(precut_dz) > c_dz) continue;
//    
//    Bool_t hpflag = pf.trackHighPurity();
//    if(!hpflag) continue;
//    if(pf.pseudoTrack().hitPattern().numberOfValidPixelHits() < 0) continue;
//    if(pf.pseudoTrack().hitPattern().numberOfValidHits() < 3) continue;
//    if(pf.pseudoTrack().normalizedChi2() > 100) continue;
//    
//    if(TMath::Abs(pf.pdgId())!=211) continue; 
//    if(TMath::Abs(pf.eta()) > 2.5) continue; 
    
    if(!aux.basicPFcut(pf)) continue;
    Float_t precut_dz = pf.vz() - closestVertex.position().z();
    if(TMath::Abs(precut_dz) > c_dz) return false;

  
    reco::TransientTrack  _track = (*builder).build(pf.pseudoTrack());
    TrajectoryStateOnSurface _tsos_pf = extrapolator.extrapolate(_track.impactPointState(), jpsi_vertex->position());

    TrajectoryStateOnSurface _tsost_pf = extrapolatort.extrapolate(_track.impactPointState(), jpsi_vertex->position());
    
    
    std::pair<bool,Measurement1D> _cur3DIP_pf = aux.signedImpactParameter3D(_tsos_pf, jpsi_vertex, closestVertex);

    std::pair<bool,Measurement1D> _cur3DIPt_pf = aux.signedTransverseImpactParameter(_tsost_pf, jpsi_vertex, closestVertex);
    
    Float_t doca3d = _cur3DIP_pf.second.value();
    Float_t doca3de = _cur3DIP_pf.second.error();
    Float_t doca3ds = _cur3DIP_pf.second.significance();

    Float_t doca2d = _cur3DIPt_pf.second.value();
    Float_t doca2de = _cur3DIPt_pf.second.error();
    Float_t doca2ds = _cur3DIPt_pf.second.significance();

    //    Float_t doca1d = TMath::Sqrt(doca3d*doca3d - doca2d*doca2d);
    //    Float_t doca1de = (doca3d/doca1d)*TMath::Sqrt(doca3de*doca3de + doca2d*doca2d*doca2de*doca2de/(doca3d*doca3d));
    //    Float_t doca1ds = doca1d/doca1de;
    

    Bool_t _isRight = false;

    if(runOnMC_){
      
      for(unsigned int mmm=0; mmm < gps.size(); mmm++){       

	std::vector<TLorentzVector> tlvs = gps[mmm];
	    
	for(unsigned int nnn=0; nnn < tlvs.size(); nnn++){

	  if(
	     reco::deltaR(pf.eta(), pf.phi(), tlvs[nnn].Eta(), tlvs[nnn].Phi()) < 0.015 &&
	     pf.pt()/tlvs[nnn].Pt() > 0.85 && 
	     pf.pt()/tlvs[nnn].Pt() < 1.15
	     ){
	    _isRight = true; 
	  }
	}
      }
    }

    
    nBranches_->JpsiTau_st_isRight.push_back(_isRight);
    nBranches_->JpsiTau_st_doca3d.push_back(doca3d);
    nBranches_->JpsiTau_st_doca2d.push_back(doca2d);
    //    nBranches_->JpsiTau_st_doca1d.push_back(doca1d);
    nBranches_->JpsiTau_st_doca3de.push_back(doca3de);
    nBranches_->JpsiTau_st_doca2de.push_back(doca2de);
    //    nBranches_->JpsiTau_st_doca1de.push_back(doca1de);
    nBranches_->JpsiTau_st_doca3ds.push_back(doca3ds);
    nBranches_->JpsiTau_st_doca2ds.push_back(doca2ds);
    //    nBranches_->JpsiTau_st_doca1ds.push_back(doca1ds);
    nBranches_->JpsiTau_st_dz.push_back(pf.vz() - closestVertex.position().z());
    nBranches_->JpsiTau_st_isAssociate.push_back((bool)(pf.vertexRef()->z()==closestVertex.position().z()));
    nBranches_->JpsiTau_st_pvAssociationQuality.push_back(pf.pvAssociationQuality());
    nBranches_->JpsiTau_st_pt.push_back(pf.pt());
    nBranches_->JpsiTau_st_eta.push_back(pf.eta());
    nBranches_->JpsiTau_st_phi.push_back(pf.phi());
    nBranches_->JpsiTau_st_charge.push_back(pf.charge());
    nBranches_->JpsiTau_st_mass.push_back(pf.mass());

    bool flag_match = false;
    bool flag_signal = false;
    int matched_pdg = -999;
    int matched_ppdg = -999;
    int matched_nprong = -999;

    if(runOnMC_){
      for(int imatch = 0; imatch < (int)match_eta.size(); imatch++){
	float _dr_ = reco::deltaR(pf.eta(), pf.phi(), match_eta[imatch], match_phi[imatch]);
	if(_dr_ < 0.015 && 
	   pf.pt()/match_pt[imatch] > 0.85 &&
	   pf.pt()/match_pt[imatch] < 1.15
	   ){
	  flag_match = true;
	  matched_pdg = match_pdg[imatch];
	  matched_ppdg = match_ppdg[imatch];
	  flag_signal = match_isSignal[imatch];
	  matched_nprong = match_nprong[imatch];
	}
      }
    }
    
    nBranches_->JpsiTau_st_isBdecay.push_back(flag_match);
    nBranches_->JpsiTau_st_isBdecaypdg.push_back(matched_pdg);
    nBranches_->JpsiTau_st_isBdecayppdg.push_back(matched_ppdg);
    nBranches_->JpsiTau_st_isSignal.push_back(flag_signal);
    nBranches_->JpsiTau_st_nprong.push_back(matched_nprong);
	

    Float_t near_dz = -999;
    for( reco::VertexCollection::const_iterator vtx = vertices_->begin(); vtx != vertices_->end(); ++vtx){
      
      if(pf.vertexRef()->z()==vtx->position().z()){
	near_dz = closestVertex.position().z() - vtx->position().z();
      }
    }


    //    nBranches_->JpsiTau_st_near_dz.push_back(near_dz);
    Bool_t isAssociate = (bool)(pf.vertexRef()->z()==closestVertex.position().z());

    attribute attr = {
      (Float_t) doca3d,
      (Float_t) doca3de,
      (Float_t) doca3ds,
      (Float_t) doca2d,
      (Float_t) doca2de,
      (Float_t) doca2ds,
      //      (Float_t) doca1d,
      //      (Float_t) doca1de,
      //      (Float_t) doca1ds,
      (Bool_t) _isRight, 
      (Float_t) precut_dz,
      (Bool_t) isAssociate,
      (Int_t) pf.pvAssociationQuality(),
      (Float_t) pf.pt(),
      (Float_t) pf.eta(),
      (Float_t) pf.phi(),
      (Int_t) pf.charge(),
      (Float_t) pf.mass(),
      (Bool_t) flag_match,
      (Int_t) matched_pdg,
      (Int_t) matched_ppdg,
      (Bool_t) flag_signal,
      (Int_t) matched_nprong,
      (Float_t) near_dz
    };

    
    //    attributes.push_back(attr);
    
    reco::TransientTrack  tt_track = (*builder).build(pf.pseudoTrack());
    //    mytracks.push_back(tt_track);

    
    pfcand_struct _cand_ = {
      (Int_t)ii,
      //      (Float_t) TMath::Abs(precut_dz),
      pf,
      tt_track,
      attr
    };
    
    reducedpfcands.push_back(_cand_);

    


//    if(_isRight){
//      doca3d_vec.push_back(doca3d);
//      doca2d_vec.push_back(doca2d);
//      doca1d_vec.push_back(doca1d);
//
//      doca3ds_vec.push_back(doca3ds);
//      doca2ds_vec.push_back(doca2ds);
//      doca1ds_vec.push_back(doca1ds);
//
//      doca3de_vec.push_back(doca3de);
//      doca2de_vec.push_back(doca2de);
//      doca1de_vec.push_back(doca1de);
//
//      dz_vec.push_back(pf.vz() - closestVertex.position().z());
//    }
//    if(_pvip_pf < 0.03){
//      npf_pvip++;
//      if(pf.vertexRef()->z()==closestVertex.position().z()) npf_pvip_dz0 ++;
//
//      Float_t precut_dz = pf.vz() - closestVertex.position().z();
//      if(TMath::Abs(precut_dz) < c_dz){
//	npf_pvip_dz2p5 ++;
//      }
//
//    }


  }

  // sorted by dz
  sort(reducedpfcands.begin(), reducedpfcands.end());

  
////  if(doca3d_vec.size()==3){
////
////    Float_t max_doca3d = -1;
////    for(int ii=0; ii < (int)doca3d_vec.size(); ii++){
////      if(max_doca3d < doca3d_vec[ii]){
////	max_doca3d = doca3d_vec[ii];
////      }
////    }
////
////    Float_t max_doca2d = -1;
////    for(int ii=0; ii < (int)doca2d_vec.size(); ii++){
////      if(max_doca2d < doca2d_vec[ii]){
////	max_doca2d = doca2d_vec[ii];
////      }
////    }
////
////    Float_t max_doca1d = -1;
////    for(int ii=0; ii < (int)doca1d_vec.size(); ii++){
////      if(max_doca1d < doca1d_vec[ii]){
////	max_doca1d = doca1d_vec[ii];
////      }
////    }
////
////
////
////    Float_t max_doca3ds = -1;
////    for(int ii=0; ii < (int)doca3ds_vec.size(); ii++){
////      if(max_doca3ds < doca3ds_vec[ii]){
////	max_doca3ds = doca3ds_vec[ii];
////      }
////    }
////
////    Float_t max_doca2ds = -1;
////    for(int ii=0; ii < (int)doca2ds_vec.size(); ii++){
////      if(max_doca2ds < doca2ds_vec[ii]){
////	max_doca2ds = doca2ds_vec[ii];
////      }
////    }
////
////    Float_t max_doca1ds = -1;
////    for(int ii=0; ii < (int)doca1ds_vec.size(); ii++){
////      if(max_doca1ds < doca1ds_vec[ii]){
////	max_doca1ds = doca1ds_vec[ii];
////      }
////    }
////
////
////    Float_t max_doca3de = -1;
////    for(int ii=0; ii < (int)doca3de_vec.size(); ii++){
////      if(max_doca3de < doca3de_vec[ii]){
////	max_doca3de = doca3de_vec[ii];
////      }
////    }
////
////    Float_t max_doca2de = -1;
////    for(int ii=0; ii < (int)doca2de_vec.size(); ii++){
////      if(max_doca2de < doca2de_vec[ii]){
////	max_doca2de = doca2de_vec[ii];
////      }
////    }
////
////    Float_t max_doca1de = -1;
////    for(int ii=0; ii < (int)doca1de_vec.size(); ii++){
////      if(max_doca1de < doca1de_vec[ii]){
////	max_doca1de = doca1de_vec[ii];
////      }
////    }
////
////
////
////    Float_t max_dz = -1;
////    Float_t max_dz_sign = -1;
////
////    for(int ii=0; ii < (int)dz_vec.size(); ii++){
////      if(max_dz < TMath::Abs(dz_vec[ii])){
////	max_dz = TMath::Abs(dz_vec[ii]);
////	max_dz_sign = dz_vec[ii];
////      }
////    }
////    //    std::vector<int>::iterator _iter = std::max_element(pvips_vec.begin(), pvips_vec.end());
////    //    size_t _index = std::distance(pvips_vec.begin(), _iter);
////    //    std::cout << *_iter << " " << pvips_vec[_index] << " (" << pvips_vec[0] << " " << pvips_vec[1] << " " << pvips_vec[2] << std::endl;
////    nBranches_->JpsiTau_st_doca3d_max.push_back(max_doca3d);
////    nBranches_->JpsiTau_st_doca2d_max.push_back(max_doca2d);
////    nBranches_->JpsiTau_st_doca1d_max.push_back(max_doca1d);
////    nBranches_->JpsiTau_st_doca3ds_max.push_back(max_doca3ds);
////    nBranches_->JpsiTau_st_doca2ds_max.push_back(max_doca2ds);
////    nBranches_->JpsiTau_st_doca1ds_max.push_back(max_doca1ds);
////    nBranches_->JpsiTau_st_doca3de_max.push_back(max_doca3de);
////    nBranches_->JpsiTau_st_doca2de_max.push_back(max_doca2de);
////    nBranches_->JpsiTau_st_doca1de_max.push_back(max_doca1de);
////    nBranches_->JpsiTau_st_dz_max.push_back(max_dz_sign);
////  }

  nBranches_->JpsiTau_st_nch = reducedpfcands.size();









  /********************************************************************
   *
   * Adding more attributes for MC!
   *
   ********************************************************************/



//  std::vector<attribute> attributes;
//
//
//  for(int iatr = 0; iatr < numOfch; iatr ++){
//
//    pat::PackedCandidate pf = pfcollection[iatr];
//
//    reco::TransientTrack  _track = (*builder).build(pf.pseudoTrack());
//    TrajectoryStateOnSurface _tsos_pf = extrapolator.extrapolate(_track.impactPointState(), jpsi_vertex->position());
//
//    TrajectoryStateOnSurface _tsost_pf = extrapolatort.extrapolate(_track.impactPointState(), jpsi_vertex->position());
//    
//    //    VertexDistance3D _a3d_pf;  
//    //    VertexDistanceXY _a3dt_pf;
//    
//    std::pair<bool,Measurement1D> _cur3DIP_pf = aux.signedImpactParameter3D(_tsos_pf, jpsi_vertex, closestVertex);
//    std::pair<bool,Measurement1D> _cur3DIPt_pf = aux.signedTransverseImpactParameter(_tsost_pf, jpsi_vertex, closestVertex);
//    
//    Float_t doca3d = _cur3DIP_pf.second.value();
//    Float_t doca3de = _cur3DIP_pf.second.error();
//    Float_t doca3ds = _cur3DIP_pf.second.significance();
//
//    Float_t doca2d = _cur3DIPt_pf.second.value();
//    Float_t doca2de = _cur3DIPt_pf.second.error();
//    Float_t doca2ds = _cur3DIPt_pf.second.significance();
//
//    Float_t doca1d = TMath::Sqrt(doca3d*doca3d - doca2d*doca2d);
//    Float_t doca1de = (doca3d/doca1d)*TMath::Sqrt(doca3de*doca3de + doca2d*doca2d*doca2de*doca2de/(doca3d*doca3d));
//    Float_t doca1ds = doca1d/doca1de;
//    
//
//    Bool_t _isRight = false;
//
//    if(runOnMC_){
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
//	    _isRight = true; 
//	  }
//	}
//      }
//    }
//
//    
//
//    bool flag_match = false;
//    bool flag_signal = false;
//    int matched_pdg = -999;
//    int matched_ppdg = -999;
//    int matched_nprong = -999;
//
//    if(runOnMC_){
//      for(int imatch = 0; imatch < (int)match_eta.size(); imatch++){
//	float _dr_ = reco::deltaR(pf.eta(), pf.phi(), match_eta[imatch], match_phi[imatch]);
//	if(_dr_ < 0.015 && 
//	   pf.pt()/match_pt[imatch] > 0.85 &&
//	   pf.pt()/match_pt[imatch] < 1.15
//	   ){
//	  flag_match = true;
//	  matched_pdg = match_pdg[imatch];
//	  matched_ppdg = match_ppdg[imatch];
//	  flag_signal = match_isSignal[imatch];
//	  matched_nprong = match_nprong[imatch];
//	}
//      }
//    }
//
//
//
//    Float_t near_dz = -999;
//    for( reco::VertexCollection::const_iterator vtx = vertices_->begin(); vtx != vertices_->end(); ++vtx){
//      
//      if(pf.vertexRef()->z()==vtx->position().z()){
//	near_dz = closestVertex.position().z() - vtx->position().z();
//      }
//    }
//
//    Float_t dz = pf.vz() - closestVertex.position().z();
//    Bool_t isAssociate = (bool)(pf.vertexRef()->z()==closestVertex.position().z());
//
//    attribute attr = {
//      (Float_t) doca3d,
//      (Float_t) doca3de,
//      (Float_t) doca3ds,
//      (Float_t) doca2d,
//      (Float_t) doca2de,
//      (Float_t) doca2ds,
//      (Float_t) doca1d,
//      (Float_t) doca1de,
//      (Float_t) doca1ds,
//      (Bool_t) _isRight, 
//      (Float_t) dz,
//      (Bool_t) isAssociate,
//      (Int_t) pf.pvAssociationQuality(),
//      (Float_t) pf.pt(),
//      (Float_t) pf.eta(),
//      (Float_t) pf.phi(),
//      (Int_t) pf.charge(),
//      (Float_t) pf.mass(),
//      (Bool_t) flag_match,
//      (Int_t) matched_pdg,
//      (Int_t) matched_ppdg,
//      (Bool_t) flag_signal,
//      (Int_t) matched_nprong,
//      (Float_t) near_dz
//    };
//
//    
//    attributes.push_back(attr);
//  }
//
//
//  std::cout << numOfch << " " << attributes.size() << " must be the same !!!!!!!!!!!!!!!!!" << std::endl;




















  /********************************************************************
   *
   * Step6: Tau selection
   *        Just select highest in pT but there might be better selection ... 
   *
   ********************************************************************/
    
    
  //  std::vector<pat::PackedCandidate> pfcollection; 
  //  std::vector<Int_t> pfcollection_id; 
  //  std::vector<Int_t> pfcollection_gen; 
  //  std::vector<reco::TransientTrack> mytracks;
  //  std::vector<Int_t> mypvassociation;
  //  std::vector<Float_t> mydoca;
  std::vector<Float_t> mydnn;
  std::vector<Float_t> mydnn_old;
    
  Int_t npf_before_dnn = 0;
  Int_t npf_qr = 0;

  if(useDNN_){

    //    std::vector<pfcand> pfcands;

/////    Int_t count_dnn = 0;
/////    Int_t count_dnn_muon = 0;
/////    
/////    for(int imu = 0; imu < (int)muoncollection.size(); imu++){
/////      
/////      if(count_dnn < numberofDNN){
/////	data.tensor<float, 3>()(0, count_dnn, 0) = mu1_fit.eta();
/////	data.tensor<float, 3>()(0, count_dnn, 1) = mu1_fit.phi();
/////	data.tensor<float, 3>()(0, count_dnn, 2) = TMath::Log(mu1_fit.pt());
/////	data.tensor<float, 3>()(0, count_dnn, 3) = muoncollection[mcidx_mu1].charge();
/////	data.tensor<float, 3>()(0, count_dnn, 4) = 0;
/////	data.tensor<float, 3>()(0, count_dnn, 5) = 0;
/////	data.tensor<float, 3>()(0, count_dnn, 6) = 0;
/////	data.tensor<float, 3>()(0, count_dnn, 7) = 0;
/////	data.tensor<float, 3>()(0, count_dnn, 8) = 0;
/////	data.tensor<float, 3>()(0, count_dnn, 9) = 0;
/////	data.tensor<float, 3>()(0, count_dnn, 10) = 0;
/////	data.tensor<float, 3>()(0, count_dnn, 11) = 0;
/////	data.tensor<float, 3>()(0, count_dnn, 12) = 1;
/////	
/////	label_perPF.matrix<int>()(0, count_dnn) = 0;
/////	label_perEVT.matrix<int>()(0, count_dnn) = 0;
/////	//	norm.matrix<float>()(0, count_dnn) = float(1);
/////	
/////
/////
/////	count_dnn_muon++;
/////	count_dnn++;
/////      }
/////    }
/////
/////
/////
/////
/////
//////////    for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){   
//////////	
//////////      pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
//////////	
//////////	
//////////      if(!pf.hasTrackDetails()) continue;
//////////      Float_t precut_dz = pf.vz() - closestVertex.position().z();
//////////      if(TMath::Abs(precut_dz) > c_dz) continue;
//////////      //      if(pf.vertexRef()->z()!=closestVertex.position().z()) continue;
//////////
//////////      npf_qr++;
//////////	
//////////      if(pf.pt() < 0.5) continue;
//////////	
//////////      Bool_t hpflag = pf.trackHighPurity();
//////////      if(!hpflag) continue;
//////////      if(pf.pseudoTrack().hitPattern().numberOfValidPixelHits() < 0) continue;
//////////      if(pf.pseudoTrack().hitPattern().numberOfValidHits() < 3) continue;
//////////      if(pf.pseudoTrack().normalizedChi2() > 100) continue;
//////////	
//////////      if(TMath::Abs(pf.pdgId())!=211) continue; 
//////////      if(TMath::Abs(pf.eta()) > 2.5) continue; 
//////////
//////////
//////////
//////////
//////////      ////////// prefiltering by the distance between J/psi and the tau candidate
//////////
//////////
/////////////      reco::TransientTrack  _track = (*builder).build(pf.pseudoTrack());
/////////////      TrajectoryStateOnSurface _tsos_pf = extrapolator.extrapolate(_track.impactPointState(), jpsi_vertex->position());
/////////////    
/////////////      //      VertexDistance3D _a3d_pf;  
/////////////      
/////////////      std::pair<bool,Measurement1D> _cur3DIP_pf = aux.absoluteImpactParameter3D(_tsos_pf, jpsi_vertex);
/////////////      
/////////////      Float_t _pvip_pf = _cur3DIP_pf.second.value();
/////////////    
/////////////      
/////////////      if(_pvip_pf > 0.03) continue;
//////////
//////////
//////////      npf_before_dnn++;	
//////////
//////////      ////////// prefiltering ////////////////////////////////////////////////////////
//////////
//////////
//////////      pfcand _cand_ = {
//////////	(Int_t)ii,
//////////	(Float_t) _pvip_pf,
//////////	(Float_t) abs(precut_dz)
//////////      };
//////////	  
//////////      pfcands.push_back(_cand_);
//////////    }
//////////
//////////    //    if(pfcands.size()>30) return false;
//////////
//////////
//////////    //sorting by distance to the vertex
//////////    sort(pfcands.begin(), pfcands.end());
/////
/////
/////    for(size_t ic = 0; ic < reducedpfcands.size(); ic++){
/////      //      Int_t idx = reducedpfcands[ic].cand_idx;
/////
/////      pat::PackedCandidate pf = reducedpfcands[ic].pfcand; //(*packedpfcandidates_)[idx];
/////      attribute attr = reducedpfcands[ic].pfaux;
/////
/////      if(count_dnn < numberofDNN){
/////	data.tensor<float, 3>()(0, count_dnn, 0) = pf.eta();
/////	data.tensor<float, 3>()(0, count_dnn, 1) = pf.phi();
/////	data.tensor<float, 3>()(0, count_dnn, 2) = TMath::Log(pf.pt());
/////	data.tensor<float, 3>()(0, count_dnn, 3) = pf.charge();
/////	data.tensor<float, 3>()(0, count_dnn, 4) = attr.pvAssociationQuality;
/////	data.tensor<float, 3>()(0, count_dnn, 5) = attr.doca3d;
/////	data.tensor<float, 3>()(0, count_dnn, 6) = attr.doca2d;
/////	data.tensor<float, 3>()(0, count_dnn, 7) = attr.doca3de;
/////	data.tensor<float, 3>()(0, count_dnn, 8) = attr.doca2de;
/////	data.tensor<float, 3>()(0, count_dnn, 9) = attr.dz;
/////	data.tensor<float, 3>()(0, count_dnn, 10) = attr.isAssociate;
/////	data.tensor<float, 3>()(0, count_dnn, 11) = attr.near_dz;
/////	data.tensor<float, 3>()(0, count_dnn, 12) = 0;
/////	
/////	label_perPF.matrix<int>()(0, count_dnn) = 0;
/////	label_perEVT.matrix<int>()(0, count_dnn) = 0;
/////	//	norm.matrix<float>()(0, count_dnn) = float(1);
/////	count_dnn++;
/////
/////
/////
/////
/////	
/////	//	pfcollection.push_back(pf);
/////	//	pfcollection_id.push_back(idx);
/////	//	reco::TransientTrack  tt_track = (*builder).build(pf.pseudoTrack());
/////	//	mytracks.push_back(tt_track);
/////	//	mypvassociation.push_back(pf.pvAssociationQuality());
/////	//	mydoca.push_back(pfcands[ic].doca);
/////      }
/////    }
/////    
/////      
/////    for(int ic=count_dnn; ic < numberofDNN; ic++){
/////      
/////      data.tensor<float, 3>()(0, ic, 0) = 0;
/////      data.tensor<float, 3>()(0, ic, 1) = 0;
/////      data.tensor<float, 3>()(0, ic, 2) = 0;
/////      data.tensor<float, 3>()(0, ic, 3) = 0;
/////      data.tensor<float, 3>()(0, ic, 4) = 0;
/////      data.tensor<float, 3>()(0, ic, 5) = 0;
/////      data.tensor<float, 3>()(0, ic, 6) = 0;
/////      data.tensor<float, 3>()(0, ic, 7) = 0;
/////      data.tensor<float, 3>()(0, ic, 8) = 0;
/////      data.tensor<float, 3>()(0, ic, 9) = 0;
/////      data.tensor<float, 3>()(0, ic, 10) = 0;
/////      data.tensor<float, 3>()(0, ic, 11) = 0;
/////      data.tensor<float, 3>()(0, ic, 12) = 0;
/////	
/////      label_perPF.matrix<int>()(0, ic) = 0;
/////      label_perEVT.matrix<int>()(0, ic) = 0;
/////      //      norm.matrix<float>()(0, ic) = float(1);
/////
/////    }
/////
/////
/////    isTraining.scalar<bool>()() = false;
/////    
/////    std::vector<tensorflow::Tensor> outputs_perPF;
/////    std::vector<tensorflow::Tensor> outputs_perEVT;
/////    //    tensorflow::run(session, {  { "Placeholder:0", data },  { "Placeholder_1:0", label }, { "Placeholder_2:0", add_global } , {"Placeholder_3:0", isTraining}, {"Placeholder_4:0", norm}}, { "Reshape_13:0" }, &outputs);
/////    tensorflow::run(session_perPF, {  { "Placeholder:0", data },  { "Placeholder_1:0", label_perPF }, { "Placeholder_2:0", isTraining}} , { "Softmax_2:0" }, &outputs_perPF);
/////		    
/////    
/////    tensorflow::run(session_perEVT, {  { "Placeholder:0", data },  { "Placeholder_1:0", label_perEVT }, { "Placeholder_2:0", isTraining}} , { "Softmax_2:0" }, &outputs_perEVT);
/////      
/////    auto finalOutputTensor_perPF = outputs_perPF[0].tensor<float, 3>();
/////    auto finalOutputTensor_perEVT = outputs_perEVT[0].tensor<float, 3>();
/////
/////    for(int ic=count_dnn_muon; ic<count_dnn; ic++){
/////      mydnn.push_back(finalOutputTensor_perPF(0, ic, 1));
/////      std::cout <<  "pefPF DNN = " << ic << " " << finalOutputTensor_perPF(0, ic, 1) << std::endl;
/////    }
/////
/////    Float_t evtDNN = finalOutputTensor_perEVT(0, 0, 1);
/////    std::cout << "evtDNN = " << evtDNN << std::endl;

    
    
    //////////////// old 

    ////////////// old 

    Int_t count_dnn_old = 0;
    Int_t count_dnn_muon_old = 0;

    for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){   
      
      pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
      
      if(TMath::Abs(pf.eta()) < 2.4 && 
	 TMath::Abs(pf.charge())==1 &&
	 pf.pt() > 4. &&
	 pf.isGlobalMuon() > 0.5 &&
	 pf.hasTrackDetails() > 0.5
	 ){
	

	if(count_dnn_old < 50){
	  data_old.tensor<float, 3>()(0, count_dnn_old, 0) = pf.eta();
	  data_old.tensor<float, 3>()(0, count_dnn_old, 1) = pf.phi();
	  data_old.tensor<float, 3>()(0, count_dnn_old, 2) = TMath::Log(pf.pt());
	  data_old.tensor<float, 3>()(0, count_dnn_old, 3) = TMath::Log(pf.energy());
	  data_old.tensor<float, 3>()(0, count_dnn_old, 4) = pf.charge();
	  data_old.tensor<float, 3>()(0, count_dnn_old, 5) = TMath::Abs(closestVertex.position().z() - pf.pseudoTrack().vz());
	  data_old.tensor<float, 3>()(0, count_dnn_old, 6) = TMath::Sqrt( TMath::Power((closestVertex.position().x() - pf.pseudoTrack().vx()), 2) + TMath::Power((closestVertex.position().y() - pf.pseudoTrack().vy()), 2));
	  data_old.tensor<float, 3>()(0, count_dnn_old, 7) = pf.isGlobalMuon();
	    
	  label_old.matrix<int>()(0, count_dnn_old) = 0;
	  norm_old.matrix<float>()(0, count_dnn_old) = float(1);

	  count_dnn_muon_old++;
	  count_dnn_old++;
	}
      }
    }

    for(size_t ic = 0; ic < reducedpfcands.size(); ic++){
      //      Int_t idx = reducedpfcands[ic].cand_idx;

      pat::PackedCandidate pf = reducedpfcands[ic].pfcand; //(*packedpfcandidates_)[idx];
      //      attribute attr = reducedpfcands[ic].pfaux;


      if(count_dnn_old < 50){
	data_old.tensor<float, 3>()(0, count_dnn_old, 0) = pf.eta();
	data_old.tensor<float, 3>()(0, count_dnn_old, 1) = pf.phi();
	data_old.tensor<float, 3>()(0, count_dnn_old, 2) = TMath::Log(pf.pt());
	data_old.tensor<float, 3>()(0, count_dnn_old, 3) = TMath::Log(pf.energy());
	data_old.tensor<float, 3>()(0, count_dnn_old, 4) = pf.charge();
	data_old.tensor<float, 3>()(0, count_dnn_old, 5) = TMath::Abs(closestVertex.position().z() - pf.pseudoTrack().vz());
	data_old.tensor<float, 3>()(0, count_dnn_old, 6) = TMath::Sqrt( TMath::Power((closestVertex.position().x() - pf.pseudoTrack().vx()), 2) + TMath::Power((closestVertex.position().y() - pf.pseudoTrack().vy()), 2));
	data_old.tensor<float, 3>()(0, count_dnn_old, 7) = pf.isGlobalMuon();
	  
	label_old.matrix<int>()(0, count_dnn_old) = 0;
	norm_old.matrix<float>()(0, count_dnn_old) = float(1);
	count_dnn_old++;

      }
    }

      
    for(int ic=count_dnn_old; ic<50; ic++){

      data_old.tensor<float, 3>()(0, ic, 0) = 0;
      data_old.tensor<float, 3>()(0, ic, 1) = 0;
      data_old.tensor<float, 3>()(0, ic, 2) = 0;
      data_old.tensor<float, 3>()(0, ic, 3) = 0;
      data_old.tensor<float, 3>()(0, ic, 4) = 0;
      data_old.tensor<float, 3>()(0, ic, 5) = 0;
      data_old.tensor<float, 3>()(0, ic, 6) = 0;
      data_old.tensor<float, 3>()(0, ic, 7) = 0;
	
      label_old.matrix<int>()(0, ic) = 0;
      norm_old.matrix<float>()(0, ic) = float(1);

    }

    add_global_old.matrix<float>()(0, 0) = float(count_dnn_muon_old); // Number of muons around 0.5 cm from PV
    add_global_old.matrix<float>()(0, 1) = float(count_dnn_old/100); //Number of charged pf candidates around 0.5 cm from PV
    isTraining_old.scalar<bool>()() = false; //Number of charged pf candidates around 0.5 cm from PV


    std::vector<tensorflow::Tensor> outputs_old;
    tensorflow::run(session_old, {  { "Placeholder:0", data_old },  { "Placeholder_1:0", label_old }, { "Placeholder_2:0", add_global_old } , {"Placeholder_3:0", isTraining_old}, {"Placeholder_4:0", norm_old}}, { "Reshape_13:0" }, &outputs_old);
      
    auto finalOutputTensor = outputs_old[0].tensor<float, 3>();

    for(int ic=count_dnn_muon_old; ic<count_dnn_old; ic++){

      mydnn.push_back(finalOutputTensor(0, ic, 1));
    }






    
    
  }
  
/////else{
/////
/////
/////    for( size_t ii = 0; ii < packedpfcandidates_->size(); ++ii ){   
/////      
/////      pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
/////	
/////      if(pf.pt() < 0.5) continue;
/////      if(!pf.hasTrackDetails()) continue;
/////	
/////      // use the PF candidates that come from closestVertex
/////      //      if(pf.vertexRef()->z()!=closestVertex.position().z()) continue;
/////	
/////      Float_t precut_dz = pf.vz() - closestVertex.position().z();
/////      if(TMath::Abs(precut_dz) > c_dz) continue;
/////	
/////      Bool_t hpflag = pf.trackHighPurity();
/////      if(!hpflag) continue;
/////      if(pf.pseudoTrack().hitPattern().numberOfValidPixelHits() < 0) continue;
/////      if(pf.pseudoTrack().hitPattern().numberOfValidHits() < 3) continue;
/////      if(pf.pseudoTrack().normalizedChi2() > 100) continue;
/////	
/////      if(TMath::Abs(pf.pdgId())!=211) continue; 
/////      if(TMath::Abs(pf.eta()) > 2.5) continue; 
/////
/////
/////      ////////// prefiltering by the distance between J/psi and the tau candidate
/////
/////      reco::TransientTrack  _track = (*builder).build(pf.pseudoTrack());
/////      TrajectoryStateOnSurface _tsos_pf = extrapolator.extrapolate(_track.impactPointState(), jpsi_vertex->position());
/////    
///////      VertexDistance3D _a3d_pf;  
/////      
/////      std::pair<bool,Measurement1D> _cur3DIP_pf = aux.absoluteImpactParameter3D(_tsos_pf, jpsi_vertex);
/////      
/////      Float_t _pvip_pf = _cur3DIP_pf.second.value();
/////    
/////      
/////      if(_pvip_pf > 0.03) continue;
/////      
/////      /////////////////////////////////////
/////
/////
/////      pfcollection.push_back(pf);
/////      //      pfcollection_id.push_back(ii);
/////      reco::TransientTrack  tt_track = (*builder).build(pf.pseudoTrack());
/////      mytracks.push_back(tt_track);
/////      mypvassociation.push_back(pf.pvAssociationQuality());
/////      mydoca.push_back(_pvip_pf);
/////    }



  

  Int_t numOfch = (size_t)reducedpfcands.size();

  nBranches_->cutflow_perevt->Fill(9);
  if(numOfch<3) return false;
  nBranches_->cutflow_perevt->Fill(10);




  //////////////////// test PVIP w.r.t J/psi ///////////////////////

//  Int_t npf_pvip = 0;
//  Int_t npf_pvip_dz0 = 0;
//  Int_t npf_pvip_dz2p5 = 0;


 

  ///////////////// check ///////////////////////////



//  for( unsigned p=0; p<genParticles_->size(); ++p ){
//    
//    bool isB( (abs((*genParticles_)[p].pdgId())>=511 && abs((*genParticles_)[p].pdgId())<=545));
//    if(!isB) continue;
//	  
//    Bool_t ismother_B = false;
//	    
//    for( unsigned int m=0; m<(*genParticles_)[p].numberOfMothers(); ++m ){
//	    //	    std::cout << (*genParticles_)[p].mother(m)->pdgId() << " ";
//
//      if( abs((*genParticles_)[p].mother(m)->pdgId())>=511 && abs((*genParticles_)[p].mother(m)->pdgId())<=545 ){
//	ismother_B = true;
//      }
//    }
//    
//    if( (*genParticles_)[p].numberOfDaughters()!=1 && 
//	(*genParticles_)[p].numberOfMothers()==1 && 
//	(*genParticles_)[p].mother(0)->pdgId()==(*genParticles_)[p].pdgId() 
//	){
//      ismother_B = false; 
//    }
//    
//
//    if(! ( (*genParticles_)[p].numberOfDaughters()==1 && (*genParticles_)[p].daughter(0)->pdgId()==(*genParticles_)[p].pdgId() ) && ismother_B==false ){
//
//      int _rank = 1;
//      std::vector<size_t> allIndices;
//      std::vector<int> pdgs;
//      std::vector<int> layers;
//      std::vector<float> ppt;
//      std::vector<float> peta;
//      std::vector<float> pphi;
//      
//      pdgs.push_back((*genParticles_)[p].pdgId());
//      layers.push_back(0);
//      ppt.push_back((*genParticles_)[p].pt());
//      peta.push_back((*genParticles_)[p].eta());
//      pphi.push_back((*genParticles_)[p].phi());
//      
//
//      std::cout << (*genParticles_)[p].pdgId() << " (" << (*genParticles_)[p].pt() << ", " << (*genParticles_)[p].eta() << ", " << (*genParticles_)[p].phi() << ")" << std::endl;
//
//      aux.recursiveDaughters(p, _rank, *genParticles_, allIndices, pdgs, layers, ppt, peta, pphi, verbose_);
//
//      for(int ivec = 0; ivec < (int)peta.size(); ivec++){
//	if(layers[ivec]==-1 && TMath::Abs(pdgs[ivec])==211){
//	  match_eta.push_back(peta[ivec]);
//	  match_phi.push_back(pphi[ivec]);
//	}
//      }
//
//    }
//  }
	



  /********************************************************************
   *
   * For Vini, training end.
   *
   ********************************************************************/
  
  //////////////////// test PVIP w.r.t J/psi ///////////////////////





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
      
    pat::PackedCandidate pf1 = reducedpfcands[iii].pfcand;
    reco::TransientTrack track1 = reducedpfcands[iii].track;

    if(useDNN_==true && mydnn[iii] < c_dnn) continue;
    npf_after_dnn++;

    for(int jjj = iii+1; jjj < numOfch; jjj ++){
	
      pat::PackedCandidate pf2 = reducedpfcands[jjj].pfcand;
      reco::TransientTrack track2 = reducedpfcands[jjj].track;

      if(useDNN_==true && mydnn[jjj] < c_dnn) continue;

      for(int kkk = jjj+1; kkk < numOfch; kkk ++){

	pat::PackedCandidate pf3 = reducedpfcands[kkk].pfcand;
	reco::TransientTrack track3 = reducedpfcands[kkk].track;

	if(useDNN_==true && mydnn[kkk] < c_dnn) continue;

	Int_t tau_charge = pf1.charge() + pf2.charge() + pf3.charge(); 

	if(TMath::Abs(tau_charge)!=(int)c_charge) continue; 


	/* reconstruct taus*/

	std::vector<RefCountedKinematicParticle> tauParticles;

	tauParticles.push_back(pFactory.particle(track1, aux.pion_mass, chi, ndf, aux.pion_sigma));
	tauParticles.push_back(pFactory.particle(track2, aux.pion_mass, chi, ndf, aux.pion_sigma));
	tauParticles.push_back(pFactory.particle(track3, aux.pion_mass, chi, ndf, aux.pion_sigma));

  
///	//reconstructing a tau decay
///	RefCountedKinematicTree tauTree = kpvFitter.fit(tauParticles);
///
///
///	if(tauTree->isEmpty() || !tauTree->isValid() || !tauTree->isConsistent()) continue;
///
///	//getting the J/Psi KinematicParticle
///	tauTree->movePointerToTheTop();
///
///	RefCountedKinematicParticle tau_part = tauTree->currentParticle();
///	if(!tau_part->currentState().isValid()) continue;
///	RefCountedKinematicVertex tau_vertex = tauTree->currentDecayVertex();
///	if(!tau_vertex->vertexIsValid()) continue; 
///
///






	RefCountedKinematicParticle tau_part;
	RefCountedKinematicVertex tau_vertex;
	RefCountedKinematicTree tauTree;
	Bool_t taufit_flag;
	std::tie(taufit_flag, tau_part, tau_vertex, tauTree) = aux.KinematicFit(tauParticles, -1, -1);
	
	if(!taufit_flag) continue;

	if(TMath::Prob(tau_vertex->chiSquared(), tau_vertex->degreesOfFreedom()) <= c_vprob) continue;
	  
	std::vector< RefCountedKinematicParticle > tau_children = tauTree->finalStateParticles();
	  
	math::PtEtaPhiMLorentzVector tau1_fit = aux.daughter_p4(tau_children, 0);
	math::PtEtaPhiMLorentzVector tau2_fit = aux.daughter_p4(tau_children, 1);
	math::PtEtaPhiMLorentzVector tau3_fit = aux.daughter_p4(tau_children, 2);


	particle_cand Taucand; 
	Taucand = aux.calculateIPvariables(extrapolator, tau_part, tau_vertex, closestVertex);

	//	std::cout << iii << " " << jjj << " " << kkk << std::endl;

	if(Taucand.fls3d < c_fsig) continue;

	std::vector<RefCountedKinematicParticle> allParticles;

	allParticles.push_back(pFactory.particle(track1, aux.pion_mass, chi, ndf, aux.pion_sigma));
	allParticles.push_back(pFactory.particle(track2, aux.pion_mass, chi, ndf, aux.pion_sigma));
	allParticles.push_back(pFactory.particle(track3, aux.pion_mass, chi, ndf, aux.pion_sigma));
	allParticles.push_back(jpsi_part);


//	RefCountedKinematicTree bcTree = kpvFitter.fit(allParticles);
//
//	if(bcTree->isEmpty() || !bcTree->isValid() || !bcTree->isConsistent()) continue;
//	  
//
//	RefCountedKinematicParticle bc_part = bcTree->currentParticle();
//	if(!bc_part->currentState().isValid()) continue;
//
//	RefCountedKinematicVertex bc_vertex = bcTree->currentDecayVertex();
//	if(!bc_vertex->vertexIsValid()) continue;


	RefCountedKinematicParticle bc_part;
	RefCountedKinematicVertex bc_vertex;
	RefCountedKinematicTree bcTree;
	Bool_t bcfit_flag;
	std::tie(bcfit_flag, bc_part, bc_vertex, bcTree) = aux.KinematicFit(allParticles, -1, -1);
	
	if(!bcfit_flag) continue;




 
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
	
	std::vector<float> iso = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
	std::vector<int> ntracks = {0,0,0,0,0,0,0,0,0, 0};
	std::vector<float> iso_mindoca = {999, 999, 999, 999, 999, 999, 999, 999, 999, 999};
	
	for(int itrk = 0; itrk < numOfch;  itrk++){
    
	  if(itrk==iii || itrk==jjj || itrk==kkk) continue;

	  Float_t iso_pt = reducedpfcands[itrk].pfcand.pt();
	  Float_t iso_eta = reducedpfcands[itrk].pfcand.eta();
	  Float_t iso_phi = reducedpfcands[itrk].pfcand.phi();

	  Float_t iso_dr = reco::deltaR(iso_eta, iso_phi, bc_part->currentState().globalMomentum().eta(), bc_part->currentState().globalMomentum().phi());

	  if(iso_dr > 0.4) continue;
	  
	  for(unsigned int iiso=0; iiso<iso.size(); iiso++){
	    Float_t iso_threshold = 0.5 + iiso*0.1;
	    if(iso_pt > iso_threshold){
	      iso[iiso] += iso_pt;
	      ntracks[iiso] += 1;

	      if(iso_mindoca[iiso] > reducedpfcands[itrk].pfaux.doca3d){
		iso_mindoca[iiso] = reducedpfcands[itrk].pfaux.doca3d;
	      }
	    }
	  }
	  

	  //	  TrajectoryStateOnSurface tsos_pf = extrapolator.extrapolate(mytracks[itrk].impactPointState(), jpsi_vertex->position());
    
    
	  //	  VertexDistance3D a3d_pf;  

	  //	  std::pair<bool,Measurement1D> cur3DIP_pf = aux.absoluteImpactParameter3D(tsos_pf, jpsi_vertex, a3d_pf);

	  //	  Float_t pvip_pf = cur3DIP_pf.second.value();
    
	  //	  if(pvip_pf < 0.03) ntracks+=1;

        }




//	Float_t iso_nocut = 0;
//	Int_t ntracks_nocut = 0;
//	Float_t iso_mindoca_nocut = 999; 
//
//	
//	for( int ii = 0; ii < (int)packedpfcandidates_->size(); ++ii ){   
//	  
//	  pat::PackedCandidate pf = (*packedpfcandidates_)[ii];
//	  
//	  //if(pf.pt() < 0.5) continue;
//	  if(!pf.hasTrackDetails()) continue;
//	  
//	  Bool_t hpflag = pf.trackHighPurity();
//	  if(!hpflag) continue;
//	  if(pf.pseudoTrack().hitPattern().numberOfValidPixelHits() < 0) continue;
//	  if(pf.pseudoTrack().hitPattern().numberOfValidHits() < 3) continue;
//	  if(pf.pseudoTrack().normalizedChi2() > 100) continue;
//	
//	  if(TMath::Abs(pf.pdgId())!=211) continue; 
//	  //	  if(TMath::Abs(pf.eta()) > 2.5) continue; 
//
//
//	  if(pfcollection_id[iii]==ii || 
//	     pfcollection_id[jjj]==ii || 
//	     pfcollection_id[kkk]==ii) continue;
//
//
//	  iso_nocut += pf.pt();
//	  
//	  ////////// prefiltering by the distance between J/psi and the tau candidate
//	  
//	  reco::TransientTrack  _track = (*builder).build(pf.pseudoTrack());
//	  TrajectoryStateOnSurface _tsos_pf = extrapolator.extrapolate(_track.impactPointState(), jpsi_vertex->position());
//    
//	  VertexDistance3D _a3d_pf;  
//      
//	  std::pair<bool,Measurement1D> _cur3DIP_pf = aux.absoluteImpactParameter(_tsos_pf, jpsi_vertex, _a3d_pf);
//	  
//	  Float_t _pvip_pf = _cur3DIP_pf.second.value();
//	  
//	  if(_pvip_pf < 0.03) ntracks_nocut+=1;
//	  if(iso_mindoca_nocut > _pvip_pf) iso_mindoca_nocut = _pvip_pf;
//
//	}





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


	

	Float_t sumofdnn_others = 0;
	for(int lll = 0; lll < numOfch; lll ++){
	  if(lll==iii || lll==jjj || lll==kkk) continue;
	  sumofdnn_others += mydnn[lll];
	}

	taucand _cand_ = {
	  //	  iii,
	  //	  jjj,
	  //	  kkk,
	  //	  (Float_t) tlv_tau_fit.Pt(),
	  //	  (Float_t) tlv_tau_fit.Eta(),
	  //	  (Float_t) tlv_tau_fit.Phi(),
	  //	  (Float_t) tlv_tau_fit.M(),
	  tlv_tau_fit,
	  Taucand,
	  //	  (Float_t) TMath::Prob(tau_vertex->chiSquared(), tau_vertex->degreesOfFreedom()),
	  //	  (Float_t) tau_vertex->vertexState().position().x(), 
	  //	  (Float_t) tau_vertex->vertexState().position().y(), 
	  //	  (Float_t) tau_vertex->vertexState().position().z(), 
	  tau_vertex,
	  (Float_t) max_dr_3prong, 
	  (Int_t) tau_charge,
	  (Bool_t) isRight,
	  (Bool_t) isRight1,
	  (Bool_t) isRight2,
	  (Bool_t) isRight3,
	  (Int_t) pid,
	  (Float_t) matched_gentaupt, 
	  (Float_t) sumofdnn,
	  (Float_t) sumofdnn_others,
	  (Float_t) dnn1,
	  (Float_t) dnn2,
	  (Float_t) dnn3,
	  //	  (Float_t) mydoca[iii],
	  //	  (Float_t) mydoca[jjj],
	  //	  (Float_t) mydoca[kkk],
	  //	  (Int_t) mypvassociation[iii],
	  //	  (Int_t) mypvassociation[jjj],
	  //	  (Int_t) mypvassociation[kkk],
	  //	  (Float_t) TMath::Prob(bc_vertex->chiSquared(), bc_vertex->degreesOfFreedom()),
	  //	  (Float_t) bc_vertex->vertexState().position().x(),
	  //	  (Float_t) bc_vertex->vertexState().position().y(),
	  //	  (Float_t) bc_vertex->vertexState().position().z(),
	  bc_vertex,
	  bc_part, 
	  //	  (Float_t) bc_part->currentState().globalMomentum().perp(),
	  //	  (Float_t) bc_part->currentState().globalMomentum().eta(),
	  //	  (Float_t) bc_part->currentState().globalMomentum().phi(),
	  //	  (Float_t) bc_part->currentState().mass(),
	  Bcand,
	  iso,
	  ntracks,
	  iso_mindoca,
	  //	  (Float_t) iso_nocut,
	  //	  (Float_t) ntracks_nocut,
	  //	  (Float_t) iso_mindoca_nocut,
	  reducedpfcands[iii],
	  reducedpfcands[jjj],
	  reducedpfcands[kkk],
	};
	  
	cands.push_back(_cand_);
      }
    }
  }
    
  sort(cands.begin(), cands.end());


  if(cands.size()==0) return false;

  if(verbose_) std::cout << "[JpsiTauNtuplizer] " << cands.size() << " tau candidates were found" << std::endl;

  nBranches_->cutflow_perevt->Fill(11);

  Int_t ncomb = 0;

  for(int ic=0; ic < (int)cands.size(); ic++){
      
    ncomb += 1;


    /********************************************************************
     *
     * Step9: Filling normal branches
     *
     ********************************************************************/


    nBranches_->JpsiTau_tau_pt.push_back(cands[ic].cand_tlv_tau_fit.Pt());
    nBranches_->JpsiTau_tau_eta.push_back(cands[ic].cand_tlv_tau_fit.Eta());
    nBranches_->JpsiTau_tau_phi.push_back(cands[ic].cand_tlv_tau_fit.Phi());
    nBranches_->JpsiTau_tau_mass.push_back(cands[ic].cand_tlv_tau_fit.M());
    nBranches_->JpsiTau_tau_q.push_back(cands[ic].cand_tau_charge);

	  //	  (Float_t) TMath::Prob(tau_vertex->chiSquared(), tau_vertex->degreesOfFreedom()),
	  //	  (Float_t) tau_vertex->vertexState().position().x(), 
	  //	  (Float_t) tau_vertex->vertexState().position().y(), 
	  //	  (Float_t) tau_vertex->vertexState().position().z(), 


    nBranches_->JpsiTau_tau_vx.push_back(cands[ic].cand_tau_vertex->vertexState().position().x());
    nBranches_->JpsiTau_tau_vy.push_back(cands[ic].cand_tau_vertex->vertexState().position().y());
    nBranches_->JpsiTau_tau_vz.push_back(cands[ic].cand_tau_vertex->vertexState().position().z());

    //    nBranches_->JpsiTau_tau_vy.push_back(cands[ic].cand_tau_vy);
    //    nBranches_->JpsiTau_tau_vz.push_back(cands[ic].cand_tau_vz);
    nBranches_->JpsiTau_tau_max_dr_3prong.push_back(cands[ic].cand_tau_max_dr_3prong);
    nBranches_->JpsiTau_tau_lip.push_back(cands[ic].cand_tau.lip);
    nBranches_->JpsiTau_tau_lips.push_back(cands[ic].cand_tau.lips);
    nBranches_->JpsiTau_tau_pvip.push_back(cands[ic].cand_tau.pvip);
    nBranches_->JpsiTau_tau_pvips.push_back(cands[ic].cand_tau.pvips);
    nBranches_->JpsiTau_tau_fl3d.push_back(cands[ic].cand_tau.fl3d);
    nBranches_->JpsiTau_tau_fls3d.push_back(cands[ic].cand_tau.fls3d);
    nBranches_->JpsiTau_tau_alpha.push_back(cands[ic].cand_tau.alpha);
    nBranches_->JpsiTau_tau_vprob.push_back(TMath::Prob(cands[ic].cand_tau_vertex->chiSquared(), cands[ic].cand_tau_vertex->degreesOfFreedom()) ); 
    nBranches_->JpsiTau_tau_isRight.push_back(cands[ic].cand_tau_isRight);
    nBranches_->JpsiTau_tau_isRight1.push_back(cands[ic].cand_tau_isRight1);
    nBranches_->JpsiTau_tau_isRight2.push_back(cands[ic].cand_tau_isRight2);
    nBranches_->JpsiTau_tau_isRight3.push_back(cands[ic].cand_tau_isRight3);
    nBranches_->JpsiTau_tau_matched_ppdgId.push_back(cands[ic].cand_tau_matched_ppdgId);
    nBranches_->JpsiTau_tau_matched_gentaupt.push_back(cands[ic].cand_tau_matched_gentaupt);
    nBranches_->JpsiTau_tau_sumofdnn.push_back(cands[ic].cand_tau_sumofdnn);
    nBranches_->JpsiTau_tau_sumofdnn_others.push_back(cands[ic].cand_tau_sumofdnn_others);
    nBranches_->JpsiTau_tau_pi1_dnn.push_back(cands[ic].cand_tau_pi1_dnn);
    nBranches_->JpsiTau_tau_pi2_dnn.push_back(cands[ic].cand_tau_pi2_dnn);
    nBranches_->JpsiTau_tau_pi3_dnn.push_back(cands[ic].cand_tau_pi3_dnn);

    //    nBranches_->JpsiTau_tau_pi1_doca.push_back(cands[ic].cand_ptau_pi1_doca);
    //    nBranches_->JpsiTau_tau_pi2_doca.push_back(cands[ic].cand_tau_pi2_doca);
    //    nBranches_->JpsiTau_tau_pi3_doca.push_back(cands[ic].cand_tau_pi3_doca);

    //    nBranches_->JpsiTau_tau_pi1_pv.push_back(cands[ic].cand_tau_pi1_pv);
    //    nBranches_->JpsiTau_tau_pi2_pv.push_back(cands[ic].cand_tau_pi2_pv);
    //    nBranches_->JpsiTau_tau_pi3_pv.push_back(cands[ic].cand_tau_pi3_pv);

    std::vector<Float_t> rhomass;
    pat::PackedCandidate pf1 = cands[ic].cand_pf1.pfcand;
    pat::PackedCandidate pf2 = cands[ic].cand_pf2.pfcand;
    pat::PackedCandidate pf3 = cands[ic].cand_pf3.pfcand;
      
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

    nBranches_->JpsiTau_tau_pi1_doca3d.push_back(cands[ic].cand_pf1.pfaux.doca3d);
    nBranches_->JpsiTau_tau_pi1_doca3de.push_back(cands[ic].cand_pf1.pfaux.doca3de);
    nBranches_->JpsiTau_tau_pi1_doca2d.push_back(cands[ic].cand_pf1.pfaux.doca2d);
    nBranches_->JpsiTau_tau_pi1_doca2de.push_back(cands[ic].cand_pf1.pfaux.doca2de);
    //    nBranches_->JpsiTau_tau_pi1_doca1d.push_back(cands[ic].cand_pf1.pfaux.doca1d);
    //    nBranches_->JpsiTau_tau_pi1_doca1de.push_back(cands[ic].cand_pf1.pfaux.doca1de);
    nBranches_->JpsiTau_tau_pi1_isRight.push_back(cands[ic].cand_pf1.pfaux.isRight);
    nBranches_->JpsiTau_tau_pi1_dz.push_back(cands[ic].cand_pf1.pfaux.dz);
    nBranches_->JpsiTau_tau_pi1_near_dz.push_back(cands[ic].cand_pf1.pfaux.near_dz);
    nBranches_->JpsiTau_tau_pi1_isAssociate.push_back(cands[ic].cand_pf1.pfaux.isAssociate);
    nBranches_->JpsiTau_tau_pi1_pvAssociationQuality.push_back(cands[ic].cand_pf1.pfaux.pvAssociationQuality);
    nBranches_->JpsiTau_tau_pi1_isBdecay.push_back(cands[ic].cand_pf1.pfaux.isBdecay);
    nBranches_->JpsiTau_tau_pi1_isBdecaypdg.push_back(cands[ic].cand_pf1.pfaux.isBdecaypdg);
    nBranches_->JpsiTau_tau_pi1_isBdecayppdg.push_back(cands[ic].cand_pf1.pfaux.isBdecayppdg);
    nBranches_->JpsiTau_tau_pi1_isSignal.push_back(cands[ic].cand_pf1.pfaux.isSignal);
    nBranches_->JpsiTau_tau_pi1_nprong.push_back(cands[ic].cand_pf1.pfaux.nprong);



      
    nBranches_->JpsiTau_tau_pi2_pt.push_back(tlv_pion2.Pt());
    nBranches_->JpsiTau_tau_pi2_eta.push_back(tlv_pion2.Eta());
    nBranches_->JpsiTau_tau_pi2_phi.push_back(tlv_pion2.Phi());
    nBranches_->JpsiTau_tau_pi2_mass.push_back(tlv_pion2.M());
    nBranches_->JpsiTau_tau_pi2_q.push_back(pf2.charge());

    nBranches_->JpsiTau_tau_pi2_doca3d.push_back(cands[ic].cand_pf2.pfaux.doca3d);
    nBranches_->JpsiTau_tau_pi2_doca3de.push_back(cands[ic].cand_pf2.pfaux.doca3de);
    nBranches_->JpsiTau_tau_pi2_doca2d.push_back(cands[ic].cand_pf2.pfaux.doca2d);
    nBranches_->JpsiTau_tau_pi2_doca2de.push_back(cands[ic].cand_pf2.pfaux.doca2de);
    //    nBranches_->JpsiTau_tau_pi2_doca1d.push_back(cands[ic].cand_pf2.pfaux.doca1d);
    //    nBranches_->JpsiTau_tau_pi2_doca1de.push_back(cands[ic].cand_pf2.pfaux.doca1de);
    nBranches_->JpsiTau_tau_pi2_isRight.push_back(cands[ic].cand_pf2.pfaux.isRight);
    nBranches_->JpsiTau_tau_pi2_dz.push_back(cands[ic].cand_pf2.pfaux.dz);
    nBranches_->JpsiTau_tau_pi2_near_dz.push_back(cands[ic].cand_pf2.pfaux.near_dz);
    nBranches_->JpsiTau_tau_pi2_isAssociate.push_back(cands[ic].cand_pf2.pfaux.isAssociate);
    nBranches_->JpsiTau_tau_pi2_pvAssociationQuality.push_back(cands[ic].cand_pf2.pfaux.pvAssociationQuality);
    nBranches_->JpsiTau_tau_pi2_isBdecay.push_back(cands[ic].cand_pf2.pfaux.isBdecay);
    nBranches_->JpsiTau_tau_pi2_isBdecaypdg.push_back(cands[ic].cand_pf2.pfaux.isBdecaypdg);
    nBranches_->JpsiTau_tau_pi2_isBdecayppdg.push_back(cands[ic].cand_pf2.pfaux.isBdecayppdg);
    nBranches_->JpsiTau_tau_pi2_isSignal.push_back(cands[ic].cand_pf2.pfaux.isSignal);
    nBranches_->JpsiTau_tau_pi2_nprong.push_back(cands[ic].cand_pf2.pfaux.nprong);



    nBranches_->JpsiTau_tau_pi3_pt.push_back(tlv_pion3.Pt());
    nBranches_->JpsiTau_tau_pi3_eta.push_back(tlv_pion3.Eta());
    nBranches_->JpsiTau_tau_pi3_phi.push_back(tlv_pion3.Phi());
    nBranches_->JpsiTau_tau_pi3_mass.push_back(tlv_pion3.M());
    nBranches_->JpsiTau_tau_pi3_q.push_back(pf3.charge());

    nBranches_->JpsiTau_tau_pi3_doca3d.push_back(cands[ic].cand_pf3.pfaux.doca3d);
    nBranches_->JpsiTau_tau_pi3_doca3de.push_back(cands[ic].cand_pf3.pfaux.doca3de);
    nBranches_->JpsiTau_tau_pi3_doca2d.push_back(cands[ic].cand_pf3.pfaux.doca2d);
    nBranches_->JpsiTau_tau_pi3_doca2de.push_back(cands[ic].cand_pf3.pfaux.doca2de);
    //    nBranches_->JpsiTau_tau_pi3_doca1d.push_back(cands[ic].cand_pf3.pfaux.doca1d);
    //    nBranches_->JpsiTau_tau_pi3_doca1de.push_back(cands[ic].cand_pf3.pfaux.doca1de);
    nBranches_->JpsiTau_tau_pi3_isRight.push_back(cands[ic].cand_pf3.pfaux.isRight);
    nBranches_->JpsiTau_tau_pi3_dz.push_back(cands[ic].cand_pf3.pfaux.dz);
    nBranches_->JpsiTau_tau_pi3_near_dz.push_back(cands[ic].cand_pf3.pfaux.near_dz);
    nBranches_->JpsiTau_tau_pi3_isAssociate.push_back(cands[ic].cand_pf3.pfaux.isAssociate);
    nBranches_->JpsiTau_tau_pi3_pvAssociationQuality.push_back(cands[ic].cand_pf3.pfaux.pvAssociationQuality);
    nBranches_->JpsiTau_tau_pi3_isBdecay.push_back(cands[ic].cand_pf3.pfaux.isBdecay);
    nBranches_->JpsiTau_tau_pi3_isBdecaypdg.push_back(cands[ic].cand_pf3.pfaux.isBdecaypdg);
    nBranches_->JpsiTau_tau_pi3_isBdecayppdg.push_back(cands[ic].cand_pf3.pfaux.isBdecayppdg);
    nBranches_->JpsiTau_tau_pi3_isSignal.push_back(cands[ic].cand_pf3.pfaux.isSignal);
    nBranches_->JpsiTau_tau_pi3_nprong.push_back(cands[ic].cand_pf3.pfaux.nprong);



    nBranches_->JpsiTau_B_pt.push_back(cands[ic].cand_b_part->currentState().globalMomentum().perp());
    nBranches_->JpsiTau_B_eta.push_back(cands[ic].cand_b_part->currentState().globalMomentum().eta());
    nBranches_->JpsiTau_B_phi.push_back(cands[ic].cand_b_part->currentState().globalMomentum().phi());
    nBranches_->JpsiTau_B_mass.push_back(cands[ic].cand_b_part->currentState().mass());
    nBranches_->JpsiTau_B_vprob.push_back( TMath::Prob(cands[ic].cand_b_vertex->chiSquared(), cands[ic].cand_b_vertex->degreesOfFreedom()) ); //cands[ic].cand_b_part->currentState().mass());
    nBranches_->JpsiTau_B_lip.push_back(cands[ic].cand_b.lip);
    nBranches_->JpsiTau_B_lips.push_back(cands[ic].cand_b.lips);
    nBranches_->JpsiTau_B_pvip.push_back(cands[ic].cand_b.pvip);
    nBranches_->JpsiTau_B_pvips.push_back(cands[ic].cand_b.pvips);
    nBranches_->JpsiTau_B_fls3d.push_back(cands[ic].cand_b.fls3d);
    nBranches_->JpsiTau_B_fl3d.push_back(cands[ic].cand_b.fl3d);
    nBranches_->JpsiTau_B_alpha.push_back(cands[ic].cand_b.alpha);

//    std::vector<RefCountedKinematicParticle> allParticles4doc;
//      
//    allParticles4doc.push_back(pFactory.particle(tt1_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));
//    allParticles4doc.push_back(pFactory.particle(tt2_muon, aux.muon_mass, chi, ndf, aux.muon_sigma));
//    allParticles4doc.push_back(pFactory.particle(mytracks[cands[ic].cand_tau_id1], aux.pion_mass, chi, ndf, aux.pion_sigma));
//    allParticles4doc.push_back(pFactory.particle(mytracks[cands[ic].cand_tau_id2], aux.pion_mass, chi, ndf, aux.pion_sigma));
//    allParticles4doc.push_back(pFactory.particle(mytracks[cands[ic].cand_tau_id3], aux.pion_mass, chi, ndf, aux.pion_sigma));
//
//
//    nBranches_->JpsiTau_B_maxdoca.push_back(aux.getMaxDoca(allParticles4doc));
//    nBranches_->JpsiTau_B_mindoca.push_back(aux.getMinDoca(allParticles4doc));
    nBranches_->JpsiTau_B_vx.push_back(cands[ic].cand_b_vertex->vertexState().position().x());
    nBranches_->JpsiTau_B_vy.push_back(cands[ic].cand_b_vertex->vertexState().position().y());
    nBranches_->JpsiTau_B_vz.push_back(cands[ic].cand_b_vertex->vertexState().position().z());

    
    for(unsigned int iiso=0; iiso<cands[ic].cand_b_iso.size(); iiso++){
      nBranches_->JpsiTau_B_iso.push_back(cands[ic].cand_b_iso[iiso]);
      nBranches_->JpsiTau_B_iso_ntracks.push_back(cands[ic].cand_b_iso_ntracks[iiso]);
      nBranches_->JpsiTau_B_iso_mindoca.push_back(cands[ic].cand_b_iso_mindoca[iiso]);
    }

    //    nBranches_->JpsiTau_B_iso_nocut.push_back(cands[ic].cand_b_iso_nocut);
    //    nBranches_->JpsiTau_B_iso_ntracks_nocut.push_back(cands[ic].cand_b_iso_ntracks_nocut);
    //    nBranches_->JpsiTau_B_iso_mindoca_nocut.push_back(cands[ic].cand_b_iso_mindoca_nocut);


      
    TLorentzVector Tlv_B;
    TLorentzVector Tlv_Jpsi;
    TLorentzVector Tlv_tau;
      
    Tlv_B.SetPtEtaPhiM(cands[ic].cand_b_part->currentState().globalMomentum().perp(),
		       cands[ic].cand_b_part->currentState().globalMomentum().eta(),
		       cands[ic].cand_b_part->currentState().globalMomentum().phi(),
		       cands[ic].cand_b_part->currentState().mass());

    // calculation of the corrected mass
    TVector3 *bvector = new TVector3(Tlv_B.Px(), Tlv_B.Py(), Tlv_B.Pz());
      
    Float_t pperp = bvector->Mag()*TMath::Sin(TMath::ACos(cands[ic].cand_b.alpha));
    
    Float_t mcorr = pperp + TMath::Sqrt(pperp*pperp + Tlv_B.M()*Tlv_B.M());
    
    nBranches_->JpsiTau_B_mcorr.push_back(mcorr);


    Tlv_B *= aux.mass_Bc/(cands[ic].cand_b_part->currentState().mass());
      
    Tlv_Jpsi.SetPtEtaPhiM(jpsi_part->currentState().globalMomentum().perp(),
			  jpsi_part->currentState().globalMomentum().eta(),
			  jpsi_part->currentState().globalMomentum().phi(),
			  jpsi_part->currentState().mass());
      
    Tlv_tau.SetPtEtaPhiM(cands[ic].cand_tlv_tau_fit.Pt(), 
			 cands[ic].cand_tlv_tau_fit.Eta(),
			 cands[ic].cand_tlv_tau_fit.Phi(),
			 cands[ic].cand_tlv_tau_fit.M() );


    
    // calculate cross products 

    TVector3 pvsv = TVector3( jpsi_vertex->vertexState().position().x() - closestVertex.position().x(), 
			      jpsi_vertex->vertexState().position().y() - closestVertex.position().y(), 
			      jpsi_vertex->vertexState().position().z() - closestVertex.position().z());


    TVector3 jpsi_vec = TVector3(Tlv_Jpsi.Px(), Tlv_Jpsi.Py(), Tlv_Jpsi.Pz());
    TVector3 tau_vec = TVector3(Tlv_tau.Px(), Tlv_tau.Py(), Tlv_tau.Pz());

    Float_t ptbal_jpsi = jpsi_vec.Cross(pvsv).Mag();
    Float_t ptbal_tau = tau_vec.Cross(pvsv).Mag();
    
    Float_t ptbal = -1;

    if(ptbal_jpsi!=0){
      ptbal = ptbal_tau/ptbal_jpsi;
    }


    nBranches_->JpsiTau_ptbal.push_back(ptbal);

    Float_t jpsi_tau_alpha = -1;

    if(jpsi_vec.Mag()!=0 && tau_vec.Mag()!=0){
      jpsi_tau_alpha = jpsi_vec.Dot(tau_vec)/(jpsi_vec.Mag()*tau_vec.Mag());
    }

    nBranches_->JpsiTau_jpsi_tau_alpha.push_back(jpsi_tau_alpha);


    //////////////////////////////////////////      

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
  //  nBranches_->JpsiTau_nch_pvipcut.push_back(npf_pvip);
  //  nBranches_->JpsiTau_nch_pvipcut_dz0.push_back(npf_pvip_dz0);
  //  nBranches_->JpsiTau_nch_pvipcut_dz2p5.push_back(npf_pvip_dz2p5);
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
    for(auto n: parName) settings["delta_" + n] = 0;
    
    //    for(auto pars: _FFErrNames) {
    //      settings[pars] = 0;
    //    }

    hammer.setFFEigenvectors("BctoJpsi", "BGLVar", settings);
	
    auto weights = hammer.getWeights("Scheme1");
    //    auto rate_orig = hammer.getRate("BcJpsiTau+Nu", "Scheme1");

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

    //    std::cout << "-----------------------" << std::endl;
    //    std::cout << "base weight = " << weight << std::endl;
    

    ///////////////////////////////////////////////////////////////////////
    // MC 
    ///////////////////////////////////////////////////////////////////////

    //    Float_t hammer_up = 0;
    //    Float_t hammer_down = 999;
    //    int ntrial = 2000;
    //    int ntrial = 500;
    //    int npass = 0;

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
	
	nBranches_->JpsiTau_hammer_ff.push_back(settings_ff);      
      }

      hweights.push_back(weight_sys);

    }

    flag_fill = true;

    nBranches_->JpsiTau_hammer_ebe_toy.push_back(hweights);


    //////////////////////////
    // ordinary method 
    //////////////////////////

//    for(int i=0; i<11; i++) { //Loop over eigenVar
//      for (int j=0; j<2; j++) { //Loop over pos/neg direction
//
//        map<string, double> settings;
//
//        for (int k=0; k<11; k++) { //Loop over parameters
//          settings["delta_" + parName[k]] = eigVar[i][k][j];
//        }
//
//        hammer.setFFEigenvectors("BctoJpsi", "BGLVar", settings);
//        auto weights = hammer.getWeights("Scheme1");
//        string var_name = "eig";
//	var_name += std::to_string(i);
//        var_name += j==0? "_Up" : "_Down";
//
//
//
//	Float_t weight_sys = -1;
//
//	if(!weights.empty()){
//	  for(auto elem: weights) {
//	    if(isnan(elem.second)) {
//	      std::cout << "[ERROR]: BGL nan weight: " << elem.second << std::endl;
//	    }else{
//	      weight_sys = elem.second;
//	    }
//	  }
//	}
//
//	if(var_name==std::string("eig0_Up")) nBranches_->JpsiTau_hammer_ebe_a0_up.push_back(weight_sys);
//	else if(var_name==std::string("eig0_Down")) nBranches_->JpsiTau_hammer_ebe_a0_down.push_back(weight_sys);
//	else if(var_name==std::string("eig1_Up")) nBranches_->JpsiTau_hammer_ebe_a1_up.push_back(weight_sys);
//	else if(var_name==std::string("eig1_Down")) nBranches_->JpsiTau_hammer_ebe_a1_down.push_back(weight_sys);
//	else if(var_name==std::string("eig2_Up")) nBranches_->JpsiTau_hammer_ebe_a2_up.push_back(weight_sys);
//	else if(var_name==std::string("eig2_Down")) nBranches_->JpsiTau_hammer_ebe_a2_down.push_back(weight_sys);
//
//	else if(var_name==std::string("eig3_Up")) nBranches_->JpsiTau_hammer_ebe_b0_up.push_back(weight_sys);
//	else if(var_name==std::string("eig3_Down")) nBranches_->JpsiTau_hammer_ebe_b0_down.push_back(weight_sys);
//	else if(var_name==std::string("eig4_Up")) nBranches_->JpsiTau_hammer_ebe_b1_up.push_back(weight_sys);
//	else if(var_name==std::string("eig4_Down")) nBranches_->JpsiTau_hammer_ebe_b1_down.push_back(weight_sys);
//	else if(var_name==std::string("eig5_Up")) nBranches_->JpsiTau_hammer_ebe_b2_up.push_back(weight_sys);
//	else if(var_name==std::string("eig5_Down")) nBranches_->JpsiTau_hammer_ebe_b2_down.push_back(weight_sys);
//
//	else if(var_name==std::string("eig6_Up")) nBranches_->JpsiTau_hammer_ebe_c1_up.push_back(weight_sys);
//	else if(var_name==std::string("eig6_Down")) nBranches_->JpsiTau_hammer_ebe_c1_down.push_back(weight_sys);
//	else if(var_name==std::string("eig7_Up")) nBranches_->JpsiTau_hammer_ebe_c2_up.push_back(weight_sys);
//	else if(var_name==std::string("eig7_Down")) nBranches_->JpsiTau_hammer_ebe_c2_down.push_back(weight_sys);
//
//	else if(var_name==std::string("eig8_Up")) nBranches_->JpsiTau_hammer_ebe_d0_up.push_back(weight_sys);
//	else if(var_name==std::string("eig8_Down")) nBranches_->JpsiTau_hammer_ebe_d0_down.push_back(weight_sys);
//	else if(var_name==std::string("eig9_Up")) nBranches_->JpsiTau_hammer_ebe_d1_up.push_back(weight_sys);
//	else if(var_name==std::string("eig9_Down")) nBranches_->JpsiTau_hammer_ebe_d1_down.push_back(weight_sys);
//	else if(var_name==std::string("eig10_Up")) nBranches_->JpsiTau_hammer_ebe_d2_up.push_back(weight_sys);
//	else if(var_name==std::string("eig10_Down")) nBranches_->JpsiTau_hammer_ebe_d2_down.push_back(weight_sys);
//	
//	//	std::cout << "test1 \t " << var_name << " " << weight_sys << std::endl;
//      }
//    }



//    for(auto pars1: _FFErrNames) {
//	  
//      for(int isUp=0; isUp < 2; isUp++) { // up, down variation    
//	    
//	std::map<std::string, double> settings;
//	    
//	int idx2 = 0;
//	    
//	for(auto pars2: _FFErrNames) {
//	      
//	  if(idx1 == idx2){
//	    if(isUp==1) settings[pars2] = _FFErr[idx2];
//	    else settings[pars2] = -_FFErr[idx2];
//	  }else{
//	    settings[pars2] = 0;
//	  }
//	  idx2 += 1;
//	}
//
//	hammer.setFFEigenvectors("BctoJpsi", "BGLVar", settings);
//	auto weights = hammer.getWeights("Scheme1");
//	std::string var_name = pars1;
//	var_name += isUp==1? "_Up" : "_Down";
//
//	Float_t weight_sys = -1;
//
//	if(!weights.empty()){
//	  for(auto elem: weights) {
//	    if(isnan(elem.second)) {
//	      std::cout << "[ERROR]: BGL nan weight: " << elem.second << std::endl;
//	    }else{
//	      weight_sys = elem.second;
//	    }
//	  }
//	}
//
//	if(isUp==0) hammer_up += TMath::Power(weight_sys, 2);
//	else if(isUp==1) hammer_down += TMath::Power(weight_sys, 2);
//
//	std::cout << "\t " << var_name << " " << weight_sys << std::endl;
//	   
//      }
//      
//      idx1 += 1;
//      
//    }






//    std::cout << "Up, down variation = " << TMath::Sqrt(hammer_up) << " " << TMath::Sqrt(hammer_down) << std::endl;
//    nBranches_->JpsiTau_hammer_ebe_up.push_back(hammer_up);
//    nBranches_->JpsiTau_hammer_ebe_down.push_back(hammer_down);

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


