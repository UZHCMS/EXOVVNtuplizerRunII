#cmsRun  config_generic_opt_skimmed.py  RunPeriod="Fall17" # for MC
#cmsRun  config_generic_opt_skimmed.py  RunPeriod="Run2017B" # for Data from 2017
#cmsRun  config_generic_opt_skimmed.py  RunPeriod="Run2018B" # for Data from 2018
#cmsRun  config_generic_opt_skimmed.py  RunPeriod="Autumn18" # for MC from 2018
#cmsRun  config_generic_opt_skimmed.py  RunPeriod="Summer16" # for MC from 2016


###### Process initialization ##########

import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuple")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')

process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string('flatTuple.root')
                                   )

#from EXOVVNtuplizerRunII.Ntuplizer.ntuplizerOptions_data_cfi import config
from EXOVVNtuplizerRunII.Ntuplizer.ntuplizerOptions_generic_cfi import config

# change from its original value
#config["DZCUT"] = 0.25
#config["FSIGCUT"] = 3
#config["VPROBCUT"] = 0.1
#config["DNNCUT"] = 0.2

				   
####### Config parser ##########

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('analysis')

options.register( 'RunPeriod',
                  '',
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "RunNumber (Default Run2017B)")

options.register( 'runUpToEarlyF',
                  'false',
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list                                                                                                                                 
                  VarParsing.VarParsing.varType.bool,          # string, int, or float                                                                                                                        
                  "false")# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAboutPythonConfigFile

####



options.maxEvents = 200
#options.maxEvents = -1

#data file

#options.inputFiles = '/store/mc/RunIIAutumn18MiniAOD/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v4/230000/322CC3C2-921E-7448-902A-8FCB1F0A2F72.root'
#options.inputFiles = '/store/user/manzoni/BcToJpsiX_TuneCP5_13TeV-pythia8/RunIISummer19UL18_MINIAODSIM_v1/201112_133556/0000/RJpsi-BcToXToJpsiMuMuSelected-RunIISummer19UL18MiniAOD_7.root'

#options.inputFiles = '/store/mc/RunIIAutumn18MiniAOD/BcToJPsiMuNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v3/100000/78A9DF86-EAC5-D242-8E7B-171AFD012CC7.root'

options.inputFiles = '/store/mc/RunIISummer19UL18MiniAOD/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1_ext1-v2/100000/02F13381-1D94-CC43-948A-2EFFB8572949.root'

#options.inputFiles = '/store/mc/RunIIAutumn18MiniAOD/OniaAndX_ToMuMu_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/00000/01325465-A815-E24E-ABB3-DAB8D4880BDE.root'

#options.inputFiles = '/store/data/Run2018D/Charmonium/MINIAOD/12Nov2019_UL2018-v1/280000/D7FD376D-30CD-AA48-8D03-E0220043BBDE.root'
#options.inputFiles = '/store/data/Run2017F/Charmonium/MINIAOD/09Aug2019_UL2017-v1/20000/00BACB48-9B0F-8F48-A68B-2F08A3E9E681.root'
#options.inputFiles = '/store/data/Run2016B/Charmonium/MINIAOD/21Feb2020_ver2_UL2016_HIPM-v1/240000/0333D5C7-28C0-7641-994A-ADE29A1EBAAD.root'
#options.inputFiles = '/store/data/Run2016B/Charmonium/MINIAOD/21Feb2020_ver2_UL2016_HIPM-v1/240000/00251310-7FD5-BE47-B127-8CFC5B8DFE6E.root'

options.parseArguments()

process.options  = cms.untracked.PSet( 
                     wantSummary = cms.untracked.bool(True),
                     SkipEvent = cms.untracked.vstring('ProductNotFound'),
                     allowUnscheduled = cms.untracked.bool(True),
                     )

#process.options.numberOfThreads=cms.untracked.uint32(2)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

# run, lumi, event ID
#                              eventsToProcess = cms.untracked.VEventRange('1:94:182460'),

if config["RUNONMC"]:
  process.source = cms.Source("PoolSource",
                              fileNames = cms.untracked.vstring(options.inputFiles),
                              # Run, lumi, Event
#                              eventsToProcess = cms.untracked.VEventRange('1:35137:14345'),
                              duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),                              
                              ) 
else:                    
  process.source = cms.Source("PoolSource",
                              fileNames = cms.untracked.vstring(options.inputFiles),
#                              skipEvents=cms.untracked.uint32(23000)
                              ) 

print " process source filenames %s" %(process.source)

#print " process source filenames ", process.source.fileNames
######## Sequence settings ##########

hltFiltersProcessName = 'RECO'
if config["RUNONMC"] or config["JSONFILE"].find('reMiniAOD') != -1:
  hltFiltersProcessName = 'PAT'

# ####### Logger ##########
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Ntuple')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(1)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

####### Define conditions ##########
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag import GlobalTag

GT = ''


jetcorr_levels=[]
jetcorr_levels_groomed=[]
if config["RUNONMC"]:
  jetcorr_levels = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])
  jetcorr_levels_groomed = cms.vstring(['L2Relative', 'L3Absolute']) # NO L1 corretion for groomed jets
else:
  jetcorr_levels = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
  jetcorr_levels_groomed = cms.vstring(['L2Relative', 'L3Absolute', 'L2L3Residual'])

   
######### read JSON file for data ##########					                                                             
if not(config["RUNONMC"]) and config["USEJSON"]:

  import FWCore.PythonUtilities.LumiList as LumiList
  import FWCore.ParameterSet.Types as CfgTypes
  process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
  myLumis = LumiList.LumiList(filename = config["JSONFILE"]).getCMSSWString().split(',')
  process.source.lumisToProcess.extend(myLumis) 

  

####### Redo Jet clustering sequence ##########
betapar = cms.double(0.0)
fatjet_ptmin = 100.0

from RecoJets.Configuration.RecoPFJets_cff import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoJets.JetProducers.PFJetParameters_cfi import *

from PhysicsTools.PatAlgos.tools.helpers import *
pattask = getPatAlgosToolsTask(process)
                                                                                                          
process.chs = cms.EDFilter("CandPtrSelector",
  src = cms.InputTag('packedPFCandidates'),
  cut = cms.string('fromPV')
)

process.ak4PFJetsCHS = ak4PFJetsCHS.clone( src = 'chs' )
process.ak4PFJetsCHS.doAreaFastjet = True


from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

#dataFormat=DataFormat.MiniAOD
#switchOnVIDElectronIdProducer(process,dataFormat,task=pattask)

#process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)
     
#my_id_modules = [
#                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff',
#                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff',
#                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
#                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',
#                 ]
           
#add them to the VID producer
#for idmod in my_id_modules:
#    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection,task=pattask)

####### Ntuplizer initialization ##########
jetsAK4 = "slimmedJets"


METS = "slimmedMETs"
METS_EGclean = "slimmedMETsEGClean"
METS_MEGclean = "slimmedMETsMuEGClean"
METS_uncorr = "slimmedMETsUncorrected"

if config["USENOHF"]: METS = "slimmedMETsNoHF"  

##___________________ MET significance and covariance matrix ______________________##


##___________________ Jets ______________________##

  

######## JEC ########
jecLevelsAK8chs = []
jecLevelsAK8Groomedchs = []
jecLevelsAK4chs = []
jecLevelsAK4 = []
jecLevelsAK8Puppi = []
jecLevelsForMET = []

print "1. options.RunPeriod ", options.RunPeriod
if options.RunPeriod=="" : options.RunPeriod=options.inputFiles[0]

if  config["RUNONMC"] :
  JECprefix = ""
  if ("Fall17" in options.RunPeriod):
    JECprefix = "Fall17_17Nov2017_V32"
    GT ='94X_mcRun2_asymptotic_v3'
  elif ("Summer16" in options.RunPeriod):
    JECprefix = "Summer16_07Aug2017_V11"
    GT = '94X_mc2017_realistic_v17'
  elif (("Autumn18" in options.RunPeriod) or ("Fall18" in options.RunPeriod)):
    JECprefix = "Autumn18_V8"
    GT = '102X_upgrade2018_realistic_v18'
  else:
    JECprefix = "Autumn18_V8"
    GT = '102X_upgrade2018_realistic_v18'

  #jecAK8chsUncFile = "JEC/%s_MC_Uncertainty_AK8PFchs.txt"%(JECprefix)
  jecAK4chsUncFile = "JEC/%s_MC_Uncertainty_AK4PFchs.txt"%(JECprefix)
 



else : #Data
   JECprefix = ""
   JEC_runDependent_suffix= ""

   if ("2017" in options.RunPeriod):
     if ("Run2017B" in  options.RunPeriod): JEC_runDependent_suffix= "B"
     elif ("Run2017C" in  options.RunPeriod): JEC_runDependent_suffix= "C"
     elif ("Run2017D" in  options.RunPeriod): JEC_runDependent_suffix= "D"
     elif ("Run2017E" in  options.RunPeriod): JEC_runDependent_suffix= "E"
     elif ("Run2017F" in  options.RunPeriod): JEC_runDependent_suffix= "F"
     
     JECprefix = "Fall17_17Nov2017"+JEC_runDependent_suffix+"_V32"
     GT = '94X_dataRun2_v11'
     
     
   elif ("2016" in options.RunPeriod):
     if ("Run2016D" in  options.RunPeriod or "Run2016B" in  options.RunPeriod  or "Run2016C" in  options.RunPeriod  ): JEC_runDependent_suffix= "ABC"
     elif ("Run2016E" in  options.RunPeriod): JEC_runDependent_suffix= "EF"
     elif ("Run2016G" in  options.RunPeriod): JEC_runDependent_suffix= "GH"
     elif ("Run2016F" in  options.RunPeriod and  not options.runUpToEarlyF): JEC_runDependent_suffix= "GH"
     elif ("Run2016F" in  options.RunPeriod and   options.runUpToEarlyF): JEC_runDependent_suffix= "EF"


     JECprefix = "Summer16_07Aug2017"+JEC_runDependent_suffix+"_V11"
     GT ='94X_dataRun2_v10'

   elif ("2018" in options.RunPeriod):
     if ("Run2018A" in  options.RunPeriod ): 
       JEC_runDependent_suffix= "A"
       GT="102X_dataRun2_Sep2018ABC_v2" 
     elif ("Run2018B" in  options.RunPeriod): 
       JEC_runDependent_suffix= "B"
       GT="102X_dataRun2_Sep2018ABC_v2"
     elif ("Run2018C" in  options.RunPeriod): 
       JEC_runDependent_suffix= "C"
       GT="102X_dataRun2_Sep2018ABC_v2"
     elif ("Run2018D" in  options.RunPeriod): 
       JEC_runDependent_suffix= "D"
       GT = '102X_dataRun2_Prompt_v16' 

     JECprefix = "Autumn18_Run"+JEC_runDependent_suffix+"_V8"
    
   #jecAK8chsUncFile = "JEC/%s_DATA_Uncertainty_AK8PFchs.txt"%(JECprefix)
   jecAK4chsUncFile = "JEC/%s_DATA_Uncertainty_AK4PFchs.txt"%(JECprefix)
 
#   GT = '106X_dataRun2_v27' 
   print "jec JEC_runDependent_suffix %s ,  prefix %s " %(JEC_runDependent_suffix,JECprefix)

print "jec prefix ", JECprefix

print "doing corrections  to met on the fly %s" ,config["CORRMETONTHEFLY"]

print "*************************************** GLOBAL TAG *************************************************" 
print GT
print "****************************************************************************************************" 
process.GlobalTag = GlobalTag(process.GlobalTag, GT)




if config["CORRMETONTHEFLY"]:  
   if config["RUNONMC"]:
     jecLevelsForMET = [				       
     	 'JEC/%s_MC_L1FastJet_AK4PFchs.txt'%(JECprefix),
     	 'JEC/%s_MC_L2Relative_AK4PFchs.txt'%(JECprefix),
     	 'JEC/%s_MC_L3Absolute_AK4PFchs.txt'%(JECprefix)
       ]
   else:       					       
     jecLevelsForMET = [
     	 'JEC/%s_DATA_L1FastJet_AK4PFchs.txt'%(JECprefix),
     	 'JEC/%s_DATA_L2Relative_AK4PFchs.txt'%(JECprefix),
     	 'JEC/%s_DATA_L3Absolute_AK4PFchs.txt'%(JECprefix),
         'JEC/%s_DATA_L2L3Residual_AK4PFchs.txt'%(JECprefix)
       ]	
      	

################# Weights for B generic background #########

IsBkgBSample = False 
fileName=str(process.source.fileNames) 
print(fileName)
if "HbToJPsiMuMu" in fileName:
   IsBkgBSample = True		    
print "IsBkgBSample" ,IsBkgBSample

################## Ntuplizer ###################
process.ntuplizer = cms.EDAnalyzer("Ntuplizer",
    runOnMC	      = cms.bool(config["RUNONMC"]),
    useDNN	      = cms.bool(config["USEDNN"]),
    useHammer	      = cms.bool(config["USEHAMMER"]),
    doGenParticles    = cms.bool(config["DOGENPARTICLES"]),
    doGenEvent	      = cms.bool(config["DOGENEVENT"]),
    doPileUp	      = cms.bool(config["DOPILEUP"]),
    doJpsiMu	      = cms.bool(config["DOJPSIMU"]),
    doJpsiTau	      = cms.bool(config["DOJPSITAU"]),
    doBsTauTau	      = cms.bool(config["DOBSTAUTAU"]),
    doBsTauTauFH      = cms.bool(config["DOBSTAUTAUFH"]),
    doBsTauTauFH_mr   = cms.bool(config["DOBSTAUTAUFH_mr"]),
    doBsDstarTauNu    = cms.bool(config["DOBSDSTARTAUNU"]),
    doVertices	      = cms.bool(config["DOVERTICES"]),
    doMissingEt       = cms.bool(config["DOMISSINGET"]),
    doGenHist         = cms.bool(config["DOGENHIST"]),
    verbose           = cms.bool(config["VERBOSE"]),
    dzcut             = cms.double(config['DZCUT']),
    fsigcut           = cms.double(config['FSIGCUT']),
    vprobcut          = cms.double(config['VPROBCUT']),
    dnncut            = cms.double(config['DNNCUT']),
    tau_charge        = cms.uint32(config['TAU_CHARGE']),
    dnnfile_old       = cms.string(config['DNNFILE_OLD']),                        
    dnnfile_perPF     = cms.string(config['DNNFILE_PERPF']),                        
    dnnfile_perEVT    = cms.string(config['DNNFILE_PEREVT']),                        
    dnnfile_perEVT_v2 = cms.string(config['DNNFILE_PEREVT_V2']),
    isBkgBSample = cms.bool(IsBkgBSample),                        
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    taus = cms.InputTag("slimmedTaus"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    ebRecHits = cms.InputTag("reducedEgamma","reducedEBRecHits"),

#    eleHEEPId51Map = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV51"),
#    eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
#    eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
#    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
#    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
#    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),

#    eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-veto"),
#    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-loose"),
#    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-medium"),
#    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-tight"),

#    eleHLTIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronHLTPreselection-Summer16-V1"), 
#    eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70"),
                                   
#    eleMVAMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp90"),
#    eleMVATightIdMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp80"),
#    mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
#    mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),
    dupCluster          = cms.InputTag("particleFlowEGammaGSFixed:dupECALClusters"),
    hitsNotReplaced     = cms.InputTag("ecalMultiAndGSGlobalRecHitEB:hitsNotReplaced"),
    mets = cms.InputTag(METS),
    mets_EGclean = cms.InputTag(METS_EGclean),
    mets_MEGclean = cms.InputTag(METS_MEGclean),
    mets_uncorr = cms.InputTag(METS_uncorr),
    mets_puppi = cms.InputTag("slimmedMETsPuppi"),
    mets_mva = cms.InputTag("MVAMET","MVAMET"),
    corrMetPx = cms.string("+0.1166 + 0.0200*Nvtx"),
    corrMetPy = cms.string("+0.2764 - 0.1280*Nvtx"),
    jecAK4forMetCorr = cms.vstring( jecLevelsForMET ),
    jetsForMetCorr = cms.InputTag(jetsAK4),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    genparticles = cms.InputTag("prunedGenParticles"),
    packedgenparticles = cms.InputTag("packedGenParticles"),
#    gentaus = cms.InputTag("tauGenJets"),
    PUInfo = cms.InputTag("slimmedAddPileupInfo"),
    genEventInfo = cms.InputTag("generator"),
    externallheProducer = cms.InputTag("externalLHEProducer"),
    HLT = cms.InputTag("TriggerResults","","HLT"),
    triggerobjects = cms.InputTag("slimmedPatTrigger"),
    triggerprescales = cms.InputTag("patTrigger"),
    noiseFilter = cms.InputTag('TriggerResults','', hltFiltersProcessName),
    jecpath = cms.string(''),
   
    
    ## Noise Filters ###################################
    # defined here: https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py
    noiseFilterSelection_HBHENoiseFilter = cms.string('Flag_HBHENoiseFilter'),   # both data and MC for 2018,
    noiseFilterSelection_HBHENoiseFilterLoose = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResultRun2Loose"),
    noiseFilterSelection_HBHENoiseFilterTight = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResultRun2Tight"),
    noiseFilterSelection_HBHENoiseIsoFilter = cms.InputTag("HBHENoiseFilterResultProducer", "HBHEIsoNoiseFilterResult"),    # both data and MC for 2018,  
    noiseFilterSelection_ecalBadCalibReducedMINIAODFilter = cms.InputTag("ecalBadCalibReducedMINIAODFilter"),  # both data and MC for 2018,
    noiseFilterSelection_CSCTightHaloFilter = cms.string('Flag_CSCTightHaloFilter'),
    noiseFilterSelection_CSCTightHalo2015Filter = cms.string('Flag_CSCTightHalo2015Filter'),
    noiseFilterSelection_hcalLaserEventFilter = cms.string('Flag_hcalLaserEventFilter'),
    noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter = cms.string('Flag_EcalDeadCellTriggerPrimitiveFilter'),  # both data and MC for 2018,
    noiseFilterSelection_goodVertices = cms.string('Flag_goodVertices'),  # both data and MC for 2018,
    noiseFilterSelection_trackingFailureFilter = cms.string('Flag_trackingFailureFilter'),
    noiseFilterSelection_eeBadScFilter = cms.string('Flag_eeBadScFilter'),
    noiseFilterSelection_ecalLaserCorrFilter = cms.string('Flag_ecalLaserCorrFilter'),
    noiseFilterSelection_trkPOGFilters = cms.string('Flag_trkPOGFilters'),
    
    #New for ICHEP 2016
    noiseFilterSelection_CSCTightHaloTrkMuUnvetoFilter = cms.string('Flag_CSCTightHaloTrkMuUnvetoFilter'),
    noiseFilterSelection_globalTightHalo2016Filter = cms.string('Flag_globalTightHalo2016Filter'),
    noiseFilterSelection_globalSuperTightHalo2016Filter = cms.string('Flag_globalSuperTightHalo2016Filter'), # both data and MC for 2018,  
    noiseFilterSelection_HcalStripHaloFilter = cms.string('Flag_HcalStripHaloFilter'),
    noiseFilterSelection_chargedHadronTrackResolutionFilter = cms.string('Flag_chargedHadronTrackResolutionFilter'),
    noiseFilterSelection_muonBadTrackFilter = cms.string('Flag_muonBadTrackFilter'),
    
    #New for Moriond
    noiseFilterSelection_badMuonsFilter = cms.string('Flag_BadPFMuonFilter'),    #('Flag_badMuons'),  # both data and MC for 2018, 
    noiseFilterSelection_duplicateMuonsFilter = cms.string('Flag_duplicateMuons'),
    noiseFilterSelection_nobadMuonsFilter = cms.string('Flag_nobadMuons'),

    # and the sub-filters
    noiseFilterSelection_trkPOG_manystripclus53X = cms.string('Flag_trkPOG_manystripclus53X'),
    noiseFilterSelection_trkPOG_toomanystripclus53X = cms.string('Flag_trkPOG_toomanystripclus53X'),
    noiseFilterSelection_trkPOG_logErrorTooManyClusters = cms.string('Flag_trkPOG_logErrorTooManyClusters'),
    # summary
    noiseFilterSelection_metFilters = cms.string('Flag_METFilters'),

    packedpfcandidates = cms.InputTag('packedPFCandidates'),
    SecondaryVertices = cms.InputTag('slimmedSecondaryVertices'),
#    losttrack = cms.InputTag('lostTracks')
)

process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

baddetEcallist = cms.vuint32(
    [872439604,872422825,872420274,872423218,
     872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,
     872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,
     872437052,872420649,872422436,872421950,
     872437185,872422564,872421566,872421695,
     872421955,872421567,872437184,872421951,
     872421694,872437056,872437057,872437313])


process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
    "EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallist, 
    taggingMode = cms.bool(True),
    debug = cms.bool(True)#False
    )




process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder')
)


####### Final path ##########
process.p = cms.Path()

process.p += process.ecalBadCalibReducedMINIAODFilter

#if config["RUNONMC"]:
#  process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
#  process.tauGenJets.GenParticles = cms.InputTag('prunedGenParticles')
#  process.p += process.tauGenJets


process.p += process.ntuplizer
process.p.associate(pattask)

print pattask

#  LocalWords:  tauIdMVAIsoDBoldDMwLT
