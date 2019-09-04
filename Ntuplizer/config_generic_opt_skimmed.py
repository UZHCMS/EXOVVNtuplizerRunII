#cmsRun  config_generic_opt_skimmed.py  RunPeriod="Fall17" # for MC
#cmsRun  config_generic_opt_skimmed.py  RunPeriod="Run2017B" # for Data


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

				   
####### Config parser ##########

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('analysis')

options.register( 'RunPeriod',
                  '',
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.string,          # string, int, or float
                  "RunNumber (Default Run2017B)")

options.register( 'runUpToEarlyF',
                  'True',
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list                                                                                                                                 
                  VarParsing.VarParsing.varType.int,          # string, int, or float                                                                                                                        
                  "1")# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAboutPythonConfigFile

####



options.maxEvents = 100

#data file
     
options.inputFiles ='file:/work/pmatorra/JpsiAnalysis/2018/CMSSW_10_2_10/src/EXOVVNtuplizerRunII/Ntuplizer/miniAOD_99.root'
#options.inputFiles ='/store/data/Run2018B/Charmonium/MINIAOD/17Sep2018-v1/10000/02CFE87F-7C17-1340-8300-FDA86C16D58C.root'

options.parseArguments()

process.options  = cms.untracked.PSet( 
                     wantSummary = cms.untracked.bool(True),
                     SkipEvent = cms.untracked.vstring('ProductNotFound'),
                     allowUnscheduled = cms.untracked.bool(True)
                     )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(options.inputFiles),
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
#                            eventsToProcess = cms.untracked.VEventRange('282917:76757818-282917:76757820'),
#                            lumisToProcess = cms.untracked.VLuminosityBlockRange('282917:126'),
                            )                     


print " process source filenames %s" %(process.source) 
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
#if config["RUNONMC"]: GT = '94X_mc2017_realistic_v12'
#if config["RUNONMC"]:#
#  if config["YEAR"] == 2018 : GT = '102X_upgrade2018_realistic_v18'

#elif config["RUNONReReco"]: GT = '94X_dataRun2_ReReco_EOY17_v2'
#elif config["RUNONPromptReco"]: GT = '92X_dataRun2_2017Prompt_v11'



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

  if config["FILTEREVENTS"]:
  
   fname = ""
   if (options.inputFiles)[0].find("SingleMuon") != -1: fname = "RunLumiEventLists/SingleMuon_csc2015_Nov14.txt"
   elif (options.inputFiles)[0].find("SingleElectron") != -1: fname = "RunLumiEventLists/SingleElectron_csc2015_Nov14.txt"
   elif (options.inputFiles)[0].find("JetHT") != -1: fname = "RunLumiEventLists/JetHT_csc2015_Nov27.txt"
   else:
    print "** WARNING: EVENT LIST NOT FOUND! exiting... "
    sys.exit()
   
   print "** FILTERING EVENT LIST: %s" %fname 
   listEventsToSkip = []
   fileEventsToSkip = open(fname,"r")

   for line in fileEventsToSkip:
     cleanLine = line.rstrip()
     listEventsToSkip.append(cleanLine+"-"+cleanLine)

   rangeEventsToSkip = cms.untracked.VEventRange(listEventsToSkip)
   process.source.eventsToSkip = rangeEventsToSkip

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






# ###### Recluster MET ##########
if config["DOMETRECLUSTERING"]:

  from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection
		  
  switchJetCollection(process,
                      jetSource = cms.InputTag('ak4PFJetsCHS'),
                      jetCorrections = ('AK4PFchs', jet_corr_levels, ''),
                      genParticles = cms.InputTag('prunedGenParticles'),
                      pvSource = cms.InputTag('offlineSlimmedPrimaryVertices')
                      )
		  		
  process.patJets.addGenJetMatch = cms.bool(False) 
  process.patJets.addGenPartonMatch = cms.bool(False) 
  process.patJets.addPartonJetMatch = cms.bool(False) 
  
  from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

  #default configuration for miniAOD reprocessing, change the isData flag to run on data
  #for a full met computation, remove the pfCandColl input
  runMetCorAndUncFromMiniAOD(process,
                           isData=not(config["RUNONMC"]),
                           )
  process.patPFMetT1T2Corr.type1JetPtThreshold = cms.double(15.0)
  process.patPFMetT2Corr.type1JetPtThreshold = cms.double(15.0)
  process.slimmedMETs.t01Variation = cms.InputTag("slimmedMETs","","RECO")
  
  if config["RUNONMC"]:
    process.patPFMetT1T2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT1T2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.patPFMetT2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
    process.shiftedPatJetEnDown.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
    process.shiftedPatJetEnUp.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
			           
####### Adding HEEP id ##########

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

dataFormat=DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process,dataFormat,task=pattask)

process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)


      
my_id_modules = [
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V1_cff',
                 ]
           

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection,task=pattask)

####### Event filters ###########

##___________________________HCAL_Noise_Filter________________________________||
if config["DOHLTFILTERS"]:
 process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
 process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
 process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
 ##___________________________BadChargedCandidate_Noise_Filter________________________________|| 
 process.load('Configuration.StandardSequences.Services_cff')
 process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
 # process.load('EXOVVNtuplizerRunII.Ntuplizer.BadChargedCandidateFilter_cfi')
 process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
 process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
 process.BadChargedCandidateFilter.debug = cms.bool(False)
 process.BadChargedCandidateSequence = cms.Sequence (process.BadChargedCandidateFilter)
 

####### Ntuplizer initialization ##########
jetsAK4 = "slimmedJets"


METS = "slimmedMETs"
METS_EGclean = "slimmedMETsEGClean"
METS_MEGclean = "slimmedMETsMuEGClean"
METS_uncorr = "slimmedMETsUncorrected"

if config["DOMETRECLUSTERING"]: jetsAK4 = "selectedPatJets"
if config["USENOHF"]: METS = "slimmedMETsNoHF"  

##___________________ MET significance and covariance matrix ______________________##

if config["DOMETSVFIT"]:
  print "Using event pfMET covariance for SVfit"
  process.load("RecoMET.METProducers.METSignificance_cfi")
  process.load("RecoMET.METProducers.METSignificanceParams_cfi")
  pattask.add(process.METSignificance)

if config["DOMVAMET"]:
  from RecoMET.METPUSubtraction.jet_recorrections import recorrectJets
  recorrectJets(process, isData=True)
  
  from RecoMET.METPUSubtraction.MVAMETConfiguration_cff import runMVAMET
  runMVAMET( process, jetCollectionPF="patJetsReapplyJEC")
  process.MVAMET.srcLeptons  = cms.VInputTag("slimmedMuons", "slimmedElectrons", "slimmedTaus")
  process.MVAMET.requireOS = cms.bool(False)

##___________________ taus ______________________##

TAUS = ""
BOOSTEDTAUS = ""

TAUS = "NewTauIDsEmbedded"

#  TAUS = "slimmedTaus"

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
    GT = '102X_mc2017_realistic_v6'
  elif ("Autumn18" in options.RunPeriod):
    JECprefix = "Autumn18_V8"
    GT = '102X_upgrade2018_realistic_v18'
  #jecAK8chsUncFile = "JEC/%s_MC_Uncertainty_AK8PFchs.txt"%(JECprefix)
  jecAK4chsUncFile = "JEC/%s_MC_Uncertainty_AK4PFchs.txt"%(JECprefix)
 



else : #Data

   JEC_runDependent_suffix= ""
   if ("2017" in options.RunPeriod):
     if ("Run2017B" in  options.RunPeriod): JEC_runDependent_suffix= "B"
     elif ("Run2017C" in  options.RunPeriod): JEC_runDependent_suffix= "C"
     elif ("Run2017D" in  options.RunPeriod): JEC_runDependent_suffix= "D"
     elif ("Run2017E" in  options.RunPeriod): JEC_runDependent_suffix= "E"
     elif ("Run2017F" in  options.RunPeriod): JEC_runDependent_suffix= "F"
     
     JECprefix = "Fall17_17Nov2017"+JEC_runDependent_suffix+"_V32"
     GT = '102X_dataRun2_v8'
     
     
   elif ("2016" in options.RunPeriod):
     if ("Run2016D" in  options.RunPeriod or "Run2016B" in  options.RunPeriod  or "Run2016C" in  options.RunPeriod  ): JEC_runDependent_suffix= "ABC"
     elif ("Run2016E" in  options.RunPeriod): JEC_runDependent_suffix= "EF"
     elif ("Run2016G" in  options.RunPeriod): JEC_runDependent_suffix= "GH"
     elif ("Run2016F" in  options.RunPeriod and  not options.runUpToEarlyF): JEC_runDependent_suffix= "GH"
     elif ("Run2016F" in  options.RunPeriod and   options.runUpToEarlyF): JEC_runDependent_suffix= "EF"


     JECprefix = "Summer16_07Aug2017"+JEC_runDependent_suffix+"_V11"
    
   elif ("2018" in options.RunPeriod):
     if ("Run2018A" in  options.RunPeriod ): JEC_runDependent_suffix= "A"
     elif ("Run2018B" in  options.RunPeriod): JEC_runDependent_suffix= "B"
     elif ("Run2018C" in  options.RunPeriod): JEC_runDependent_suffix= "C"
     elif ("Run2018D" in  options.RunPeriod): JEC_runDependent_suffix= "D"
    

     JECprefix = "Autumn18_Run"+JEC_runDependent_suffix+"_V8"
     GT ='94X_dataRun2_v10'
   #jecAK8chsUncFile = "JEC/%s_DATA_Uncertainty_AK8PFchs.txt"%(JECprefix)
   jecAK4chsUncFile = "JEC/%s_DATA_Uncertainty_AK4PFchs.txt"%(JECprefix)
   GT = '102X_dataRun2_v10' 
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
      			    

                                                                       
################## Ntuplizer ###################
process.ntuplizer = cms.EDAnalyzer("Ntuplizer",
    runOnMC	      = cms.bool(config["RUNONMC"]),
    doGenParticles    = cms.bool(config["DOGENPARTICLES"]),
    doGenEvent	      = cms.bool(config["DOGENEVENT"]),
    doPileUp	      = cms.bool(config["DOPILEUP"]),
    doElectrons       = cms.bool(config["DOELECTRONS"]),
    doMuons	      = cms.bool(config["DOMUONS"]),
    doTaus	      = cms.bool(config["DOTAUS"]),
    doJpsiMu	      = cms.bool(config["DOJPSIMU"]),
    doJpsiEle	      = cms.bool(config["DOJPSIELE"]),
    doVertices	      = cms.bool(config["DOVERTICES"]),
    doTriggerDecisions= cms.bool(config["DOTRIGGERDECISIONS"]),
    doTriggerObjects  = cms.bool(config["DOTRIGGEROBJECTS"]),
    doHltFilters      = cms.bool(config["DOHLTFILTERS"]),
    doMissingEt       = cms.bool(config["DOMISSINGET"]),
    doMETSVFIT        = cms.bool(config["DOMETSVFIT"]),
    doMVAMET          = cms.bool(config["DOMVAMET"]),
    doMultipleTauMVAversions = cms.bool(config["DOMULTIPLETAUMVAVERSIONS"]),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    ebRecHits = cms.InputTag("reducedEgamma","reducedEBRecHits"),

#    eleHEEPId51Map = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV51"),
#    eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
#    eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
#    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
#    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
#    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),

    eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-veto"),
    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-loose"),
    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-medium"),
    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Fall17-94X-V1-tight"),

    eleHLTIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronHLTPreselection-Summer16-V1"), 
    eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70"),
                                   
    eleMVAMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp90"),
    eleMVATightIdMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Fall17-noIso-V1-wp80"),
    mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
    mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),
    dupCluster          = cms.InputTag("particleFlowEGammaGSFixed:dupECALClusters"),
    hitsNotReplaced     = cms.InputTag("ecalMultiAndGSGlobalRecHitEB:hitsNotReplaced"),
    taus = cms.InputTag(TAUS),
    #tausBoostedTau = cms.InputTag(BOOSTEDTAUS),
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

    packedpfcandidates = cms.InputTag('packedPFCandidates')
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



#### Needed to load the new tau ID payloads ####

from RecoTauTag.RecoTau.TauDiscriminatorTools import noPrediscriminants
from RecoTauTag.RecoTau.PATTauDiscriminationByMVAIsolationRun2_cff import patDiscriminationByIsolationMVArun2v1raw, patDiscriminationByIsolationMVArun2v1VLoose



tauIdDiscrMVA_trainings_run2_2017 = {
  'tauIdMVAIsoDBoldDMwLT2017' : "tauIdMVAIsoDBoldDMwLT2017",
  }
tauIdDiscrMVA_WPs_run2_2017 = {
  'tauIdMVAIsoDBoldDMwLT2017' : {
    'Eff95' : "DBoldDMwLTEff95",
    'Eff90' : "DBoldDMwLTEff90",
    'Eff80' : "DBoldDMwLTEff80",
    'Eff70' : "DBoldDMwLTEff70",
    'Eff60' : "DBoldDMwLTEff60",
    'Eff50' : "DBoldDMwLTEff50",
    'Eff40' : "DBoldDMwLTEff40"
    }
  }
tauIdDiscrMVA_2017_version = "v1"


def loadMVA_WPs_run2_2017(process):
                print "LoadMVA_WPs_run2_2017: performed::::"
                #global cms



		for training, gbrForestName in tauIdDiscrMVA_trainings_run2_2017.items():

			process.loadRecoTauTagMVAsFromPrepDB.toGet.append(
				cms.PSet(
					record = cms.string('GBRWrapperRcd'),
					tag = cms.string("RecoTauTag_%s%s" % (gbrForestName, tauIdDiscrMVA_2017_version)),
					label = cms.untracked.string("RecoTauTag_%s%s" % (gbrForestName, tauIdDiscrMVA_2017_version))
				)
			)

			for WP in tauIdDiscrMVA_WPs_run2_2017[training].keys():
				process.loadRecoTauTagMVAsFromPrepDB.toGet.append(
					cms.PSet(
						record = cms.string('PhysicsTGraphPayloadRcd'),
						tag = cms.string("RecoTauTag_%s%s_WP%s" % (gbrForestName, tauIdDiscrMVA_2017_version, WP)),
						label = cms.untracked.string("RecoTauTag_%s%s_WP%s" % (gbrForestName, tauIdDiscrMVA_2017_version, WP))
					)
				)

			process.loadRecoTauTagMVAsFromPrepDB.toGet.append(
				cms.PSet(
					record = cms.string('PhysicsTFormulaPayloadRcd'),
					tag = cms.string("RecoTauTag_%s%s_mvaOutput_normalization" % (gbrForestName, tauIdDiscrMVA_2017_version)),
					label = cms.untracked.string("RecoTauTag_%s%s_mvaOutput_normalization" % (gbrForestName, tauIdDiscrMVA_2017_version))
				)
)


####### Tau new MVA ##########

from RecoTauTag.RecoTau.TauDiscriminatorTools import noPrediscriminants
process.load('RecoTauTag.Configuration.loadRecoTauTagMVAsFromPrepDB_cfi')
from RecoTauTag.RecoTau.PATTauDiscriminationByMVAIsolationRun2_cff import *
loadMVA_WPs_run2_2017(process)

process.rerunDiscriminationByIsolationMVArun2v1raw = patDiscriminationByIsolationMVArun2v1raw.clone(
   PATTauProducer = cms.InputTag('slimmedTaus'),
   Prediscriminants = noPrediscriminants,
   loadMVAfromDB = cms.bool(True),
   mvaName = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1"),# training with 2017 MC_v1 for oldDM # name of the training you want to use
   mvaOpt = cms.string("DBoldDMwLTwGJ"), # option you want to use for your training (i.e., which variables are used to compute the BDT score)
   requireDecayMode = cms.bool(True),
   verbosity = cms.int32(0)
)


process.rerunDiscriminationByIsolationMVArun2v1VLoose = patDiscriminationByIsolationMVArun2v1VLoose.clone(
   PATTauProducer = cms.InputTag('slimmedTaus'),    
   Prediscriminants = noPrediscriminants,
   toMultiplex = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1raw'),
   key = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1raw:category'),
   loadMVAfromDB = cms.bool(True),
   mvaOutput_normalization = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_mvaOutput_normalization"), # normalization fo the training you want to use
   mapping = cms.VPSet(
      cms.PSet(
         category = cms.uint32(0),
         cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff90"), # this is the name of the working point you want to use
         variable = cms.string("pt"),
      )
   )
)

# here we produce all the other working points for the training
process.rerunDiscriminationByIsolationMVArun2v1VVLoose = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
process.rerunDiscriminationByIsolationMVArun2v1VVLoose.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff95")
process.rerunDiscriminationByIsolationMVArun2v1Loose = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
process.rerunDiscriminationByIsolationMVArun2v1Loose.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff80")
process.rerunDiscriminationByIsolationMVArun2v1Medium = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
process.rerunDiscriminationByIsolationMVArun2v1Medium.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff70")
process.rerunDiscriminationByIsolationMVArun2v1Tight = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
process.rerunDiscriminationByIsolationMVArun2v1Tight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff60")
process.rerunDiscriminationByIsolationMVArun2v1VTight = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
process.rerunDiscriminationByIsolationMVArun2v1VTight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff50")
process.rerunDiscriminationByIsolationMVArun2v1VVTight = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
process.rerunDiscriminationByIsolationMVArun2v1VVTight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff40")

# this sequence has to be included in your cms.Path() before your analyzer which accesses the new variables is called.
process.rerunMvaIsolation2SeqRun2 = cms.Sequence(
   process.rerunDiscriminationByIsolationMVArun2v1raw
   *process.rerunDiscriminationByIsolationMVArun2v1VLoose
   *process.rerunDiscriminationByIsolationMVArun2v1VVLoose
   *process.rerunDiscriminationByIsolationMVArun2v1Loose
   *process.rerunDiscriminationByIsolationMVArun2v1Medium
   *process.rerunDiscriminationByIsolationMVArun2v1Tight
   *process.rerunDiscriminationByIsolationMVArun2v1VTight
   *process.rerunDiscriminationByIsolationMVArun2v1VVTight
)


if  not config["DOMULTIPLETAUMVAVERSIONS"]: ## skip this partial inclusion, if you are also using a *db file and just do it all together at the end.
# embed new id's into new tau collection
  embedID = cms.EDProducer("PATTauIDEmbedder",
     src = cms.InputTag('slimmedTaus'),
     tauIDSources = cms.PSet(
        byIsolationMVArun2v1DBoldDMwLTrawNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1raw'),
        byVLooseIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1VLoose'),
        byVVLooseIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1VVLoose'),
        byLooseIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1Loose'),
        byMediumIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1Medium'),
        byTightIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1Tight'),
        byVTightIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1VTight'),
        byVVTightIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1VVTight'),
        ),
     )

  setattr(process, "NewTauIDsEmbedded", embedID)

  ##===== V2 of the MVA isolation from qslite file:
  #Do_Tau_isolation_FromDB =True
if  config["DOMULTIPLETAUMVAVERSIONS"]:


  from CondCore.DBCommon.CondDBSetup_cfi import *


  process.loadRecoTauTagMVAsFromPrepDB2 = cms.ESSource("PoolDBESSource",

                                                       CondDBSetup,
                                                       DumpStat = cms.untracked.bool(False),
                                                       toGet = cms.VPSet(),                                             
                                                       #  connect = cms.string("frontier://FrontierPrep/CMS_COND_PHYSICSTOOLS") # prep database
                                                       #connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS') # prod database
                                                       connect = cms.string('sqlite_file:RecoTauTag_MVAs_2018Mar15.db')
  )

  #process.loadRecoTauTagMVAsFromPrepDB2

  tauIdDiscrMVA_trainings_run2_2017_v2 = {
    'tauIdMVAIsoDBoldDMwLT2017' : "tauIdMVAIsoDBoldDMwLT2017",
    'tauIdMVAIsoDBnewDMwLT2017' : "tauIdMVAIsoDBnewDMwLT2017",
    'tauIdMVAIsoDBoldDMdR0p3wLT2017' : "tauIdMVAIsoDBoldDMdR0p3wLT2017",
    }
  tauIdDiscrMVA_WPs_run2_2017_v2 = {
    'tauIdMVAIsoDBoldDMwLT2017' : {
      'Eff95' : "DBoldDMwLTEff95",
      'Eff90' : "DBoldDMwLTEff90",
      'Eff80' : "DBoldDMwLTEff80",
      'Eff70' : "DBoldDMwLTEff70",
      'Eff60' : "DBoldDMwLTEff60",
      'Eff50' : "DBoldDMwLTEff50",
      'Eff40' : "DBoldDMwLTEff40"
      } ,
    'tauIdMVAIsoDBnewDMwLT2017' : {
      'Eff95' : "DBnewDMwLTEff95",
      'Eff90' : "DBnewDMwLTEff90",
      'Eff80' : "DBnewDMwLTEff80",
      'Eff70' : "DBnewDMwLTEff70",
      'Eff60' : "DBnewDMwLTEff60",
      'Eff50' : "DBnewDMwLTEff50",
      'Eff40' : "DBnewDMwLTEff40"
      } ,
    'tauIdMVAIsoDBoldDMdR0p3wLT2017' : {
      'Eff95' : "DBoldDMdR0p3wLTEff95",
      'Eff90' : "DBoldDMdR0p3wLTEff90",
      'Eff80' : "DBoldDMdR0p3wLTEff80",
      'Eff70' : "DBoldDMdR0p3wLTEff70",
      'Eff60' : "DBoldDMdR0p3wLTEff60",
      'Eff50' : "DBoldDMdR0p3wLTEff50",
      'Eff40' : "DBoldDMdR0p3wLTEff40"
     }
    }
  tauIdDiscrMVA_2017_v2_version = "v2"

  tauIdDiscrMVA_mvaOutput_normalizations_v2 = {
    'tauIdMVAIsoDBoldDMwLT2017' : "tauIdMVAIsoDBoldDMwLT_2017v2_dR0p5",
    'tauIdMVAIsoDBnewDMwLT2017' : "tauIdMVAIsoDBnewDMwLT_2017v2_dR0p5",
    'tauIdMVAIsoDBoldDMdR0p3wLT2017' : "tauIdMVAIsoDBoldDMwLT_2017v2_dR0p3",
    }
  def loadMVA_WPs_run2_2017_v2(process):
                  print "LoadMVA_WPs_run2_2017_v2: performed::::"
                  #global cms


                  for training, gbrForestName in tauIdDiscrMVA_trainings_run2_2017_v2.items():
                          print " printing tauIdDiscrMVA_trainings_run2_2017_v2.items(),  training %s , gbrForestName %s" %(training, gbrForestName)
                          print " tauIdDiscrMVA_mvaOutput_normalizations_v2[training] %s"% tauIdDiscrMVA_mvaOutput_normalizations_v2[training]
                          print "printato"
 
                          process.loadRecoTauTagMVAsFromPrepDB2.toGet.append(
                                  cms.PSet(
                                          record = cms.string('GBRWrapperRcd'),
                                          tag = cms.string("RecoTauTag_%s%s" % (gbrForestName, tauIdDiscrMVA_2017_v2_version)),
                                          label = cms.untracked.string("RecoTauTag_%s%s" % (gbrForestName, tauIdDiscrMVA_2017_v2_version))
                                  )
                          )

                          for WP in tauIdDiscrMVA_WPs_run2_2017_v2[training].keys():
                                  process.loadRecoTauTagMVAsFromPrepDB2.toGet.append(
                                          cms.PSet(
                                                  record = cms.string('PhysicsTGraphPayloadRcd'),
                                                  tag = cms.string("RecoTauTag_%s%s_WP%s" % (gbrForestName, tauIdDiscrMVA_2017_v2_version, WP)),
                                                  label = cms.untracked.string("RecoTauTag_%s%s_WP%s" % (gbrForestName, tauIdDiscrMVA_2017_v2_version, WP))
                                          )
                                  )

                          process.loadRecoTauTagMVAsFromPrepDB2.toGet.append(
                                  cms.PSet(
                                          record = cms.string('PhysicsTFormulaPayloadRcd'),
                                          tag = cms.string("RecoTauTag_%s_mvaOutput_normalization" % (tauIdDiscrMVA_mvaOutput_normalizations_v2[training])),
                                          label = cms.untracked.string("RecoTauTag_%s_mvaOutput_normalization" % (tauIdDiscrMVA_mvaOutput_normalizations_v2[training]))### https://github.com/cgalloni/cmssw/blob/CMSSW_8_0_X/RecoTauTag/Configuration/python/loadRecoTauTagMVAsFromPrepDB_cfi.py

                                  )
  )







  #======= for v2


  loadMVA_WPs_run2_2017_v2(process)

  process.rerunDiscriminationByIsolationMVArun2v2raw = patDiscriminationByIsolationMVArun2v1raw.clone(
     PATTauProducer = cms.InputTag('slimmedTaus'),
     Prediscriminants = noPrediscriminants,
     loadMVAfromDB = cms.bool(True),
     mvaName = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2"),# training with 2017 MC_v1 for oldDM # name of the training you want to use
     mvaOpt = cms.string("DBoldDMwLTwGJ"), # option you want to use for your training (i.e., which variables are used to compute the BDT score)
     requireDecayMode = cms.bool(True),
     verbosity = cms.int32(0)
  )


  process.rerunDiscriminationByIsolationMVArun2v2VLoose = patDiscriminationByIsolationMVArun2v1VLoose.clone(
     PATTauProducer = cms.InputTag('slimmedTaus'),    
     Prediscriminants = noPrediscriminants,
     toMultiplex = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2raw'),
     key = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2raw:category'),
     loadMVAfromDB = cms.bool(True),
     mvaOutput_normalization = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT_2017v2_dR0p5_mvaOutput_normalization"), # normalization fo the training you want to use
     mapping = cms.VPSet(
        cms.PSet(
           category = cms.uint32(0),
           cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff90"), # this is the name of the working point you want to use
           variable = cms.string("pt"),
        )
     )
  )

  # here we produce all the other working points for the training
  process.rerunDiscriminationByIsolationMVArun2v2VVLoose = process.rerunDiscriminationByIsolationMVArun2v2VLoose.clone()
  process.rerunDiscriminationByIsolationMVArun2v2VVLoose.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff95")
  process.rerunDiscriminationByIsolationMVArun2v2Loose = process.rerunDiscriminationByIsolationMVArun2v2VLoose.clone()
  process.rerunDiscriminationByIsolationMVArun2v2Loose.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff80")
  process.rerunDiscriminationByIsolationMVArun2v2Medium = process.rerunDiscriminationByIsolationMVArun2v2VLoose.clone()
  process.rerunDiscriminationByIsolationMVArun2v2Medium.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff70")
  process.rerunDiscriminationByIsolationMVArun2v2Tight = process.rerunDiscriminationByIsolationMVArun2v2VLoose.clone()
  process.rerunDiscriminationByIsolationMVArun2v2Tight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff60")
  process.rerunDiscriminationByIsolationMVArun2v2VTight = process.rerunDiscriminationByIsolationMVArun2v2VLoose.clone()
  process.rerunDiscriminationByIsolationMVArun2v2VTight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff50")
  process.rerunDiscriminationByIsolationMVArun2v2VVTight = process.rerunDiscriminationByIsolationMVArun2v2VLoose.clone()
  process.rerunDiscriminationByIsolationMVArun2v2VVTight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v2_WPEff40")

  # this sequence has to be included in your cms.Path() before your analyzer which accesses the new variables is called.
  process.rerunMvaIsolation2SeqRun2_2 = cms.Sequence(
     process.rerunDiscriminationByIsolationMVArun2v2raw
     *process.rerunDiscriminationByIsolationMVArun2v2VLoose
     *process.rerunDiscriminationByIsolationMVArun2v2VVLoose
     *process.rerunDiscriminationByIsolationMVArun2v2Loose
     *process.rerunDiscriminationByIsolationMVArun2v2Medium
     *process.rerunDiscriminationByIsolationMVArun2v2Tight
     *process.rerunDiscriminationByIsolationMVArun2v2VTight
     *process.rerunDiscriminationByIsolationMVArun2v2VVTight
  )

  #========== MVA newDM
  process.rerunDiscriminationByIsolationMVArun2v2newDMraw = patDiscriminationByIsolationMVArun2v1raw.clone(
     PATTauProducer = cms.InputTag('slimmedTaus'),
     Prediscriminants = noPrediscriminants,
     loadMVAfromDB = cms.bool(True),
     mvaName = cms.string("RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2"),# training with 2017 MC_v1 for newDM # name of the training you want to use                                                                                                 
     mvaOpt = cms.string("DBnewDMwLTwGJ"), # option you want to use for your training (i.e., which variables are used to compute the BDT score)                                                                                                 
     requireDecayMode = cms.bool(True),
     verbosity = cms.int32(0)
  )


  process.rerunDiscriminationByIsolationMVArun2v2newDMVLoose = patDiscriminationByIsolationMVArun2v1VLoose.clone(
     PATTauProducer = cms.InputTag('slimmedTaus'),
     Prediscriminants = noPrediscriminants,
     toMultiplex = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2newDMraw'),
     key = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2newDMraw:category'),
     loadMVAfromDB = cms.bool(True),
     mvaOutput_normalization = cms.string("RecoTauTag_tauIdMVAIsoDBnewDMwLT_2017v2_dR0p5_mvaOutput_normalization"), # normalization fo the training you want to use                                                                             
     mapping = cms.VPSet(
        cms.PSet(
           category = cms.uint32(0),
           cut = cms.string("RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff90"), # this is the name of the working point you want to use                                                                                                          
           variable = cms.string("pt"),
        )
     )
  )

  # here we produce all the other working points for the training                                                                                                                                                                               
  process.rerunDiscriminationByIsolationMVArun2v2newDMVVLoose = process.rerunDiscriminationByIsolationMVArun2v2newDMVLoose.clone()
  process.rerunDiscriminationByIsolationMVArun2v2newDMVVLoose.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff95")
  process.rerunDiscriminationByIsolationMVArun2v2newDMLoose = process.rerunDiscriminationByIsolationMVArun2v2newDMVLoose.clone()
  process.rerunDiscriminationByIsolationMVArun2v2newDMLoose.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff80")
  process.rerunDiscriminationByIsolationMVArun2v2newDMMedium = process.rerunDiscriminationByIsolationMVArun2v2newDMVLoose.clone()
  process.rerunDiscriminationByIsolationMVArun2v2newDMMedium.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff70")
  process.rerunDiscriminationByIsolationMVArun2v2newDMTight = process.rerunDiscriminationByIsolationMVArun2v2newDMVLoose.clone()
  process.rerunDiscriminationByIsolationMVArun2v2newDMTight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff60")
  process.rerunDiscriminationByIsolationMVArun2v2newDMVTight = process.rerunDiscriminationByIsolationMVArun2v2newDMVLoose.clone()
  process.rerunDiscriminationByIsolationMVArun2v2newDMVTight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff50")
  process.rerunDiscriminationByIsolationMVArun2v2newDMVVTight = process.rerunDiscriminationByIsolationMVArun2v2newDMVLoose.clone()
  process.rerunDiscriminationByIsolationMVArun2v2newDMVVTight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBnewDMwLT2017v2_WPEff40")

  # this sequence has to be included in your cms.Path() before your analyzer which accesses the new variables is called.                                                                                                                        
  process.rerunMvaIsolation2SeqRun2_2_newDM = cms.Sequence(
     process.rerunDiscriminationByIsolationMVArun2v2newDMraw
     *process.rerunDiscriminationByIsolationMVArun2v2newDMVLoose
     *process.rerunDiscriminationByIsolationMVArun2v2newDMVVLoose
     *process.rerunDiscriminationByIsolationMVArun2v2newDMLoose
     *process.rerunDiscriminationByIsolationMVArun2v2newDMMedium
     *process.rerunDiscriminationByIsolationMVArun2v2newDMTight
     *process.rerunDiscriminationByIsolationMVArun2v2newDMVTight
     *process.rerunDiscriminationByIsolationMVArun2v2newDMVVTight
  )


  #========== MVA oldDM _dR0p3

  process.rerunDiscriminationByIsolationMVArun2v2dR0p3raw = patDiscriminationByIsolationMVArun2v1raw.clone(
     PATTauProducer = cms.InputTag('slimmedTaus'),
     Prediscriminants = noPrediscriminants,
     loadMVAfromDB = cms.bool(True),
     mvaName = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2"),#RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2 training with 2017 MC_v1 for oldDM # name of the training you want to use
     mvaOpt = cms.string("DBoldDMwLTwGJ"), # option you want to use for your training (i.e., which variables are used to compute the BDT score)
     requireDecayMode = cms.bool(True),
     verbosity = cms.int32(0)
  )


  process.rerunDiscriminationByIsolationMVArun2v2dR0p3VLoose = patDiscriminationByIsolationMVArun2v1VLoose.clone(
     PATTauProducer = cms.InputTag('slimmedTaus'),    
     Prediscriminants = noPrediscriminants,
     toMultiplex = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2dR0p3raw'),
     key = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2dR0p3raw:category'),
     loadMVAfromDB = cms.bool(True),
     mvaOutput_normalization = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT_2017v2_dR0p3_mvaOutput_normalization"), #RecoTauTag_tauIdMVAIsoDBoldDMwLT_2017v2_dR0p3_mvaOutput_normalization normalization fo the training you want to use
     mapping = cms.VPSet(
        cms.PSet(
           category = cms.uint32(0),
           cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff90"), # this is the name of the working point you want to use
           variable = cms.string("pt"),
        )
     )
  )

  # here we produce all the other working points for the training
  process.rerunDiscriminationByIsolationMVArun2v2dR0p3VVLoose = process.rerunDiscriminationByIsolationMVArun2v2dR0p3VLoose.clone()
  process.rerunDiscriminationByIsolationMVArun2v2dR0p3VVLoose.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff95")
  process.rerunDiscriminationByIsolationMVArun2v2dR0p3Loose = process.rerunDiscriminationByIsolationMVArun2v2dR0p3VLoose.clone()
  process.rerunDiscriminationByIsolationMVArun2v2dR0p3Loose.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff80")
  process.rerunDiscriminationByIsolationMVArun2v2dR0p3Medium = process.rerunDiscriminationByIsolationMVArun2v2dR0p3VLoose.clone()
  process.rerunDiscriminationByIsolationMVArun2v2dR0p3Medium.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff70")
  process.rerunDiscriminationByIsolationMVArun2v2dR0p3Tight = process.rerunDiscriminationByIsolationMVArun2v2dR0p3VLoose.clone()
  process.rerunDiscriminationByIsolationMVArun2v2dR0p3Tight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff60")
  process.rerunDiscriminationByIsolationMVArun2v2dR0p3VTight = process.rerunDiscriminationByIsolationMVArun2v2dR0p3VLoose.clone()
  process.rerunDiscriminationByIsolationMVArun2v2dR0p3VTight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff50")
  process.rerunDiscriminationByIsolationMVArun2v2dR0p3VVTight = process.rerunDiscriminationByIsolationMVArun2v2dR0p3VLoose.clone()
  process.rerunDiscriminationByIsolationMVArun2v2dR0p3VVTight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMdR0p3wLT2017v2_WPEff40")

  process.rerunMvaIsolation2SeqRun2_2_dR0p3 = cms.Sequence(
    process.rerunDiscriminationByIsolationMVArun2v2dR0p3raw
    *process.rerunDiscriminationByIsolationMVArun2v2dR0p3VLoose
    *process.rerunDiscriminationByIsolationMVArun2v2dR0p3VVLoose
    *process.rerunDiscriminationByIsolationMVArun2v2dR0p3Loose
    *process.rerunDiscriminationByIsolationMVArun2v2dR0p3Medium
    *process.rerunDiscriminationByIsolationMVArun2v2dR0p3Tight
    *process.rerunDiscriminationByIsolationMVArun2v2dR0p3VTight
    *process.rerunDiscriminationByIsolationMVArun2v2dR0p3VVTight
    )

  # this se






  # embed new id's into new tau collection
  embedID = cms.EDProducer("PATTauIDEmbedder",
     src = cms.InputTag('slimmedTaus'),
     tauIDSources = cms.PSet(
        byIsolationMVArun2v1DBoldDMwLTrawNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1raw'),
        byVLooseIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1VLoose'),
        byVVLooseIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1VVLoose'),
        byLooseIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1Loose'),
        byMediumIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1Medium'),
        byTightIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1Tight'),
        byVTightIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1VTight'),
        byVVTightIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1VVTight'),
        byIsolationMVArun2v2DBoldDMwLTrawNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2raw'),
        byVLooseIsolationMVArun2v2DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2VLoose'),
        byVVLooseIsolationMVArun2v2DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2VVLoose'),
        byLooseIsolationMVArun2v2DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2Loose'),
        byMediumIsolationMVArun2v2DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2Medium'),
        byTightIsolationMVArun2v2DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2Tight'),
        byVTightIsolationMVArun2v2DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2VTight'),
        byVVTightIsolationMVArun2v2DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2VVTight'), 

        byIsolationMVArun2v2DBnewDMwLTrawNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2newDMraw'),
        byVLooseIsolationMVArun2v2DBnewDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2newDMVLoose'),
        byVVLooseIsolationMVArun2v2DBnewDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2newDMVVLoose'),
        byLooseIsolationMVArun2v2DBnewDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2newDMLoose'),
        byMediumIsolationMVArun2v2DBnewDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2newDMMedium'),
        byTightIsolationMVArun2v2DBnewDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2newDMTight'),
        byVTightIsolationMVArun2v2DBnewDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2newDMVTight'),
        byVVTightIsolationMVArun2v2DBnewDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2newDMVVTight'),

        byIsolationMVArun2v2DBoldDMdR0p3wLTrawNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2dR0p3raw'),
        byVLooseIsolationMVArun2v2DBoldDMdR0p3wLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2dR0p3VLoose'),
        byVVLooseIsolationMVArun2v2DBoldDMdR0p3wLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2dR0p3VVLoose'),
        byLooseIsolationMVArun2v2DBoldDMdR0p3wLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2dR0p3Loose'),
        byMediumIsolationMVArun2v2DBoldDMdR0p3wLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2dR0p3Medium'),
        byTightIsolationMVArun2v2DBoldDMdR0p3wLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2dR0p3Tight'),
        byVTightIsolationMVArun2v2DBoldDMdR0p3wLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2dR0p3VTight'),
        byVVTightIsolationMVArun2v2DBoldDMdR0p3wLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v2dR0p3VVTight'),
        ),
     )

  setattr(process, "NewTauIDsEmbedded", embedID)


process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder')
)


####### Final path ##########
process.p = cms.Path()
if config["DOHLTFILTERS"]:
 process.p += process.HBHENoiseFilterResultProducer
 process.p += process.BadChargedCandidateSequence

# For new MVA ID !
process.p += process.rerunMvaIsolation2SeqRun2 
if config["DOMULTIPLETAUMVAVERSIONS"]:
  process.p += process.rerunMvaIsolation2SeqRun2_2
  process.p += process.rerunMvaIsolation2SeqRun2_2_dR0p3
  process.p += process.rerunMvaIsolation2SeqRun2_2_newDM

process.p += getattr(process, "NewTauIDsEmbedded")
# For new MVA ID END!
process.p += process.ecalBadCalibReducedMINIAODFilter
process.p += process.ntuplizer
process.p.associate(pattask)

print pattask

#  LocalWords:  tauIdMVAIsoDBoldDMwLT
