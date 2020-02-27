#cmsRun  config_generic_opt_skimmed.py  RunPeriod="Fall17" # for MC
#cmsRun  config_generic_opt_skimmed.py  RunPeriod="Run2017B" # for Data
#cmsRun  config_generic_opt_skimmed.py  RunPeriod="Autumn18" # for MC from 2018


###### Process initialization ##########

import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuple")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')


#from EXOVVNtuplizerRunII.Ntuplizer.ntuplizerOptions_data_cfi import config
from EXOVVNtuplizerRunII.Ntuplizer.ntuplizerOptions_generic_cfi import config

# change from its original value
#config["DZCUT"] = 0.25
#config["FSIGCUT"] = 3
#config["VPROBCUT"] = 0.1
#config["DNNCUT"] = 0.2

import sys
args = sys.argv

if len(args)!=4:
  print 'input file name missing'
  sys.exit(0)
else:
  inputFile = args[2]
  outputFile = args[3]
  
RunPeriod='Autumn18'

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outputFile)
                                   )


####### Config parser ##########

#import FWCore.ParameterSet.VarParsing as VarParsing
#
#options = VarParsing.VarParsing ('analysis')
#
#options.register( 'RunPeriod',
#                  '',
#                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
#                  VarParsing.VarParsing.varType.string,          # string, int, or float
#                  "RunNumber (Default Run2017B)")
#
#options.register( 'runUpToEarlyF',
#                  'false',
#                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list                                                                                                                                 
#                  VarParsing.VarParsing.varType.bool,          # string, int, or float                                                                                                                        
#                  "false")# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAboutPythonConfigFile
#


#options.maxEvents = 1000
maxEvents = -1

#data file
     
#options.inputFiles ='file:/work/pmatorra/JpsiAnalysis/2018/CMSSW_10_2_10/src/EXOVVNtuplizerRunII/Ntuplizer/miniAOD_99.root'
#options.inputFiles = '/store/user/cgalloni/BJpsiX_MuMu_230819/Autumn18_10_2_9_miniAOD/190823_131752/0000/miniAOD_57.root'
#options.inputFiles = '/store/user/cgalloni/BcJpsiMuNu_020519/Fall18_10_2_9-MINIAODSIM_noDuplCheck/190506_100026/0000/miniAOD_99.root'
#options.inputFiles = '/store/mc/RunIIAutumn18MiniAOD/BuToKJpsi_ToMuMu_probefilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/120000/FDD42175-87DC-D648-860B-F240C5E2CB91.root'
#options.inputFiles ='/store/data/Run2018B/Charmonium/MINIAOD/17Sep2018-v1/10000/02CFE87F-7C17-1340-8300-FDA86C16D58C.root'
#options.inputFiles ='/store/user/cgalloni/BJpsiX_MuMu_270819/Autumn18_10_2_9_miniAOD/190827_143312/0005/miniAOD_5000.root'
#options.inputFiles = '/store/user/cgalloni/BcJpsiTauNu_020519/Fall18_10_2_9-MINIAODSIM_noDuplCheck_020519/190505_141436/0000/miniAOD_99.root'
#options.inputFiles = '/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
#options.inputFiles = '/store/mc/RunIIAutumn18MiniAOD/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/260000/E342DC18-5142-1545-B077-D4405CE0BF05.root'
#options.inputFiles = '/store/mc/RunIIAutumn18MiniAOD/BsToTauTau_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v1/280000/FD3C4F1A-A2ED-AE43-A68E-76E59838E891.root'
#options.inputFiles = 'file:/scratch/ytakahas/E342DC18-5142-1545-B077-D4405CE0BF05.root'
#options.inputFiles = '/store/data/Run2018A/ParkingBPH1/MINIAOD/05May2019-v1/280001/1CE00F5E-9F52-4045-BC7C-C9178E71DB9E.root'
#options.inputFiles = 'file:/work/ytakahas/work/NtuplizerProd/CMSSW_10_6_8/src/EXOVVNtuplizerRunII/Ntuplizer/E30B9F61-DB46-0446-9B6A-2E21B806D4CE.root'
#options.inputFiles = 'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/USER/v1/A2581C35-B56A-5046-A02D-8C8C1562DEEE.root'
#options.inputFiles = 'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/UpsilonToTauTau_inclusive_5f_Pythia_LO/USER/v1/UpsilonToTauTau_3prong_miniaod_part0.root'

#print 'File =', options.Files

#inputFiles = 'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/' + inputfilename

#options.inputFiles = [
#  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/UpsilonToTauTau_inclusive_5f_Pythia_LO/USER/v1/UpsilonToTauTau_3prong_miniaod_part0.root',
#  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/UpsilonToTauTau_inclusive_5f_Pythia_LO/USER/v1/UpsilonToTauTau_3prong_miniaod_part1.root',
#  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/UpsilonToTauTau_inclusive_5f_Pythia_LO/USER/v1/UpsilonToTauTau_3prong_miniaod_part2.root',
#  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/UpsilonToTauTau_inclusive_5f_Pythia_LO/USER/v1/UpsilonToTauTau_3prong_miniaod_part3.root',
#  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/UpsilonToTauTau_inclusive_5f_Pythia_LO/USER/v1/UpsilonToTauTau_3prong_miniaod_part4.root',
#  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/UpsilonToTauTau_inclusive_5f_Pythia_LO/USER/v1/UpsilonToTauTau_3prong_miniaod_part5.root',
#  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/UpsilonToTauTau_inclusive_5f_Pythia_LO/USER/v1/UpsilonToTauTau_3prong_miniaod_part6.root',
#  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/UpsilonToTauTau_inclusive_5f_Pythia_LO/USER/v1/UpsilonToTauTau_3prong_miniaod_part7.root',
#  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/UpsilonToTauTau_inclusive_5f_Pythia_LO/USER/v1/UpsilonToTauTau_3prong_miniaod_part8.root',
#  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/UpsilonToTauTau_inclusive_5f_Pythia_LO/USER/v1/UpsilonToTauTau_3prong_miniaod_part9.root',
#]


#options.inputFiles = [
#  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/USER/v1/E342DC18-5142-1545-B077-D4405CE0BF05.root',
#  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/USER/v1/E30B9F61-DB46-0446-9B6A-2E21B806D4CE.root',
#  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/USER/v1/C7AA83E5-B4AC-2A41-AF7F-D5EBB5CFD8B0.root',
#  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/USER/v1/A2581C35-B56A-5046-A02D-8C8C1562DEEE.root',
#  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/USER/v1/7BF4CF93-B45C-7A4F-BBAF-D732781383EC.root',
#  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/USER/v1/770EA2B8-56DC-DE4B-9D93-31827AF593F8.root',
#  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/USER/v1/759D2273-789E-7949-B52E-AAE06B8D04FF.root',
#  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/USER/v1/58EB9039-F398-A64A-AD21-A6CC073C841B.root',
#  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/USER/v1/4808A71A-BA10-CD4E-B519-11A3468B6A93.root',
#  'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/ytakahas/mc/ytakahas/BcToJPsiTauNu_TuneCP5_13TeV-bcvegpy2-pythia8-evtgen/USER/v1/2377586D-0618-DB49-87EC-79FFF1DA2AC3.root'
#]

#options.parseArguments()

process.options  = cms.untracked.PSet( 
                     wantSummary = cms.untracked.bool(True),
                     SkipEvent = cms.untracked.vstring('ProductNotFound'),
                     allowUnscheduled = cms.untracked.bool(True)
                     )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(maxEvents) )
if config["RUNONMC"]:
  process.source = cms.Source("PoolSource",
                              fileNames = cms.untracked.vstring(inputFile),
                              duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                              
                              ) 
else:                    
  process.source = cms.Source("PoolSource",
                              fileNames = cms.untracked.vstring(inputFile),
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



print "1. RunPeriod ", RunPeriod
if RunPeriod=="" : RunPeriod=inputFile

if  config["RUNONMC"] :
  JECprefix = ""
  if ("Fall17" in RunPeriod):
    JECprefix = "Fall17_17Nov2017_V32"
    GT ='94X_mcRun2_asymptotic_v3'
  elif ("Summer16" in RunPeriod):
    JECprefix = "Summer16_07Aug2017_V11"
    GT = '94X_mc2017_realistic_v17'
  elif (("Autumn18" in RunPeriod) or ("Fall18" in RunPeriod)):
    JECprefix = "Autumn18_V8"
    GT = '102X_upgrade2018_realistic_v18'
  else:
    JECprefix = "Autumn18_V8"
    GT = '102X_upgrade2018_realistic_v18'

  #jecAK8chsUncFile = "JEC/%s_MC_Uncertainty_AK8PFchs.txt"%(JECprefix)
  jecAK4chsUncFile = "JEC/%s_MC_Uncertainty_AK4PFchs.txt"%(JECprefix)
 



else : #Data

   JEC_runDependent_suffix= ""
   if ("2017" in RunPeriod):
     if ("Run2017B" in  RunPeriod): JEC_runDependent_suffix= "B"
     elif ("Run2017C" in  RunPeriod): JEC_runDependent_suffix= "C"
     elif ("Run2017D" in  RunPeriod): JEC_runDependent_suffix= "D"
     elif ("Run2017E" in  RunPeriod): JEC_runDependent_suffix= "E"
     elif ("Run2017F" in  RunPeriod): JEC_runDependent_suffix= "F"
     
     JECprefix = "Fall17_17Nov2017"+JEC_runDependent_suffix+"_V32"
     GT = '94X_dataRun2_v11'
     
     
   elif ("2016" in RunPeriod):
     if ("Run2016D" in  RunPeriod or "Run2016B" in  RunPeriod  or "Run2016C" in  RunPeriod  ): JEC_runDependent_suffix= "ABC"
     elif ("Run2016E" in  RunPeriod): JEC_runDependent_suffix= "EF"
     elif ("Run2016G" in  RunPeriod): JEC_runDependent_suffix= "GH"
     elif ("Run2016F" in  RunPeriod): JEC_runDependent_suffix= "GH"
#     elif ("Run2016F" in  RunPeriod and   options.runUpToEarlyF): JEC_runDependent_suffix= "EF"


     JECprefix = "Summer16_07Aug2017"+JEC_runDependent_suffix+"_V11"
     GT ='94X_dataRun2_v10'

   elif ("2018" in RunPeriod):
     if ("Run2018A" in  RunPeriod ): 
       JEC_runDependent_suffix= "A"
       GT="102X_dataRun2_Sep2018ABC_v2" 
     elif ("Run2018B" in  RunPeriod): 
       JEC_runDependent_suffix= "B"
       GT="102X_dataRun2_Sep2018ABC_v2"
     elif ("Run2018C" in  RunPeriod): 
       JEC_runDependent_suffix= "C"
       GT="102X_dataRun2_Sep2018ABC_v2"
     elif ("Run2018D" in  RunPeriod): 
       JEC_runDependent_suffix= "D"
       GT = '102X_dataRun2_Prompt_v13' 

     JECprefix = "Autumn18_Run"+JEC_runDependent_suffix+"_V8"
    
   #jecAK8chsUncFile = "JEC/%s_DATA_Uncertainty_AK8PFchs.txt"%(JECprefix)
   jecAK4chsUncFile = "JEC/%s_DATA_Uncertainty_AK4PFchs.txt"%(JECprefix)
 
   print "jec JEC_runDependent_suffix %s ,  prefix %s " %(JEC_runDependent_suffix,JECprefix)

#print "jec prefix ", JECprefix

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
    useDNN	      = cms.bool(config["USEDNN"]),
    doGenParticles    = cms.bool(config["DOGENPARTICLES"]),
    doGenEvent	      = cms.bool(config["DOGENEVENT"]),
    doPileUp	      = cms.bool(config["DOPILEUP"]),
    doJpsiMu	      = cms.bool(config["DOJPSIMU"]),
    doJpsiTau	      = cms.bool(config["DOJPSITAU"]),
    doBsTauTau	      = cms.bool(config["DOBSTAUTAU"]),
    doVertices	      = cms.bool(config["DOVERTICES"]),
    doMissingEt       = cms.bool(config["DOMISSINGET"]),
    doGenHist         = cms.bool(config["DOGENHIST"]),
    dzcut             = cms.double(config['DZCUT']),
    fsigcut           = cms.double(config['FSIGCUT']),
    vprobcut          = cms.double(config['VPROBCUT']),
    dnncut            = cms.double(config['DNNCUT']),
    dnnfile           = cms.string(config['DNNFILE']),                        
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
    gentaus = cms.InputTag("tauGenJets"),
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
    losttrack = cms.InputTag('lostTracks')
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

if config["RUNONMC"]:
  process.load("PhysicsTools.JetMCAlgos.TauGenJets_cfi")
  process.tauGenJets.GenParticles = cms.InputTag('prunedGenParticles')
  process.p += process.tauGenJets


process.p += process.ntuplizer
process.p.associate(pattask)

print pattask

#  LocalWords:  tauIdMVAIsoDBoldDMwLT
