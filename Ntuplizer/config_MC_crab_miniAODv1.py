####### Process initialization ##########
#### To use with MC MiniAOD_v1 ####
import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuple")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')

process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string('flatTuple.root')
                                   )

from EXOVVNtuplizerRunII.Ntuplizer.ntuplizerOptions_MC_crab_miniAODv1_cfi import config

				   
####### Config parser ##########

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('analysis')

options.maxEvents = 100

#data file

#options.inputFiles = '/store/mc/RunIISpring15MiniAODv2/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/003F1529-D36D-E511-9E33-001E6724816F.root'
#options.inputFiles = '/store/mc/RunIISpring16MiniAODv1/RadionToWWToWlepWhad_width0p3_M-3000_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_v3-v1/50000/20BB570E-2F16-E611-BF25-02163E015F8E.root'
options.inputFiles = '/store/mc/RunIISpring16MiniAODv2/WJetsToQQ_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/30000/0EA1D6CA-931A-E611-BFCD-BCEE7B2FE01D.root'
#options.inputFiles = 'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/mc/RunIIFall15MiniAODv2/WprimeToWhToWhadhbb_narrow_M-1000_13TeV-madgraph/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/8EE77064-42B8-E511-8BA1-003048895D40.root'
#options.inputFiles = 'xroot://t3dcachedb.psi.ch:1094//pnfs/psi.ch/cms/trivcat/store/user/cgalloni/RunII/RadionTohhTohtatahbb_narrow_M-1000_13TeV-madgraph/MiniAOD_TauBoosted_v0667_maxDepth100_jetPt100/miniAOD_4.root'
#options.inputFiles = 'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/cgalloni/RunII/RadionTohhTohtatahbb_narrow_M-1000_13TeV-madgraph/MiniAOD_TauBoosted_v0667_maxDepth100_jetPt100/miniAOD_4.root'
#options.inputFiles = '/store/user/cgalloni/MiniAOD_191115/RadionTohhTohtatahbb_narrow_M-1000_13TeV-madgraph/BoostedTaus_RadionTohhTohtatahbb_narrow_M-1000_13TeV-madgraph_v0/151119_130555/0000/miniAOD_1.root'
options.parseArguments()

process.options  = cms.untracked.PSet( 
                     wantSummary = cms.untracked.bool(False),
                     SkipEvent = cms.untracked.vstring('ProductNotFound'),
                     allowUnscheduled = cms.untracked.bool(True)
                     )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(options.inputFiles)
                            )                     

######## Sequence settings ##########

# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015#ETmiss_filters
# For the RunIISpring15DR74 MC campaing, the process name in PAT.
# For Run2015B PromptReco Data, the process name is RECO.
# For Run2015B re-MiniAOD Data 17Jul2015, the process name is PAT.
hltFiltersProcessName = 'RECO'
if config["RUNONMC"] or config["JSONFILE"].find('reMiniAOD') != -1:
  hltFiltersProcessName = 'PAT'
reclusterPuppi=config["DOAK8PUPPIRECLUSTERING"]
if reclusterPuppi:
  print "RECLUSTERING PUPPI with latest tune from CMSSW_8_0_20"
  
#! To recluster and add AK8 Higgs tagging and softdrop subjet b-tagging (both need to be simoultaneously true or false, if not you will have issues with your softdrop subjets!)
#If you use the softdrop subjets from the slimmedJetsAK8 collection, only CSV seems to be available?
doAK8softdropReclustering = False
if config["DOAK8RECLUSTERING"] == True: doAK8softdropReclustering = True

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
if config["RUNONMC"]: GT = '80X_mcRun2_asymptotic_2016_miniAODv2'
elif not(config["RUNONMC"]): GT = '80X_dataRun2_Prompt_v8'


print "*************************************** GLOBAL TAG *************************************************" 
print GT
print "****************************************************************************************************" 
process.GlobalTag = GlobalTag(process.GlobalTag, GT)
   
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
                                                                                                          
process.chs = cms.EDFilter("CandPtrSelector",
  src = cms.InputTag('packedPFCandidates'),
  cut = cms.string('fromPV')
)

process.ak4PFJetsCHS = ak4PFJetsCHS.clone( src = 'chs' )
process.ak4PFJetsCHS.doAreaFastjet = True
process.ak8CHSJets = ak8PFJetsCHS.clone( src = 'chs', jetPtMin = fatjet_ptmin )

process.NjettinessAK8 = cms.EDProducer("NjettinessAdder",
             src = cms.InputTag("ak8CHSJets"),
             Njets = cms.vuint32(1, 2, 3, 4),
             # variables for measure definition :
             measureDefinition = cms.uint32( 0 ), # CMS default is normalized measure
             beta = cms.double(1.0),        # CMS default is 1
             R0 = cms.double( 0.8 ),        # CMS default is jet cone size
             Rcutoff = cms.double( 999.0),      # not used by default
             # variables for axes definition :
             axesDefinition = cms.uint32( 6 ),    # CMS default is 1-pass KT axes
             nPass = cms.int32(999),       # not used by default
             akAxesR0 = cms.double(-999.0)      # not used by default
             )

if config["DOAK10TRIMMEDRECLUSTERING"]:			       
  process.ECFAK10 = cms.EDProducer("ECFAdder",
             src = cms.InputTag("ak10CHSJetsTrimmed"),
             Njets = cms.vuint32(1, 2, 3),
             beta = cms.double(1.0),        # CMS default is 1
             )
			       

process.ak8CHSJetsPruned = ak8PFJetsCHSPruned.clone( src = 'chs', jetPtMin = fatjet_ptmin )
process.ak8CHSJetsSoftDrop = ak8PFJetsCHSSoftDrop.clone( src = 'chs', jetPtMin = fatjet_ptmin, beta = betapar  )


if config["DOAK10TRIMMEDRECLUSTERING"]:			       
  process.ak10CHSJetsTrimmed = ak8PFJetsCHSTrimmed.clone( src = 'chs', jetPtMin = fatjet_ptmin, rParam = 1.0, rFilt = 0.2, trimPtFracMin = 0.05 )

if reclusterPuppi:
  process.load('CommonTools/PileupAlgos/Puppi_cff')
  process.puppi.useExistingWeights = False
  process.puppi.candName = cms.InputTag('packedPFCandidates')
  process.puppi.vertexName = cms.InputTag('offlineSlimmedPrimaryVertices')  
  process.ak8PuppiJets = ak8PFJetsCHS.clone( src = 'puppi', jetPtMin = fatjet_ptmin )
  process.ak8PuppiJetsPruned = ak8PFJetsCHSPruned.clone( src = 'puppi', jetPtMin = fatjet_ptmin )
  process.ak8PuppiJetsSoftDrop = ak8PFJetsCHSSoftDrop.clone( src = 'puppi', jetPtMin = fatjet_ptmin, beta = betapar  )
  process.NjettinessAK8Puppi = process.NjettinessAK8.clone( src = 'ak8PuppiJets' )

if config["GETJECFROMDBFILE"]:
  process.load("CondCore.DBCommon.CondDBCommon_cfi")
  process.jec = cms.ESSource("PoolDBESSource",
            DBParameters = cms.PSet(
                messageLevel = cms.untracked.int32(5)
                ),
            timetype = cms.string('runnumber'),
            toGet = cms.VPSet(
            cms.PSet(
                 record = cms.string('JetCorrectionsRecord'),
                 tag    = cms.string('JetCorrectorParametersCollection_Summer15_50nsV5_MC_AK4PFchs'),
                 label  = cms.untracked.string('AK4PFchs')
                 ),
            cms.PSet(
                 record = cms.string('JetCorrectionsRecord'),
                 tag    = cms.string('JetCorrectorParametersCollection_Summer15_50nsV5_MC_AK8PFchs'),
                 label  = cms.untracked.string('AK8PFchs')
                 ),
            cms.PSet(
                 record = cms.string('JetCorrectionsRecord'),
                 tag    = cms.string('JetCorrectorParametersCollection_Summer15_50nsV5_MC_AK8PFPuppi'),
                 label  = cms.untracked.string('AK8PFPuppi')
                 ),
            ),
            connect = cms.string('sqlite:Summer15_50nsV5_MC.db')
            )
  if not config["RUNONMC"]:
    process.jec.toGet[0].tag =  cms.string('JetCorrectorParametersCollection_Summer15_50nsV5_DATA_AK4PFchs')
    process.jec.toGet[1].tag =  cms.string('JetCorrectorParametersCollection_Summer15_50nsV5_DATA_AK8PFchs')
    process.jec.toGet[2].tag =  cms.string('JetCorrectorParametersCollection_Summer15_50nsV5_DATA_AK8PFPuppi')
    process.jec.connect = cms.string('sqlite:Summer15_50nsV5_DATA.db')
  process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')



####### Add AK8 GenJets ##########
if config["ADDAK8GENJETS"]:

  from RecoJets.Configuration.RecoGenJets_cff import ak8GenJets
  process.ak8GenJets = ak8GenJets.clone(src = 'packedGenParticles')

  process.NjettinessGenAK8 = cms.EDProducer("NjettinessAdder",
                              src=cms.InputTag("ak8GenJets"),
                              Njets=cms.vuint32(1,2,3,4),          # compute 1-, 2-, 3-, 4- subjettiness
                              # variables for measure definition : 
                              measureDefinition = cms.uint32( 0 ), # CMS default is normalized measure
                              beta = cms.double(1.0),              # CMS default is 1
                              R0 = cms.double( 0.8 ),              # CMS default is jet cone size
                              Rcutoff = cms.double( 999.0),       # not used by default
                              # variables for axes definition :
                              axesDefinition = cms.uint32( 6 ),    # CMS default is 1-pass KT axes
                              nPass = cms.int32(999),             # not used by default
                              akAxesR0 = cms.double(-999.0)        # not used by default
                              )

  process.genParticlesForJets = cms.EDProducer("InputGenJetsParticleSelector",
                                               src = cms.InputTag("packedGenParticles"),
                                               ignoreParticleIDs = cms.vuint32(
                                                                  1000022,
                                                                  1000012, 1000014, 1000016,
                                                                  2000012, 2000014, 2000016,
                                                                  1000039, 5100039,
                                                                  4000012, 4000014, 4000016,
                                                                  9900012, 9900014, 9900016,
                                                                  39),
                                              partonicFinalState = cms.bool(False),
                                              excludeResonances = cms.bool(False),
                                              excludeFromResonancePids = cms.vuint32(12, 13, 14, 16),
                                              tausAsJets = cms.bool(False)
                                              )

  from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters

  process.ak8GenJetsPruned = ak8GenJets.clone(
              SubJetParameters,
              usePruning = cms.bool(True),
              writeCompound = cms.bool(True),
              jetCollInstanceName=cms.string("SubJets")
              )
            
  process.ak8GenJetsSoftDrop = ak8GenJets.clone(
              SubJetParameters,
              useSoftDrop = cms.bool(True),
              R0 = cms.double(0.8),
              beta = betapar,
              writeCompound = cms.bool(True),
              jetCollInstanceName=cms.string("SubJets")
              )

  process.ak8GenJetsPrunedMass = cms.EDProducer("RecoJetDeltaRValueMapProducer",
                                            src = cms.InputTag("ak8GenJets"),
                                            matched = cms.InputTag("ak8GenJetsPruned"),
                                            distMax = cms.double(0.8),
                                            value = cms.string('mass')
                                            )

  process.ak8GenJetsSoftDropMass = cms.EDProducer("RecoJetDeltaRValueMapProducer",
                                            src = cms.InputTag("ak8GenJets"),
                                            matched = cms.InputTag("ak8GenJetsSoftDrop"),                                         
                                            distMax = cms.double(0.8),
                                            value = cms.string('mass') 
                                            )
                                                      

  # process.ak8GenJetsPrunedMass = ak8PFJetsCHSPrunedMass.clone(
  #             matched = cms.InputTag("ak8GenJetsPruned"),
  #             src = cms.InputTag("ak8GenJets")
  #             )
  #
  # process.ak8GenJetsSoftDropMass = ak8PFJetsCHSSoftDropMass.clone(
  #             matched = cms.InputTag("ak8GenJetsSoftDrop"),
  #             beta = betapar,
  #             src = cms.InputTag("ak8GenJets")
  #             )

  # process.substructureSequenceGen+=process.ak8GenJets
  # process.substructureSequenceGen+=process.NjettinessGenAK8
  #
  # process.substructureSequenceGen += process.ak8GenJetsSoftDrop + process.ak8GenJetsSoftDropMass
  # process.substructureSequenceGen += process.ak8GenJetsPruned + process.ak8GenJetsPrunedMass

  from EXOVVNtuplizerRunII.Ntuplizer.redoPatJets_cff import patJetCorrFactorsAK8, patJetsAK8, selectedPatJetsAK8

  # Redo pat jets from gen AK8

  process.genJetsAK8 = patJetsAK8.clone( jetSource = 'ak8GenJets' )
  process.genJetsAK8.userData.userFloats.src = [ cms.InputTag("ak8GenJetsPrunedMass"), cms.InputTag("ak8GenJetsSoftDropMass"), cms.InputTag("NjettinessGenAK8:tau1"), cms.InputTag("NjettinessGenAK8:tau2"), cms.InputTag("NjettinessGenAK8:tau3")]
  process.genJetsAK8.addJetCorrFactors = cms.bool(False)
  process.genJetsAK8.jetCorrFactorsSource = cms.VInputTag( cms.InputTag("") )
  process.selectedGenJetsAK8 = selectedPatJetsAK8.clone( src = 'genJetsAK8', cut = cms.string('pt > 20') )

################# Prepare recluster jets with b-tagging ######################
bTagDiscriminators = [
    # 'pfJetProbabilityBJetTags',
    # 'pfJetBProbabilityBJetTags',
    # 'pfSimpleSecondaryVertexHighEffBJetTags',
    # 'pfSimpleSecondaryVertexHighPurBJetTags',
    'pfCombinedInclusiveSecondaryVertexV2BJetTags',
    # 'pfTrackCountingHighPurBJetTags',
    # 'pfTrackCountingHighEffBJetTags',
    'pfBoostedDoubleSecondaryVertexAK8BJetTags'    
]

#Needed in 80X to get the latest Hbb training
if config["UpdateJetCollection"]:
  from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
## Update the slimmedJets in miniAOD: corrections from the chosen Global Tag are applied and the b-tag discriminators are re-evaluated
  updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJetsAK8'),
    jetCorrections = ('AK8PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    btagDiscriminators = bTagDiscriminators
  )
## Update to latest PU jet ID training
  process.load("RecoJets.JetProducers.PileupJetID_cfi")
  process.pileupJetIdUpdated = process.pileupJetId.clone(
    jets=cms.InputTag("slimmedJets"),
    inputIsCorrected=True,
    applyJec=True,
    vertexes=cms.InputTag("offlineSlimmedPrimaryVertices")
  )
  from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors, updatedPatJets
  process.patJetCorrFactorsReapplyJEC = updatedPatJetCorrFactors.clone(
    src = cms.InputTag("slimmedJets"),
    levels = ['L1FastJet', 'L2Relative', 'L3Absolute'] )
  process.updatedJets = updatedPatJets.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
  )
  process.updatedJets.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
  process.updatedJets.userData.userInts.src += ['pileupJetIdUpdated:fullId']

def cap(s): return s[0].upper() + s[1:]

from PhysicsTools.PatAlgos.tools.jetTools import *
#process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

bTagParameters = dict(
    #trackSource = cms.InputTag('unpackedTracksAndVertices'),
    pfCandidates = cms.InputTag('packedPFCandidates'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    elSource = cms.InputTag('slimmedElectrons'),
    muSource = cms.InputTag('slimmedMuons'),
    btagDiscriminators = bTagDiscriminators
) 

def recluster_addBtagging(process, fatjets_name, groomed_jets_name, jetcorr_label = 'AK8PFchs', jetcorr_label_subjets = 'AK4PFchs', genjets_name = None, verbose = False, btagging = True, subjets = True):
    rParam = getattr(process, fatjets_name).rParam.value()
    algo = None
    if 'ca' in fatjets_name.lower():
	algo = 'ca'
	assert getattr(process, fatjets_name).jetAlgorithm.value() == 'CambridgeAachen'
    elif 'ak' in fatjets_name.lower():
	algo = 'ak'
	assert getattr(process, fatjets_name).jetAlgorithm.value() == 'AntiKt'
    else:
	raise RuntimeError, "Unknown jet algorithm for fatjets name %s" % fatjets_name
  
    subjets_name = groomed_jets_name + 'Subjets' # e.g. AK8CHSPruned + Subjets
  
    # add genjet producers, if requested:
    groomed_genjets_name = 'INVALID'
    ungroomed_genjets_name = 'INVALID'
  
    if genjets_name is not None:
	    groomed_jetproducer = getattr(process, groomed_jets_name)
	    assert groomed_jetproducer.type_() in ('FastjetJetProducer', 'CATopJetProducer'), "do not know how to construct genjet collection for %s" % repr(groomed_jetproducer)
	    groomed_genjets_name = genjets_name(groomed_jets_name)
	    if verbose: print "Adding groomed genjets ", groomed_genjets_name
	    setattr(process, groomed_genjets_name, groomed_jetproducer.clone(src = cms.InputTag('packedGenParticles'), jetType = 'GenJet'))
	    # add for ungroomed jets if not done yet (maybe never used in case ungroomed are not added, but that's ok ..)
	    ungroomed_jetproducer = getattr(process, fatjets_name)
	    assert ungroomed_jetproducer.type_() == 'FastjetJetProducer'
	    ungroomed_genjets_name = genjets_name(fatjets_name)
	    if verbose: print "Adding ungroomed genjets ", ungroomed_genjets_name
	    setattr(process, ungroomed_genjets_name, ungroomed_jetproducer.clone(src = cms.InputTag('packedGenParticles'), jetType = 'GenJet'))
      

    # patify ungroomed jets, if not already done:
    if config["RUNONMC"]:
      jetcorr_levels = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute'])
    else:
      jetcorr_levels = cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'])
    add_ungroomed = not hasattr(process, 'patJets' + cap(fatjets_name))
    addJetCollection(process, labelName = fatjets_name, jetSource = cms.InputTag(fatjets_name), algo = algo, rParam = rParam,
	    jetCorrections = (jetcorr_label, jetcorr_levels, 'None'),
	    genJetCollection = cms.InputTag(ungroomed_genjets_name),
	    **bTagParameters
	)

    # patify groomed fat jets, with b-tagging:
    if config["RUNONMC"]:
      jetcorr_levels_groomed = cms.vstring(['L2Relative', 'L3Absolute']) # NO L1 corretion for groomed jets
    else:
      jetcorr_levels_groomed = cms.vstring(['L2Relative', 'L3Absolute', 'L2L3Residual'])
    addJetCollection(process, labelName = groomed_jets_name, jetSource = cms.InputTag(groomed_jets_name), algo = algo, rParam = rParam,
       jetCorrections = (jetcorr_label, jetcorr_levels_groomed, 'None'),
       **bTagParameters)
    # patify subjets, with subjet b-tagging:
    if subjets:
      addJetCollection(process, labelName = subjets_name, jetSource = cms.InputTag(groomed_jets_name, 'SubJets'), algo = algo, rParam = rParam,
	jetCorrections = (jetcorr_label_subjets, jetcorr_levels, 'None'),
	explicitJTA = True,
	svClustering = True,
	fatJets = cms.InputTag(fatjets_name), groomedFatJets = cms.InputTag(groomed_jets_name),
	genJetCollection = cms.InputTag(groomed_genjets_name, 'SubJets'),
	**bTagParameters)
  
      # add the merged jet collection which contains the links from fat jets to subjets:
      setattr(process, 'patJets' + cap(groomed_jets_name) + 'Packed',cms.EDProducer("BoostedJetMerger",
	jetSrc=cms.InputTag("patJets" + cap(groomed_jets_name)),
	subjetSrc=cms.InputTag("patJets" + cap(subjets_name)))
	)
  
    # adapt all for b-tagging, and switch off some PAT features not supported in miniAOD:
    if subjets:
      module_names = [subjets_name, groomed_jets_name]
    else:
      module_names = [groomed_jets_name]
    if add_ungroomed: module_names += [fatjets_name]
    for name in module_names:
	if hasattr(process,'pfInclusiveSecondaryVertexFinderTagInfos' + cap(name)):
	    getattr(process,'pfInclusiveSecondaryVertexFinderTagInfos' + cap(name)).extSVCollection = cms.InputTag('slimmedSecondaryVertices')
	getattr(process, 'patJetPartonMatch' + cap(name)).matched = 'prunedGenParticles'
	producer = getattr(process, 'patJets' + cap(name))
	producer.addJetCharge = False
	producer.addAssociatedTracks = False
	if not config["DOHBBTAG"]:
	    producer.addDiscriminators = True
	    producer.addBTagInfo = True
	producer.addGenJetMatch = genjets_name is not None
	# for fat groomed jets, gen jet match and jet flavor is not working, so switch it off:
	if name == groomed_jets_name:
	    producer.addGenJetMatch = False
	    producer.getJetMCFlavour = False

################# Recluster jets with b-tagging ######################
if config["DOAK8RECLUSTERING"]: 
    recluster_addBtagging(process, 'ak8CHSJets', 'ak8CHSJetsSoftDrop', genjets_name = lambda s: s.replace('CHS', 'Gen'))
    recluster_addBtagging(process, 'ak8CHSJets', 'ak8CHSJetsPruned', genjets_name = lambda s: s.replace('CHS', 'Gen'))
    process.ak8PFJetsCHSPrunedMass = cms.EDProducer("RecoJetDeltaRValueMapProducer",
                                            src = cms.InputTag("ak8CHSJets"),
                                            matched = cms.InputTag("ak8CHSJetsPruned"),
                                            distMax = cms.double(0.8),
                                            value = cms.string('mass')
                                            )

    process.ak8PFJetsCHSSoftDropMass = cms.EDProducer("RecoJetDeltaRValueMapProducer",
                                            src = cms.InputTag("ak8CHSJets"),
                                            matched = cms.InputTag("ak8CHSJetsSoftDrop"),                                         
                                            distMax = cms.double(0.8),
                                            value = cms.string('mass') 
                                            )         
    process.ak8PFJetsCHSPrunedMassCorrected = cms.EDProducer("RecoJetDeltaRValueMapProducer",
                                            src = cms.InputTag("ak8CHSJets"),
                                            matched = cms.InputTag("patJetsAk8CHSJetsPrunedPacked"),
                                            distMax = cms.double(0.8),
                                            value = cms.string('mass')
                                            )

    process.ak8PFJetsCHSSoftDropMassCorrected = cms.EDProducer("RecoJetDeltaRValueMapProducer",
                                            src = cms.InputTag("ak8CHSJets"),
                                            matched = cms.InputTag("patJetsAk8CHSJetsSoftDropPacked"),                                         
                                            distMax = cms.double(0.8),
                                            value = cms.string('mass') 
                                            )    
					         
    process.patJetsAk8CHSJets.userData.userFloats.src += ['ak8PFJetsCHSPrunedMass','ak8PFJetsCHSSoftDropMass','ak8PFJetsCHSPrunedMassCorrected','ak8PFJetsCHSSoftDropMassCorrected']
    process.patJetsAk8CHSJets.userData.userFloats.src += ['NjettinessAK8:tau1','NjettinessAK8:tau2','NjettinessAK8:tau3']
    process.patJetsAk8CHSJets.addTagInfos = True
    #process.patJetsAk8CHSJetsSoftDropSubjets.addBTagInfo = True

################# Recluster trimmed jets ######################
if config["DOAK10TRIMMEDRECLUSTERING"]:	
    recluster_addBtagging(process, 'ak8CHSJets', 'ak10CHSJetsTrimmed', genjets_name = lambda s: s.replace('CHS', 'Gen'), verbose = False, btagging = False, subjets = False)
    process.patJetsAk10CHSJetsTrimmed.userData.userFloats.src += ['ECFAK10:ecf1','ECFAK10:ecf2','ECFAK10:ecf3']
    
################# Recluster puppi jets ######################
if reclusterPuppi:
    recluster_addBtagging(process, 'ak8PuppiJets', 'ak8PuppiJetsSoftDrop', jetcorr_label = 'AK8PFPuppi', genjets_name = lambda s: s.replace('Puppi', 'Gen'))
    recluster_addBtagging(process, 'ak8PuppiJets', 'ak8PuppiJetsPruned', jetcorr_label = 'AK8PFPuppi', genjets_name = lambda s: s.replace('Puppi', 'Gen'))
  
    process.ak8PFJetsPuppiPrunedMass = cms.EDProducer("RecoJetDeltaRValueMapProducer",
    					  src = cms.InputTag("ak8PuppiJets"),
    					  matched = cms.InputTag("ak8PuppiJetsPruned"),
    					  distMax = cms.double(0.8),
    					  value = cms.string('mass')
    					  )

    process.ak8PFJetsPuppiSoftDropMass = cms.EDProducer("RecoJetDeltaRValueMapProducer",
    					  src = cms.InputTag("ak8PuppiJets"),
    					  matched = cms.InputTag("ak8PuppiJetsSoftDrop"),					  
    					  distMax = cms.double(0.8),
    					  value = cms.string('mass') 
    					  )	    
    process.ak8PFJetsPuppiPrunedMassCorrected = cms.EDProducer("RecoJetDeltaRValueMapProducer",
    					  src = cms.InputTag("ak8PuppiJets"),
    					  matched = cms.InputTag("patJetsAk8PuppiJetsPrunedPacked"),
    					  distMax = cms.double(0.8),
    					  value = cms.string('mass')
    					  )

    process.ak8PFJetsPuppiSoftDropMassCorrected = cms.EDProducer("RecoJetDeltaRValueMapProducer",
    					  src = cms.InputTag("ak8PuppiJets"),
    					  matched = cms.InputTag("patJetsAk8PuppiJetsSoftDropPacked"),					 
    					  distMax = cms.double(0.8),
    					  value = cms.string('mass') 
    					  )	    

    process.patJetsAk8PuppiJets.userData.userFloats.src += ['ak8PFJetsPuppiSoftDropMass','ak8PFJetsPuppiSoftDropMassCorrected']
    #process.patJetsAk8PuppiJets.userData.userFloats.src += ['ak8PFJetsPuppiPrunedMass','ak8PFJetsPuppiPrunedMassCorrected']
    process.patJetsAk8PuppiJets.userData.userFloats.src += ['NjettinessAK8Puppi:tau1','NjettinessAK8Puppi:tau2','NjettinessAK8Puppi:tau3']
    process.patJetsAk8PuppiJets.addTagInfos = True

    process.packedJetsAk8PuppiJets = cms.EDProducer("JetSubstructurePacker",
            jetSrc = cms.InputTag("patJetsAk8PuppiJets"),
            distMax = cms.double(0.8),
            algoTags = cms.VInputTag(
                cms.InputTag("patJetsAk8PuppiJetsSoftDropPacked")
            ),
            algoLabels = cms.vstring(
                'SoftDropPuppi'
                ),
            fixDaughters = cms.bool(False),
            packedPFCandidates = cms.InputTag("packedPFCandidates"),
    )

# ###### Recluster MET ##########
if config["DOMETRECLUSTERING"]:

  from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection
		  
  if config["RUNONMC"]:
     switchJetCollection(process,
                         jetSource = cms.InputTag('ak4PFJetsCHS'),
		         jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], ''),
		         genParticles = cms.InputTag('prunedGenParticles'),
		         pvSource = cms.InputTag('offlineSlimmedPrimaryVertices')
     )
  else:
     switchJetCollection(process,
                         jetSource = cms.InputTag('ak4PFJetsCHS'),
		         jetCorrections = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute','L2L3Residual'], ''),
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
switchOnVIDElectronIdProducer(process,dataFormat)

process.egmGsfElectronIDSequence = cms.Sequence(process.egmGsfElectronIDs)

# define which IDs we want to produce
#my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff',
#		 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff']
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
                 'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

####### Event filters ###########

##___________________________HCAL_Noise_Filter________________________________||
if config["DOHLTFILTERS"]:
 process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
 process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
 process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 

####### Ntuplizer initialization ##########
jetsAK4 = "slimmedJets"
jetsAK8 = "slimmedJetsAK8"
jetsAK8pruned = ""
# jetsAK8softdrop = "slimmedJetsAK8PFCHSSoftDropPacked" (if you want to add this subjet collection, changes need to be made in plugins/JetsNtuplizer.cc! Not needed to obtain subjets)
jetsAK8softdrop = ""
jetsAK10trimmed = ""
jetsAK8Puppi = ""  

METS = "slimmedMETs"
if config["DOMETRECLUSTERING"]: jetsAK4 = "selectedPatJets"
if config["USENOHF"]: METS = "slimmedMETsNoHF"  

##___________________ MET significance and covariance matrix ______________________##

if config["DOMETSVFIT"]:
  print "Using event pfMET covariance for SVfit"
  process.load("RecoMET.METProducers.METSignificance_cfi")
  process.load("RecoMET.METProducers.METSignificanceParams_cfi")
  process.METSequence = cms.Sequence (process.METSignificance)

if config["DOMVAMET"]:
  from RecoMET.METPUSubtraction.jet_recorrections import recorrectJets
  recorrectJets(process, isData=False)

  from RecoMET.METPUSubtraction.MVAMETConfiguration_cff import runMVAMET
  runMVAMET( process, jetCollectionPF="patJetsReapplyJEC")
  process.MVAMET.srcLeptons  = cms.VInputTag("slimmedMuons", "slimmedElectrons", "slimmedTaus")
  process.MVAMET.requireOS = cms.bool(False)

##___________________ taus ______________________##
TAUS = ""
BOOSTEDTAUS = ""
genAK8 = ""

if config["ADDAK8GENJETS"]:
  genAK8 = 'selectedGenJetsAK8'
    
if config["DOAK8RECLUSTERING"]:
  jetsAK8 = "patJetsAk8CHSJets"
if config["UpdateJetCollection"]:
  jetsAK4 = "updatedJets"
  jetsAK8 = "updatedPatJetsTransientCorrected"
if doAK8softdropReclustering:  
  jetsAK8softdrop = "patJetsAk8CHSJetsSoftDropPacked"  
if config["DOAK8PRUNEDRECLUSTERING"]:  
  jetsAK8pruned = "patJetsAk8CHSJetsPrunedPacked"
if config["DOAK10TRIMMEDRECLUSTERING"]:  
  jetsAK10trimmed = "patJetsAk10CHSJetsTrimmed"
if reclusterPuppi:  
  jetsAK8Puppi = "packedJetsAk8PuppiJets"  

if config["DOTAUSBOOSTED"]:
  TAUS = "slimmedTaus"
  BOOSTEDTAUS = "slimmedTausBoosted"     
else:
  TAUS = "slimmedTaus"
  BOOSTEDTAUS = "slimmedTaus"
 
######## JEC ########
jecLevelsAK8chs = []
jecLevelsAK8Groomedchs = []
jecLevelsAK4chs = []
jecLevelsAK4 = []
jecLevelsAK8Puppi = []
jecLevelsForMET = []

if config["BUNCHSPACING"] == 25 and config["RUNONMC"] and config["SPRING16"]:
  JECprefix = "Spring16_25nsV3"
elif config["BUNCHSPACING"] == 25 and not(config["RUNONMC"]) and config["SPRING16"]:
  JECprefix = "Spring16_25nsV3"

jecAK8chsUncFile = "%s_MC_Uncertainty_AK8PFchs.txt"%(JECprefix)
jecAK4chsUncFile = "%s_MC_Uncertainty_AK4PFchs.txt"%(JECprefix)


if config["CORRJETSONTHEFLY"]:
   if config["RUNONMC"]:
     jecLevelsAK8chs = [
     	 '%s_MC_L1FastJet_AK8PFchs.txt'%(JECprefix),
     	 '%s_MC_L2Relative_AK8PFchs.txt'%(JECprefix),
     	 '%s_MC_L3Absolute_AK8PFchs.txt'%(JECprefix)
       ]
     jecLevelsAK8Groomedchs = [
     	 '%s_MC_L2Relative_AK8PFchs.txt'%(JECprefix),
     	 '%s_MC_L3Absolute_AK8PFchs.txt'%(JECprefix)
       ]
     jecLevelsAK8Puppi = [
     	 '%s_MC_L2Relative_AK8PFPuppi.txt'%(JECprefix),
     	 '%s_MC_L3Absolute_AK8PFPuppi.txt'%(JECprefix)
       ]
     jecLevelsAK4chs = [
     	 '%s_MC_L1FastJet_AK4PFchs.txt'%(JECprefix),
     	 '%s_MC_L2Relative_AK4PFchs.txt'%(JECprefix),
     	 '%s_MC_L3Absolute_AK4PFchs.txt'%(JECprefix)
       ]
   else:
     jecLevelsAK8chs = [
     	 '%s_DATA_L1FastJet_AK8PFchs.txt'%(JECprefix),
     	 '%s_DATA_L2Relative_AK8PFchs.txt'%(JECprefix),
     	 '%s_DATA_L3Absolute_AK8PFchs.txt'%(JECprefix),
	     '%s_DATA_L2L3Residual_AK8PFchs.txt'%(JECprefix)
       ]
     jecLevelsAK8Groomedchs = [
     	 '%s_DATA_L2Relative_AK8PFchs.txt'%(JECprefix),
     	 '%s_DATA_L3Absolute_AK8PFchs.txt'%(JECprefix),
	     '%s_DATA_L2L3Residual_AK8PFchs.txt'%(JECprefix)
       ]
     jecLevelsAK8Puppi = [
     	 '%s_DATA_L2Relative_AK8PFPuppi.txt'%(JECprefix),
     	 '%s_DATA_L3Absolute_AK8PFPuppi.txt'%(JECprefix),
	     '%s_DATA_L2L3Residual_AK8PFPuppi.txt'%(JECprefix)
       ]
     jecLevelsAK4chs = [
     	 '%s_DATA_L1FastJet_AK4PFchs.txt'%(JECprefix),
     	 '%s_DATA_L2Relative_AK4PFchs.txt'%(JECprefix),
     	 '%s_DATA_L3Absolute_AK4PFchs.txt'%(JECprefix),
	     '%s_DATA_L2L3Residual_AK4PFchs.txt'%(JECprefix)
       ]   
if config["CORRMETONTHEFLY"]:  
   if config["RUNONMC"]:
     jecLevelsForMET = [				       
     	 '%s_MC_L1FastJet_AK4PFchs.txt'%(JECprefix),
     	 '%s_MC_L2Relative_AK4PFchs.txt'%(JECprefix),
     	 '%s_MC_L3Absolute_AK4PFchs.txt'%(JECprefix)
       ]
   else:       					       
     jecLevelsForMET = [
     	 '%s_DATA_L1FastJet_AK4PFchs.txt'%(JECprefix),
     	 '%s_DATA_L2Relative_AK4PFchs.txt'%(JECprefix),
     	 '%s_DATA_L3Absolute_AK4PFchs.txt'%(JECprefix),
       '%s_DATA_L2L3Residual_AK4PFchs.txt'%(JECprefix)
       ]	
      			    
#from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
#process.goodSlimmedJets = cms.EDFilter("PFJetIDSelectionFunctorFilter",
#                        filterParams = pfJetIDSelector.clone(),
#                        src = cms.InputTag("slimmedJets")
#                        )
#process.goodFatJets = cms.EDFilter("PFJetIDSelectionFunctorFilter",
#                        filterParams = pfJetIDSelector.clone(),
#                        src = cms.InputTag(jetsAK8)
#                        )
                                                                                      
################## Ntuplizer ###################
process.ntuplizer = cms.EDAnalyzer("Ntuplizer",
    runOnMC	      = cms.bool(config["RUNONMC"]),
    doGenParticles    = cms.bool(config["DOGENPARTICLES"]),
    doGenJets	      = cms.bool(config["DOGENJETS"]),
    doGenEvent	      = cms.bool(config["DOGENEVENT"]),
    doPileUp	      = cms.bool(config["DOPILEUP"]),
    doElectrons       = cms.bool(config["DOELECTRONS"]),
    doMuons	      = cms.bool(config["DOMUONS"]),
    doTaus	      = cms.bool(config["DOTAUS"]),
    doAK8Jets	      = cms.bool(config["DOAK8JETS"]),
    doAK4Jets	      = cms.bool(config["DOAK4JETS"]),
    doVertices	      = cms.bool(config["DOVERTICES"]),
    doTriggerDecisions= cms.bool(config["DOTRIGGERDECISIONS"]),
    doTriggerObjects  = cms.bool(config["DOTRIGGEROBJECTS"]),
    doHltFilters      = cms.bool(config["DOHLTFILTERS"]),
    doMissingEt       = cms.bool(config["DOMISSINGET"]),
    doHbbTag	      = cms.bool(config["DOHBBTAG"]),
    doPrunedSubjets   = cms.bool(config["DOAK8PRUNEDRECLUSTERING"]),
    doTrimming        = cms.bool(config["DOAK10TRIMMEDRECLUSTERING"]),
    doPuppi           = cms.bool(config["DOAK8PUPPI"]),
    doBoostedTaus     = cms.bool(config["DOTAUSBOOSTED"]),
    doMETSVFIT        = cms.bool(config["DOMETSVFIT"]),
    doMVAMET        = cms.bool(config["DOMVAMET"]),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    eleHEEPId51Map = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV51"),
    eleHEEPIdMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
    eleVetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
    eleLooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
    taus = cms.InputTag(TAUS),
    tausBoostedTau = cms.InputTag(BOOSTEDTAUS),
    jets = cms.InputTag(jetsAK4),
    fatjets = cms.InputTag(jetsAK8),
    prunedjets = cms.InputTag(jetsAK8pruned),
    softdropjets = cms.InputTag(jetsAK8softdrop),
    trimmedjets = cms.InputTag(jetsAK10trimmed),
    puppijets = cms.InputTag(jetsAK8Puppi),
    genJets = cms.InputTag("slimmedGenJets"),
    genJetsAK8 = cms.InputTag(genAK8),
    subjetflavour = cms.InputTag("AK8byValAlgo"),
    mets = cms.InputTag(METS),
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
    HLT = cms.InputTag("TriggerResults","","HLT2"),
    triggerobjects = cms.InputTag("selectedPatTrigger"),
    triggerprescales = cms.InputTag("patTrigger"),
    noiseFilter = cms.InputTag('TriggerResults','', hltFiltersProcessName),
    jecAK8chsPayloadNames = cms.vstring( jecLevelsAK8chs ),
    jecAK8chsUnc = cms.string( jecAK8chsUncFile ),
    jecAK8GroomedchsPayloadNames = cms.vstring( jecLevelsAK8Groomedchs ),
    jecAK8PuppiPayloadNames = cms.vstring( jecLevelsAK8Puppi ),
    jecAK4chsPayloadNames = cms.vstring( jecLevelsAK4chs ),
    jecAK4chsUnc = cms.string( jecAK4chsUncFile ),
    jecpath = cms.string(''),
    
    ## Noise Filters ###################################
    # defined here: https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/PhysicsTools/PatAlgos/python/slimming/metFilterPaths_cff.py
    noiseFilterSelection_HBHENoiseFilter = cms.string('Flag_HBHENoiseFilter'),
    noiseFilterSelection_HBHENoiseFilterLoose = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResultRun2Loose"),
    noiseFilterSelection_HBHENoiseFilterTight = cms.InputTag("HBHENoiseFilterResultProducer", "HBHENoiseFilterResultRun2Tight"),
    noiseFilterSelection_HBHENoiseIsoFilter = cms.InputTag("HBHENoiseFilterResultProducer", "HBHEIsoNoiseFilterResult"),
    noiseFilterSelection_CSCTightHaloFilter = cms.string('Flag_CSCTightHaloFilter'),
    noiseFilterSelection_CSCTightHalo2015Filter = cms.string('Flag_CSCTightHalo2015Filter'),
    noiseFilterSelection_hcalLaserEventFilter = cms.string('Flag_hcalLaserEventFilter'),
    noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter = cms.string('Flag_EcalDeadCellTriggerPrimitiveFilter'),
    noiseFilterSelection_goodVertices = cms.string('Flag_goodVertices'),
    noiseFilterSelection_trackingFailureFilter = cms.string('Flag_trackingFailureFilter'),
    noiseFilterSelection_eeBadScFilter = cms.string('Flag_eeBadScFilter'),
    noiseFilterSelection_ecalLaserCorrFilter = cms.string('Flag_ecalLaserCorrFilter'),
    noiseFilterSelection_trkPOGFilters = cms.string('Flag_trkPOGFilters'),
    # and the sub-filters
    noiseFilterSelection_trkPOG_manystripclus53X = cms.string('Flag_trkPOG_manystripclus53X'),
    noiseFilterSelection_trkPOG_toomanystripclus53X = cms.string('Flag_trkPOG_toomanystripclus53X'),
    noiseFilterSelection_trkPOG_logErrorTooManyClusters = cms.string('Flag_trkPOG_logErrorTooManyClusters'),
    # summary
    noiseFilterSelection_metFilters = cms.string('Flag_METFilters'),

)

####### Final path ##########
process.p = cms.Path()
if config["DOHLTFILTERS"]:
 process.p += process.HBHENoiseFilterResultProducer
process.p += process.ntuplizer
