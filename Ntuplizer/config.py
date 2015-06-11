####### Process initialization ##########

import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuple")

# process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
# process.load('Configuration.StandardSequences.Geometry_cff')
# process.load('Configuration.StandardSequences.MagneticField_38T_cff')
# #process.GlobalTag.globaltag = 'PHYS14_25_V2'
# process.GlobalTag.globaltag = 'MCRUN2_74_V7'

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')


process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('flatTuple.root')
                                   )
				   
####### Config parser ##########

import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing ('analysis')

options.maxEvents = -1
options.inputFiles = 'file:/afs/cern.ch/work/c/cgalloni/RunII/CMSSW_7_4_1/src/ZZ_MuTauEleTau_2000.root '#B2G-RunIISpring15DR74-00001_MuTauEleTau.root'
options.parseArguments()

process.options  = cms.untracked.PSet( 
                     wantSummary = cms.untracked.bool(True),
                     SkipEvent = cms.untracked.vstring('ProductNotFound'),
                     allowUnscheduled = cms.untracked.bool(True)
                     )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(options.inputFiles)
                            )                     

######## Sequence settings ##########


#! To recluster and add AK8 Higgs tagging and softdrop subjet b-tagging (both need to be simoultaneously true or false, if not you will have issues with your softdrop subjets!)
#If you use the softdrop subjets from the slimmedJetsAK8 collection, only CSV seems to be available?
doAK8reclustering = True
doAK8softdropReclustering = False
doBtagging = False

#! To add pruned jet and pruned subjet collection (not in MINIAOD)
doAK8prunedReclustering = False

#! To recluster MET with new corrections
doMETReclustering = False
doSemileptonicTausBoosted = False
####### Logger ##########

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Ntuple')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(5)
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


####### Redo Jet clustering sequence ##########

fatjet_ptmin = 100.0

from RecoJets.Configuration.RecoPFJets_cff import *
                                                                                                          
process.chs = cms.EDFilter("CandPtrSelector",
  src = cms.InputTag('packedPFCandidates'),
  cut = cms.string('fromPV')
)

process.ak4PFJetsCHS = ak4PFJetsCHS.clone( src = 'chs' )
process.ak8CHSJets = ak8PFJetsCHS.clone( src = 'chs', jetPtMin = fatjet_ptmin )


process.NjettinessAK8 = cms.EDProducer("NjettinessAdder",
             src = cms.InputTag("ak8CHSJets"),
             Njets = cms.vuint32(1, 2, 3, 4),
             # variables for measure definition :
             measureDefinition = cms.uint32( 0 ), # CMS default is normalized measure
             beta = cms.double(1.0),        # CMS default is 1
             R0 = cms.double( 0.8 ),        # CMS default is jet cone size
             Rcutoff = cms.double( -999.0),      # not used by default
             # variables for axes definition :
             axesDefinition = cms.uint32( 6 ),    # CMS default is 1-pass KT axes
             nPass = cms.int32(-999),       # not used by default
             akAxesR0 = cms.double(-999.0)      # not used by default
             )
			       

process.ak8CHSJetsPruned = ak8PFJetsCHSPruned.clone( src = 'chs', jetPtMin = fatjet_ptmin  )
process.ak8CHSJetsSoftDrop = ak8PFJetsCHSSoftDrop.clone( src = 'chs', jetPtMin = fatjet_ptmin  )


################# Recluster jets with b-tagging ######################

bTagDiscriminators = [
    'pfJetProbabilityBJetTags',
    'pfJetBProbabilityBJetTags',
    'pfSimpleSecondaryVertexHighEffBJetTags',
    'pfSimpleSecondaryVertexHighPurBJetTags',
    'pfCombinedInclusiveSecondaryVertexV2BJetTags',
    'pfTrackCountingHighPurBJetTags',
    'pfTrackCountingHighEffBJetTags',
    'pfBoostedDoubleSecondaryVertexAK8BJetTags'
    
]

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

def recluster_addBtagging(process, fatjets_name, groomed_jets_name, jetcorr_label = 'AK7PFchs', jetcorr_label_subjets = 'AK4PFchs', genjets_name = None, verbose = False, btagging = True):
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
    add_ungroomed = not hasattr(process, 'patJets' + cap(fatjets_name))
    addJetCollection(process, labelName = fatjets_name, jetSource = cms.InputTag(fatjets_name), algo = algo, rParam = rParam,
            jetCorrections = (jetcorr_label, cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
            genJetCollection = cms.InputTag(ungroomed_genjets_name),
            **bTagParameters
        )

    # patify groomed fat jets, with b-tagging:
    addJetCollection(process, labelName = groomed_jets_name, jetSource = cms.InputTag(groomed_jets_name), algo = algo, rParam = rParam,
       jetCorrections = (jetcorr_label, cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
       **bTagParameters)
    # patify subjets, with subjet b-tagging:
    addJetCollection(process, labelName = subjets_name, jetSource = cms.InputTag(groomed_jets_name, 'SubJets'), algo = algo, rParam = rParam,
        jetCorrections = (jetcorr_label_subjets, cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
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
    module_names = [subjets_name, groomed_jets_name]
    if add_ungroomed: module_names += [fatjets_name]
    for name in module_names:
        if hasattr(process,'pfInclusiveSecondaryVertexFinderTagInfos' + cap(name)):
            getattr(process,'pfInclusiveSecondaryVertexFinderTagInfos' + cap(name)).extSVCollection = cms.InputTag('slimmedSecondaryVertices')
        getattr(process, 'patJetPartonMatch' + cap(name)).matched = 'prunedGenParticles'
        producer = getattr(process, 'patJets' + cap(name))
        producer.addJetCharge = False
        producer.addAssociatedTracks = False
        if not doBtagging:
            producer.addDiscriminators = False
            producer.addBTagInfo = False
        producer.addGenJetMatch = genjets_name is not None
        # for fat groomed jets, gen jet match and jet flavor is not working, so switch it off:
        if name == groomed_jets_name:
            producer.addGenJetMatch = False
            producer.getJetMCFlavour = False
    
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

process.patJetsAk8CHSJets.userData.userFloats.src += ['ak8PFJetsCHSPrunedMass','ak8PFJetsCHSSoftDropMass']

process.patJetsAk8CHSJets.userData.userFloats.src += ['NjettinessAK8:tau1','NjettinessAK8:tau2','NjettinessAK8:tau3']
process.patJetsAk8CHSJets.addTagInfos = True

# ###### Recluster MET ##########

process.redoPatMET = cms.Sequence()

# process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")

from RecoMET.METProducers.PFMET_cfi import pfMet
from JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff import corrPfMetType1

from JetMETCorrections.Configuration.JetCorrectors_cff import ak4PFCHSL1FastL2L3Corrector,ak4PFCHSL1FastjetCorrector,ak4PFCHSL2RelativeCorrector,ak4PFCHSL3AbsoluteCorrector
# process.ak4PFCHSL1FastL2L3Corrector = ak4PFCHSL1FastL2L3Corrector.clone()

from JetMETCorrections.Type1MET.correctedMet_cff import pfMetT1
from PhysicsTools.PatAlgos.producersLayer1.metProducer_cfi import patMETs

process.pfMet = pfMet.clone(src = "packedPFCandidates")
process.pfMet.calculateSignificance = False

corrPfMetType1.jetCorrLabel = cms.InputTag('ak4PFCHSL1FastL2L3Corrector')
corrPfMetType1.src = cms.InputTag('ak4PFJetsCHS')
corrPfMetType1.offsetCorrLabel = cms.InputTag("ak4PFCHSL1FastjetCorrector")

process.pfMetT1 = pfMetT1.clone()
process.patMETs = patMETs.clone()
process.patMETs.addGenMET = False # There's no point in recalculating this, and we can't remake it since we don't have genParticles beyond |eta|=5

if doMETReclustering:
    process.redoPatMET+=process.ak4PFJetsCHS
    process.redoPatMET+=process.pfMet
    process.redoPatMET+=process.ak4PFCHSL1FastjetCorrector
    process.redoPatMET+=process.ak4PFCHSL2RelativeCorrector
    process.redoPatMET+=process.ak4PFCHSL3AbsoluteCorrector
    process.redoPatMET+=process.ak4PFCHSL1FastL2L3Corrector
    process.redoPatMET+=process.corrPfMetType1
    process.redoPatMET+=process.pfMetT1
    process.redoPatMET+=process.patMETs

####### Ntuplizer initialization ##########

runOnMC = True

jetsAK8 = "slimmedJetsAK8"
jetsAK8pruned = ""
# jetsAK8softdrop = "slimmedJetsAK8PFCHSSoftDropPacked" (if you want to add this subjet collection, changes need to be made in plugins/JetsNtuplizer.cc! Not needed to obtain subjets)
jetsAK8softdrop = ""

METS = "slimmedMETs"
reclusteredMETS = ""
pfMETS = ""

if doAK8reclustering:
  jetsAK8 = "patJetsAk8CHSJets"
if doAK8softdropReclustering:  
  jetsAK8softdrop = "patJetsAk8CHSJetsSoftDropPacked"  
if doAK8prunedReclustering:  
  jetsAK8pruned = "patJetsAk8CHSJetsPrunedPacked"
   
if doMETReclustering:
    METS = "" 
    reclusteredMETS = "patMETs"
    pfMETS = "pfMet" 
     

######## JEC ########
jecLevelsAK8chs = [
    'JEC/MCRUN2_74_V7::All_L1FastJet_AK8PFchs.txt', #JEC for 74X
    'JEC/MCRUN2_74_V7::All_L2Relative_AK8PFchs.txt',
    'JEC/MCRUN2_74_V7::All_L3Absolute_AK8PFchs.txt'
  ]
jecLevelsAK4chs = [
    'JEC/MCRUN2_74_V7::All_L1FastJet_AK4PFchs.txt',
    'JEC/MCRUN2_74_V7::All_L2Relative_AK4PFchs.txt',
    'JEC/MCRUN2_74_V7::All_L3Absolute_AK4PFchs.txt'
  ]
   

process.ntuplizer = cms.EDAnalyzer("Ntuplizer",
    runOnMC = cms.bool(runOnMC),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    tausMuTau = cms.InputTag("slimmedTausMuTau"),
    tausEleTau = cms.InputTag("slimmedTausEleTau"),
    jets = cms.InputTag("slimmedJets"),
    fatjets = cms.InputTag(jetsAK8),
    prunedjets = cms.InputTag(jetsAK8pruned),
    softdropjets = cms.InputTag(jetsAK8softdrop),
    genJets = cms.InputTag("slimmedGenJets"),
    subjetflavour = cms.InputTag("AK8byValAlgo"),
    mets = cms.InputTag(METS),
    reclusteredmets = cms.InputTag(reclusteredMETS),
    pfmets = cms.InputTag(pfMETS),
    corrMetPx = cms.string("+0.1166 + 0.0200*Nvtx"),
    corrMetPy = cms.string("+0.2764 - 0.1280*Nvtx"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    genparticles = cms.InputTag("prunedGenParticles"),
    PUInfo = cms.InputTag("addPileupInfo"),
    genEventInfo = cms.InputTag("generator"),
    HLT = cms.InputTag("TriggerResults","","HLT"),
    triggerobjects = cms.InputTag("selectedPatTrigger"),
    triggerprescales = cms.InputTag("patTrigger"),
    jecAK8chsPayloadNames = cms.vstring( jecLevelsAK8chs ),
    jecAK4chsPayloadNames = cms.vstring( jecLevelsAK4chs ),
    jecpath = cms.string(''),
)

####### Final path ##########
process.p = cms.Path(process.ntuplizer) 
