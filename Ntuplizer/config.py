####### Process initialization ##########

import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuple")

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(options.inputFiles)
                            )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('flatTuple.root')
                                   )
				   
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.GlobalTag.globaltag = 'PLS170_V7AN1::All'

####### Logger ##########

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Ntuple')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(5)
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

####### Config parser ##########

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis')
options.maxEvents = 1
options.inputFiles ='file:testZZ.root'

options.parseArguments()

process.options  = cms.untracked.PSet( 
                     wantSummary = cms.untracked.bool(True),
                     SkipEvent = cms.untracked.vstring('ProductNotFound')
                     )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )


####### Redo Jet clustering sequence ##########

from RecoJets.Configuration.RecoPFJets_cff import ak8PFJetsCHS, ak8PFJetsCHSPruned, ak8PFJetsCHSSoftDrop, ak8PFJetsCHSPrunedLinks, ak8PFJetsCHSSoftDropLinks

process.chs = cms.EDFilter("CandPtrSelector",
  src = cms.InputTag('packedPFCandidates'),
  cut = cms.string('fromPV')
)

process.ak8PFJetsCHS = ak8PFJetsCHS.clone( src = cms.InputTag('chs') )
process.ak8PFJetsCHSPruned = ak8PFJetsCHSPruned.clone( src = cms.InputTag('chs') )
process.ak8PFJetsCHSSoftDrop = ak8PFJetsCHSPruned.clone( src = cms.InputTag('chs') )
process.ak8PFJetsCHSPrunedLinks = ak8PFJetsCHSPrunedLinks.clone()
process.ak8PFJetsCHSSoftDropLinks = ak8PFJetsCHSSoftDropLinks.clone()

process.NjettinessAK8 = cms.EDProducer("NjettinessAdder",
                               src = cms.InputTag("ak8PFJetsCHS"),
                               Njets = cms.vuint32(1, 2, 3, 4),
                               # variables for measure definition : 
                               measureDefinition = cms.uint32( 0 ), # CMS default is normalized measure
                               beta = cms.double(1.0),              # CMS default is 1
                               R0 = cms.double( 0.8 ),              # CMS default is jet cone size
                               Rcutoff = cms.double( -999.0),       # not used by default
                               # variables for axes definition :
                               axesDefinition = cms.uint32( 6 ),    # CMS default is 1-pass KT axes
                               nPass = cms.int32(-999),             # not used by default
                               akAxesR0 = cms.double(-999.0)        # not used by default
                               )

process.substructureSequence = cms.Sequence(process.chs + 
                                    	    process.ak8PFJetsCHS +
                                    	    process.ak8PFJetsCHSPruned +
                                    	    process.ak8PFJetsCHSSoftDrop +
                                    	    process.NjettinessAK8 +
                                    	    process.ak8PFJetsCHSPrunedLinks +
                                    	    process.ak8PFJetsCHSSoftDropLinks)

####### Redo pat jets sequence ##########

from ExoDiBosonResonances.EDBRJets.redoPatJets_cff import patJetCorrFactorsAK8, patJetsAK8, selectedPatJetsAK8

process.patJetCorrFactorsAK8 = patJetCorrFactorsAK8.clone( src = 'ak8PFJetsCHS' )
process.patJetsAK8 = patJetsAK8.clone( jetSource = 'ak8PFJetsCHS' )
process.selectedPatJetsAK8 = selectedPatJetsAK8.clone(cut = cms.string('pt > 20'))

process.redoPatJets = cms.Sequence( process.patJetCorrFactorsAK8 + process.patJetsAK8 + process.selectedPatJetsAK8 )

####### ExoDiBosonResonances objects ##########

import ExoDiBosonResonances.EDBRCommon.goodElectrons_cff
import ExoDiBosonResonances.EDBRCommon.goodMuons_cff
import ExoDiBosonResonances.EDBRCommon.goodJets_cff

process.load("ExoDiBosonResonances.EDBRCommon.goodElectrons_cff")
process.load("ExoDiBosonResonances.EDBRCommon.goodMuons_cff")
process.load("ExoDiBosonResonances.EDBRCommon.goodJets_cff")

process.leptonSequence = cms.Sequence(process.muSequence + process.eleSequence)
                                    
process.goodJets.src = 'selectedPatJetsAK8'
process.jetSequence = cms.Sequence(process.fatJetsSequence)

####### Ntuplizer initialization ##########

runOnMC = True

process.ntuplizer = cms.EDAnalyzer("Ntuplizer",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    taus = cms.InputTag("slimmedTaus"),
    photons = cms.InputTag("slimmedPhotons"),
    jets = cms.InputTag("slimmedJets"),
    fatjets = cms.InputTag("goodJets"),
    mets = cms.InputTag("slimmedMETs"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
)

####### Final path ##########

process.p = cms.Path(process.substructureSequence*process.redoPatJets*process.leptonSequence*process.jetSequence*process.ntuplizer)
