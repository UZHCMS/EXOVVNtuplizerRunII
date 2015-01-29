import FWCore.ParameterSet.Config as cms

process = cms.Process("Ntuple")

runOnMC = True

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Ntuple')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(5)
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis')
options.maxEvents = -1
#options.inputFiles = 'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/user/jngadiub/Wprime-M1000_WH_lvqq/patTuple/patTuple_Wprime_WH_lvqq_1.root'
#options.inputFiles = 'dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/pattuple/mc/v1/background/TT_CT10_TuneZ2star_8TeV-powheg-tauola/EDBR_PATtuple_edbr_TTBAR_03012013/d697403925b5703dc02456d20f4b4874/TT_CT10_TuneZ2star_8TeV-powheg-tauola__Summer12_DR53X-PU_S10_START53_V7A-v2__AODSIM_3503_1_LRP.root'
options.inputFiles ='root://xrootd.unl.edu//store/mc/Phys14DR/RSGravitonToWW_kMpl01_M_1000_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/10000/CE92595C-9676-E411-A785-00266CF2E2C8.root'

options.parseArguments()

process.options  = cms.untracked.PSet( 
                     wantSummary = cms.untracked.bool(True),
                     SkipEvent = cms.untracked.vstring('ProductNotFound')
                     )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(options.inputFiles)
                            )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('flatTuple.root')
                                   )

process.ntuplizer = cms.EDAnalyzer('Ntuplizer')

process.p = cms.Path(process.ntuplizer)
