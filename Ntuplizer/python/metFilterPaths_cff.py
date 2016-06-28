import FWCore.ParameterSet.Config as cms

## We don't use "import *" because the cff contains some modules for which the C++ class doesn't exist
## and this triggers an error under unscheduled mode
from RecoMET.METFilters.metFilters_cff import HBHENoiseFilter, CSCTightHaloFilter, hcalLaserEventFilter, EcalDeadCellTriggerPrimitiveFilter, eeBadScFilter, ecalLaserCorrFilter
from RecoMET.METFilters.metFilters_cff import goodVertices, trackingFailureFilter, trkPOGFilters, manystripclus53X, toomanystripclus53X, logErrorTooManyClusters
from RecoMET.METFilters.metFilters_cff import metFilters

# individual filters
Flag_HBHENoiseFilter = cms.Path(HBHENoiseFilter)
Flag_CSCTightHaloFilter = cms.Path(CSCTightHaloFilter)
Flag_hcalLaserEventFilter = cms.Path(hcalLaserEventFilter)
Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(EcalDeadCellTriggerPrimitiveFilter)
Flag_goodVertices = cms.Path(goodVertices)
Flag_trackingFailureFilter = cms.Path(goodVertices + trackingFailureFilter)
Flag_eeBadScFilter = cms.Path(eeBadScFilter)
Flag_ecalLaserCorrFilter = cms.Path(ecalLaserCorrFilter)
Flag_trkPOGFilters = cms.Path(trkPOGFilters)

# NEW ICHEP!! *Not added to All_METFilters
    
Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(CSCTightHaloTrkMuUnvetoFilter)
Flag_CSCTightHalo2015Filter = cms.Path(CSCTightHalo2015Filter)
Flag_globalTightHalo2016Filter = cms.Path(globalTightHalo2016Filter)
Flag_globalSuperTightHalo2016Filter = cms.Path(globalSuperTightHalo2016Filter)
Flag_HcalStripHaloFilter = cms.Path(HcalStripHaloFilter)
Flag_chargedHadronTrackResolutionFilter = cms.Path(chargedHadronTrackResolutionFilter)
Flag_muonBadTrackFilter = cms.Path(muonBadTrackFilter)


# and the sub-filters
Flag_trkPOG_manystripclus53X = cms.Path(~manystripclus53X)
Flag_trkPOG_toomanystripclus53X = cms.Path(~toomanystripclus53X)
Flag_trkPOG_logErrorTooManyClusters = cms.Path(~logErrorTooManyClusters)


# and the summary
Flag_METFilters = cms.Path(metFilters)

#jen addition
def miniAOD_insertMETFiltersInSchedule(process):
    All_METFilters = [ 'Flag_HBHENoiseFilter', 'Flag_CSCTightHaloFilter', 'Flag_hcalLaserEventFilter', 'Flag_EcalDeadCellTriggerPrimitiveFilter', 'Flag_goodVertices', 'Flag_trackingFailureFilter', 'Flag_eeBadScFilter', 'Flag_ecalLaserCorrFilter', 'Flag_trkPOGFilters', 'Flag_trkPOG_manystripclus53X', 'Flag_trkPOG_toomanystripclus53X', 'Flag_trkPOG_logErrorTooManyClusters', 'Flag_METFilters'] 
    if hasattr(process,'schedule'):
        found = -1
        for i,p in enumerate(process.schedule):
            if isinstance(p, cms.EndPath):
                found = i
                break
        for pname in All_METFilters:
            path = getattr(process,pname)
            if not isinstance(path, cms.Path): raise RuntimeError, "Path %s is not a cms.Path?" % pathname
            if found == -1: process.schedule.append(path)
            else:           process.schedule.insert(found, path)
	    
def miniAOD_customizeMETFiltersFastSim(process):
    """Replace some MET filters that don't work in FastSim with trivial bools"""
    for X in 'CSCTightHaloFilter', 'HBHENoiseFilter':
        process.globalReplace(X, cms.EDFilter("HLTBool", result=cms.bool(True)))
    for X in 'manystripclus53X', 'toomanystripclus53X', 'logErrorTooManyClusters':
        process.globalReplace(X, cms.EDFilter("HLTBool", result=cms.bool(False)))
    return process
