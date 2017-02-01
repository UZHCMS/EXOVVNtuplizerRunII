#!/usr/bin/env python
"""
This is a small script that submits a config over many datasets
"""
import os
from optparse import OptionParser

def getOptions() :
    """
    Parse and return the arguments provided by the user.
    """
    usage = ('usage: python submit_all.py -c CONFIG -d DIR ')

    parser = OptionParser(usage=usage)
    parser.add_option("-c", "--config", dest="config",
        help=("The crab script you want to submit "),
        metavar="CONFIG")
    parser.add_option("-d", "--dir", dest="dir",
        help=("The crab directory you want to use "),
        metavar="DIR")
    parser.add_option("-f", "--datasets", dest="datasets",
        help=("File listing datasets to run over"),
        metavar="FILE")

    parser.add_option("-s", "--string", dest="string_to_add",
        help=("Splitting job by lumi sections or by files"),
        metavar="STRING")
    (options, args) = parser.parse_args()
    if options.config == None or options.dir == None:
        parser.error(usage)

    return options


def main():

    options = getOptions()

    #from WMCore.Configuration import Configuration
    from CRABClient.UserUtilities import config
    config = config()

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.section_("General")
    config.General.workArea = options.dir
    config.General.transferLogs = True

    config.section_("JobType")
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = options.config
    config.JobType.allowUndistributedCMSSW = True
    config.JobType.sendExternalFolder = True
   # config.JobType.pyCfgParams = ['DataProcessing=MC25ns_MiniAODv2','lheLabel=externalLHEProducer']
    config.JobType.inputFiles = [
        './JEC/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4Calo.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4JPT.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PF.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PFPuppi.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PFchs.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L1FastJet_AK8PF.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L1FastJet_AK8PFPuppi.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L1FastJet_AK8PFchs.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L1JPTOffset_AK4JPT.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L1RC_AK4PF.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L1RC_AK4PFchs.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L1RC_AK8PF.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L1RC_AK8PFchs.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2L3Residual_AK4Calo.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2L3Residual_AK4JPT.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2L3Residual_AK4PF.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2L3Residual_AK4PFPuppi.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2L3Residual_AK4PFchs.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2L3Residual_AK8PF.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2L3Residual_AK8PFPuppi.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2L3Residual_AK8PFchs.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2Relative_AK4Calo.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2Relative_AK4JPT.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2Relative_AK4PF.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2Relative_AK4PFPuppi.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2Relative_AK4PFchs.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2Relative_AK8PF.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2Relative_AK8PFPuppi.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2Relative_AK8PFchs.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2Residual_AK4Calo.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2Residual_AK4JPT.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2Residual_AK4PF.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2Residual_AK4PFPuppi.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2Residual_AK4PFchs.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2Residual_AK8PF.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2Residual_AK8PFPuppi.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L2Residual_AK8PFchs.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L3Absolute_AK4Calo.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L3Absolute_AK4JPT.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L3Absolute_AK4PF.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L3Absolute_AK4PFPuppi.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L3Absolute_AK4PFchs.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L3Absolute_AK8PF.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L3Absolute_AK8PFPuppi.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_L3Absolute_AK8PFchs.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_Uncertainty_AK4Calo.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_Uncertainty_AK4JPT.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_Uncertainty_AK4PF.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_Uncertainty_AK4PFPuppi.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_Uncertainty_AK4PFchs.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_Uncertainty_AK8PF.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_Uncertainty_AK8PFPuppi.txt',
        './JEC/Summer16_23Sep2016GV3_DATA_Uncertainty_AK8PFchs.txt',
        './JER/Spring16_25nsV6_MC_PtResolution_AK8PFchs.txt',
        './JER/Spring16_25nsV6_MC_PtResolution_AK4PFchs.txt',
        './JER/Spring16_25nsV6_MC_PtResolution_AK8PFPuppi.txt',
        './JER/Spring16_25nsV6_MC_PtResolution_AK4PFPuppi.txt',
        './JER/Spring16_25nsV6_MC_SF_AK8PFchs.txt',
        './JER/Spring16_25nsV6_MC_SF_AK4PFchs.txt',
        './JER/Spring16_25nsV6_MC_SF_AK8PFPuppi.txt',
        './JER/Spring16_25nsV6_MC_SF_AK4PFPuppi.txt',
        './JSON/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt',
]


    config.section_("Data")
    config.Data.inputDataset = None
    # config.Data.inputDBS = 'phys03' #to be commented in case of global#
    config.Data.splitting = 'LumiBased'#'LumiBased'#
    config.Data.unitsPerJob = 20
    config.Data.runRange= "278802-280385"
    config.Data.ignoreLocality = True
    config.Data.publication = False
    config.Data.lumiMask      = '/mnt/t3nfs01/data01/shome/cgalloni/RunII/CMSSW_8_0_21/src/EXOVVNtuplizerRunII/Ntuplizer/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
    config.section_("Site")
    config.Site.storageSite = 'T3_CH_PSI'
    config.Site.blacklist=['T2_US_Nebraska','T2_US_Wisconsin','T2_FR_IPHC','T2_EE_Estonia','T2_DE_RWTH']
    #config.Site.whitelist=['T2_US_Nebraska','T2_US_Wisconsin','T2_FR_IPHC','T2_EE_Estonia',
#    config.Data.outLFNDirBase = '/store/user/cgalloni/Ntuple_80_190916'
    config.Data.outLFNDirBase = '/store/user/cgalloni/Ntuple_80_20170131'


    print 'Using config ' + options.config
    print 'Writing to directory ' + options.dir


    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print 'Cannot execute command'
            print hte.headers

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    datasetsFile = open( options.datasets )
    jobsLines = datasetsFile.readlines()
    jobs = []
    for ijob in jobsLines :
        s = ijob.rstrip()
        jobs.append( s )
        print '  --> added ' + s


    for ijob, job in enumerate(jobs) :

        ptbin = job.split('/')[1]
        cond = job.split('/')[2]

        config.General.requestName = ptbin +'_'+ cond + options.string_to_add
        config.Data.inputDataset = job
        config.Data.outputDatasetTag =  ptbin +'_'+ cond + options.string_to_add
        print "ptbin :%s and cond: %s " %(ptbin, cond)
        print 'Submitting ' + config.General.requestName + ', dataset = ' + job
        print 'Configuration :'
        print config

        try :
            from multiprocessing import Process
            p = Process(target=submit, args=(config,))
            p.start()
            p.join()
           # submit(config)
        except :
            print 'Not submitted.'






if __name__ == '__main__':
    main()
