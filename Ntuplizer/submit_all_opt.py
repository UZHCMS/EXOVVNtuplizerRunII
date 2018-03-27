
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
    parser.add_option("-l", "--luminosity", dest="luminosity",
        help=("Splitting job by lumi sections or by files"),
        metavar="LUMI")
    parser.add_option("-D", "--isData", dest="isData",
        help=("If is data saving the run period in the name"),
        metavar="isData")
    parser.add_option("-s", "--string", dest="string_to_add",
        help=("Splitting job by lumi sections or by files"),
        metavar="STRING")


    (options, args) = parser.parse_args()


    if options.config == None or options.dir == None:
        parser.error(usage)
    if options.luminosity == None:
        options.luminosity = False
    else:
        options.luminosity = True
    if options.isData == None:
        options.isData = False
    else:
        options.isData = True
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
    #config.JobType.pyCfgParams = ['RunPeriod']
    config.JobType.inputFiles = [
        'RecoTauTag_MVAs_2018Mar15.db',
        './JSON/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
        ]

    config.section_("Data")
    config.Data.inputDataset = None
    config.Data.allowNonValidInputDataset = True #To allow to run on non valid dataset
    # config.Data.inputDBS = 'phys03' #to be commented in case of global#
    if options.luminosity == True :
        config.Data.splitting = 'LumiBased'
        config.Data.unitsPerJob = 25
    else:
        config.Data.splitting = 'FileBased'
        config.Data.unitsPerJob = 1
    #config.Data.ignoreLocality = True
    config.Data.publication = False
    #config.Data.outLFNDirBase = '/store/user/cgalloni/Ntuple_2017_94v2_preliminary'
    config.Data.outLFNDirBase = '/store/user/cgalloni/Ntuple_2017_94v2_preliminary_JEC_V6_tau_iso_JECUNC'
    #config.Data.outLFNDirBase = '/pnfs/psi.ch/cms/trivcat/store/t3groups/uniz-higgs/Fall17'

    config.section_("Site")
    #config.Site.storageSite = 'T2_CH_CSCS'
    config.Site.storageSite = 'T3_CH_PSI'
    #config.Site.blacklist=['T1_US_FNAL','T2_US_Wisconsin','T2_FR_IPHC','T2_EE_Estonia','T2_DE_RWTH','T2_KR_KNU','T2_KR_KISTI','T2_BR_SPRACE']
    #config.Site.whitelist=['T2_US_Nebraska','T2_US_Purdue','T2_CH_CSCS', 'T2_CH_CERN', 'T2_IT_Pisa','T2_US_MIT', 'T2_US_Florida', 'T2_US_UCSD', 'T2_IT_Bari','T2_IT_Legnaro']
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

        config.General.requestName =  ptbin + (("_"+cond)if options.isData else "")  + options.string_to_add
        config.Data.inputDataset = job
        config.JobType.pyCfgParams = ["RunPeriod="+job]  
        config.Data.outputDatasetTag = ptbin  + (("_"+cond)if options.isData else "" ) +  options.string_to_add
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
