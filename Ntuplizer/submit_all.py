#!/usr/bin/env python
"""
This is a small script that submits a config over many datasets
"""
import os
from optparse import OptionParser

import datetime

#d_today = datetime.datetime.now()
#
##d_today = datetime.date.today() 
#
#d_today = str(d_today).split('.')[0]
#d_today = d_today.replace(' ', '-').replace(':','')
#
#print d_today

username = os.environ['USER']

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
    
    parser.add_option('-l', '--luminosity', action="store_true", default=False, dest='luminosity')

#    parser.add_option("-l", "--luminosity", dest="luminosity",
#                      help=("Splitting job by lumi sections or by files"),
#                      metavar="LUMI")
    parser.add_option('-D', '--isData', action="store_true", default=False, dest='isData')
#    parser.add_option("-D", "--isData", dest="isData",
#                      help=("If is data saving the run period in the name"),
#                      metavar="isData")
    parser.add_option("-s", "--string", dest="string_to_add",
                      help=("Splitting job by lumi sections or by files"),
                      metavar="STRING")
    parser.add_option("-n", "--numOfFiles", dest="numOfFiles", default=None,
                      help=("Number of files per job"))

    parser.add_option('-g', '--isGlobal', action="store_true", default=False, dest='isGlobal', metavar="STRING")


    (options, args) = parser.parse_args()

#    print 'options.luminosity = ', options.luminosity

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
    config.General.transferLogs = True
    
    config.section_("JobType")
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = options.config

    if options.isData:
        config.JobType.maxMemoryMB = 5000
        config.JobType.numCores = 2
    else:
        config.JobType.maxMemoryMB = 4000
#        config.JobType.numCores = 2

    config.JobType.allowUndistributedCMSSW = True
    config.JobType.sendExternalFolder = True
    # config.JobType.pyCfgParams = ['DataProcessing=MC25ns_MiniAODv2','lheLabel=externalLHEProducer']
    #config.JobType.pyCfgParams = ['RunPeriod']
    config.JobType.inputFiles = [
#        'RecoTauTag_MVAs_2018Mar15.db',
        './JSON/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt',
        './data/DNN/tau_10_small.root'
        ]

    config.section_("Data")
    config.Data.inputDataset = None
    config.Data.allowNonValidInputDataset = True #To allow to run on non valid dataset

    if not options.isGlobal:
        config.Data.inputDBS = 'phys03' #to be commented in case of global#
    if options.luminosity == True :
#        config.Data.splitting = 'Automatic'
#        config.Data.splitting = 'LumiBased'
        config.Data.splitting = 'FileBased'
        config.Data.unitsPerJob = 2
    else:
        config.Data.splitting = 'FileBased'
        config.Data.unitsPerJob = 10

    if options.numOfFiles!=None:
        config.Data.unitsPerJob = int(options.numOfFiles)
        
#    config.Data.ignoreLocality = True
    config.Data.publication = False

    #config.Data.outLFNDirBase = '/store/user/cgalloni/Ntuple_2017_94v2_preliminary'
#    config.Data.outLFNDirBase = '/store/user/ytakahas/RJpsi_20191002_BJpsiX_020519'



    config.section_("Site")
    config.Site.storageSite = 'T2_CH_CSCS'
#    config.Site.storageSite = 'T3_CH_PSI'
#    config.Site.blacklist=['T2_US_Florida', 'T2_UA_KIPT']
#    config.Site.blacklist=['T2_US_Florida', 'T2_UA_KIPT']
#    config.Site.whitelist=['T2_US_Nebraska','T2_FR_GRIF_IRFU', 'T2_IT_Legnaro', 'T2_US_Purdue','T2_CH_CSCS', 'T2_CH_CERN', 'T2_US_Florida', 'T2_ES_CIEMAT', 'T2_US_Wisconsin']
    print 'Using config ' + options.config
    print 'Writing to directory ' + options.dir


    config.section_('Debug')
    config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']



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
        if s.find('#')!=-1: continue
        jobs.append( s )
        print '  --> added ' + s


    for ijob, job in enumerate(jobs) :

#        if job.find('#')!=-1: continue

        ptbin = job.split('/')[1]
        cond = job.split('/')[2]

#        print 'ptbin =', ptbin
#        print '2nd = ', (("_"+cond) if options.isData else "")
#        print '3rd = ', options.string_to_add
        config.General.requestName =  ptbin + (("_"+cond) if options.isData else "")  + "_" + options.string_to_add
        config.Data.inputDataset = job
#        config.Data.outputDatasetTag = ptbin  + (("_"+cond) if options.isData else "") + "_" + options.string_to_add
        config.Data.outputDatasetTag = cond #.split('-')[0] + "_" + options.string_to_add
        config.Data.outLFNDirBase = '/store/user/' + username + '/' + ptbin.split('_')[0] + '_' + options.string_to_add #+ '/' + cond + '/' + options.string_to_add
        config.General.workArea = options.dir + '_' + ptbin + '_' + cond
        print "ptbin :%s and cond: %s " %(ptbin, cond)
        print 'Submitting ' + config.General.requestName + ', dataset = ' + job
        print 'Configuration :'
        print config

        try :
            from multiprocessing import Process
            p = Process(target=submit, args=(config,))
            p.start()
            p.join()
#            submit(config)
        except :
            print 'Not submitted.'






if __name__ == '__main__':
    main()
