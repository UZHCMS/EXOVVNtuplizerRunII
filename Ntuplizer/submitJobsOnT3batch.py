import os, commands
import random
import sys
from optparse import OptionParser
import ConfigParser
import time

#-----------------------------------------------------------------------------------------
def getLocalJobsDir(localdir):

   return commands.getoutput("pwd")+"/"+localdir

#-----------------------------------------------------------------------------------------
def getJobsDirs(outdir,jobname):
   user = commands.getoutput("whoami")
   
   path = "/pnfs/psi.ch/cms/trivcat/store/user/%s"%(user+"/"+outdir)
   cmd = "uberftp t3se01.psi.ch 'ls %s'" %path
   ls_la = commands.getoutput(cmd)
   
   list_ = []
   list_.extend(ls_la.split(os.linesep))
      
   jobsdir = []
   for a in list_:
      b = a.split(" ")
      if b[-1:][0].find(jobname) != -1:
         c = b[-1:][0].split()
         jobsdir.append(path+"/"+c[0])

   return jobsdir	 
	 
#-----------------------------------------------------------------------------------------
def getFileListT3(src):

   cmd = "uberftp t3se01.psi.ch 'ls %s'" %src
   ls_la = commands.getoutput(cmd)
   
   list_ = []
   list_.extend(ls_la.split(os.linesep))   
   files = []
   for a in list_:
      b = a.split(" ")
      if b[-1:][0].find("root") != -1:
         c = b[-1:][0].split(".")
         filename = c[0] + ".root"
         files.append( filename )

   if len(files) == 0 :
      print ("No files found in the directory %s" %(src))
      sys.exit()
   else:
      print ("Found %i files in the directory %s" %(len(files),src))  
   
   return files   

#-----------------------------------------------------------------------------------------
def getFileListFromSampleT3(sample):

   path = ""
   if sample == "WJetsToLNu_HT-100to200":
      path = "/pnfs/psi.ch/cms/trivcat/store/user/jngadiub/Wprime-M1000_WH_lvqq/patTuple/"
   else:
      print "Sample %s unknown! Exiting..." %sample
      sys.exit()
   
   return getFileListT3(path) 

#-----------------------------------------------------------------------------------------
def getFileListFromSampleDAS(sample):

   dataset = ""
   if sample == "WJetsToLNu_HT-100to200":
      dataset="/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM"
   elif sample == "WJetsToLNu_HT-200to400":
      dataset = "/WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM"
   elif sample == "WJetsToLNu_HT-400to600":
      dataset = "/WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM"  
   elif sample == "WJetsToLNu_HT-600toInf":
      dataset = "/WJetsToLNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM" 
   elif sample == "TBarToLeptons_s-channel":
      dataset = "/TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM" 
   elif sample == "TBarToLeptons_t-channel":
      dataset = "/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM"
   elif sample == "TTJets":
      dataset = "/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM"  
   elif sample == "TToLeptons_s-channel":
      dataset = "/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM"
   elif sample == "TToLeptons_t-channel":
      dataset = "/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM"
   elif sample == "T_tW-channel":
      dataset = "/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM" 
   elif sample == "Tbar_tW-channel":
      dataset = "/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM"
   elif sample == "QCD_HT_500To1000_ext":
      dataset = "/QCD_HT-500To1000_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1/MINIAODSIM"
   elif sample == "QCD_HT_250To500_ext":
      dataset = "/QCD_HT_250To500_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v2/MINIAODSIM"
   elif sample == "QCD_HT_1000ToInf_ext":
      dataset ="/QCD_HT_1000ToInf_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1_ext1-v1/MINIAODSIM"
   elif sample == "QCD_HT-500To1000":
      dataset = "/QCD_HT-500To1000_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM"       
   elif sample == "QCD_HT_1000ToInf":
      dataset = "/QCD_HT_1000ToInf_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM"      
   elif sample == "QCD_HT_250To500":
      dataset = "/QCD_HT_250To500_13TeV-madgraph/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM"      
   elif sample == "QCD_Pt-1000to1400":
      dataset = "/QCD_Pt-1000to1400_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_trkalmb_castor_PHYS14_25_V1-v1/MINIAODSIM"      
   elif sample == "QCD_Pt-1400to1800":
      dataset = "/QCD_Pt-1400to1800_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_trkalmb_castor_PHYS14_25_V1-v1/MINIAODSIM"      
   elif sample == "QCD_Pt-15to3000":
      dataset = "/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/Phys14DR-PU20bx25_trkalmb_PHYS14_25_V1-v1/MINIAODSIM"      
   elif sample == "QCD_Pt-1800to2400":
      dataset = "/QCD_Pt-1800to2400_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_trkalmb_PHYS14_25_V1-v2/MINIAODSIM"      
   elif sample == "QCD_Pt-2400to3200":
      dataset = "/QCD_Pt-2400to3200_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_trkalmb_PHYS14_25_V1-v1/MINIAODSIM"      
   elif sample == "QCD_Pt-300to470":
      dataset = "/QCD_Pt-300to470_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_trkalmb_castor_PHYS14_25_V1-v2/MINIAODSIM"      
   elif sample == "QCD_Pt-3200":
      dataset = "/QCD_Pt-3200_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_trkalmb_PHYS14_25_V1-v1/MINIAODSIM"      
   elif sample == "QCD_Pt-470to600":
      dataset = "/QCD_Pt-470to600_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_trkalmb_castor_PHYS14_25_V1-v2/MINIAODSIM"      
   elif sample == "QCD_Pt_600to800":
      dataset = "/QCD_Pt-600to800_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_trkalmb_castor_PHYS14_25_V1-v1/MINIAODSIM"      
   elif sample == "QCD_Pt_800to1000":
      dataset = "/QCD_Pt-800to1000_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_trkalmb_castor_PHYS14_25_V1-v2/MINIAODSIM"      
   elif sample == "RSGravitonToWW_M_1000":
      dataset = "/RSGravitonToWW_kMpl01_M_1000_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v2/MINIAODSIM"      
   elif sample == "RSGravitonToWW_M_2000":
      dataset = "/RSGravitonToWW_kMpl01_M_2000_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM"      
   elif sample == "RSGravitonToWW_M_3000":
      dataset = "/RSGravitonToWW_kMpl01_M_3000_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM"
   elif sample == "RSGravitonToWW_M_4000":
      dataset = "/RSGravitonToWW_kMpl01_M_4000_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM"
   else:
      print "Sample %s unknown! Exiting..." %sample
      sys.exit()
   
   return getFileListDAS(dataset) 


    
#-----------------------------------------------------------------------------------------
def getFileListDAS(dataset):

   cmd = './das_client.py --query="file dataset=%s" --limit=1000' %(dataset)
   cmd_out = commands.getoutput( cmd )
   tmpList = cmd_out.split(os.linesep)
   files = []
   for l in tmpList:
      if l.find(".root") != -1:
         files.append(l)
            
   return files 
   
#-----------------------------------------------------------------------------------------   
argv = sys.argv
parser = OptionParser()

parser.add_option("-C", "--config", dest="config", default=[], action="append",
                              help="file defining the job configuration")
parser.add_option("-t", "--test", dest="test", default=False, action="store_true",
                              help="test job config")			      
parser.add_option("-c", "--copy", dest="copyfiles", default=False, action="store_true",
                              help="copy job output")
parser.add_option("--clean", "--clean", dest="clean", default=False, action="store_true",
                              help="clean job dir")
parser.add_option("--useDAS", "--useDAS", dest="useDAS", default=False, action="store_true",
                              help="use das query")
parser.add_option("-H", "-H", dest="help", default=False, action="store_true",
                              help="usage")
			      			      
(opts, args) = parser.parse_args(argv)

if opts.config =="":
    opts.config = "config"

if opts.help:
   print ""
   print "Usage for testing: "
   print "          python submitJobsOnT3batch.py -C submitJobsOnT3batch.cfg --test --useDAS"
   print "Or: "
   print "          python submitJobsOnT3batch.py -C submitJobsOnT3batch.cfg --useDAS"
   print ""
   sys.exit()
   
print opts.config
config = ConfigParser.ConfigParser()
config.read(opts.config)
nfiles = config.getint('JobsConfig','nfiles')
src = config.get('JobsConfig','src')
sample = config.get('JobsConfig','sample')
outdir = config.get('JobsConfig','outdir')
localjobdir = config.get('JobsConfig','localjobdir')
jobname = config.get('JobsConfig','jobname')
outfile = config.get('JobsConfig','outfile')
queue = config.get('JobsConfig','queue')
cmsswdir = config.get('JobsConfig','cmsswdir')
cfg = config.get('JobsConfig','cfg')
bashfile = config.get('JobsConfig','bashfile')
maxevents = config.getint('JobsConfig','maxevents')
prefix = config.get('JobsConfig','prefix')
newdir = config.get('JobsConfig','newdir')
   
#-----------------------------------------------------------------------------------------
if opts.copyfiles:

   jobsdir = getJobsDirs(outdir,jobname)
   user = commands.getoutput("whoami")	
    
   for j in jobsdir:
      a = j.split("-")
      jobid = a[1]
      inputpath = "srm://t3se01.psi.ch:8443/srm/managerv2?SFN="+j+"/"+outfile
      outputpath = "srm://t3se01.psi.ch:8443/srm/managerv2?SFN=%s"%(newdir+"/"+prefix+"_"+jobid+".root")
      checkfile = "ls -l %s"%(newdir+"/"+prefix+"_"+jobid+".root") 
      status,cmd_out = commands.getstatusoutput(checkfile)
      if not(status): 
         cmd = "gfal-rm %s"%(outputpath)
	 print cmd
	 os.system(cmd)
      cmd = "gfal-copy %s %s" %(inputpath,outputpath)
      print cmd
      os.system(cmd)
          	 
   sys.exit()	 

if opts.clean:

   jobsdir = getJobsDirs(outdir,jobname)
   
   for j in jobsdir:
      status,cmd_out = commands.getstatusoutput( 'ls '+j )
      if not(status):
         a = j.split("-")
         jobid = a[1]
         inputpath = "srm://t3se01.psi.ch:8443/srm/managerv2?SFN="+j+"/"+outfile
	 cmd = "gfal-rm "+inputpath
	 print cmd
	 os.system(cmd)
      jdir = "srm://t3se01.psi.ch:8443/srm/managerv2?SFN="+j
      cmd = "lcg-del -d %s" %jdir
      print cmd
      os.system(cmd)
    
   print "" 
   cmd = "rm -rf "+getLocalJobsDir(localjobdir) 
   print cmd 
   os.system(cmd)   
   user = commands.getoutput("whoami")
   dir = '/pnfs/psi.ch/cms/trivcat/store/user/'+user+"/"+outdir
   cmd = "lcg-del -d srm://t3se01.psi.ch:8443/srm/managerv2?SFN=" + dir
   print cmd
   os.system(cmd)
         
   sys.exit()   
	 
#-----------------------------------------------------------------------------------------	 
files = []

if not(opts.useDAS):
   if sample != "": files = getFileListFromSampleT3(sample)
   else: files = getFileListT3(src)   
else:
   if sample != "": files = getFileListFromSampleDAS(sample)
   else: files = getFileListDAS(src)   

print len(files)

status,cmd_out = commands.getstatusoutput( 'ls '+getLocalJobsDir(localjobdir) )
if status: os.system('mkdir '+getLocalJobsDir(localjobdir))
   
it = 0
jobIndex = 1
cmds = []
while it < len(files):
  newlist = files[it:it+nfiles]
  inputFiles = ""
  print ""
  print "Submitting job %i : %i files" %(jobIndex,len(newlist))
  for f in newlist:
    if f != "":
       print "    * %s" %(f)
       if not(opts.useDAS): inputFiles+="dcap://t3se01.psi.ch:22125/%s," %(src+"/"+f) 
       else: inputFiles+="root://xrootd.unl.edu/%s," %(f)	  
  tmpList = list(inputFiles)
  tmpList.pop(len(tmpList)-1)
  inputFiles = "".join(tmpList)
  cmd = "qsub -o %s -e %s -q %s %s %s %s %s %s %s %s %s %s %s" %(getLocalJobsDir(localjobdir),getLocalJobsDir(localjobdir),queue,bashfile,outfile,outdir,cmsswdir,cfg,inputFiles,maxevents,jobIndex,jobname,getLocalJobsDir(localjobdir))
  cmds.append(cmd)
  it+=nfiles
  jobIndex+=1

print "==============================================="

for c in cmds:
   print ""
   print c
   print ""
   if not(opts.test):
      os.system(c)
   
#python submitJobsOnT3batch.py -C submitJobsOnT3batch.cfg --useDAS

#crash M1000:
#job 1 (file 63) --> 6902352
#job 6 (file 49) --> 6902353
#job 9 (file 17) --> 6902354

#crash M2000:
#job 1
#job 5
#job 8
#job 9
