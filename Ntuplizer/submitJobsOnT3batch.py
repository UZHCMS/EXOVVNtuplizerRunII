import os, commands
import random
import sys
from optparse import OptionParser
import ConfigParser
import time
import xml.etree.ElementTree as ET
import ROOT

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

  # find all .root files in src directory and add to fileList
  fileList = []
  for subdir, dirs, files in os.walk(src):
    for file in files:
      if (file.find(".root") >= 0):
        print os.path.join(src, file)
        fileList.append(os.path.join(src, subdir, file).replace("/pnfs/psi.ch/cms/trivcat",""))
   
   
  if len(fileList) == 0:
    print ("No files found in the directory %s" %(src))
    sys.exit(1)
  else:
    print ("Found %i files in the directory %s" %(len(fileList),src))  
   
  return fileList   


#-----------------------------------------------------------------------------------------
def writeXMLfile(xmlfile,cmds):

  print "Writing XML configuration file: "
  print "  * local job output dir = " + cmds[0].split(" ")[2]
  print "  * queue = " + cmds[0].split(" ")[6]
  print "  * bash file = " + cmds[0].split(" ")[7]
  print "  * output file = " + cmds[0].split(" ")[8]
  print "  * SE output dir = " + cmds[0].split(" ")[9]
  print "  * CMSSW dir = " + cmds[0].split(" ")[10]
  print "  * CMSSW config file = " + cmds[0].split(" ")[11]
  print "  * jobs name = " + cmds[0].split(" ")[15]
  
  root = ET.Element("JobsConfiguration")
  tree = ET.ElementTree(root)																	      

  generaltag = ET.SubElement(root,'General')
  tag = ET.SubElement(generaltag,'LocalOutputDir')
  tag.text = cmds[0].split(" ")[2]
  tag = ET.SubElement(generaltag,'Queue')
  tag.text = cmds[0].split(" ")[6]
  tag = ET.SubElement(generaltag,'BashFile')
  tag.text = cmds[0].split(" ")[7]    
  tag = ET.SubElement(generaltag,'OutFile')
  tag.text = cmds[0].split(" ")[8]
  tag = ET.SubElement(generaltag,'SEOutDir')
  tag.text = cmds[0].split(" ")[9]
  tag = ET.SubElement(generaltag,'CMSSWDir')
  tag.text = cmds[0].split(" ")[10]
  tag = ET.SubElement(generaltag,'CMSSWConfigFile')
  tag.text = cmds[0].split(" ")[11]
  tag = ET.SubElement(generaltag,'JobName')
  tag.text = cmds[0].split(" ")[15]    
    
  for c in cmds:
     jobtag = ET.SubElement(root,'Jobs')
     jobtag.attrib['InputFiles'] = c.split(" ")[12]
     jobtag.attrib['MaxEvents'] = c.split(" ")[13]
     jobtag.attrib['ID'] = c.split(" ")[14]
     
  indent(root)       
  tree.write(xmlfile)

#-----------------------------------------------------------------------------------------
def indent(elem, level=0):
  i = "\n" + level*"  "
  if len(elem):
    if not elem.text or not elem.text.strip():
      elem.text = i + "  "
    if not elem.tail or not elem.tail.strip():
      elem.tail = i
    for elem in elem:
      indent(elem, level+1)
    if not elem.tail or not elem.tail.strip():
      elem.tail = i
  else:
    if level and (not elem.tail or not elem.tail.strip()):
      elem.tail = i

#-----------------------------------------------------------------------------------------
def submitJobsFromXML(xmlfile,jobslist):

   tree = ET.parse(xmlfile)																		      
   root = tree.getroot()
   
   cmds = []
   for j in root.findall('Jobs'):
      for id in jobslist:
         if j.get('ID') == id:
            cmd = "qsub -o %s -e %s -q %s %s %s %s %s %s %s %s %s %s %s" %(root.find('General').find('LocalOutputDir').text,root.find('General').find('LocalOutputDir').text,root.find('General').find('Queue').text,root.find('General').find('BashFile').text,root.find('General').find('OutFile').text,root.find('General').find('SEOutDir').text,root.find('General').find('CMSSWDir').text,root.find('General').find('CMSSWConfigFile').text,j.get('InputFiles'),j.get('MaxEvents'),j.get('ID'),root.find('General').find('JobName').text,root.find('General').find('LocalOutputDir').text)
            cmds.append(cmd) 

   status,cmd_out = commands.getstatusoutput( 'ls '+getLocalJobsDir(localjobdir) )
   if status: os.system('mkdir '+getLocalJobsDir(localjobdir))

   for c in cmds:
      print ""
      print c
      if not(opts.test):
         os.system(c)

#----------------------------------------------------------------------------------------- 
def checkJobsOutputFromXML(xmlfile):

   tree = ET.parse(xmlfile)																		      
   root = tree.getroot()
   ROOT.gROOT.ProcessLine('.x rootlogon.C')

   jobsdir = getJobsDirs(outdir,jobname)
   user = commands.getoutput("whoami")	

   jobsevents = {}
   print "=========================================="          
   for j in jobsdir:
      a = j.rsplit("-",1)
      jobid = a[1]
      inputpath = j+"/"+outfile
      checkfile = "ls -l %s"%(inputpath) 
      status,cmd_out = commands.getstatusoutput(checkfile)
      if status: continue
      tfile = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125/"+inputpath)
      if not tfile.Get("ntuplizer/tree"):
         print "WARNING: tree not found! Job %s : found 0 events" %(jobid)
	 jobsevents[jobid] = 0
	 continue
      ttree = ROOT.TTree()
      tfile.GetObject("ntuplizer/tree",ttree)
      print "Job %s : found %i events" %(jobid,ttree.GetEntries())
      jobsevents[jobid] = ttree.GetEntries()
    
   expevents = {} 
   print "=========================================="      
   for j in root.findall('Jobs'):
      inputFiles = j.get('InputFiles')
      filelist = inputFiles.split(",")
      count = 0
      for f in filelist:
         if f.find('xrootd') == -1: tfile = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat"+f)
         else: tfile = ROOT.TFile.Open(f)
	 ttree = ROOT.TTree()
         tfile.GetObject("Events",ttree)
         count+=ttree.GetEntries()
      print "Job %s : expected %i events" %(j.get('ID'),count) 
      expevents[j.get('ID')] = count

   print "==========================================" 
   print "Expected %i jobs, found %i" %(len(expevents),len(jobsevents))           
   for j,e in expevents.iteritems():
      for jj,ee in jobsevents.iteritems():
         if j == jj and e != ee: print "Job "+str(j)+": found " + str(jobsevents[j]) + " expected " + str(e)    
            
#-----------------------------------------------------------------------------------------
def getFileListDAS(dataset,instance="prod/global",run=-1):

   cmd = 'das_client.py --query="file dataset=%s instance=%s" --limit=10000' %(dataset,instance)
   if run != -1: cmd = 'das_client.py --query="file run=%i dataset=%s instance=%s" --limit=10000' %(run,dataset,instance)
   print cmd
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
parser.add_option("--check", "--checkoutput", dest="checkoutput", default=False, action="store_true",
                              help="check jobs output")
parser.add_option("-r", "--resubmit", dest="resubmit", default="-1", action="store", type="string",
                              help="resubmit jobs from xml file")
parser.add_option("-H", "-H", dest="help", default=False, action="store_true",
                              help="usage")	      			      
(opts, args) = parser.parse_args(argv)

if opts.config =="":
    opts.config = "config"

if opts.help:
   print ""
   print "Usage for testing: "
   print "          python submitJobsOnT3batch.py -C submitJobsOnT3batch.cfg --test"
   print "Or: "
   print "          python submitJobsOnT3batch.py -C submitJobsOnT3batch.cfg"
   print ""
   sys.exit()
   
print opts.config
config = ConfigParser.ConfigParser()
config.read(opts.config)
nfiles = config.getint('JobsConfig','nfiles')
src = config.get('JobsConfig','src')
outdir = config.get('JobsConfig','outdir')
localjobdir = config.get('JobsConfig','localjobdir')
jobname = config.get('JobsConfig','jobname')
outfile = config.get('JobsConfig','outfile')
queue = config.get('JobsConfig','queue')
cmsswdir = config.get('JobsConfig','cmsswdir')
cfg = config.get('JobsConfig','cfg')
bashfile = config.get('JobsConfig','bashfile')
xmlfile = config.get('JobsConfig','xmlfile')
maxevents = config.getint('JobsConfig','maxevents')
prefix = config.get('JobsConfig','prefix')
newdir = config.get('JobsConfig','newdir')
instance = config.get('JobsConfig','instance')
run = int(config.get('JobsConfig','run'))

# check if input data are on local SE (path starting with/containing /pnfs) or if CMS dataset name has been given
# in case CMS dataset name has been given, need to use DAS script later to determine file names
useDAS = True
if (src.find("/pnfs") >= 0):
  useDAS = False
 
#-----------------------------------------------------------------------------------------
if opts.copyfiles:

   jobsdir = getJobsDirs(outdir,jobname)
   user = commands.getoutput("whoami")	
    
   print ""
   print "" 
   print "..... Copying configuration ......" 
   checkdir = "ls -l " + newdir + "/config"
   status,cmd_out = commands.getstatusoutput(checkdir)
   if status: 
     cmd = "srmmkdir srm://t3se01.psi.ch/%s"%(newdir) #No gFAL!
     print cmd
     os.system(cmd)
     # cmd = "gfal-mkdir srm://t3se01.psi.ch:8443/srm/managerv2?SFN=%s"%(newdir+"/config")
     cmd = "srmmkdir srm://t3se01.psi.ch/%s"%(newdir+"/config") #No gFAL!
     # cmd = "env -i X509_USER_PROXY=~/.x509up_u`id -u` gfal-mkdir -p gsiftp://t3se01.psi.ch%s" %(newdir)
     print cmd
     os.system(cmd)
     
   else:
     cmd = "uberftp t3se01.psi.ch 'ls %s'" %(newdir+"/config")
     ls_la = commands.getoutput(cmd)

     list_ = []
     list_.extend(ls_la.split(os.linesep))   
     dirs = []
     for a in list_:
       b = a.split(" ")
       status = os.path.exists(newdir + "/config/" + b[-1:][0].strip('\r'))
       if status:
         # cmd = "gfal-rm srm://t3se01.psi.ch:8443/srm/managerv2?SFN=%s"%(newdir + "/config/" + b[-1:][0].strip('\r'))
         cmd = "srm-rm srm://t3se01.psi.ch%s"%(newdir + "/config/" + b[-1:][0].strip('\r')) #No gFAL!
         os.system(cmd)
   
   # cmd = "lcg-cp -b -D srmv2 " + xmlfile + " srm://t3se01.psi.ch:8443/srm/managerv2?SFN=%s"%(newdir+"/config/"+os.path.basename(xmlfile))
   cmd = "xrdcp -d 1 " + xmlfile + " root://t3dcachedb.psi.ch:1094///%s"%(newdir+"/config/"+os.path.basename(xmlfile)) #No gFAL!
   print cmd
   os.system(cmd)
   # cmd = "lcg-cp -b -D srmv2 " + cmsswdir + "/src/" + cfg + " srm://t3se01.psi.ch:8443/srm/managerv2?SFN=%s -f"%(newdir+"/config/"+os.path.basename(cfg))
   cmd = "xrdcp -d 1 " + cmsswdir + "/src/" + cfg + " root://t3dcachedb.psi.ch:1094///%s -f"%(newdir+"/config/"+os.path.basename(cfg)) #No gFAL!
   print cmd
   os.system(cmd)
   f = open(cmsswdir + "/src/" + cfg, 'r')
   ismc = False
   for l in f:
      if l.find("ntuplizerOptions_MC_cfi") != -1: ismc = True
   f.close()
   if ismc:
      # cmd = "lcg-cp -b -D srmv2 python/ntuplizerOptions_MC_cfi.py srm://t3se01.psi.ch:8443/srm/managerv2?SFN=%s"%(newdir+"/config/ntuplizerOptions_MC_cfi.py")
      cmd = "xrdcp -d 1 python/ntuplizerOptions_MC_cfi.py root://t3dcachedb.psi.ch:1094///%s"%(newdir+"/config/ntuplizerOptions_MC_cfi.py")
   else:
      # cmd = "lcg-cp -b -D srmv2 python/ntuplizerOptions_data_cfi.py srm://t3se01.psi.ch:8443/srm/managerv2?SFN=%s"%(newdir+"/config/ntuplizerOptions_data_cfi.py")
      cmd = "xrdcp -d 1 python/ntuplizerOptions_data_cfi.py root://t3dcachedb.psi.ch:1094///%s"%(newdir+"/config/ntuplizerOptions_data_cfi.py")
   print cmd
   os.system(cmd)
   # cmd = "lcg-cp -b -D srmv2 " + opts.config[0] + " srm://t3se01.psi.ch:8443/srm/managerv2?SFN=%s"%(newdir+"/config/"+os.path.basename(opts.config[0]))
   cmd = "xrdcp -d 1 " + opts.config[0] + " root://t3dcachedb.psi.ch:1094///%s"%(newdir+"/config/"+os.path.basename(opts.config[0]))
   print cmd
   os.system(cmd)

   print ""
   print ""
   print "..... Copying job files ......" 
    
   for j in jobsdir:
      jobid = j.rsplit("-",1)[1]
      # inputpath = "srm://t3se01.psi.ch:8443/srm/managerv2?SFN="+j+"/"+outfile
      # outputpath = "srm://t3se01.psi.ch:8443/srm/managerv2?SFN=%s"%(newdir+"/"+prefix+"_"+jobid+".root")
      inputpath = "root://t3dcachedb.psi.ch:1094//"+j+"/"+outfile #USE XROOTD!!
      outputpath = "root://t3dcachedb.psi.ch:1094//%s"%(newdir+"/"+prefix+"_"+jobid+".root") #USE XROOTD!!
      checkfile = "ls -l %s"%(newdir+"/"+prefix+"_"+jobid+".root") 
      status,cmd_out = commands.getstatusoutput(checkfile)
      if not(status): 
         # cmd = "gfal-rm %s"%(outputpath)
         cmd = "srmrm %s"%(outputpath)
         print cmd
         os.system(cmd)
      # cmd = "gfal-copy %s %s" %(inputpath,outputpath)
      cmd = "xrdcp -d 1 %s %s -f" %(inputpath,outputpath)#TEST
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
	 cmd = "srmrm "+inputpath
	 print cmd
	 os.system(cmd)
      #jdir = "srm://t3se01.psi.ch:8443/srm/managerv2?SFN="+j
      cmd = "uberftp t3se01.psi.ch 'rm -r %s'" %j
      print cmd
      os.system(cmd)
    
   print ""
   cmd = "rm -rf "+getLocalJobsDir(localjobdir)
   print cmd
   os.system(cmd)
   user = commands.getoutput("whoami")
   dir = '/pnfs/psi.ch/cms/trivcat/store/user/'+user+"/"+outdir
   cmd = "uberftp t3se01.psi.ch 'rm -r " + dir + "'"
   print cmd
   os.system(cmd)
      
   sys.exit()   

if opts.checkoutput:
   checkJobsOutputFromXML(xmlfile)
   sys.exit()

#-----------------------------------------------------------------------------------------	 
if opts.resubmit != "-1": 
  submitJobsFromXML(xmlfile,opts.resubmit.split(','))
  sys.exit()	
 
#-----------------------------------------------------------------------------------------	 
files = []

if not(useDAS):
   files = getFileListT3(src)   
else:
   files = getFileListDAS(src,instance,run)   

print "****************************************"
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
       if not(useDAS): inputFiles+=f+"," 
       else: inputFiles+="root://xrootd.unl.edu/%s," %(f)	  
  tmpList = list(inputFiles)
  tmpList.pop(len(tmpList)-1)
  inputFiles = "".join(tmpList)
  cmd = "qsub -o %s -e %s -q %s %s %s %s %s %s %s %s %s %s %s" %(getLocalJobsDir(localjobdir),getLocalJobsDir(localjobdir),queue,bashfile,outfile,outdir,cmsswdir,cfg,inputFiles,maxevents,jobIndex,jobname,getLocalJobsDir(localjobdir))
  cmds.append(cmd)
  it+=nfiles
  jobIndex+=1

print "==============================================="
writeXMLfile(xmlfile,cmds)

for c in cmds:
   print ""
   print c
   if not(opts.test):
      os.system(c)
