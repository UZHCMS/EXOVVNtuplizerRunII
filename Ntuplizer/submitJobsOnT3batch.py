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

   cmd = "uberftp t3se01.psi.ch 'ls %s'" %src
   ls_la = commands.getoutput(cmd)

   list_ = []
   list_.extend(ls_la.split(os.linesep))   
   dirs = []
   for a in list_:
      b = a.split(" ")
      status = os.path.isdir(src + b[-1:][0].strip('\r'))
      if status: dirs.append(src + b[-1:][0].strip('\r'))

   if len(dirs) == 0: dirs.append(src)
   
   files = []
   for d in dirs:
      ls_la = commands.getoutput("ls " + d)
      list_ = []
      list_.extend(ls_la.split(os.linesep))
      for a in list_:
         b = a.split(" ")
         if b[-1:][0].find("root") != -1:
            files.append(d.replace('/pnfs/psi.ch/cms/trivcat','')+"/"+b[-1:][0])
   
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
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/"
   elif sample == "WJetsToLNu_HT-200to400":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/WJetsToLNu_HT-200to400_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/"
   elif sample == "WJetsToLNu_HT-400to600":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/WJetsToLNu_HT-400to600_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/"
   elif sample == "WJetsToLNu_HT-600toInf":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/WJetsToLNu_HT-600toInf_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/"         
   elif sample == "TBarToLeptons_s-channel":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/"    
   elif sample == "TBarToLeptons_t-channel":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/"   
   elif sample == "TTJets":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/"     
   elif sample == "TToLeptons_s-channel":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/"   
   elif sample == "TToLeptons_t-channel":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/"   
   elif sample == "T_tW-channel":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/"    
   elif sample == "Tbar_tW-channel":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/"   
   elif sample == "VV":
      path = "/pnfs/psi.ch/cms/trivcat/store/user/jngadiub/Phys14/VV/"
   elif sample == "QCD_HT_500To1000_ext":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/QCD_HT-500To1000_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1_ext1-v1/"   
   elif sample == "QCD_HT_250To500_ext":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/QCD_HT_250To500_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1_ext1-v2/"   
   elif sample == "QCD_HT_1000ToInf_ext":
      dataset ="/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/QCD_HT_1000ToInf_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1_ext1-v1/"   
   elif sample == "QCD_HT-500To1000":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/QCD_HT-500To1000_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/"          
   elif sample == "QCD_HT_1000ToInf":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/QCD_HT_1000ToInf_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/"            
   elif sample == "QCD_HT_250To500":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/QCD_HT_250To500_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/"         
   elif sample == "QCD_Pt-1000to1400":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/QCD_Pt-1000to1400_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_castor_PHYS14_25_V1-v1/"         
   elif sample == "QCD_Pt-1400to1800":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/QCD_Pt-1400to1800_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_castor_PHYS14_25_V1-v1/"         
   elif sample == "QCD_Pt-15to3000":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/"      
   elif sample == "QCD_Pt-1800to2400":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/QCD_Pt-1800to2400_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v2/"         
   elif sample == "QCD_Pt-2400to3200":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/QCD_Pt-2400to3200_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/"         
   elif sample == "QCD_Pt-300to470":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/QCD_Pt-300to470_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_castor_PHYS14_25_V1-v2/"         
   elif sample == "QCD_Pt-3200":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/QCD_Pt-3200_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/"         
   elif sample == "QCD_Pt-470to600":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/QCD_Pt-470to600_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_castor_PHYS14_25_V1-v2/"         
   elif sample == "QCD_Pt_600to800":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/QCD_Pt-600to800_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_castor_PHYS14_25_V1-v1/"         
   elif sample == "QCD_Pt_800to1000":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/QCD_Pt-800to1000_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_castor_PHYS14_25_V1-v2/"         
   elif sample == "RSGravitonToWW_M_1000":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/RSGravitonToWW_kMpl01_M_1000_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/"         
   elif sample == "RSGravitonToWW_M_2000":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/RSGravitonToWW_kMpl01_M_2000_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/"        
   elif sample == "RSGravitonToWW_M_3000":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/RSGravitonToWW_kMpl01_M_3000_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/"   
   elif sample == "RSGravitonToWW_M_4000":
      path = "/pnfs/psi.ch/cms/trivcat/store/mc/Phys14DR/RSGravitonToWW_kMpl01_M_4000_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/"
   else:
      print "Sample %s unknown! Exiting..." %sample
      sys.exit()
   
   return getFileListT3(path) 

#-----------------------------------------------------------------------------------------
def getFileListFromSampleDAS(sample):

   dataset = ""
   instance = "prod/global"
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
   elif sample == "VV":
      dataset = "/LOWW-lvjj-PTWgt180/qili-Q-Test-v3-7d492cb64f2cdaff326f939f96e45c96/USER"
      instance = "prod/phys03"
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
   
   return getFileListDAS(dataset,instance) 

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

   cmd = './das_client.py --query="file dataset=%s instance=%s" --limit=2000' %(dataset,instance)
   if run != -1: cmd = './das_client.py --query="file run=%i dataset=%s instance=%s" --limit=2000' %(run,dataset,instance)
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
parser.add_option("--useDAS", "--useDAS", dest="useDAS", default=False, action="store_true",
                              help="use das query")
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
xmlfile = config.get('JobsConfig','xmlfile')
maxevents = config.getint('JobsConfig','maxevents')
prefix = config.get('JobsConfig','prefix')
newdir = config.get('JobsConfig','newdir')
instance = config.get('JobsConfig','instance')
run = int(config.get('JobsConfig','run'))
 
#-----------------------------------------------------------------------------------------
if opts.copyfiles:

   jobsdir = getJobsDirs(outdir,jobname)
   user = commands.getoutput("whoami")	
    
   for j in jobsdir:
      jobid = j.rsplit("-",1)[1]
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

if opts.checkoutput:
   checkJobsOutputFromXML(xmlfile)
   sys.exit()

#-----------------------------------------------------------------------------------------	 
if opts.resubmit != "-1": 
  submitJobsFromXML(xmlfile,opts.resubmit.split(','))
  sys.exit()	
 
#-----------------------------------------------------------------------------------------	 
files = []

if not(opts.useDAS):
   if sample != "": files = getFileListFromSampleT3(sample)
   else: files = getFileListT3(src)   
else:
   if sample != "": files = getFileListFromSampleDAS(sample)
   else: files = getFileListDAS(src,instance,run)   

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
       if not(opts.useDAS): inputFiles+=f+"," 
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
