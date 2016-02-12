#!/usr/bin/env python

import os
import sys
from optparse import OptionParser

class SampleObject():
  def __init__(self, sample, isData):
    self.fullName = sample.strip("\n")
    splitSample = self.fullName.split("/")
    self.base = splitSample[1]
    self.campaign = splitSample[2].split("-",1)[0]
    self.scenario = splitSample[2].split("-",1)[1]
    self.type = splitSample[3]
    self.isData = isData
    self.name = "%s_%s_%s" %(self.base, self.campaign, self.scenario)
    
# a few global variables
cmsswBase = ""
nfiles = -1


def main():
  
  global cmsswBase
  global nfiles

  parser = OptionParser(usage="usage: %prog [options] filelist")
  parser.add_option("-d", "--data", dest="data", default=False, action="store_true",
                    help="data sets are collision data and not MC [default: %default]")
  parser.add_option("-t", "--template", dest="template", default="submitJobsOnT3batch.cfg", action="store",
                    help="template file for config [default: %default]")
  parser.add_option("-p", "--pnfs", dest="pnfs", default="/pnfs/psi.ch/cms/trivcat/store", action="store",
                    help="/pnfs base location for data sets (assumes next subdirectory is either mc or data) [default: %default]")
  parser.add_option("-n", "--nfiles", dest="nfiles", default="10", action="store",
                    help="number of files per job [default: %default]")

  
  (options, args) = parser.parse_args()
  
  if len(args) != 1:
    parser.error("Please provide at file list name")

  fileList = args[0]
  isData = options.data
  templateCfg = options.template
  pnfs = options.pnfs
  nfiles = options.nfiles
  try:
    int(nfiles)
  except ValueError:
    print "nfiles = %s - is not an integer value"
  cmsswBase = os.getenv("CMSSW_BASE")
  if not cmsswBase:
    print "CMSSW not set up (CMSSW_BASE does not exist)."
    sys.exit(1)

  print "Using file list:", fileList
  print "data:", isData
  print "Using template config:", templateCfg
  print "Using /pnfs base directory:", pnfs
  print "Using %s files per job" %nfiles
  print "Using CMSSW_BASE:", cmsswBase
  print "-"*80
  
  with open(fileList) as fileListFile:
    for sample in fileListFile:
      if sample.startswith("#") or sample.isspace():
        continue
      sample = SampleObject(sample, isData)
      print "Working on sample:", sample.fullName
      print "Sample name:", sample.name
      createConfig(sample, templateCfg, pnfs)

    
def createConfig(sample, templateCfg, pnfs):
  cfgFileName = "%s.cfg" %(sample.name)
  cfgFile = open(cfgFileName, "w")
  sampleLocation = findSampleLocation(sample, pnfs)
  with open(templateCfg) as templateCfgFile:
    for line in templateCfgFile:
      if line.startswith("#") or line.isspace():
        continue
      if (line.find("src") < 0):
        line = line.replace("FULLSAMPLE", sample.fullName)
      else:
        line = line.replace("FULLSAMPLE", sampleLocation)
      line = line.replace("SAMPLE", sample.name)
      line = line.replace("NFILES", nfiles)
      line = line.replace("CMSSW_BASE", cmsswBase)
      if sample.isData:
        # change this if the config names should ever change
        line = line.replace("config_MC.py", "config_data.py")
      cfgFile.write(line)
  cfgFile.close()


def findSampleLocation(sample, pnfs):
  subDir = "%s/%s/%s/%s" %(sample.campaign, sample.base, sample.type, sample.scenario)
  if sample.isData:
    pnfs = pnfs+"/data"
  else:
    pnfs = pnfs+"/mc"
  fullPath = "%s/%s/" %(pnfs, subDir)
  # print fullPath
  for subdir, dirs, files in os.walk(fullPath):
    for file in files:
      if (file.find(".root") >= 0):
        return fullPath
  print "No ROOT file found for sample", sample.fullName, " in ", fullPath,"- using global data set name."
  return sample.fullName


if __name__ == "__main__":
  main()
