#!/usr/bin/env python

import os
import sys
from optparse import OptionParser

class SampleObject():
  def __init__(self, sample):
    self.fullName = sample.strip("\n")
    splitSample = self.fullName.split("/")
    self.base = splitSample[1]
    self.campaign = splitSample[2].split("-",1)[0]
    self.scenario = splitSample[2].split("-",1)[1]
    self.type = splitSample[3]

def main():

  parser = OptionParser(usage="usage: %prog [options] filelist")
  parser.add_option("-d", "--data", dest="data", default=False, action="store_true",
                    help="data sets are collision data and not MC")
  parser.add_option("-t", "--template", dest="template", default="submitJobsOnT3batch.cfg", action="store",
                    help="template file for config")
  parser.add_option("-p", "--pnfs", dest="pnfs", default="/pnfs/psi.ch/cms/trivcat/store", action="store",
                    help="/pnfs base location for data sets (assumes next subdirectory is either mc or data)")
  
  (options, args) = parser.parse_args()
  
  if len(args) != 1:
    parser.error("Please provide at file list name")

  fileList = args[0]
  isData = options.data
  templateCfg = options.template
  pnfs = options.pnfs

  print "Using file list:", fileList
  print "data:", isData
  print "Using template config:", templateCfg
  print "Using /pnfs base directory:", pnfs
  
  with open(fileList) as fileListFile:
    for sample in fileListFile:
      if sample.startswith("#") or sample.isspace():
        continue
      sample = SampleObject(sample)
      print "Working on sample:", sample.fullName
      createConfig(sample, templateCfg, pnfs, isData)

    
def createConfig(sample, templateCfg, pnfs, isData):
  cfgFileName = "%s_%s_%s.cfg" %(sample.base, sample.campaign, sample.scenario)
  cfgFile = open(cfgFileName, "w")
  with open(templateCfg) as templateCfgFile:
    for line in templateCfgFile:
      if line.startswith("#") or line.isspace():
        continue
      if (line.find("src") < 0):
        # print line.replace("FULLSAMPLE", sample.fullName).replace("SAMPLE", sample.base).strip("\n")
        cfgFile.write(line.replace("FULLSAMPLE", sample.fullName).replace("SAMPLE", sample.base))
      else:
        sampleLocation = findSampleLocation(sample, pnfs, isData)
        # print line.replace("FULLSAMPLE", sampleLocation).replace("SAMPLE", sample.base).strip("\n")
        cfgFile.write(line.replace("FULLSAMPLE", sampleLocation).replace("SAMPLE", sample.base))
  cfgFile.close()


def findSampleLocation(sample, pnfs, isData):
  subDir = "%s/%s/%s/%s" %(sample.campaign, sample.base, sample.type, sample.scenario)
  if isData:
    pnfs = pnfs+"/data"
  else:
    pnfs = pnfs+"/mc"
  fullPath = "%s/%s/" %(pnfs, subDir)
  # print fullPath
  for subdir, dirs, files in os.walk(fullPath):
    for file in files:
      if (file.find(".root") >= 0):
        return fullPath
  print "No ROOT file found for sample", sample.fullName, "- using global data set name."
  return sample.fullName


if __name__ == "__main__":
  main()