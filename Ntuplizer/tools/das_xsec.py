#!/usr/bin/env python

import subprocess
import sys
import os
import json
import optparse

def queryDAS(command):
  subProcess = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
  output = subProcess.communicate()
  jsonOutput = json.loads(output[0])
  if (jsonOutput["status"] != "ok"):
    print "Failure in query: %s" %command
    sys.exit()
  return jsonOutput


def main():

  dataset = "/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/Fall13-POSTLS162_V1-v1/GEN-SIM"
  dataset = "/WJetsToLNu_HT-100to200_Tune4C_13TeV-madgraph-tauola/*/*"
  dataset = "/QCD*Pt*13TeV*/*/*"
  dataset = "/QCDddbar_Pt-15to3000_TuneZ2star_Flat_13TeV_pythia6/RunIISpring15DR74-AsymptNoPU_MCRUN2_74_V9A-v1/MINIAODSIM"
  
  parser=optparse.OptionParser(usage="%prog -j jobOption.py")
  parser.add_option("-d", "--dataset", dest="dataset",
                      action="store", default=dataset,
                      help="dataset pattern to look for")
  (options, args)=parser.parse_args()
  dataset=options.dataset
  
  print "Looking for datasets of pattern: %s" %dataset

  command = 'das_client.py --query=\"dataset=%s\" --format=json --limit=0' %(dataset)
  # print command
  jsonOutput = queryDAS(command)
  
  

  # print jsonOutput
  # print "----"
  # print jsonOutput.keys()
  # print "----"
  # print jsonOutput["nresults"]
  # print "----"
  # print jsonOutput["data"][0]
  # print "----"
  # print jsonOutput["data"][0]["dataset"][0].keys()
  # print "----"
  # print jsonOutput["data"][0]["dataset"][0]["name"]
  #
  # print "---- loop ----"
  
  print "Found %i datasets" %(jsonOutput["nresults"])
  
  print "Will loop over datasets found to find MINIAODs"
  for i in range(jsonOutput["nresults"]):
    thisDataset = jsonOutput["data"][i]["dataset"][0]["name"]
    # print thisDataset
    if (thisDataset.find("MINIAOD") >= 0):
      print "Investigating %s" %thisDataset
      eventsCommand = 'das_client.py --query=\"dataset=%s\" --format=json --limit=0' %(thisDataset)
      eventsOutput = queryDAS(eventsCommand)
      # print eventsOutput["data"][0]["dataset"]
      print "Number of events: %i" %eventsOutput["data"][0]["dataset"][0]["nevents"]
      print "Finding parent:"
      parentDataset = thisDataset
      maxSteps = 5
      currentStep = 0
      while ((parentDataset.find("/GEN") < 0) and (currentStep < maxSteps)):
        parentCommand = 'das_client.py --query=\"parent dataset=%s\" --format=json --limit=0' %(parentDataset)
        parentOutput = queryDAS(parentCommand)
        print parentOutput["data"][0]["parent"][0]["name"]
        parentDataset = parentOutput["data"][0]["parent"][0]["name"]
        currentStep += 1
      if (currentStep == maxSteps):
        print "ERROR: Couldn't find parent of type GEN for %s" %thisDataset
        continue
      xsecCommand = 'das_client.py --query=\"mcm dataset=%s | grep mcm.generator_parameters\" --format=json --limit=0' %(parentDataset)
      xsecOutput = queryDAS(xsecCommand)
      # print xsecOutput["data"][0]["mcm"][0]["generator_parameters"]
      print "Cross-section [pb]: %f" %xsecOutput["data"][0]["mcm"][0]["generator_parameters"][0]["cross_section"]
      print "Filter efficiency: %f +/- %f" %(xsecOutput["data"][0]["mcm"][0]["generator_parameters"][0]["filter_efficiency"], xsecOutput["data"][0]["mcm"][0]["generator_parameters"][0]["filter_efficiency_error"])
      print "Matching efficiency: %f +/- %f" %(xsecOutput["data"][0]["mcm"][0]["generator_parameters"][0]["match_efficiency"], xsecOutput["data"][0]["mcm"][0]["generator_parameters"][0]["match_efficiency_error"])
      print "%s %i %f pb" %(dataset,eventsOutput["data"][0]["dataset"][0]["nevents"],xsecOutput["data"][0]["mcm"][0]["generator_parameters"][0]["cross_section"])
      
  
  
  
  # number of events from MINIAOD
  # xsec for GEN-SIM

  # os.system('python das_client.py --query=\"parent dataset=%s\" --format=json' %(dataset))

if __name__ == "__main__":
  main()