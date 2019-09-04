import FWCore.ParameterSet.Config as cms

config = dict()

#--------- general ----------#

#--------- Set Just one to true ----------#
config["RUNONMC"] = True
#-----------------------------------------#

config["USEJSON"] = not (config["RUNONMC"])
config["JSONFILE"] = "JSON/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt" #"Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt"

config["USENOHF"] = False
config["FILTEREVENTS"] = False

#--------- basic sequences ----------#
config["DOGENPARTICLES"] = (True and config["RUNONMC"])
config["DOGENEVENT"] = (True and config["RUNONMC"])
config["DOPILEUP"] = (True and config["RUNONMC"])
config["DOELECTRONS"] = False #True
config["DOMUONS"] = False #True
config["DOTAUS"] = False #True
config["DOVERTICES"] = True #True
config["DOTRIGGERDECISIONS"] = True
config["DOTRIGGEROBJECTS"] = True
config["DOHLTFILTERS"] = True
config["DOMISSINGET"] = True
config["DOMVAMET"] = False
config["DOJPSIMU"] = True
config["DOJPSIELE"] = False

#--------- MET reclustering ----------#
config["DOMETRECLUSTERING"] = False
config["DOMETSVFIT"] = True

#--------- JEC ----------#

config["CORRMETONTHEFLY"] = True  # at the moment JEC available just for MC Fall17

#--------- TAU ----------#
config["DOMULTIPLETAUMVAVERSIONS"] = False#True #This flag eneables the possibility to access a sqlite *db file and save the latest training of the tau MVA isolation "v2" in parellel as the one of "v1" taken from the CMSSW database.
