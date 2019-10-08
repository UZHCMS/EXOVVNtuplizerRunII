import FWCore.ParameterSet.Config as cms

config = dict()

#--------- general ----------#

#--------- Set Just one to true ----------#
config["RUNONMC"] = True
#-----------------------------------------#

config["USEJSON"] = not (config["RUNONMC"])
#config["JSONFILE"] = "JSON/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt" #data 2017
#config["JSONFILE"] = "JSON/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt" # data 2016
config["JSONFILE"] = "JSON/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt" # data 2018
config["USENOHF"] = False


#--------- basic sequences ----------#
config["DOGENPARTICLES"] = (True and config["RUNONMC"])
config["DOGENEVENT"] = (True and config["RUNONMC"])
config["DOPILEUP"] = (True and config["RUNONMC"])
config["DOVERTICES"] = True #True
config["DOTRIGGERDECISIONS"] = True
config["DOTRIGGEROBJECTS"] = False
config["DOHLTFILTERS"] = True
config["DOMISSINGET"] = True
config["DOMVAMET"] = False
config["DOJPSIMU"] = True
config["DOJPSIELE"] = False
config["DOGENHIST"] = (True and config["RUNONMC"]);
config["DOCUTFLOW"] = True;


#--------- MET reclustering ----------#
config["DOMETRECLUSTERING"] = False
config["DOMETSVFIT"] = True

#--------- JEC ----------#

config["CORRMETONTHEFLY"] = True  # at the moment JEC available just for MC Fall17

