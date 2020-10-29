import FWCore.ParameterSet.Config as cms

config = dict()

#--------- general ----------#

#--------- Set Just one to true ----------#
config["RUNONMC"] = True
#-----------------------------------------#
#config["USEHAMMER"] = (True and config["RUNONMC"])
config["USEHAMMER"] = False
config["VERBOSE"] = False

#--------- For taus ----------#
config["USEDNN"] = True
config["DZCUT"] = 0.25 # this is fixed !!
config["FSIGCUT"] = 3
config["VPROBCUT"] = 0.05
config["DNNCUT"] = 0.1443
#config["DNNCUT"] = 0.0012
config["TAU_CHARGE"] = 1

config["USEJSON"] = not (config["RUNONMC"])
#config["USEJSON"] = False
#config["JSONFILE"] = "JSON/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt" #data 2017
#config["JSONFILE"] = "JSON/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt" # data 2016
#config["JSONFILE"] = "JSON/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt" # data 2017UL

#config["JSONFILE"] = "JSON/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt" # data 2018
config["JSONFILE"] = "JSON/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
#config["JSONFILE"] = "JSON/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
#config["JSONFILE"] = "JSON/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"

config["USENOHF"] = False


#--------- basic sequences ----------#
config["DOGENPARTICLES"] = (True and config["RUNONMC"])
config["DOGENEVENT"] = (True and config["RUNONMC"])
config["DOPILEUP"] = (True and config["RUNONMC"])
config["DOVERTICES"] = True
config["DOMISSINGET"] = True

config["DOJPSIMU"] = False
config["DOJPSITAU"] = True
config["DOBSTAUTAU"] = False
config["DOBSTAUTAUFH"] = False
config["DOBSTAUTAUFH_mr"] = False # mass regression
config["DOBSDSTARTAUNU"] = False
config["ISTRUTH"] = False

config["DOGENHIST"] = (True and config["RUNONMC"]);

if config["DOJPSIMU"]:
    config["USEDNN"] = False

if config["DOJPSITAU"]:
    config["DNNFILE"] = "data/DNN/BcJPsi/DUMMY"

elif config["DOBSTAUTAU"] or config["DOBSDSTARTAUNU"]:
    config["DNNFILE"] = "data/DNN/BsTauTau_semilep/DUMMY"    

elif config["DOBSTAUTAUFH"] or config["DOBSTAUTAUFH_mr"]:
    config["DNNFILE"] = "data/DNN/BsTauTau_fullhad/DUMMY"
    config["DNNCUT"] = 0.47
    config["FSIGCUT"] = 2
    config["VPROBCUT"] = 0.05
else:
    config["DNNFILE"] = "data/DNN/BcJPsi/DUMMY"


#--------- JEC ----------#

config["CORRMETONTHEFLY"] = False  # at the moment JEC available just for MC Fall17

