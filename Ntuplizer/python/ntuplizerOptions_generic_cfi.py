import FWCore.ParameterSet.Config as cms

config = dict()

# only common flags are described here !!!
# the rest will be dynamically switched later in the config


#--------- Verbose setting ----------#
config["VERBOSE"] = False

#--------- For filtering ----------#
config["DZCUT"] = 0.12
config["FSIGCUT"] = 3
config["VPROBCUT"] = 0.00001
config["TAU_CHARGE"] = 1

#config["JSONFILE"] = "JSON/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt" # data 2018
config["JSONFILE"] = "JSON/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"
#config["JSONFILE"] = "JSON/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
#config["JSONFILE"] = "JSON/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"


#config["JSONFILE"] = "JSON/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt"
#config["JSONFILE"] = "JSON/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt"
#config["JSONFILE"] = "JSON/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt"


#--------- basic sequences ----------#
config["DOVERTICES"] = True
config["DOMISSINGET"] = False

config["DOJPSIMU"] = False
config["DOJPSITAU"] = True
config["DOJPSIK"] = False
config["DOJPSIKE"] = False
config["DOJPSIKE"] = False
config["DOBSTAUTAUFH"] = False

config["DNNFILE_PERPF"] = "data/DNN/BcJPsi/TAU_UL18/DUMMY"
config["DNNFILE_PEREVT_MC"] = "data/DNN/BcJPsi/TAU_MCBKG/DUMMY"
config["DNNFILE_PEREVT_DATA"] = "data/DNN/BcJPsi/TAU_DATABKG/DUMMY"
config["BWEIGHTFILE"] = "data/BWEIGHT/decay_weight.root"

#--------- JEC ----------#

config["CORRMETONTHEFLY"] = False  # at the moment JEC available just for MC Fall17

