import FWCore.ParameterSet.Config as cms

config = dict()

#--------- general ----------#

#--------- Set Just one to true ----------#
config["RUNONMC"] = False
config["RUNONPromptReco"] = True
config["RUNONReReco"] = False
#-----------------------------------------#

config["USEJSON"] = False
config["JSONFILE"] = "JSON/Cert_294927-305636_13TeV_PromptReco_Collisions17_JSON.txt"
config["BUNCHSPACING"] = 25
config["USENOHF"] = False
config["FILTEREVENTS"] = False

#--------- basic sequences ----------#
config["DOGENPARTICLES"] = (True and config["RUNONMC"])
config["DOGENJETS"] = (True and config["RUNONMC"])
config["DOGENEVENT"] = (True and config["RUNONMC"])
config["DOPILEUP"] = (True and config["RUNONMC"])
config["DOELECTRONS"] = True
config["DOMUONS"] = True
config["DOTAUS"] = True
config["DOAK8JETS"] = True
config["DOAK4JETS"] = True
config["DOVERTICES"] = True
config["DOTRIGGERDECISIONS"] = True
config["DOTRIGGEROBJECTS"] = False
config["DOHLTFILTERS"] = True
config["DOMISSINGET"] = True
config["DOTAUSBOOSTED"] = True
config["DOMETSVFIT"] = True
config["DOMVAMET"] = False

#--------- AK8 jets reclustering ----------#
config["ADDAK8GENJETS"] = True #! Add AK8 gen jet collection with pruned and softdrop mass
config["DOAK8RECLUSTERING"] = False
config["DOAK8PRUNEDRECLUSTERING"] = False #! To add pruned jet and pruned subjet collection (not in MINIAOD)
config["DOAK8PUPPI"] = False
config["DOAK10TRIMMEDRECLUSTERING"] = False #ATLAS sequence
config["DOHBBTAG"] = False #Higgs-tagger
config["DOAK8PUPPIRECLUSTERING"] = False
config["UpdateJetCollection"] = False #needed for Higgs-tagger in 80X

#--------- MET reclustering ----------#
config["DOMETRECLUSTERING"] = False

#--------- JEC ----------#
config["CORRJETSONTHEFLY"] = False
config["CORRMETONTHEFLY"] = False
config["GETJECFROMDBFILE"] = False # If not yet in global tag, but db file available
