import FWCore.ParameterSet.Config as cms

config = dict()

#--------- general ----------#
config["FALL15"] = True
config["RUNONMC"] = False
config["USEJSON"] = True
config["JSONFILE"] = "JSON/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_Silver_v2.txt"
config["FILTEREVENTS"] = True
config["BUNCHSPACING"] = 25
config["USENOHF"] = False

#--------- basic sequences ----------#
config["DOGENPARTICLES"] = False
config["DOGENJETS"] = False
config["DOGENEVENT"] = False
config["DOPILEUP"] = False
config["DOELECTRONS"] = True
config["DOMUONS"] = True
config["DOTAUS"] = False
config["DOAK8JETS"] = True
config["DOAK4JETS"] = True
config["DOVERTICES"] = True
config["DOTRIGGERDECISIONS"] = True
config["DOTRIGGEROBJECTS"] = True
config["DOHLTFILTERS"] = True
config["DOMISSINGET"] = True
config["DOTAUSBOOSTED"] = False
config["DOMETSVFIT"] = False


#--------- AK8 jets reclustering ----------#
config["ADDAK8GENJETS"] = False #! Add AK8 gen jet collection with pruned and softdrop mass
config["DOAK8RECLUSTERING"] =False
config["DOAK8PRUNEDRECLUSTERING"] = False #! To add pruned jet and pruned subjet collection (not in MINIAOD)
config["DOAK8PUPPIRECLUSTERING"] = True
config["DOAK10TRIMMEDRECLUSTERING"] = True #ATLAS sequence
config["DOHBBTAG"] = False #Higgs-tagger

#--------- MET reclustering ----------#
config["DOMETRECLUSTERING"] = False

#--------- JEC ----------#
config["CORRJETSONTHEFLY"] = True
config["CORRMETONTHEFLY"] = True
config["GETJECFROMDBFILE"] = False # If not yet in global tag, but db file available
