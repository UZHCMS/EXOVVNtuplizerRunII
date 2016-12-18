import FWCore.ParameterSet.Config as cms

config = dict()

#--------- general ----------#
config["SPRING16"] = True
config["RUNONMC"] = True
config["USEJSON"] = False
config["JSONFILE"] = "Cert_271036-279931_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt"
config["BUNCHSPACING"] = 25
config["USENOHF"] = False
config["FILTEREVENTS"] = False

#--------- basic sequences ----------#
config["DOGENPARTICLES"] = True
config["DOGENJETS"] = True
config["DOGENEVENT"] = True
config["DOPILEUP"] = True
config["DOELECTRONS"] = True
config["DOMUONS"] = True
config["DOTAUS"] = True
config["DOAK8JETS"] = True
config["DOAK4JETS"] = True
config["DOVERTICES"] = True
config["DOTRIGGERDECISIONS"] = True
config["DOTRIGGEROBJECTS"] = True
config["DOHLTFILTERS"] = True
config["DOMISSINGET"] = True
config["DOTAUSBOOSTED"] = True
config["DOMETSVFIT"] = True
config["DOMVAMET"] = False

#--------- AK8 jets reclustering ----------#
config["ADDAK8GENJETS"] = True #! Add AK8 gen jet collection with pruned and softdrop mass
config["DOAK8RECLUSTERING"] = False
config["DOAK8PRUNEDRECLUSTERING"] = False #! To add pruned jet and pruned subjet collection (not in MINIAOD)
config["DOAK8PUPPI"] = True
config["DOAK10TRIMMEDRECLUSTERING"] = False #ATLAS sequence
config["DOHBBTAG"] = True #Higgs-tagger
config["DOAK8PUPPIRECLUSTERING"] = False
config["UpdateJetCollection"] = True #needed for Higgs-tagger in 80X

#--------- MET reclustering ----------#
config["DOMETRECLUSTERING"] = False

#--------- JEC ----------#
config["CORRJETSONTHEFLY"] = True
config["CORRMETONTHEFLY"] = True
config["GETJECFROMDBFILE"] = False # If not yet in global tag, but db file available
