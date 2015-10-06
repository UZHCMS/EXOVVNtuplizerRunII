import FWCore.ParameterSet.Config as cms

config = dict()

#--------- general ----------#
config["RUNONMC"] = True
config["USEJSON"] = False
config["JSONFILE"] = 'goldenJSON_PromptReco.txt'
config["BUNCHSPACING"] = 25

#--------- basic sequences ----------#
config["DOGENPARTICLES"] = True
config["DOGENJETS"] = False
config["DOGENEVENT"] = True
config["DOPILEUP"] = True
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
config["DOSEMILEPTONICTAUSBOOSTED"] = False

#--------- AK8 jets reclustering ----------#
config["ADDAK8GENJETS"] = False #! Add AK8 gen jet collection with pruned and softdrop mass
config["DOAK8RECLUSTERING"] = True
config["DOAK8PRUNEDRECLUSTERING"] = True #! To add pruned jet and pruned subjet collection (not in MINIAOD)
config["DOAK8PUPPIRECLUSTERING"] = False #ATLAS sequence
config["DOAK10TRIMMEDRECLUSTERING"] = False
config["DOHBBTAG"] = False #Higgs-tagger

#--------- MET reclustering ----------#
config["DOMETRECLUSTERING"] = False

#--------- JEC ----------#
config["CORRJETSONTHEFLY"] = True
config["CORRMETONTHEFLY"] = True
config["GETJECFROMDBFILE"] = False # If not yet in global tag, but db file available
