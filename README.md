	# EXOVVNtuplizerRunII

Ntuplizer for searches for heavy resonances decaying to dibosons

## installation instructions

Setting up CMSSW (for september reprocessing):

```
cmsrel CMSSW_9_2_13
cd CMSSW_9_2_13/src
cmsenv
git init
```

### getting the latest b-tagger
(https://twiki.cern.ch/twiki/bin/viewauth/CMS/Hbbtagging#V4_training)

=> included by default in 92X

### update MET filter. Skip this!

=> no obvious branch exists (seems to run fine 92X baseline, can be reomved?)
```
git cms-merge-topic -u cms-met:CMSSW_8_0_X-METFilterUpdate
```


### getting the ntuplizer codes
```
cd $CMSSW_BASE/src
export GITUSER=`git config user.github`
git clone https://github.com/${GITUSER}/EXOVVNtuplizerRunII 
cd EXOVVNtuplizerRunII
git remote add UZHCMS https://github.com/UZHCMS/EXOVVNtuplizerRunII
git fetch UZHCMS
git checkout -b DevelopmentBranch_9_2_13 UZHCMS/92X_legacy
cd $CMSSW_BASE/src
```

The flags for running on Spring15(74) or Fall15(76) or Spring16(80) samples have to be changed with config["FALL15"]=False/True and config["SPRING16"]=False/True in python/ntuplizerOptions_*_cfi.py


### updates of latest cut-based electron ID
(https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_for_8_0)

=> cut based ele id included by default in 92X

### updates of latest HEEP electron ID
(https://twiki.cern.ch/twiki/bin/viewauth/CMS/HEEPElectronIdentificationRun2#Recipe_for_regular_users)

=> HEEP ID from here https://twiki.cern.ch/twiki/bin/viewauth/CMS/HEEPElectronIdentificationRun2

```
cmsenv
git cms-addpkg RecoEgamma/EgammaIsolationAlgos 
git cms-merge-topic rgoldouz:TrkIsoFix -u 
scramv1 b -j 16
```




### compile first, before adding MVA-based electron ID 
(I don't konw why, but otherwise, $CMSSW_BASE/external/slc6_amd64_gcc530/data area will be deleted)
```
cd $CMSSW_BASE/src
scram b -j8
```

### updates of latest MVA-based electron ID. Skip this!
(https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentificationRun2#Recipes_and_implementation)

=> I didn't look into this one yet, but seems to be included in 92X
```
mkdir -p $CMSSW_BASE/external/slc6_amd64_gcc530/data/RecoEgamma/ElectronIdentification/
cd $CMSSW_BASE/external/slc6_amd64_gcc530/
git clone https://github.com/ikrav/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
cd data/RecoEgamma/ElectronIdentification/data
git checkout egm_id_80X_v1
cd $CMSSW_BASE/src
scram b -j8
```






### running for data and MC! just set the proper flag in python/ntuplizerOptions_generic_cfi.py

```
cmsRun config_generic.py 

```


to recluster jets and MET, or to add the Higgs-tagger the following flags can be changed:
```
config["DOAK8RECLUSTERING"] = False
config["DOHBBTAG"] = False
config["DOAK8PRUNEDRECLUSTERING"] = False
config["DOMETRECLUSTERING"] = False
```
If you want to use Higgs tagger the first two flags must all be set to True.

### Batch submission

#### Config file creation

Config file creation can be done via the [createConfig.py](Ntuplizer/tools/createConfig.py) script. It requires a text file with a list of input data sets, see e.g. [samples/QCD_HT_RunIISpring15MiniAODv2.txt](Ntuplizer/samples/QCD_HT_RunIISpring15MiniAODv2.txt). To run:
```
python tools/createConfig.py samples/QCD_HT_RunIISpring15MiniAODv2.txt
```
When running over *data*, this requires the ```-d``` flag. The script will automatically determine if the data sets are available on the T3 storage element. Also, ```--help``` will provide more information (e.g. allows changing the default number of jobs per event). If you run the script from a different directory, you need to provide the location of the [template file](Ntuplizer/submitJobsOnT3batch.cfg).

#### Job submission

Submit your jobs using the [submitJobsOnT3batch.py](Ntuplizer/submitJobsOnT3batch.py) script with the generated config files like this:
```
python submitJobsOnT3batch.py -C myconfig.cfg
```
Once the jobs are done, they can be checked for completeness like this:
```
python submitJobsOnT3batch.py -C myconfig.cfg --check
```
Resubmit jobs like this:
```
python submitJobsOnT3batch.py -C myconfig.cfg --resubmit 1,4,7
```
And eventually copied to the SE (path given in the config file):
```
python submitJobsOnT3batch.py -C myconfig.cfg --copy
```

Finally, note that, when you run on crab, you have to enable 
```
config.JobType.sendExternalFolder = True
```
as described at https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentificationRun2#Recipes_and_implementation
