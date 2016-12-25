# EXOVVNtuplizerRunII

Ntuplizer for searches for heavy resonances decaying to dibosons

## installation instructions


For Spring16(80):

```
cmsrel CMSSW_8_0_21
cd CMSSW_8_0_21/src
cmsenv

```

### getting the latest b-tagger
### (see, https://twiki.cern.ch/twiki/bin/viewauth/CMS/Hbbtagging#V4_training)

```
export CMSSW_GIT_REFERENCE="/cvmfs/cms.cern.ch/cmssw.git.daily"
git cms-init
git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git
git fetch --tags btv-cmssw
git cms-merge-topic -u cms-btv-pog:BoostedDoubleSVTaggerV4-WithWeightFiles-v1_from-CMSSW_8_0_21
```

### update MET filter

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
git checkout -b MoriondDevelopmentBranch UZHCMS/80X_ntuplizer
cd $CMSSW_BASE/src
```

The flags for running on Spring15(74) or Fall15(76) or Spring16(80) samples have to be changed with config["FALL15"]=False/True and config["SPRING16"]=False/True in python/ntuplizerOptions_*_cfi.py


### updates with latest electron ID

```
git cms-merge-topic ikrav:egm_id_80X_v2
git cms-merge-topic Sam-Harper:HEEPV70VID
git cms-merge-topic ikrav:egm_id_80X_v3

git cms-merge-topic Sam-Harper:PackedCandNoPuppi
cd $CMSSW_BASE/src
scram b -j8

```

### add MVA based electron ID

```
mkdir -p $CMSSW_BASE/external/slc6_amd64_gcc530/data/RecoEgamma/ElectronIdentification/
cd $CMSSW_BASE/external/slc6_amd64_gcc530/
git clone https://github.com/ikrav/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
cd data/RecoEgamma/ElectronIdentification/data
git checkout egm_id_80X_v1

cd $CMSSW_BASE/src
scram b -j8
```


### running

```
cmsRun config_data.py (for data)
cmsRun config_MC.py (for MC)
```

the flags for running on data can be changed in python/ntuplizerOptions_data_cfi.py
the flags for running on MC can be changed in python/ntuplizerOptions_MC_cfi.py

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
