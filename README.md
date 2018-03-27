	# EXOVVNtuplizerRunII

Ntuplizer for searches for heavy resonances decaying to dibosons

## installation instructions

Setting up CMSSW (for september reprocessing):

```
cmsrel CMSSW_9_4_0
cd CMSSW_9_4_0/src
cmsenv
git init
```



### getting the ntuplizer codes
```
cd $CMSSW_BASE/src
export GITUSER=`git config user.github`
git clone https://github.com/${GITUSER}/EXOVVNtuplizerRunII 
cd EXOVVNtuplizerRunII
git remote add UZHCMS https://github.com/UZHCMS/EXOVVNtuplizerRunII
git fetch UZHCMS
git checkout -b DevelopmentBranch_9_4_0 UZHCMS/94X_ntuplizer
cd $CMSSW_BASE/src
```


### updates of latest cut-based electron ID and MVAID. HEEP are default in 94X
(https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2 and 
```
cd $CMSSW_BASE/src
cmsenv
git cms-init
git cms-merge-topic lsoffi:CMSSW_9_4_0_pre3_TnP
git cms-merge-topic guitargeek:ElectronID_MVA2017_940pre3
scram b -j 8
cd $CMSSW_BASE/external
cd slc6_amd64_gcc630/
git clone https://github.com/lsoffi/RecoEgamma-PhotonIdentification.git data/RecoEgamma/PhotonIdentification/data
cd data/RecoEgamma/PhotonIdentification/data
git checkout CMSSW_9_4_0_pre3_TnP
cd $CMSSW_BASE/external
cd slc6_amd64_gcc630/
git clone https://github.com/lsoffi/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
cd data/RecoEgamma/ElectronIdentification/data
git checkout CMSSW_9_4_0_pre3_TnP
cd $CMSSW_BASE/src
scram b distclean
scramv1 b -j 16
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
