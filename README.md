# EXOVVNtuplizerRunII

Ntuplizer for searches for heavy resonances decaying to dibosons

## installation instructions

Setting up CMSSW (for september reprocessing):

```
cmsrel CMSSW_10_2_10
cd CMSSW_10_2_10/src
cmsenv
git cms-init
```



### getting the ntuplizer code
```
cd $CMSSW_BASE/src
export GITUSER=`git config user.github`
git clone https://github.com/${GITUSER}/EXOVVNtuplizerRunII 
cd EXOVVNtuplizerRunII
git remote add UZHCMS https://github.com/UZHCMS/EXOVVNtuplizerRunII
git fetch UZHCMS
git checkout -b Development_BcMu_10210_NoLepProducers UZHCMS/BcMu_10210_NoLepProducers
cd $CMSSW_BASE/src
scram b -j 8
```

### JETMET filters to remove noise events
```
cd $CMSSW_BASE/src
git cms-addpkg RecoMET/METFilters 
scram b -j 8
```


### Update the electronIDs: not really used, but for consistency
(https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2)
```
git cms-merge-topic cms-egamma:EgammaPostRecoTools 
scram b -j 8

git cms-addpkg EgammaAnalysis/ElectronTools
rm EgammaAnalysis/ElectronTools/data -rf
git clone https://github.com/cms-egamma/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data
cd EgammaAnalysis/ElectronTools/data
git checkout ScalesSmearing2018_Dev
cd -
git cms-merge-topic cms-egamma:EgammaPostRecoTools_dev
scram b clean 
scam b distclean 
scram b -j 8 
```


### running for data and MC
Just set the proper flags in python/ntuplizerOptions_generic_cfi.py
and then run your config, as for example so :

```
cmsRun config_generic_opt_skimmed.py RunPeriod="Fall17" (on 2017 MC)
```
or
```
cmsRun config_generic_opt_skimmed.py RunPeriod="Run2017B" (for 2017 data)
```


### CRAB submission 
List the dataset list you would like to process in sample/sample.txt
Modify in the submission script the info relative to the user, the storage element (T2_CH_CSCS, T3_CH_PSI).
Examples of submissions are available in files like `commands_ForPeople.txt` and `commands_sample_submission_CRAB.txt`, and here:
you need to specify the directory where you want to store the crab project directory, the configaruion file that has to be run with `cmsRun`, the txt file with the sample list and any additional string you want to attach to the dataset name. Do a print of the command, before sending the jobs, to be sure that all options are well taken.

```
cmsenv
source /cvmfs/cms.cern.ch/crab3/crab.sh
python submit_all_opt.py -d CRAB_dir -c config_generic_opt_skimmed.py -f samples/sample.txt  -s ""

```


to recluster the MET the following flags can be changed:
```
config["DOMETRECLUSTERING"] = False
```


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
