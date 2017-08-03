# EXOVVNtuplizerRunII

Ntuplizer for searches for heavy resonances decaying to dibosons

## installation instructions

```
# You can see my setup script here
/mnt/t3nfs01/data01/shome/ytakahas/public/setup_80X_legacy.sh

or in case you cannot access above, 

/afs/cern.ch/user/y/ytakahas/public/forUZH/setup_80X_legacy.sh
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

Finally, note that, when you run on crab, you have to enable 
```
config.JobType.sendExternalFolder = True
```
as described at https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariateElectronIdentificationRun2#Recipes_and_implementation
