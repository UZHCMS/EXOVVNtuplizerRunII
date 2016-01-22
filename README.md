# EXOVVNtuplizerRunII

Ntuplizer for searches for heavy resonances decaying to dibosons

## installation instructions

```
export SCRAM_ARCH=slc6_amd64_gcc491
cmsrel CMSSW_7_4_12_patch2
cd CMSSW_7_4_12_patch2/src
cmsenv
git cms-init
```

### optional packages

For the boosted Hbb tagger (will add a lot of packages that will take a long time to compile):
```
git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git
git fetch btv-cmssw
git cms-merge-topic -u cms-btv-pog:BoostedDoubleSVTagger-WithWeightFiles-v2_from-CMSSW_7_4_1
```


### getting the code

```
git cms-addpkg RecoJets/Configuration
export GITUSER=`git config user.github`
echo "Your github username has been set to \"$GITUSER\""
git clone git@github.com:${GITUSER}/EXOVVNtuplizerRunII.git
cd EXOVVNtuplizerRunII
git remote add UZHCMS git@github.com:UZHCMS/EXOVVNtuplizerRunII.git
git fetch UZHCMS
git checkout -b DevelopmentBranch UZHCMS/master
cd $CMSSW_BASE/src
scram b distclean
scram b -j8
cd EXOVVNtuplizerRunII/Ntuplizer
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



