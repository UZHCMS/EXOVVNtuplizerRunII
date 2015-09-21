# EXOVVNtuplizerRunII
Ntuplizer for searches for heavy resonances decaying to dibosons

## installation instructions

```
export SCRAM_ARCH=slc6_amd64_gcc491
cmsrel CMSSW_7_4_7_patch2
cd CMSSW_7_4_7_patch2/src
cmsenv
git cms-init
```

For the boosted Hbb tagger (will add a lot of packages that will take a long time to compile):
```
git remote add btv-cmssw https://github.com/cms-btv-pog/cmssw.git
git fetch btv-cmssw
git cms-merge-topic -u cms-btv-pog:BoostedDoubleSVTagger-WithWeightFiles-v2_from-CMSSW_7_4_1
```

For energy correlation function and D2 (ATLAS W-tagger):
```
git cms-merge-topic ahinzmann:ecf_74
```

For PUPPI v8:
```
git cms-addpkg CommonTools/PileupAlgos
git cms-merge-topic nhanvtran:puppi-etadep-746p2-v8
```

getting the code:
```
git cms-addpkg RecoJets/Configuration
export GITUSER=`git config user.github`
echo "Your github username has been set to \"$GITUSER\""
git clone https://github.com/$GITUSER/EXOVVNtuplizerRunII
cd EXOVVNtuplizerRunII
git remote add UZHCMS https://github.com/UZHCMS/EXOVVNtuplizerRunII
git fetch UZHCMS
git checkout -b DevelopmentBranch UZHCMS/master
cd $CMSSW_BASE/src
scram b distclean
scram b -j8
cd EXOVVNtuplizerRunII/Ntuplizer
```

to run:
```
cmsRun config.py
```

to recluster jets and MET, the following flags can be changed:
```
doAK8reclustering = False
doAK8softdropReclustering = False
doBtagging = False
doAK8prunedReclustering = False
doMETReclustering = False
```
If you want to use Higgs tagger the first three flags must all be set to True.
