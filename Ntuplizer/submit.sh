today=`date "+%Y%m%d"`
today2=`date "+%Y%m%d-%H%M%S"`


sed -i -e "s/isData = True/isData = False/" config_generic_opt_skimmed.py


#for year in 2016 2017 2018
#for year in 2018
for year in 2018
#for year in 2016 2017 2018
do
    echo $year

    postfix="MC_legacy_mc_${year}_${today}"
    python submit_all.py -d crab/Crab_${today}_${postfix} -c config_generic_opt_skimmed.py -f samples/BcJpsiTauNu_official_${year}.txt -s "_v1" --string "${postfix}" --isGlobal --numOfFiles 1


#    postfix="denominator_mc_${year}"
##    echo $psotfix
#    python submit_all.py -d crab/Crab_${today}_${postfix} -c config_generic_opt_skimmed.py -f samples/BcJpsiTauNu_official_${year}.txt -s "_v1" --string "${postfix}" --isGlobal --numOfFiles 1
#    python submit_all.py -d crab/Crab_${today}_${postfix} -c config_generic_opt_skimmed.py -f samples/BcJpsiMuNu_official_${year}.txt -s "_v1" --string "${postfix}" --isGlobal --numOfFiles 1

    if [ $year = "2018" ]; then
#	postfix="den_mc_${year}"
	postfix="MC_legacy_${year}_${today}"
	echo $postfix
	python submit_all.py -d crab/Crab_${today2}_${postfix} -c config_generic_opt_skimmed.py -f samples/BcJpsiX.txt -s "_v1" --string "${postfix}" --isGlobal --numOfFiles 5
#	python submit_all.py -d crab/Crab_${today2}_${postfix} -c config_generic_opt_skimmed.py -f samples/BcJpsiX_inclusive.txt -s "_v1" --string "${postfix}" --numOfFiles 3
    fi


done


