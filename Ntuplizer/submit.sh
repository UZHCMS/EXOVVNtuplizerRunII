today=`date "+%Y%m%d"`
today2=`date "+%Y%m%d-%H%M%S"`


sed -i -e "s/isData = True/isData = False/" config_generic_opt_skimmed.py

for year in 2018

do
    echo $year

#    postfix="Legacy_${year}_${today}"
    postfix="JpsiK_multiple_${year}_${today}"

    sed -i -e "s/config\[\"ISBKG\"\] = True/config\[\"ISBKG\"\] = False/" config_generic_opt_skimmed.py
    
#    python submit_all.py -d crab/Crab_${today2}_${postfix} -c config_generic_opt_skimmed.py -f samples/BcJpsiTauNu_official_${year}.txt -s "_v1" --string "${postfix}" --isGlobal --numOfFiles 1


    if [ $year = "2018" ]; then

	sed -i -e "s/config\[\"ISBKG\"\] = False/config\[\"ISBKG\"\] = True/" config_generic_opt_skimmed.py

	echo $postfix
	python submit_all.py -d crab/Crab_${today2}_${postfix} -c config_generic_opt_skimmed.py -f samples/BcJpsiX.txt -s "_v1" --string "${postfix}" --isGlobal --numOfFiles 15


#	python submit_all.py -d crab/Crab_${today2}_${postfix} -c config_generic_opt_skimmed.py -f samples/JpsiX.txt -s "_v1" --string "${postfix}" --numOfFiles 20

    fi


done


