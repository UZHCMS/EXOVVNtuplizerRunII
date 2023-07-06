today=`date "+%Y%m%d"`
today2=`date "+%Y%m%d-%H%M%S"`


sed -i -e "s/isData = True/isData = False/" config_generic_opt_skimmed.py

for year in 2018

do
    echo $year

    postfix="TauValidation_${year}_${today}"
#    postfix="Approval_${year}_${today}"
#    postfix="Freeze_${year}_${today}"
#    postfix="JpsiPi_${year}_${today}"
#    postfix="JpsiK_${year}_${today}"


    sed -i -e "s/config\[\"ISBKG\"\] = True/config\[\"ISBKG\"\] = False/" config_generic_opt_skimmed.py
    
    python submit_all.py -d crab/Crab_${today2}_${postfix} -c config_generic_opt_skimmed.py -f samples/TauValidation_${year}.txt -s "_v1" --string "${postfix}" --isGlobal --numOfFiles 2



#    sed -i -e "s/config\[\"ISBKG\"\] = False/config\[\"ISBKG\"\] = True/" config_generic_opt_skimmed.py

#    python submit_all.py -d crab/Crab_${today2}_${postfix} -c config_generic_opt_skimmed.py -f samples/BJpsiX_${year}.txt -s "_v1" --string "${postfix}" --isGlobal --numOfFiles 15

done


