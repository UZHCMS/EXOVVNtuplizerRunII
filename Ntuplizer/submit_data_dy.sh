today=`date "+%Y%m%d"`


sed -i -e "s/isData = False/isData = True/" config_generic_opt_skimmed.py
sed -i -e "s/config\[\"ISBKG\"\] = True/config\[\"ISBKG\"\] = False/" config_generic_opt_skimmed.py

for year in 2018
#for year in 2017
do
    echo $year
    postfix="TauValidation_${year}_${today}"
#    postfix="Approval_GJSON_${year}_${today}"
#    postfix="den_data_${year}"
#    postfix="JpsiPi_GJSON_${year}_${today}"
#    postfix="JpsiK_multiple_${year}_${today}"

    python submit_all.py -d crab/Crab_${today}_${postfix} -c config_generic_opt_skimmed.py -f samples/data_singleMuon_${year}.txt -s "_v1" --string "${postfix}"  --luminosity --isGlobal --isData
done


