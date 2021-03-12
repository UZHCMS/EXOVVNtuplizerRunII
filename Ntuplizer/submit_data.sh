today=`date "+%Y%m%d"`


#for year in 2016
for year in 2018
do
    echo $year
#    postfix="den_data_${year}"
    postfix="Legacy_v2_data_${year}_${today}"

    python submit_all.py -d crab/Crab_${today}_${postfix} -c config_generic_opt_skimmed.py -f samples/data_charmonium_${year}.txt -s "_v1" --string "${postfix}"  --luminosity --isGlobal --isData
done


