for dir in /pnfs/psi.ch/cms/trivcat/store/user/cgalloni/MINIAOD/*
do
    dir=${dir%*/}      # remove the trailing "/"
    name="${dir##*/}"
#    echo     # print everything after the final "/"
    echo ${dir}, ${name}

    python getDataset.py --path ${dir} --analysis ${name}

done



for dir in /pnfs/psi.ch/cms/trivcat/store/user/sleontsi/MINIAOD/*
do
    dir=${dir%*/}      # remove the trailing "/"
    name="${dir##*/}"
#    echo     # print everything after the final "/"
    echo ${dir}, ${name}

    python getDataset.py --path ${dir} --analysis ${name}

done



