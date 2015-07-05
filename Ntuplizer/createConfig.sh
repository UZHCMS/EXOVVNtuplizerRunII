#!/bin/zsh

if [ "$#" -ne 1 ]; then
  echo "Please provide text file containing sample names as argument"
  exit 1
fi

for SAMPLE in `cat $1`; do
  BASENAME=`echo $SAMPLE | cut -d"/" -f 2`
  echo $SAMPLE
  echo $BASENAME
  sed -e "s|FULLSAMPLE|${SAMPLE}|g" submitJobsOnT3batch.cfg > tmp.cfg
  sed -e "s|SAMPLE|${BASENAME}|g" tmp.cfg > $BASENAME.cfg
  rm tmp.cfg
  
done