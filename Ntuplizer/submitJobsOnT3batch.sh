#!/bin/bash
#################################
# PSI Tier-3 example batch Job  #
#################################

##### CONFIGURATION ##############################################
# Output files to be copied back to the User Interface
# (the file path must be given relative to the working directory)
OUTFILES="myout.txt myerr.txt"

# Output files to be copied to the SE
# (as above the file path must be given relative to the working directory)
SEOUTFILES=$1

#
# By default, the files will be copied to $USER_SRM_HOME/$username/$JOBDIR,
# but here you can define the subdirectory under your SE storage path
# to which files will be copied (uncomment line)
SEUSERSUBDIR=$2

#
# User's CMS hypernews name (needed for user's SE storage home path
# USER_SRM_HOME below)
HN_NAME=`whoami`

# set DBG=1 for additional debug output in the job report files
# DBG=2 will also give detailed output on SRM operations
DBG=2

#### The following configurations you should not need to change
# The SE's user home area (SRMv2 URL)
USER_SRM_HOME="srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user"
USER_XRD_HOME="root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/"

# Top working directory on worker node's local disk. The batch
# job working directory will be created below this
TOPWORKDIR=/scratch/`whoami`

JOB_ID=$7

# Basename of job sandbox (job workdir will be $TOPWORKDIR/$JOBDIR)
JOBDIR=$8-$JOB_ID
##################################################################


############ BATCH QUEUE DIRECTIVES ##############################
# Lines beginning with #$ are used to set options for the SGE
# queueing system (same as specifying these options to the qsub
# command

# Job name (defines name seen in monitoring by qstat and the
#     job script's stderr/stdout names)
#$ -N flatTuple_job

### Specify the queue on which to run
#$ -q all.q

# Change to the current working directory from which the job got
# submitted (will also result in the job report stdout/stderr being
# written to this directory)
#$ -cwd

# here you could change location of the job report stdout/stderr files
#  if you did not want them in the submission directory

##$ -o $8
##$ -e $8

##################################################################

##### MONITORING/DEBUG INFORMATION ###############################
DATE_START=`date +%s`
echo "Job started at " `date`
cat <<EOF
################################################################
## QUEUEING SYSTEM SETTINGS:
HOME=$HOME
USER=$USER
JOB_ID=$JOB_ID
JOB_NAME=$JOB_NAME
HOSTNAME=$HOSTNAME
TASK_ID=$TASK_ID
QUEUE=$QUEUE

EOF

echo "################################################################"

if test 0"$DBG" -gt 0; then
   echo "######## Environment Variables ##########"
   env
   echo "################################################################"
fi


##### SET UP WORKDIR AND ENVIRONMENT ######################################
STARTDIR=`pwd`
WORKDIR=$TOPWORKDIR/$JOBDIR
RESULTDIR=$STARTDIR/$JOBDIR

gfal-mkdir $USER_SRM_HOME/$HN_NAME/$SEUSERSUBDIR/$JOBDIR
if test x"$SEUSERSUBDIR" = x; then
   SERESULTDIR=$USER_SRM_HOME/$HN_NAME/$JOBDIR
   SERESULTDIRXRD=$USER_XRD_HOME/$HN_NAME/$JOBDIR
else
   SERESULTDIR=$USER_SRM_HOME/$HN_NAME/$SEUSERSUBDIR/$JOBDIR  
   SERESULTDIRXRD=$USER_XRD_HOME/$HN_NAME/$SEUSERSUBDIR/$JOBDIR    

fi
if test -e "$WORKDIR"; then
   echo "ERROR: WORKDIR ($WORKDIR) already exists! Aborting..." >&2
   rm -rf $WORKDIR
   #exit 1
fi
mkdir -p $WORKDIR
if test ! -d "$WORKDIR"; then
   echo "ERROR: Failed to create workdir ($WORKDIR)! Aborting..." >&2
   exit 1
fi

cd $WORKDIR
cat <<EOF
################################################################
## JOB SETTINGS:
STARTDIR=$STARTDIR
WORKDIR=$WORKDIR
RESULTDIR=$RESULTDIR
SERESULTDIR=$SERESULTDIR
EOF

###########################################################################
## YOUR FUNCTIONALITY CODE GOES HERE
# set up CMS environment
source $VO_CMS_SW_DIR/cmsset_default.sh

# Here we produce some output for the files that are copied back to
#  our shared home
scramv1 list > myout.txt 2>myerr.txt

gfal-ls srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/ >> myout.txt 2>>myerr.txt
#lcg-ls -b -D srmv2 -l srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/ >> myout.txt 2>>myerr.txt

# create a dummy file for copying back to the SE
dd if=/dev/urandom of=mybigfile count=100 &>/dev/null

CMSSW_DIR=$3
CMSSW_CONFIG_FILE=$CMSSW_DIR/src/$4


cd $CMSSW_DIR/src
eval `scramv1 runtime -sh`
if test $? -ne 0; then
   echo "ERROR: Failed to source scram environment" >&2
   exit 1
fi

cd $WORKDIR

cp -r $CMSSW_DIR/src/EXOVVNtuplizerRunII/Ntuplizer/JEC .
cp -r $CMSSW_DIR/src/EXOVVNtuplizerRunII/Ntuplizer/JER .
cp -r $CMSSW_DIR/src/EXOVVNtuplizerRunII/Ntuplizer/JSON .
cp -r $CMSSW_DIR/src/EXOVVNtuplizerRunII/Ntuplizer/RunLumiEventLists .

cmsRun $CMSSW_CONFIG_FILE inputFiles=$5 maxEvents=$6>> myout.txt 2>>myerr.txt


#### RETRIEVAL OF OUTPUT FILES AND CLEANING UP ############################
cd $WORKDIR
if test 0"$DBG" -gt 0; then
    echo "########################################################"
    echo "############# Working directory contents ###############"
    echo "pwd: " `pwd`
    ls -Rl
    echo "########################################################"
    echo "YOUR OUTPUT WILL BE MOVED TO $RESULTDIR"
    echo "########################################################"
fi

if test x"$OUTFILES" != x; then
   mkdir -p $RESULTDIR
   if test ! -e "$RESULTDIR"; then
          echo "ERROR: Failed to create $RESULTDIR ...Aborting..." >&2
          exit 1
   fi
   for n in $OUTFILES; do
       if test ! -e $WORKDIR/$n; then
          echo "WARNING: Cannot find output file $WORKDIR/$n. Ignoring it" >&2
       else
          cp -a $WORKDIR/$n $RESULTDIR/$n
          if test $? -ne 0; then
             echo "ERROR: Failed to copy $WORKDIR/$n to $RESULTDIR/$n" >&2
          fi
   fi
   done
fi

if test x"$SEOUTFILES" != x; then
   if test 0"$DBG" -ge 2; then
      srmdebug="-v"
   fi
   for n in $SEOUTFILES; do
       if test ! -e $WORKDIR/$n; then
          echo "WARNING: Cannot find output file $WORKDIR/$n. Ignoring it" >&2
       else
          #lcg-cp $srmdebug -b -D srmv2 file:$WORKDIR/$n $SERESULTDIR/$n
          # gfal-copy $srmdebug file:$WORKDIR/$n $SERESULTDIR/$n
          xrdcp -d 1 $WORKDIR/$n $SERESULTDIRXRD/$n -f
          if test $? -ne 0; then
             echo "ERROR: Failed to copy $WORKDIR/$n to $SERESULTDIR/$n" >&2
          fi
   fi
   done
fi

echo "Cleaning up $WORKDIR"
rm -rf $WORKDIR

cd $STARTDIR
mv $JOBDIR $9/


###########################################################################
DATE_END=`date +%s`
RUNTIME=$((DATE_END-DATE_START))
echo "################################################################"
echo "Job finished at " `date`
echo "Wallclock running time: $RUNTIME s"
exit 0





