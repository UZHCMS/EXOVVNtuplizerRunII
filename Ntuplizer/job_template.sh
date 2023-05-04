#!/bin/bash
# 
#SBATCH -p wn
#SBATCH --account=t3
#SBATCH --time 01:00:00
#SBATCH -e cn-test.err 
#SBATCH -o cn-test.out  # replace default slurm-SLURM_JOB_ID.out

echo HOME: $HOME 
echo USER: $USER 
echo SLURM_JOB_ID: $SLURM_JOB_ID
echo HOSTNAME: $HOSTNAME

# each worker node has local /scratch space to be used during job run
mkdir -p /scratch/$USER/${SLURM_JOB_ID}
export TMPDIR=/scratch/$USER/${SLURM_JOB_ID}

cmsRun config_generic_opt_skimmed_batch.py INFILE $TMPDIR/tmp_ONAME_IDJ.root
xrdcp -f $TMPDIR/tmp_ONAME_IDJ.root root://t3dcachedb03.psi.ch/OUTFILE

# cleaning of temporal working dir when job was completed:
rm -rf /scratch/$USER/${SLURM_JOB_ID}
