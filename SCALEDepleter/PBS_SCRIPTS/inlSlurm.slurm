#!/bin/bash

# Submit this script with: sbatch thefilename

#SBATCH --time=1:00:00 # walltime
#SBATCH --ntasks-per-node=48 # number of processor cores (i.e. tasks)
#SBATCH --nodes=1 # number of nodes
#SBATCH --wckey edu_res # Project Code
#SBATCH -J "variableTemperature" # job name
#SBATCH -e slurm-%j.err # %j will be replaced by the SLURM_JOB_ID
#SBATCH -o slurm-%j.out # %j will be replaced by the SLURM_JOB_ID
#SBATCH --mail-user=jfaulkner31@gatech.edu # email address
#SBATCH --mail-type=BEGIN

## some echo stuff
echo omp_num_threads = $OMP_NUM_THREADS

## loading scale and data for scale 6.3
module load use.exp_ctl
module load scale/6.3.1
export DATA=/hpc-common/data/scale/6.3

# assign tmpdir
export TMPDIR=tmp
echo $TMPDIR


echo $SLURM_JOB_NODELIST
echo
echo
echo ALL DONE!






# cd into working dir
#cd $PBS_O_WORKDIR

# use this for python -- python -u (unbuffered mode) inputfilename.py nodefile tmpdir numcores
#python -u _INPUT_FILE_.py ${PBS_NODEFILE} $TMPDIR 47

# old run for scalerte on inl
### scalerte -vvv -N 47 -M ${PBS_NODEFILE} homogenization.inp -m -T /tmp
