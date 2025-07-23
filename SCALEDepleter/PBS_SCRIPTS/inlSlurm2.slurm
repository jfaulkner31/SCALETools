#!/bin/bash

# Submit this script with: sbatch thefilename
################# #SBATCH --mem-per-cpu=4G
################# #SBATCH -p short # queue (short queue use -p short after sbatch


#SBATCH --time=120:00:00 # walltime
#SBATCH --ntasks-per-node=35 # number of processor cores (i.e. tasks)
#SBATCH --nodes=2 # number of nodes
#SBATCH --mem-per-cpu=6G
#SBATCH --wckey edu_res # Project Code
#SBATCH -J "homgnzd_case1" # job name
#SBATCH -e slurm-%j.err # %j will be replaced by the SLURM_JOB_ID
#SBATCH -o slurm-%j.out # %j will be replaced by the SLURM_JOB_ID
#SBATCH --mail-user=jfaulkner31@gatech.edu # email address
#SBATCH --mail-type=BEGIN


## some echo stuff
echo omp_num_threads = $OMP_NUM_THREADS

## loading scale and data for scale 6.3
module load use.exp_ctl
module load scale/6.3.2
export DATA=/hpc-common/data/scale/6.3

# assign tmpdir
export TMPDIR=tmp
echo $TMPDIR

# Make the machine file for slurm
echo NODELIST is: $SLURM_JOB_NODELIST
MACHINEFILE="nodes.$SLURM_JOB_ID"
srun -l /bin/hostname | sort -n |awk 'BEGIN{FS="[:.]"};{print $2}' > "$MACHINEFILE"

# Load python that has what we need (np and  mpl)
ml py-openmc/0.15.0-gcc-13.3.0-ompi5-mylq

# cd into working dir
cd $SLURM_SUBMIT_DIR

# Set PML for OpenMPI -- use ob1 and ofi since this is Infiniband network using OmniPath (bitterroot)
export OMPI_MCA_pml=ob1
export OMPI_MCA_mtl=psm2

# use this for python -- python -u (unbuffered mode) inputfilename.py nodefile tmpdir numcores
python -u _INPUT_FILE_.py ${MACHINEFILE} $TMPDIR 69


##################
##### NOTES #####
##################


# old run for scalerte on inl
### scalerte -vvv -N 47 -M ${PBS_NODEFILE} homogenization.inp -m -T /tmp

# more hints:
# https://www.nrel.gov/hpc/assets/pdfs/pbs-to-slurm-translation-sheet.pdf



########## NOTES #################
### when i ssh into thje node i see that every core is using about 4.6g
## what i did to get this to run basic FHR depletion problem ----
# run on bitterroot with 256 gb per node.
# do: #SBATCH --mem-per-cpu=6G
# do 38 cores
# 375 iterations
# 0.32 min per iteration
# 2 hours for 1 MC run
# 17*4 MC runs = 68
# = 136 hours for MC
# = 6 days runtime on bitterroot
