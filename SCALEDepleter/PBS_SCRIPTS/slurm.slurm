#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=24
#SBATCH --time=60:00:00
#SBATCH --mem-per-cpu=11GB
#SBATCH --account=cgts-bpetrovic3-montecarlo-CODA20
#SBATCH --job-name=SCALE
#SBATCH --output=%x.%j.out
#SBATCH --nodelist=atl1-1-02-001-26-1,atl1-1-02-001-26-2,atl1-1-02-001-10-1,atl1-1-02-001-10-2

# good nodes:
### --nodelist=
#atl1-1-02-001-25-2
#atl1-1-02-001-26-1
#atl1-1-02-001-26-2
#atl1-1-02-001-3-1
#atl1-1-02-001-3-2
#atl1-1-02-001-4-1
#atl1-1-02-001-4-2
#atl1-1-02-001-5-1
#atl1-1-02-001-5-
######
#run1 through run10 node selections:
## atl1-1-02-001-1-1,atl1-1-02-001-1-2
## atl1-1-02-001-10-1,atl1-1-02-001-10-2
## atl1-1-02-001-11-1,atl1-1-02-001-11-2
## atl1-1-02-001-12-1,atl1-1-02-001-12-2
## atl1-1-02-001-13-1,atl1-1-02-001-14-1
## atl1-1-02-001-15-1,atl1-1-02-001-15-2
## atl1-1-02-001-16-1,atl1-1-02-001-16-2
## atl1-1-02-001-17-1,atl1-1-02-001-17-2
## atl1-1-02-001-18-1,atl1-1-02-001-18-2
## atl1-1-02-001-19-1,atl1-1-02-001-19-2

# create the machine file for slurm to run on
MACHINEFILE="nodes.$SLURM_JOB_ID" 
srun -l /bin/hostname | sort -n |awk 'BEGIN{FS="[:.]"};{print $2}' > "$MACHINEFILE"


# SLURM_NTASKS is automatically set by slrum
cd $SLURM_SUBMIT_DIR
module purge
module load scale/6.3.1
module load openblas
export TMPDIR=tmp
echo $TMPDIR

# final 
echo MACHINEFILE IS 
echo $MACHINEFILE
module rm python
source $HOME/miniconda3/bin/activate
conda activate scale
python --version
which python

# run
python -u celiCalculationInput.py $MACHINEFILE $TMPDIR 95
