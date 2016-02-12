#!/bin/sh
#SBATCH --partition=debug
#SBATCH --time=00:50:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name="orcatest"
#SBATCH --output=slurm_orca.out
##SBATCH --mail-user=wevangel@buffalo.edu
##SBATCH --mail-type=END
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR

echo "working directory = "$SLURM_SUBMIT_DIR
tmp=(${SLURM_SUBMIT_DIR//\// })
job=${tmp[${#tmp[@]} - 1]}
ulimit -s unlimited
module load orca
module list
which orca

echo "Launch job"
orca $job".inp" >& $job".out"
#
echo "All Done!"
