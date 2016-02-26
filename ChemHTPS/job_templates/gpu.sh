#!/usr/bin/env bash
#SBATCH --partition=gpu
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name="orcatest"
#SBATCH --output=slurm_orca.out

echo "SLURM job ID         = "$SLURM_JOB_ID
echo "Working Dir          = "$SLURM_SUBMIT_DIR
echo "Temporary scratch    = "$SLURMTMPDIR
echo "Compute Nodes        = "$SLURM_JOB_NODELIST
echo "Number of Processors = "$SLURM_NPROCS
echo "Number of Nodes      = "$SLURM_JOB_NUM_NODES
echo "Tasks per Node       = "$SLURM_TASKS_PER_NODE
echo "Memory per Node      = "$SLURM_MEM_PER_NODE

tmp=(${SLURM_SUBMIT_DIR//\// })
job=${tmp[${#tmp[@]} - 1]}
ulimit -s unlimited
module load orca
module list
which orca

echo "Launch job"
srun python job_runscript.py --scratch_dir $SLURMTMPDIR --submit_dir $SLURM_SUBMIT_DIR
#
echo "All Done!"
