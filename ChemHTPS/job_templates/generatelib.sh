#!/usr/bin/env bash
Timehere
Nodeshere
#SBATCH --job-name="generatelib"
#SBATCH --output=generatelib.out
Cpushere

# ====================================================
# For 16-core nodes
# ====================================================
##SBATCH --constraint=CPU-E5-2630v3
##SBATCH --cpus-per-task=16
##SBATCH --mem=64000

echo "SLURM job ID         = "$SLURM_JOB_ID
echo "Working Dir          = "$SLURM_SUBMIT_DIR
echo "Temporary scratch    = "$SLURMTMPDIR
echo "Compute Nodes        = "$SLURM_JOB_NODELIST
echo "Number of Processors = "$SLURM_NPROCS
echo "Number of Nodes      = "$SLURM_JOB_NUM_NODES
echo "Tasks per Node       = "$SLURM_TASKS_PER_NODE
echo "Memory per Node      = "$SLURM_MEM_PER_NODE

ulimit -s unlimited
module load python
module load intel-mpi

export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi.so

echo "Launch job"
Runlinehere
echo "All Done!"
