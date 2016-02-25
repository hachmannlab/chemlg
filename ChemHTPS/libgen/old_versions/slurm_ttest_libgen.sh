#!/bin/sh
#SBATCH --partition=general-compute
#SBATCH --time=100:00:00
#SBATCH --job-name="ttest_libgen2"
#SBATCH --output=ttest_libgen2.out
#SBATCH --error=ttest_libgen2.err
#SBATCH --clusters=chemistry
#SBATCH --partition=beta
#SBATCH --account=pi-hachmann

#SBATCH --nodes=6
#SBATCH --cpus-per-task=16

tic=`date +%s`
echo "Start Time = "`date`
echo "Loading modules ..."

module load python
module load intel intel-mpi

ulimit -s unlimited
module list



echo "SLURM job ID         = "$SLURM_JOB_ID
echo "Working Dir          = "`pwd`

echo "Compute Nodes        = "`nodeset -e $SLURM_NODELIST`
echo "Number of Processors = "$SLURM_NPROCS
echo "Number of Nodes      = "$SLURM_NNODES 


# for i in $(seq 12 12 240); do  
#     echo $i
#     mpirun -np $i python readpara.py 
# done 



#mpirun -np 240 python lib_generator_v09.py
#mpirun -np 216 python lib_generator_v09.py
#mpirun -np 192 python lib_generator_v09.py
#mpirun -np 168 python lib_generator_v09.py
#mpirun -np 144 python lib_generator_v09.py
#mpirun -np 120 python lib_generator_v09.py
mpirun -np 96 python lib_generator_v09.py
mpirun -np 80 python lib_generator_v09.py
mpirun -np 72 python lib_generator_v09.py
mpirun -np 60 python lib_generator_v09.py
mpirun -np 48 python lib_generator_v09.py
mpirun -np 36 python lib_generator_v09.py
mpirun -np 24 python lib_generator_v09.py
mpirun -np 12 python lib_generator_v09.py
mpirun -np 10 python lib_generator_v09.py
mpirun -np 8 python lib_generator_v09.py
mpirun -np 6 python lib_generator_v09.py
mpirun -np 4 python lib_generator_v09.py
mpirun -np 2 python lib_generator_v09.py
mpirun -np 1 python lib_generator_v09.py


echo "All Done!"

echo "End Time = "`date`
toc=`date +%s`

elapsedTime=`expr $toc - $tic`
echo "Elapsed Time = $elapsedTime seconds"
