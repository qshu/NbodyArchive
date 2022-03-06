#!/bin/bash -x

#SBATCH --account=chhd28 

#SBATCH --job-name=albr1.1m

###SBATCH --partition=gpus
#SBATCH --partition=develgpus

#SBATCH --gres=gpu:1
###SBATCH --gres=gpu:1

###Number of nodes requested

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
###SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12

#SBATCH --output=job-%j.out
#SBATCH --error=job-%j.err

#SBATCH --mail-user=spurzem@ari.uni-heidelberg.de
#SBATCH --mail-type=END

###SBATCH --time=24:00:00
#SBATCH --time=00:30:00

module list

###By default sets to OpenMP threads to cpus-per-task; Change or remove
###this to fit your needs, has to match with ntasks-per-node and 
###cpus-per-task omp_threads != cpus-per-task, cpus-per-task * ntasks-per-node <= 16

if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  omp_threads=$SLURM_CPUS_PER_TASK
else
  omp_threads=1
fi

export OMP_NUM_THREADS=$omp_threads

echo $OMP_NUM_THREADS

mpiprocs=$(($SLURM_NTASKS_PER_NODE * $SLURM_NNODES))

echo " OMP_NUM_THREADS=" $OMP_NUM_THREADS " mpiprocs=" $mpiprocs

#EXEC=./nbody6++.sse.gpu.mpi.b1m
EXEC=./nbody6++.avx.gpu.mpi.b1m

srun -n $mpiprocs $EXEC < NGC3201_Sana.inp 1> NGC3201_Sana.$SLURM_JOB_ID.out 2> NGC3201_Sana.$SLURM_JOB_ID.err

exit
