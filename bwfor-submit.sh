#!/bin/bash
#MSUB -l nodes=2:ppn=16:gpu
#MSUB -l walltime=0:00:30:00
#MSUB -v EXECUTABLE=nbody6++.avx.gpu.mpi
#MSUB -N NB_test

submitdir=$MOAB_SUBMITDIR
cd $submitdir

module load compiler/gnu/5.2
module load mpi/openmpi/1.10-gnu-5.2
module load devel/cuda/8.0
export OMP_NUM_THREADS=$((${MOAB_PROCCOUNT}/${MOAB_NODECOUNT}))
# Use when loading OpenMPI in version 1.8.x
mpirun --bind-to core --map-by core -report-bindings $EXECUTABLE < bwfor.inp 1> test.dat 2> test.err
