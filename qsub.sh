#!/bin/bash
#PBS -N job_field_random
#PBS -q intermediate
#PBS -l nodes=4:ppn=8
#PBS -j oe

cd $PBS_O_WORKDIR
module load anaconda2 mpich-3.3.1 gcc-8.2.0

mpirun -np 32 python mpi_weaklens_aroundsource.py --config redmapper_randoms_radnom_rotation2.config --random True --start $start


#to run a job  qsub  qsub.sh -v start=0


