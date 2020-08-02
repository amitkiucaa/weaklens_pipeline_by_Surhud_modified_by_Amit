#!/bin/bash
#PBS -N job_field_random
#PBS -q intermediate
#PBS -l nodes=1:ppn=6
#PBS -j oe
#PBS -l walltime=120:00:00
cd $PBS_O_WORKDIR
module load anaconda2 mpich-3.3.1 gcc-8.2.0

mpirun -np 320 python mpi_weaklens_aroundsource.py --config rr3.config --random True --start $start


#to run a job  qsub  qsub.sh -v start=0


