#!/bin/bash
#PBS -N job_field_random
#PBS -q intermediate
#PBS -l nodes=1:ppn=32
#PBS -j oe
#PBS -l walltime=01:00:00
 
cd $PBS_O_WORKDIR
module load anaconda2 mpich-3.3.1 gcc-8.2.0
python weaklens_aroundsource.py --config arindam_lenses.config
