#!/bin/sh --login
#
# -- Request 24 cores
#PBS -l procs=8
#
# -- Specifya maximum wallclock of 4 hours
#PBS -l walltime=0:10:00
#
# -- Specify under which account a job should run
#PBS -A gsd-fv3 
#
# -- Set the name of the job, or moab will default to STDIN
##PBS -N 

# change directory to the working directory of the job
# Use the if clause so that this script stays portable
#

if [ x$PBS_O_WORKDIR != x ]; then
   cd $PBS_O_WORKDIR
fi

np=$PBS_NP

# run Fortran version
module load intel/14.0.2 impi/5.1.2.150
# add esmf if needed
# module load esmf/7.0.0

mpirun -np $np ./src/chemdrv

#### OR ####
# run NUOPC version 0.3
# module use /scratch4/BMC/esmf/Raffaele.Montuoro/dev/sw/modulefiles
# module load intel/14.0.2 impi/5.1.2.150 netcdf/4.3.0 esmf/7.1.0.dev
#
# mpirun -np $np ./src/nuopcdrv
