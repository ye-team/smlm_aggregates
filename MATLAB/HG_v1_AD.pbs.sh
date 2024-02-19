#!/bin/bash
#PBS -l walltime=12:00:00,select=1:ncpus=32:mem=124gb
#PBS -N HG_Matlab

module load matlab/R2020a
cd $PBS_O_WORKDIR
matlab -batch new_matlab_script_v3_AD
