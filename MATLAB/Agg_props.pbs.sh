#!/bin/bash
#PBS -l walltime=12:00:00,select=1:ncpus=32:mem=124gb
#PBS -N Agg_props

module load matlab/R2020a
cd $PBS_O_WORKDIR
matlab -batch Agg_properties_v1