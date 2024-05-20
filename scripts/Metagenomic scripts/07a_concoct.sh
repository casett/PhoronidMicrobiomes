#!/bin/bash -l
#
#SBATCH --ntasks 24 #number cores
#SBATCH --mem=250G #memory
#SBATCH -p intel
#SBATCH -D /rhome/cassande/bigdata/eisenlab/sg_phor/
#SBATCH -o logs/07a_concoct.log
#SBATCH -e logs/07a_concoct.log
#SBATCH -J sgphor_concoct


COVDIR=coverage
CPU=24
MIN=1000
DIR=data/filtered
SAMPFILE=data/raw_data/File_names.txt
PREFIX=PHOR
ASSEMDIR=Phor_CoAssem
PROFDIR=anvio_profiles_1000minlen

source activate anvio-7.1
#must also install concoct to conda anvio env

anvi-cluster-contigs -p ${PROFDIR}/$PREFIX'_SAMPLES_MERGED_1000min'/PROFILE.db -c $PREFIX.db -C CONCOCT --driver concoct -T $CPU --just-do-it

