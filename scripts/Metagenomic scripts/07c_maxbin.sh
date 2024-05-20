#!/bin/bash -l
#
#SBATCH --ntasks 24 #number cores
#SBATCH --mem=400G #memory
#SBATCH -p intel
#SBATCH -D /rhome/cassande/bigdata/eisenlab/sg_phor/
#SBATCH -o logs/07c_maxbin.log
#SBATCH -e logs/07c_maxbin.log
#SBATCH -J sgphor_metabat


COVDIR=coverage
CPU=24
MIN=1000
DIR=data/filtered
SAMPFILE=data/raw_data/File_names.txt
PREFIX=PHOR
ASSEMDIR=Phor_CoAssem
PROFDIR=anvio_profiles_1000minlen

source activate anvio-7.1
#must also install maxbin to conda anvio env

anvi-cluster-contigs -p ${PROFDIR}/$PREFIX'_SAMPLES_MERGED_1000min'/PROFILE.db -c $PREFIX.db -C MAXBIN --driver maxbin2 -T $CPU --just-do-it

