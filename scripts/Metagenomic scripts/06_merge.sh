#!/bin/bash -l
#
#SBATCH --ntasks 24 #number cores
#SBATCH --mem=350G #memory
#SBATCH -p intel
#SBATCH -D /rhome/cassande/bigdata/eisenlab/sg_phor/
#SBATCH -o logs/06_merge.log
#SBATCH -e logs/06_merge.log
#SBATCH -J phor_merge_anvio


COVDIR=coverage
CPU=24
MIN=1000
DIR=data/filtered
SAMPFILE=data/raw_data/File_names.txt
PREFIX=PHOR
ASSEMDIR=Phor_CoAssem


PROFDIR=anvio_profiles_1000minlen

source activate anvio-7.1

anvi-merge ${PROFDIR}/*profile/PROFILE.db -o ${PROFDIR}/$PREFIX'_SAMPLES_MERGED_1000min' -c $PREFIX.db

