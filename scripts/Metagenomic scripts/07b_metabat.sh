#!/bin/bash -l
#
#SBATCH --ntasks 24 #number cores
#SBATCH --mem=100G #memory
#SBATCH -p intel
#SBATCH -D /rhome/cassande/bigdata/eisenlab/sg_phor/
#SBATCH -o logs/07b_metabat.log
#SBATCH -e logs/07b_metabat.log
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
#must also install metabat2 to conda anvio env

anvi-cluster-contigs -p ${PROFDIR}/$PREFIX'_SAMPLES_MERGED_1000min'/PROFILE.db -c $PREFIX.db -C METABAT --driver metabat2 -T $CPU --just-do-it

