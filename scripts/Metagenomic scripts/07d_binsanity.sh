#!/bin/bash -l
#
#SBATCH --ntasks 24 #number cores
#SBATCH --mem=100G #memory
#SBATCH -p intel
#SBATCH -D /rhome/cassande/bigdata/eisenlab/sg_phor/
#SBATCH -o logs/07d_binsanity.log
#SBATCH -e logs/07d_binsanity.log
#SBATCH -J sgphor_binsanity


COVDIR=coverage
CPU=24
MIN=1000
DIR=data/filtered
SAMPFILE=data/raw_data/File_names.txt
PREFIX=PHOR
ASSEMDIR=Phor_CoAssem
PROFDIR=anvio_profiles_1000minlen

source activate anvio-7.1
#must also install binsanity to conda anvio env

anvi-cluster-contigs -p ${PROFDIR}/$PREFIX'_SAMPLES_MERGED_1000min'/PROFILE.db -c $PREFIX.db -C BINSANITY --driver binsanity -T $CPU --just-do-it

