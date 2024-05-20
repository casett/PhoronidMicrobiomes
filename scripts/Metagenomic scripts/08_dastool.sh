#!/bin/bash -l
#
#SBATCH --ntasks 24 #number cores
#SBATCH --mem=98G #memory
#SBATCH -p intel
#SBATCH -D /rhome/cassande/bigdata/eisenlab/sg_phor/
#SBATCH -o logs/08_dastool.log
#SBATCH -e logs/08_dastool.log
#SBATCH -J sgphor_das_anvio


COVDIR=coverage
CPU=24
MIN=1000
DIR=data/filtered
SAMPFILE=data/raw_data/File_names.txt
PREFIX=PHOR
ASSEMDIR=Phor_CoAssem
PROFDIR=anvio_profiles_1000minlen


source activate anvio-7.1
#made sure to install dastool via conda to anvio-7 env


anvi-cluster-contigs -p ${PROFDIR}/$PREFIX'_SAMPLES_MERGED_1000min'/PROFILE.db -c $PREFIX.db  -C DASTOOL --driver dastool -S CONCOCT,METABAT,MAXBIN,BINSANITY -T $CPU --search-engine diamond --just-do-it --score-threshold 0
