#!/bin/bash -l
#
#SBATCH --ntasks 24 #number cores
#SBATCH --mem=98G #memory
#SBATCH -p intel
#SBATCH -D /rhome/cassande/bigdata/eisenlab/sg_phor/
#SBATCH -o logs/09_summary.log
#SBATCH -e logs/09_summary.log
#SBATCH -J sgphor_summary

COVDIR=coverage
CPU=24
MIN=1000
DIR=data/filtered
SAMPFILE=data/raw_data/File_names.txt
PREFIX=PHOR
ASSEMDIR=Phor_CoAssem
PROFDIR=anvio_profiles_1000minlen


source activate anvio-7.1

#anvi-summarize -p ${PROFDIR}/$PREFIX'_SAMPLES_MERGED_1000min'/PROFILE.db -c $PREFIX.db -o ${PROFDIR}/$PREFIX'_SAMPLES_MERGED'/sample_summary_DASTOOL -C DASTOOL

#anvi-summarize -p ${PROFDIR}/$PREFIX'_SAMPLES_MERGED_1000min'/PROFILE.db -c $PREFIX.db -o ${PROFDIR}/$PREFIX'_SAMPLES_MERGED'/sample_summary_METABAT -C METABAT

#anvi-summarize -p ${PROFDIR}/$PREFIX'_SAMPLES_MERGED_1000min'/PROFILE.db -c $PREFIX.db -o ${PROFDIR}/$PREFIX'_SAMPLES_MERGED'/sample_summary_CONCOCT -C CONCOCT

#anvi-summarize -p ${PROFDIR}/$PREFIX'_SAMPLES_MERGED_1000min'/PROFILE.db -c $PREFIX.db -o ${PROFDIR}/$PREFIX'_SAMPLES_MERGED'/sample_summary_MAXBIN -C MAXBIN

#anvi-summarize -p ${PROFDIR}/$PREFIX'_SAMPLES_MERGED_1000min'/PROFILE.db -c $PREFIX.db -o ${PROFDIR}/$PREFIX'_SAMPLES_MERGED'/sample_summary_BINSANITY -C BINSANITY

anvi-summarize -p ${PROFDIR}/$PREFIX'_SAMPLES_MERGED_1000min'/PROFILE.db -c $PREFIX.db -o ${PROFDIR}/$PREFIX'_SAMPLES_MERGED'/sample_summary_DASTOOL_manual_curation -C DASTOOL


