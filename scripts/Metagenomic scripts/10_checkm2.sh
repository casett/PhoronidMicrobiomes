#!/bin/bash -l
#
#SBATCH -n 24 #number cores
#SBATCH -p intel
#SBATCH -D /rhome/cassande/bigdata/eisenlab/sg_phor/
#SBATCH -o logs/10_checkm2.log
#SBATCH -e logs/10_checkm2.log
#SBATCH -J sgphor_checkM
#SBATCH --mem 96G #memory in Gb
#SBATCH -t 96:00:00 #time in hours:min:sec


conda activate checkm2

BINFOLDER=anvio_profiles_1000minlen/PHOR_SAMPLES_MERGED/sample_summary_DASTOOL_manual_curation/all_Bins
OUTPUT=DASTOOL_checkM2
CPU=24

#checkm2 database --download

checkm2 predict --threads $CPU --input $BINFOLDER --output-directory $OUTPUT -x .fa
