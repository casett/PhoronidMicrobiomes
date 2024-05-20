#!/bin/bash -l
#
#SBATCH -e logs/03_metabolic.log
#SBATCH -o logs/03_metabolic.log
#SBATCH --nodes 1 --ntasks 24 -p intel
#SBATCH -J phor_metabolic --mem 100G
#SBATCH -D /rhome/cassande/bigdata/eisenlab/sg_phor


conda activate /rhome/cassande/bigdata/.conda/envs/METABOLIC_v4.0

ASSEM=anvio_profiles_1000minlen/PHOR_SAMPLES_MERGED/sample_summary_DASTOOL_manual_curation/all_bins_metabolic
CPU=24
OUT=metabolic_input/metabolic_output_mags
META=metabolic_input/METABOLIC

mkdir $OUT

perl $META/METABOLIC-G.pl -t $CPU -in-gn $ASSEM -p meta -o $OUT
