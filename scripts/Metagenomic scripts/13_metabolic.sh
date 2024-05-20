#!/bin/bash -l
#
#SBATCH -e logs/02_metabolic.log
#SBATCH -o logs/02_metabolic.log
#SBATCH --nodes 1 --ntasks 24 -p intel
#SBATCH -J phor_metabolic --mem 100G
#SBATCH -D /rhome/cassande/bigdata/eisenlab/sg_phor


conda activate METABOLIC_v4.0

ASSEM=metabolic_input/data
CPU=24
OUT=metabolic_input/metabolic_output
META=metabolic_input/METABOLIC

perl $META/METABOLIC-G.pl -t $CPU -in-gn $ASSEM -p meta -o $OUT
