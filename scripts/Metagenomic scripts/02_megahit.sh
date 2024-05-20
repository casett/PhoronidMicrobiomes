#!/bin/bash -l
#
#SBATCH -e logs/02_megahit.log
#SBATCH -o logs/02_megahit.log
#SBATCH --nodes 1 --ntasks 24 -p highmem
#SBATCH -J phor_megahit --mem 600G
#SBATCH -D /rhome/cassande/bigdata/eisenlab/sg_phor

module load megahit
CPU=24
DIR=data/filtered
OUT=Phor_CoAssem

megahit -1 $DIR/PHOR-511_S2/PHOR-511_S2_1_filtered.fastq,$DIR/PHOR-512_S3/PHOR-512_S3_1_filtered.fastq -2 $DIR/PHOR-511_S2/PHOR-511_S2_2_filtered.fastq,$DIR/PHOR-512_S3/PHOR-512_S3_2_filtered.fastq -o $OUT -t $CPU

pigz $DIR/*/*.fastq
