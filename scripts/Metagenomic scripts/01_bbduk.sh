#!/bin/bash -l
#
#SBATCH -n 16 #number cores
#SBATCH -D /rhome/cassande/bigdata/eisenlab/sg_phor/
#SBATCH -e logs/01_bbduk.log
#SBATCH -o logs/01_bbduk.log
#SBATCH -J phor_qc
#SBATCH --mem 50G #memory per node in Gb
#SBATCH -p intel

IN=data/raw_data
SAMPFILE=$IN/File_names.txt
OUT=data/filtered

mkdir $OUT

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read PREFIX 
do
	mkdir $OUT/$PREFIX

	/rhome/cassande/bigdata/software/bbmap/bbduk.sh in1=$IN/$PREFIX'_L001_R1_001_comb.fastq.gz' in2=$IN/$PREFIX'_L001_R2_001_comb.fastq.gz' out1=$OUT/$PREFIX/$PREFIX'_1_filtered.fastq' out2=$OUT/$PREFIX/$PREFIX'_2_filtered.fastq' ref=/rhome/cassande/bigdata/software/bbmap/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=10 maq=10


done

