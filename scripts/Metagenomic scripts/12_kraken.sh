#!/bin/bash -l
#
#SBATCH --ntasks 32 #number cores
#SBATCH -J sgphor_kraken
#SBATCH --mem=480G #memory
#SBATCH -p batch
#SBATCH -D /rhome/cassande/bigdata/eisenlab/sg_phor/
#SBATCH -o logs/01_krak.log
#SBATCH -e logs/01_krak.log

#module load bracken

conda activate /rhome/cassande/bigdata/.conda/envs/kraken2

CPU=32

export KRAKEN2_DB_PATH="/rhome/cassande/bigdata/eisenlab/sg_phor/"

#kraken2-build --download-taxonomy --db nt
#kraken2-build --download-library nt --db nt
#kraken2-build --build --threads $CPU --db nt --fast-build --max-db-size 475000000000
#bracken-build -d nt -t $CPU -k 35 -l 250

#echo "built"

for prefix in $(ls data/filtered);
do
	echo $prefix
	read1=( data/filtered/${prefix}/${prefix}_1_filtered.fastq.gz ) #the parentheses assign the globbed filename to an array (of length 1)
	read2=( data/filtered/${prefix}/${prefix}_2_filtered.fastq.gz )

	kraken2 --db nt --threads $CPU --report ${prefix}_kraken2_report_paired.tsv --gzip-compressed --paired ${read1} ${read2} > ${prefix}.log

	bracken -d nt -i ${prefix}_kraken2_report_paired.tsv -o ${prefix}.S.bracken -r 250 -l 'S' -t 10
	bracken -d nt -i ${prefix}_kraken2_report_paired.tsv -o ${prefix}.F.bracken -r 250 -l 'F' -t 10
	bracken -d nt -i ${prefix}_kraken2_report_paired.tsv -o ${prefix}.G.bracken -r 250 -l 'G' -t 10
	bracken -d nt -i ${prefix}_kraken2_report_paired.tsv -o ${prefix}.K.bracken -r 250 -l 'K' -t 10
	bracken -d nt -i ${prefix}_kraken2_report_paired.tsv -o ${prefix}.P.bracken -r 250 -l 'P' -t 10
	bracken -d nt -i ${prefix}_kraken2_report_paired.tsv -o ${prefix}.C.bracken -r 250 -l 'C' -t 10
	bracken -d nt -i ${prefix}_kraken2_report_paired.tsv -o ${prefix}.O.bracken -r 250 -l 'O' -t 10
done
