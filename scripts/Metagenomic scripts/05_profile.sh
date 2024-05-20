#!/bin/bash -l
#
#SBATCH --ntasks 16 #number cores
#SBATCH --mem=250G #memory
#SBATCH -p intel
#SBATCH -D /rhome/cassande/bigdata/eisenlab/sg_phor/
#SBATCH -o logs/05_profile.log
#SBATCH -e logs/05_profile.log
#SBATCH -J phor_profile_anvio


COVDIR=coverage
CPU=16
MIN=1000
DIR=data/filtered
SAMPFILE=data/raw_data/File_names.txt
PREFIX=PHOR
ASSEMDIR=Phor_CoAssem


OUT=anvio_profiles_1000minlen

mkdir $OUT

source activate anvio-7.1

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read SAMPLE
do 
	anvi-profile -i ${COVDIR}/$SAMPLE'.bam' -c $PREFIX.db --num-threads $CPU --min-contig-length $MIN --cluster-contigs --output-dir ${OUT}/$SAMPLE'_profile'
done	


