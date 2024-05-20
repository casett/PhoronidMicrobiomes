#!/bin/bash -l
#
#SBATCH --ntasks 16 #number cores
#SBATCH --mem=250G #memory
#SBATCH -p intel
#SBATCH -o logs/04_kaiju.log
#SBATCH -e logs/04_kaiju.log
#SBATCH -J sgphor_kaiju_anvio
#SBATCH -D /rhome/cassande/bigdata/eisenlab/sg_phor/


DB=/rhome/cassande/bigdata/software/databases/kaiju
CPU=16
DIR=data/filtered
SAMPFILE=raw_data/File_names.txt
PREFIX=PHOR
ASSEMDIR=Phor_CoAssem


module load kaiju
source activate anvio-7.1


OUT=kaiju_output
mkdir $OUT

kaiju -z $CPU -t $DB/nodes.dmp -f $DB/kaiju_db_nr_euk.fmi -i  $PREFIX.gene.calls.fa -o ${OUT}/$PREFIX.kaiju.out -v

kaiju-addTaxonNames -t $DB/nodes.dmp -n $DB/names.dmp -i ${OUT}/$PREFIX.kaiju.out -o ${OUT}/$PREFIX.kaiju.names.out -r superkingdom,phylum,class,order,family,genus,species

anvi-import-taxonomy-for-genes -i ${OUT}/$PREFIX.kaiju.names.out -c $PREFIX.db -p kaiju --just-do-it
	



