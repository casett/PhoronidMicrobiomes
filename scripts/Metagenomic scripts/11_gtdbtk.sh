#!/bin/bash -l
#SBATCH --ntasks=16 # Number of cores
#SBATCH --mem=400G # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -p intel
#SBATCH -o logs/11_gtdbtk.log
#SBATCH -e logs/11_gtdbtk.log
#SBATCH -J sgphor_gtdbtk
#SBATCH -D /rhome/cassande/bigdata/eisenlab/sg_phor/



conda activate gtdbtk-2.2.6

#download-db.sh


INPUT=anvio_profiles_1000minlen/PHOR_SAMPLES_MERGED/sample_summary_DASTOOL_manual_curation/all_Bins
OUTPUT=DASTOOL_gtbdk_results_v2.2.6
CPU=16
PREFIX=PHOR

GTDBTK_DATA_PATH=/rhome/cassande/bigdata/.conda/envs/gtdbtk-2.2.6/share/gtdbtk-2.2.6/db

gtdbtk classify_wf --genome_dir $INPUT --out_dir $OUTPUT -x .fa --cpus $CPU --prefix $PREFIX.gtbdk --mash_db $OUTPUT/mash_db
