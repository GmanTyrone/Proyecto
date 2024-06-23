#!/bin/bash

picrust2_pipeline.py \
  -s phyloseq/dna-sequences.fna \
  -i phyloseq/feature-table.biom \
  -o picrust2_output_K &&

cd ./picrust2_output_K &&
hsp.py \
  -i 16S \
  -t out.tre -o marker_predicted_and_nsti.tsv.gz \
  -p 1 \
  -n &&

hsp.py \
  -i KO -t out.tre \
  -o KO_predicted.tsv.gz \
  -p 1 &&

metagenome_pipeline.py \
  -i ../phyloseq/feature-table.biom \
  -m marker_predicted_and_nsti.tsv.gz \
  -f KO_predicted.tsv.gz \
  -o KO_metagenome_out \
  --strat_out &&

pathway_pipeline.py \
  -i KO_metagenome_out/pred_metagenome_contrib.tsv.gz \
  -o KEGG-Pathways \
  --no_regroup \
  --map /tmp/picrust2-2.4.2/picrust2/default_files/pathway_mapfiles/KEGG_pathways_to_KO.tsv &&

add_descriptions.py \
  -i KEGG-Pathways/path_abun_unstrat.tsv.gz \
  --custom_map_table /tmp/picrust2-2.4.2/picrust2/default_files/description_mapfiles/KEGG_pathways_info.tsv.gz \
  -o KEGG-Pathways/path_abun_unstrat_descrip.tsv.gz &&

add_descriptions.py \
  -i KEGG-Pathways/path_abun_unstrat.tsv.gz \
  --custom_map_table /tmp/KO_LEVEL23_2022_1214.tsv \
  -o KEGG-Pathways/L2_3_path_abun_unstrat_descrip.tsv.gz &&

gzip -d KEGG-Pathways/L2_3_path_abun_unstrat_descrip.tsv.gz
