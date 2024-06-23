#!/bin/bash

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path manifest.tsv \
  --input-format SingleEndFastqManifestPhred33V2 \
  --output-path demultiplexed-seqs.qza
qiime demux summarize \
  --i-data demultiplexed-seqs.qza \
  --o-visualization demultiplexed-seqs.qzv
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demultiplexed-seqs.qza \
  --p-trim-left 20 \
  --p-trunc-len 150 \
  --p-n-threads 0 \
  --o-representative-sequences rep-seqs.qza \
  --o-table table.qza \
  --o-denoising-stats stats-dada2.qza
qiime metadata tabulate \
  --m-input-file stats-dada2.qza \
  --o-visualization stats-dada2.qzv
qiime feature-table summarize \
  --i-table table.qza \
  --m-sample-metadata-file hillburns_metadata_Clean.tsv \
  --o-visualization table.qzv
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
qiime feature-table rarefy \
  --i-table table.qza \
  --p-sampling-depth 1412 \
  --p-with-replacement \
  --o-rarefied-table rarefied_table.qza
qiime feature-table relative-frequency \
  --i-table table.qza \
  --o-relative-frequency-table rel-feature-table.qza
mkdir metrics-results
for METRIC in observed_features chao1 shannon simpson
do
  qiime diversity alpha \
    --i-table rarefied_table.qza \
    --p-metric ${METRIC} \
    --o-alpha-diversity metrics-results/${METRIC}_vector.qza
done
mkdir alpha-diversity-visualization
for METRIC in observed_features chao1 shannon simpson
do
  qiime diversity alpha-group-significance \
    --i-alpha-diversity metrics-results/${METRIC}_vector.qza \
    --m-metadata-file hillburns_metadata_Clean.tsv \
    --o-visualization alpha-diversity-visualization/${METRIC}-group-significance.qzv
done
for METRIC in unweighted_unifrac weighted_unifrac
do
  qiime diversity beta-phylogenetic \
    --i-table rarefied_table.qza \
    --i-phylogeny rooted-tree.qza \
    --p-metric ${METRIC} \
    --o-distance-matrix metrics-results/${METRIC}_distance_matrix.qza
done
for METRIC in canberra
do
  qiime diversity beta \
    --i-table rarefied_table.qza \
    --p-metric ${METRIC} \
    --o-distance-matrix metrics-results/${METRIC}_distance_matrix.qza
done
mkdir beta-diversity-visualization
 for METRIC in unweighted_unifrac weighted_unifrac canberra
do
  qiime diversity beta-group-significance \
    --i-distance-matrix metrics-results/${METRIC}_distance_matrix.qza \
    --m-metadata-file hillburns_metadata_Clean.tsv \
    --m-metadata-column alcohol_how_much \
    --p-method permanova \
    --p-pairwise \
    --o-visualization beta-diversity-visualization/${METRIC}-disease-status-significance.qzv
done
 for METRIC in unweighted_unifrac weighted_unifrac canberra
do
  qiime diversity pcoa \
    --i-distance-matrix metrics-results/${METRIC}_distance_matrix.qza \
    --o-pcoa  metrics-results/${METRIC}_pcoa_matrix.qza
done
 for METRIC in unweighted_unifrac weighted_unifrac canberra
do
  qiime emperor plot \
    --i-pcoa metrics-results/${METRIC}_pcoa_matrix.qza \
    --m-metadata-file hillburns_metadata_Clean.tsv \
    --o-visualization beta-diversity-visualization/${METRIC}_emperor.qzv
done
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-min-depth 1 \
  --p-max-depth 45500 \
  --p-steps 50 \
  --p-metrics observed_features chao1 shannon simpson \
  --m-metadata-file hillburns_metadata_Clean.tsv \
  --o-visualization alpha-rarefaction.qzv
wget 'https://data.qiime2.org/2022.2/common/silva-138-99-seqs.qza'
wget 'https://data.qiime2.org/2022.2/common/silva-138-99-tax.qza'
qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --o-reads ref-seqs-520f-926r.qza
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs-520f-926r.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --o-classifier silva-138-99-520f-926r-classifier.qza
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-520f-926r-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file hillburns_metadata_Clean.tsv \
  --o-visualization taxa-bar-plots.qzv
mkdir phyloseq
taxon=( domain-kingdom phylum class order family genus species )
for i in {1..7}
do
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level ${i} \
  --o-collapsed-table phyloseq/${taxon[${i}-1]}-table.qza
done
for taxon in domain-kingdom phylum class order family genus species
do
  qiime feature-table relative-frequency \
    --i-table phyloseq/${taxon}-table.qza \
    --o-relative-frequency-table phyloseq/rel-${taxon}-table.qza
done
for taxon in order family genus
do
  qiime feature-table heatmap \
    --i-table phyloseq/${taxon}-table.qza \
    --o-visualization heatmap_${taxon}.qzv \
    --p-color-scheme YlOrRd \
    --p-cluster both
done
qiime tools export \
  --input-path taxonomy.qza \
  --output-path phyloseq/taxonomy
mv phyloseq/taxonomy/taxonomy.tsv phyloseq/taxonomy.tsv
rm -r phyloseq/taxonomy
qiime tools export \
  --input-path table.qza \
  --output-path phyloseq/feature-table
biom convert -i phyloseq/feature-table/feature-table.biom -o phyloseq/feature-table.tsv --to-tsv
cp phyloseq/feature-table/feature-table.biom phyloseq/feature-table.biom
rm -r phyloseq/feature-table
qiime tools export \
  --input-path rel-feature-table.qza \
  --output-path phyloseq/rel-feature-table
biom convert -i phyloseq/rel-feature-table/feature-table.biom -o phyloseq/rel-feature-table.tsv --to-tsv
rm -r phyloseq/rel-feature-table
for taxon in domain-kingdom phylum class order family genus species
do
  qiime tools export \
    --input-path phyloseq/${taxon}-table.qza \
    --output-path phyloseq/${taxon}-table
  biom convert -i phyloseq/${taxon}-table/feature-table.biom -o phyloseq/${taxon}-table.tsv --to-tsv
  rm -r phyloseq/${taxon}-table
done
for taxon in domain-kingdom phylum class order family genus species
do
  qiime tools export \
    --input-path phyloseq/rel-${taxon}-table.qza \
    --output-path phyloseq/rel-${taxon}-table
  biom convert -i phyloseq/rel-${taxon}-table/feature-table.biom -o phyloseq/rel-${taxon}-table.tsv --to-tsv
  rm -r phyloseq/rel-${taxon}-table
done
qiime tools export \
 --input-path rep-seqs.qza \
 --output-path phyloseq/rep-seqs
mv phyloseq/rep-seqs/dna-sequences.fasta phyloseq/dna-sequences.fna
rm -r phyloseq/rep-seqs
qiime tools export \
  --input-path unrooted-tree.qza \
  --output-path phyloseq
mv phyloseq/tree.nwk phyloseq/unrooted_tree.nwk
qiime tools export \
  --input-path rooted-tree.qza \
  --output-path phyloseq
mv phyloseq/tree.nwk phyloseq/rooted_tree.nwk
cp hillburns_metadata_Clean.tsv phyloseq/hillburns_metadata_Clean.tsv
cp rep-seqs.qza phyloseq/rep-seqs.qza
cp table.qza phyloseq/feature-table.qza
cp rel-feature-table.qza phyloseq/rel-feature-table.qza
cp unrooted-tree.qza phyloseq/unrooted-tree.qza
cp rooted-tree.qza phyloseq/rooted-tree.qza
cp taxonomy.qza phyloseq/taxonomy.qza
cp rarefied_table.qza phyloseq/rarefied_table.qza
