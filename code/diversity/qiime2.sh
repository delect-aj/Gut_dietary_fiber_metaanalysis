#!/bin/bash
qiime tools import \
  --input-path table_6721_2378.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path ../diversity/table.qza
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table-summary.qzv
qiime tools import \
  --input-path /beegfs/dongbiao/greengene2/2024.09.phylogeny.id.nwk \
  --type 'Phylogeny[Rooted]' \
  --output-path greengene_rooted-tree.qza
qiime phylogeny filter-tree \
  --i-tree greengene_rooted-tree.qza \
  --i-table table.qza \
  --o-filtered-tree rooted-tree.qza
#核心多样性分析
time qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 1000 \
  --m-metadata-file metadata.tsv \
  --output-dir core-metrics-results
~

~

