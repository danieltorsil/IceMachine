#!/bin/bash

# Activate the conda environment named "qiime2-2021.2"
conda activate qiime2-2021.2

# Run qiime2input.py script to create a Manifest file for input
qiime2input.py -i /media/darwin/Proyectos/2022/MP_02_2019/ITS/QIIME2/fastq/ > Manifest

# Import paired-end sequence data using the Manifest file
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path Manifest --output-path ./demux-paired-end.qza --input-format PairedEndFastqManifestPhred33

# Summarize the imported data
qiime demux summarize --i-data ./demux-paired-end.qza --o-visualization ./demux-paired-end.qzv

# Perform denoising and quality filtering using DADA2
qiime dada2 denoise-paired \
--i-demultiplexed-seqs ./demux-paired-end.qza \
--verbose \
--p-trim-left-f 17 \
--p-trim-left-r 22 \
--p-trunc-len-f 280 \
--p-trunc-len-r 246 \
--p-n-threads 16 \
--o-table ./table.qza \
--o-representative-sequences ./rep-seqs.qza \
--o-denoising-stats ./denoising-stats.qza

# Tabulate denoising statistics
qiime metadata tabulate --m-input-file ./denoising-stats.qza --o-visualization ./stats-dada2.qzv

# Classify sequences against a reference database (UNITE classifier)
qiime feature-classifier classify-sklearn   --i-classifier /home/darwin/Escritorio/DaniTS/BBDD/unite-ver8-99-classifier-04.02.2020_entrenado.qza   --i-reads ./rep-seqs.qza   --verbose   --p-n-jobs 12   --o-classification ./taxonomy.qza

# Export the feature table
qiime tools export --input-path ./table.qza --output-path ./table

# List files in the "table" directory
ls table

# Export the taxonomy data
qiime tools export --input-path ./taxonomy.qza --output-path ./taxonomy

# List files in the "taxonomy" directory
ls taxonomy

# Modify the column names in the taxonomy file
sed -i -e 's/Feature ID/#OTUID/g' ./taxonomy/taxonomy.tsv
sed -i -e 's/Taxon/taxonomy/g' ./taxonomy/taxonomy.tsv
sed -i -e 's/Confidence/confidence/g' ./taxonomy/taxonomy.tsv

# Display the first few lines of the modified taxonomy file
head taxonomy/taxonomy.tsv

# Add taxonomy metadata to the feature table
biom add-metadata -i ./table/feature-table.biom -o ./table-with-taxonomy.biom --observation-metadata-fp ./taxonomy/taxonomy.tsv --sc-separated taxonomy

# Convert the feature table to JSON format
biom convert -i ./table-with-taxonomy.biom -o ./table-with-taxonomy-json2.biom --table-type="OTU table" --to-json
