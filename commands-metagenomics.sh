# Run FastQC on input files IM_R1.fq.gz and IM_R2.fq.gz with 16 threads,
# and save the results in the 'fastQCin' directory.
fastqc -t 16 -o fastQCin IM_R1.fq.gz IM_R2.fq.gz

# Use cutadapt to trim adapter sequences from the input files, and save the
# trimmed reads in IM_R1_cutadapt.fq.gz and IM_R2_cutadapt.fq.gz.
cutadapt -j 16 -g AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
-a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG \
-G AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
-A GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG \
-o IM_R1_cutadapt.fq.gz -p IM_R2_cutadapt.fq.gz IM_R1_cutadapt.fq.gz IM_R2_cutadapt.fq.gz

# Use bbduk.sh to further trim and filter the reads, and save them in
# IM_R1_bbduk.fq.gz and IM_R2_bbduk.fq.gz.
bbduk.sh threads=16 in=IM_R1_cutadapt.fq.gz in2=IM_R2_cutadapt.fq.gz \
out=IM_R1_bbduk.fq.gz out2=IM_R2_bbduk.fq.gz qtrim=lr trimq=20 maq=20 minlen=75

# Align the trimmed reads to a reference genome using Bowtie2 with 16 threads,
# and save the aligned and unaligned reads in IM_mapped_and_unmapped.sam.
bowtie2 -p 16 -x /media/darwin/Externo/DanielTS/BBDD/GRCh38_noalt_as/GRCh38_noalt_as \
-1 IM_R1_bbduk.fq.gz -2 IM_R2_bbduk.fq.gz --un-conc-gz filtered_ > IM_mapped_and_unmapped.sam

# Rename the filtered output files from Bowtie2.
mv filtered_IM_R1_bbduk.fq.gz.1 IM_filtered_R1.fq.gz
mv filtered_IM_R2_bbduk.fq.gz.2 IM_filtered_R2.fq.gz

# Remove the intermediate SAM file.
rm -r IM_mapped_and_unmapped.sam

# Run FastQC on the filtered reads IM_filtered_R1.fq.gz and IM_filtered_R2.fq.gz
# with 16 threads, and save the results in the 'fastQCour' directory.
fastqc -t 16 -o fastQCout IM_filtered_R1.fq.gz IM_filtered_R2.fq.gz

# Assemble the filtered reads using Megahit with 16 threads and save the output
# in the 'megahit/' directory with the prefix 'a_'.
megahit -t 16 -1 IM_filtered_R1.fq.gz -2 IM_filtered_R2.fq.gz -o megahit/ --out-prefix a_

# Run Quast to assess the quality of the assembly and save the results in the 'QUAST/' directory
# using 16 threads and the contigs from the Megahit assembly.
quast.py -o QUAST/ -t 16 megahit/IM.contigs.fa

# Create directories for Bowtie2 output.
mkdir bw2
mkdir bw2/bw2_index

# Build Bowtie2 index for the Megahit assembly.
bowtie2-build megahit/IM.contigs.fa --threads 16 bw2/bw2_index/IM_assembly_index

# Align the filtered reads to the Bowtie2 index and save the output as a SAM file.
bowtie2 --threads 16 -x  bw2/bw2_index/IM_assembly_index-1 IM_filtered_R1.fq.gz -2 IM_filtered_R2.fq.gz -S bw2/bw2_index/IM_assembly_coverage.sam

# Convert SAM to BAM format.
samtools view -b bw2/bw2_index/IM_assembly_coverage.sam --threads 16 > bw2/bw2_index/IM_assembly_coverage.bam

# Sort the BAM file.
samtools sort --threads 16 bw2/bw2_index/IM_assembly_coverage.bam > bw2/bw2_index/IM_assembly_coverage_sorted.bam

# Generate depth information for the assembly.
jgi_summarize_bam_contig_depths --outputDepth bw2/depth.txt bw2/bw2_index/IM_assembly_coverage_sorted.bam

# Remove intermediate SAM and BAM files.
rm bw2/bw2_index/IM_assembly_coverage.sam
rm bw2/bw2_index/IM_assembly_coverage.bam
rm bw2/bw2_index/IM_assembly_coverage_sorted.bam

# Create directories for Metabat2 and MaxBin2 output.
mkdir metabat2_depth
mkdir maxbin2_depth

# Run Metabat2 to bin contigs.
metabat2 -i megahit/IM.contigs.fa -o metabat2_depth/bin -t 16 -a bw2/depth.txt

# Prepare depth information for MaxBin2.
cut -f1,3 bw2/depth.txt | tail -n+2 > bw2/depth_maxbin.txt

# Run MaxBin2 to bin contigs.
run_MaxBin.pl -contig megahit/IM.contigs.fa -out maxbin2_depth/bin -abund bw2/depth_maxbin.txt

# Create a directory for DAS_Tool output.
mkdir dastool

# Use DAS_Tool to integrate Metabat2 and MaxBin2 binning results.
Fasta_to_Scaffolds2Bin.sh -i metabat2_depth -e fa > dastool/metabat2.scaffolds2bin.tsv
Fasta_to_Scaffolds2Bin.sh -i maxbin2_depth -e fasta > dastool/maxbin2.scaffolds2bin.tsv
DAS_Tool -i dastool/metabat2.scaffolds2bin.tsv,dastool/maxbin2.scaffolds2bin.tsv -l metabat,maxbin -c megahit/IM.contigs.fa -o dastool/ --search_engine diamond --threads 16 --write_bins

# Create directories for CheckM output.
mkdir checkm_metabat2_depth
mkdir checkm_maxbin2_depth
mkdir dastool/checkm

# Run CheckM to assess bin quality.
checkm lineage_wf -t 16 -f checkm_metabat2_depth/quality_output.txt -x fa metabat2_depth
checkm lineage_wf -t 16 -f checkm_maxbin2_depth/quality_output.txt -x fasta maxbin2_depth
checkm lineage_wf -t 16 -f dastool/quality_dastool_output.txt -x fa dastool/*"DASTool_bins" dastool/checkm


# Calculate genome coverage using 'coverm' tool.
coverm genome --coupled IM_filtered_R1.fq.gz IM_filtered_R2.fq.gz --genome-fasta-files dastool/*fa -o output.tsv

# Classify genomes with 'gtdbtk'.
gtdbtk classify_wf --genome_dir dastool/ --out_dir MAGS_classification/ -x fa --cpus 16


# Annotate genomes with 'prokka'.
for MAG in dastool/*fa
	do
    		prokka --prefix ${MAG} MAG
    	done

# Run the Prodigal tool to predict protein-coding genes and save the output as IM_protein_sequences.faa
prodigal -i megahit/IM.contigs.fa -a IM_protein_sequences.faa

# Iterate over a list of functions: EPS, fimbriae, oligotrophic, and psychrotrophic.
for function in EPS fimbriae oligotrophic psychrotrophic
	do
		# Create a BLAST database from the protein sequences file for the current function.
		makeblastdb -in ${function}.faa -dbtype prot -parse_seqids -out ${function}_database -title "${function}_database"
		
		# Perform a BLAST search using the created database and save the results in a log file.
		blastp -db "${function}_database" -out ./${function}_out_blast.log -outfmt '6 qacc qseqid stitle staxids sscinames sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore' -query IM_protein_sequences.faa -max_target_seqs 1 -evalue 1e-5 -num_threads 16
		
		# Filter the BLAST results based on a bitscore threshold of 50 and save the filtered results to a new log file.
		filter_by_bitscore.py --value 50 ${function}_out_blast.log -o ${function}_out_blast_filtered.log
	done

