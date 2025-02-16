# preprocessing.sh
# Single-Cell RNA-seq Quality Control and Preprocessing Pipeline

# Step 1: Quality Control on FASTQ Files
mkdir -p results/fastqc_reports
echo "Running FastQC..."
fastqc data/fastq/*.fastq.gz -o results/fastqc_reports/

# Step 2: Summarize QC Reports
mkdir -p results/multiqc_reports
echo "Generating MultiQC report..."
multiqc results/fastqc_reports/ -o results/multiqc_reports/

# Step 3: Adapter Trimming
mkdir -p results/trimmed_fastq
echo "Running Trim Galore..."
trim_galore --paired data/fastq/*.fastq.gz -o results/trimmed_fastq/

# Step 4: Alignment with STAR
mkdir -p results/aligned/
echo "Aligning with STAR..."
STAR --runThreadN 8 \
     --genomeDir reference/star_index/ \
     --readFilesIn results/trimmed_fastq/sample_R1.fastq.gz results/trimmed_fastq/sample_R2.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix results/aligned/sample_ \
     --outSAMtype BAM SortedByCoordinate

# Step 5: Process 10X Data with Cell Ranger
mkdir -p results/cellranger/
echo "Processing 10X data with Cell Ranger..."
cellranger count --id=sample \
    --transcriptome=reference/cellranger_index/ \
    --fastqs=data/fastq/ \
    --sample=sample \
    --localcores=8 --localmem=32

# Step 6: Deduplication and UMI Correction
mkdir -p results/umi_corrected/
echo "Correcting UMIs..."
umi_tools dedup -I results/aligned/sample_Aligned.sortedByCoord.out.bam --output-stats=results/umi_corrected/sample_umi_stats.tsv \
    -S results/umi_corrected/sample_dedup.bam
