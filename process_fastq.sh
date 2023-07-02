#!/bin/bash

#SBATCH -c 20                              # Request cores
#SBATCH -t 0-05:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=100G                         # Memory total (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID (%j)
                                           # You can change the filenames given with -o and -e to any filenames you'd like

# Load modules here!
# module load python/3.7.4 picard/2.27.5

# Merge fastq files across lanes for R1
if [ -s $1_merge_R1_001.fastq.gz ]; then
    echo "$1_merge_R1_001.fastq.gz already exists."
else 
    cat $1_L001_R1_001.fastq.gz \
        $1_L002_R1_001.fastq.gz \
        $1_L003_R1_001.fastq.gz \
        $1_L004_R1_001.fastq.gz > $1_merge_R1_001.fastq.gz
fi

# Merge fastq files across lanes for R2
if [ -s $1_merge_R2_001.fastq.gz ]; then
    echo "$1_merge_R2_001.fastq.gz already exists."
else 
	cat $1_L001_R2_001.fastq.gz \
        $1_L002_R2_001.fastq.gz \
        $1_L003_R2_001.fastq.gz \
        $1_L004_R2_001.fastq.gz > $1_merge_R2_001.fastq.gz
fi

# Trim fastq files for adapters and remove low-quality base-calls
if [[ -s $1_merge_R1_001.trimmed.fastq.gz && $1_merge_R2_001.trimmed.fastq.gz ]]; then
    echo "$1_merge_R1_001.trimmed.fastq.gz and $1_merge_R2_001.trimmed.fastq.gz already exist."
else 
	cutadapt \
    	-j 20 \
    	-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    	-q 30 \
    	-o $1_merge_R1_001.trimmed.fastq.gz -p $1_merge_R2_001.trimmed.fastq.gz \
    	$1_merge_R1_001.fastq.gz $1_merge_R2_001.fastq.gz
fi

# Align trimmed fastq files with bowtie2 (change -x parameter to directory containing bowtie2 index!!!)
if [ -s $1.sam ] || [ -s $1.sort.bam ] && [ -s $1_alignment_file.log ]; then
    echo "$1 alignment is already complete."
else 
	(bowtie2 \
	    --threads 20 \
        --local \
        --dovetail \
        --very-sensitive-local \
        --no-unal \
        --no-mixed \
        --no-discordant \
        --phred33 \
        -x ~/bowtie2_mm10/mm10 \
        -1 $1_merge_R1_001.trimmed.fastq.gz \
        -2 $1_merge_R2_001.trimmed.fastq.gz \
        -S $1.sam) 2>$1_alignment_file.log
fi

# Align to spike-in genome (E. coli strain K12, MG1655) (change -x parameter to directory containing bowtie2 index!!!)
if [ -s $1_spike_in.sam ]; then
    echo "$1_spike_in.sam already exists."
else
	(bowtie2 \
		--threads 20 \
		--local \
		--very-sensitive-local \
		--no-unal \
		--no-mixed \
		--no-discordant \
		--phred33 \
		--no-overlap \
		-x ~/bowtie2_MG1655/MG1655 \
		-1 $1_merge_R1_001.trimmed.fastq.gz \
		-2 $1_merge_R2_001.trimmed.fastq.gz \
		-S $1_spike_in.sam) 2>$1_spike_in_alignment_file.log
fi

# Convert sam file to bam file
if [ -s $1.bam ]; then
    echo "$1.bam already exists."
else 
	samtools view -@ 20 -S -b $1.sam > $1.bam
fi

# Sort the bam file
if [ -s $1.sort.bam ]; then
    echo "$1.sort.bam already exists."
else 
	samtools sort -@ 20 $1.bam > $1.sort.bam
fi

# Index the bam file
samtools index -@ 20 $1.sort.bam

# Check if bam file is indexed, and if so, remove the sam file
if [ -s $1.sort.bam ]; then
    rm -rf $1.sam
fi

# Remove duplicate reads
if [ -s $1.sort.rmdup.bam ]; then
    echo "$1.sort.rmdup.bam already exists."
else 
	java -jar /n/app/picard/2.8.0/bin/picard-2.8.0.jar MarkDuplicates REMOVE_DUPLICATES=true I=$1.sort.bam O=$1.sort.rmdup.bam M=$1.picard.output.txt
fi

# Index the bam file after removing duplicates
samtools index -@ 20 $1.sort.rmdup.bam

# Filter reads by chromosome and mapping quality
if [ -s $1.sort.rmdup.chrom.bam ]; then
    echo "$1.sort.rmdup.chrom.bam already exists."
else 
	samtools view \
		-@ 20 \
		-q 40 \
		-b $1.sort.rmdup.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY > $1.sort.rmdup.chrom.bam 
fi

# Index the filtered bam file
samtools index -@ 20 $1.sort.rmdup.chrom.bam

# Filter reads by fragment size (keep reads <= 120 bp)
if [ -s $1.sort.rmdup.chrom.120bp.bam ]; then
    echo "$1.sort.rmdup.chrom.120bp.bam already exists."
else 
	alignmentSieve -p 20 -b $1.sort.rmdup.chrom.bam -o $1.sort.rmdup.chrom.120bp.bam --maxFragmentLength 120
fi

# Index the size-filtered bam file
samtools index -@ 20 $1.sort.rmdup.chrom.120bp.bam
 
# Generate a spike-in normalized coverage file
yeast_reads=`awk '/aligned concordantly exactly 1 time/ {print $1}' $1_spike_in_alignment_file.log`
scale_factor=`echo "scale=5 ; 10000 / $yeast_reads" | bc -l`
#scale_factor = samtools flagstat -@ 20 $1.sort.rmdup.chrom.120bp.bam | awk '/read1/ {print $1}' * S

if [ -s $1.sort.rmdup.chrom.120bp.spikein.bw ]; then
    echo "$1.sort.rmdup.chrom.120bp.spikein.bw already exists."
else 
	bamCoverage -bs 1 -p 20 --scaleFactor $scale_factor -b $1.sort.rmdup.chrom.120bp.bam -o $1.sort.rmdup.chrom.120bp.spikein.bw
fi

# Call peaks on the processed bam file
if [ -s $1_macs2.narrowPeak ]; then
    echo "$1_macs2.narrowPeak already exists."
else 
    macs2 callpeak -t $1.sort.rmdup.chrom.120bp.bam -f BAMPE -q 0.05 --keep-dup all -g mm -n $1_macs2
fi
