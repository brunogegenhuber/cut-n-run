# CUT&RUN processing script for Gegenhuber and Tollkuhn (2023) Methods in Molecular Biology

This script (process_fastq.sh) processes Illumina fastq files with a SLURM job scheduler to generate: 1) duplicate-filtered BAM files, 2) spike-in normalized bigwig tracks, and 3) MACS2 peaks. This script requires: cutadapt, bowtie2, samtools, picard, deepTools, and MACS2.
