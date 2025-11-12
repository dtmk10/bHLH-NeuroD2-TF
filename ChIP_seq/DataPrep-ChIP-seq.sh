#!/bin/bash


eval "$(conda shell.bash hook)"
mkdir ChIP-seq-Analysis
cd $PWD/ChIP-seq-Analysis

conda create -y -n chip-seq python=3.7.7
conda activate chip-seq
 conda install -y bioconda::sra-tools
 conda install -y bioconda::cutadapt
 conda install -y bioconda::bwa
 conda install -y bioconda::samtools
 conda install -y bioconda::macs2

fastq-dump --split-files -v --gzip --outdir $PWD/01_SRR3955800_RawData SRR3955800
fastq-dump --split-files -v --gzip --outdir $PWD/01_SRR3955796_RawData SRR3955796

cutadapt -m 20 -q 30 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGC -o SRR3955800_1.trim.fastq 01_SRR3955800_RawData/SRR3955800_1.fastq.gz
cutadapt -m 20 -q 30 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGC -o SRR3955796_1.trim.fastq 01_SRR3955796_RawData/SRR3955796_1.fastq.gz

mkdir mapped_reads
bwa mem ../mouse_genome/mm10.fa ./SRR3955800_1.trim.fastq > mapped_reads/Input_SRR3955800.sam
bwa mem ../mouse_genome/mm10.fa ./SRR3955796_1.trim.fastq > mapped_reads/Peaks_SRR3955796.sam

cd mapped_reads
samtools sort -O bam Input_SRR3955800.sam > Input_SRR3955800.bam
samtools sort -O bam Peaks_SRR3955796.sam > Peaks_SRR3955796.bam

samtools index Input_SRR3955800.bam
samtools index Peaks_SRR3955796.bam

cd ..
mkdir PeakCalling
macs2 callpeak -g mm -f BAM -t mapped_reads/Peaks_SRR3955796.bam -c mapped_reads/Input_SRR3955800.bam --bw 200 --outdir PeakCalling/ -n Neuro_D2

meme-chip -oc IP_meme --db JASPAR_db/JASPAR2024_CORE_non-redundant_pfms_meme.txt IP_peaks.resized.bed.fa