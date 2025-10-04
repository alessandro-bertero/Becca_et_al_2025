#!/bin/bash

fastq_files=/home/shared_folder/fastq
fastq_dir=/home/shared_folder/fastq/trimmed
fasta=/home/shared_folder/genome_files/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
gtf=/home/shared_folder/Homo_sapiens.GRCh38.111.gtf
STAR --runThreadN 3 --runMode geomeGenerate --limitGenomeGenerateRAM=150000000000 --genomeDir /home/shared_folder/hg38 --genomeFastaFiles $fasta --sjdbGTFfile $gtf
index=/home/shared_folder/hg38
SBe=("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20" "21" "22" "23" "24")
for i in ${SBe[@]}
do
	file1=$fastq_files/sample"$i"_1.fq.gz
	file2=$fastq_files/sample"$i"_2.fq.gz
	trim_galore --adapter CTGTCTCTTATACACATCT --paired --trim1 --clip_R1 4 --gzip --fastqc --output_dir $fastq_dir $file1 $file2
	
	base=sample"$i"
	fq1=$fastq_dir/sample"$i"_1_val_1.fq.gz
	fq2=$fastq_dir/sample"$i"_2_val_2.fq.gz
	STAR --runThreadN 3 --genomeDir $index --readFilesIn $fq1 $fq2 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --readFilesCommand zcat --outFileNamePrefix $base"_"
done
