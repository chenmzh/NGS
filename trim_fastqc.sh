#!/bin/bash

for f in /local0/students/Bubble8/data/*; do
#	java -jar trimmomatic-0.39.jar SE -phred33 -threads 16 $f
        java -jar /local0/students/Bubble8/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 -threads 16 $f \
  "/local0/students/Bubble8/trimmed_recom/$(basename $f)" \
  ILLUMINACLIP:/local0/students/Bubble8/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10\
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 && \
  /local0/scratch/sysgen_2021/tools/FastQC/fastqc -o\
  /local0/students/Bubble8/trimmed_recom/\
  "/local0/students/Bubble8/trimmed_recom/$(basename $f)";
done
