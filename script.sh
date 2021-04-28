#!/bin/bash


for f in /local0/students/Bubble8/data/*; do
  fastqc -o /local0/students/Bubble8/raw_data_fastqc $f;
done 

# Get all the file name
# Iterate over them to do fastqc

 
# comment
