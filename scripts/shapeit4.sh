#!/bin/bash
#SBATCH -N1 -c16 --mem=320G -t 60:00:00
module load samtools
module load boost
bcftools index $1
shapeit4.2 --sequencing --effective-size 1000000 --use-PS 0.0001 --input $1 --region $2 --output $(basename $1 .bcf).phased.bcf --map $3 --thread 10 --log $1.shapeit.log
