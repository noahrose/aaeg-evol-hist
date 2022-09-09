#!/bin/bash
module load samtools

SAMP=$1
REF=$2
CHR=$3
BAM=/tigress/noahr/afrI3.bams/${SAMP}.*bam
PHASE=$4
echo $SAMP
echo $BAM

DEPTH=$(samtools depth -r $CHR $BAM | awk '{sum += $3} END {print sum / NR}')

echo $DEPTH

bcftools view -Oz -s $SAMP $4 > $SAMP.$CHR.phased.vcf.gz

samtools mpileup -B -q20 -Q20 -C50 -g -r $CHR -f $REF $BAM | bcftools call -c -V indels |\
    bamCaller.py $DEPTH $SAMP.mask.$CHR.bed.gz | bgzip -c > $SAMP.$CHR.vcf.gz

bcftools index $SAMP.$CHR.vcf.gz
bcftools index $SAMP.$CHR.phased.vcf.gz

bcftools merge --force-samples $SAMP.$CHR.vcf.gz $SAMP.$CHR.phased.vcf.gz | awk 'BEGIN {OFS="\t"}
    $0 ~ /^##/ {print $0}
    $0 ~ /^#CHROM/ {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}
    $0 !~ /^#/ {
        if(substr($11, 1, 3) != "./.")
            $10 = $11
        print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10
    }' |  bcftools view -O z > $SAMP.$CHR.final.vcf.gz



