#!/bin/bash
module load samtools

date

echo "getting sample-specific vcfs..."
bcftools view -R norepeat_5_30x_coverage.bed -Ou -s $1 $2.bcf | bcftools filter -S . -Ou -i "FMT/GQ >= 20 & FMT/DP >= 8" |  > temp_$2/$1_temp.bcf
echo "...hets"
bcftools view -Ov -g het temp_$2/$1_temp.bcf > temp_$2/$1_het_temp.vcf
echo "...non-hets"
bcftools view -Ob -g ^het temp_$2/$1_temp.bcf > temp_$2/$1_nohet_temp.bcf

date

echo "getting informative reads..."
extractHAIRS --bam /tigress/noahr/afrI3.bams/$1.*.bam --VCF temp_$2/$1_het_temp.vcf --out temp_$2/$1.hairs
echo "prephasing..."
HAPCUT2 --fragments temp_$2/$1.hairs --vcf temp_$2/$1_het_temp.vcf  --out temp_$2/$1_het --outvcf 1

date

echo "adding back prephased variants..."
mkdir temp_$2/$1_sort
bcftools concat -Ou temp_$2/$1_nohet_temp.bcf temp_$2/$1_het.phased.VCF | bcftools annotate -Ou -x INFO,^FMT/GT,FMT/PS | bcftools sort -T temp_$2/$1_sort -Ob > prephase_$2/$1_prephased.bcf

date

rmdir temp_$2/$1_sort
rm temp_$2/$1_temp.bcf
rm temp_$2/$1_het_temp.vcf
rm temp_$2/$1_nohet_temp.bcf
rm temp_$2/$1.hairs
rm temp_$2/$1_het.phased.VCF
rm temp_$2/$1_het
