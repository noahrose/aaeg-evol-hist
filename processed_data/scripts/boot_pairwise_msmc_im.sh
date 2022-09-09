#!/bin/bash
#SBATCH -N1 -c3 -t 2:00:00

mkdir -p $6

generate_multihetsep.py --chr NC_035107.1 --mask $1.mask.NC_035107.1.bed.gz \
--mask $2.mask.NC_035107.1.bed.gz \
--mask $3.mask.NC_035107.1.bed.gz \
--mask $4.mask.NC_035107.1.bed.gz \
--mask $5.chr1.bed.gz \
$1.NC_035107.1.final.vcf.gz $2.NC_035107.1.final.vcf.gz $3.NC_035107.1.final.vcf.gz $4.NC_035107.1.final.vcf.gz\
> $6/$6.NC_035107.1.multihetsep.txt &

generate_multihetsep.py --chr NC_035108.1 --mask $1.mask.NC_035108.1.bed.gz \
--mask $2.mask.NC_035108.1.bed.gz \
--mask $3.mask.NC_035108.1.bed.gz \
--mask $4.mask.NC_035108.1.bed.gz \
--mask $5.chr2.bed.gz \
$1.NC_035108.1.final.vcf.gz $2.NC_035108.1.final.vcf.gz $3.NC_035108.1.final.vcf.gz $4.NC_035108.1.final.vcf.gz\
> $6/$6.NC_035108.1.multihetsep.txt &

generate_multihetsep.py --chr NC_035109.1 --mask $1.mask.NC_035109.1.bed.gz \
--mask $2.mask.NC_035109.1.bed.gz \
--mask $3.mask.NC_035109.1.bed.gz \
--mask $4.mask.NC_035109.1.bed.gz \
--mask $5.chr3.bed.gz \
$1.NC_035109.1.final.vcf.gz $2.NC_035109.1.final.vcf.gz $3.NC_035109.1.final.vcf.gz $4.NC_035109.1.final.vcf.gz\
> $6/$6.NC_035109.1.multihetsep.txt &

wait

sbatch msmc_dir.sh $6

multihetsep_bootstrap.py -n $7 -s 20000000 --chunks_per_chromosome 20 --nr_chromosomes 3  $6/boot $6/$6.NC_035107.1.multihetsep.txt $6/$6.NC_035108.1.multihetsep.txt $6/$6.NC_035109.1.multihetsep.txt

for i in $6/boot*; do
sbatch msmc_dir.sh $i
done

