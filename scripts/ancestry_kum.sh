#!/bin/bash
#SBATCH -N1 -c1 -t 2:00:00

cat chr1_mask.vcf | perl vcf2ahmm_kum.pl > kum_chr1_mask.inputfile
cat chr2_mask.vcf | perl vcf2ahmm_kum.pl > kum_chr2_mask.inputfile
cat chr3_mask.vcf | perl vcf2ahmm_kum.pl > kum_chr3_mask.inputfile
cat chr1.vcf | perl vcf2ahmm_kum.pl > kum_chr1.inputfile
cat chr2.vcf | perl vcf2ahmm_kum.pl > kum_chr2.inputfile
cat chr3.vcf | perl vcf2ahmm_kum.pl > kum_chr3.inputfile

Rscript --vanilla fixAhmm.R chr1.gmap.gz kum_chr1.inputfile kum_chr1.gmap.inputfile
Rscript --vanilla fixAhmm.R chr2.gmap.gz kum_chr2.inputfile kum_chr2.gmap.inputfile
Rscript --vanilla fixAhmm.R chr3.gmap.gz kum_chr3.inputfile kum_chr3.gmap.inputfile
Rscript --vanilla fixAhmm.R chr1.gmap.gz kum_chr1_mask.inputfile kum_chr1_mask.gmap.inputfile
Rscript --vanilla fixAhmm.R chr2.gmap.gz kum_chr2_mask.inputfile kum_chr2_mask.gmap.inputfile
Rscript --vanilla fixAhmm.R chr3.gmap.gz kum_chr3_mask.inputfile kum_chr3_mask.gmap.inputfile

ancestry_hmm -v -i kum_chr1.gmap.inputfile -s kum.txt -a 2 0.94 0.06 -p 0 10000 0.94 -p 1 -500 0.06 > kum_posterior/kum_chr1.log
mv *Kumasi*.viterbi kum_posterior/chr1
ancestry_hmm -v -i kum_chr2.gmap.inputfile -s kum.txt -a 2 0.91 0.09 -p 0 10000 0.91 -p 1 -500 0.09 > kum_posterior/kum_chr2.log
mv *Kumasi*.viterbi kum_posterior/chr2
ancestry_hmm -v -i kum_chr3.gmap.inputfile -s kum.txt -a 2 0.95 0.05 -p 0 10000 0.95 -p 1 -500 0.05 > kum_posterior/kum_chr3.log
mv *Kumasi*.viterbi kum_posterior/chr3
ancestry_hmm -v -i kum_chr1_mask.gmap.inputfile -s kum.txt -b 80 1000 -a 2 0.94 0.06 -p 0 10000 0.94 -p 1 -500 0.06 > kum_posterior/kum_chr1_mask.log
mv *Kumasi*.viterbi kum_posterior/chr1_mask
ancestry_hmm -v -i kum_chr2_mask.gmap.inputfile -s kum.txt -b 80 1000 -a 2 0.94 0.06 -p 0 10000 0.94 -p 1 -500 0.06 > kum_posterior/kum_chr2_mask.log
mv *Kumasi*.viterbi kum_posterior/chr2_mask
ancestry_hmm -v -i kum_chr3_mask.gmap.inputfile -s kum.txt -b 80 1000 -a 2 0.96 0.04 -p 0 10000 0.96 -p 1 -500 0.04 > kum_posterior/kum_chr3_mask.log
mv *Kumasi*.viterbi kum_posterior/chr3_mask

cat kum_posterior/kum_chr1_mask.log | grep -P "\t1\t" | cut -f3 > kum_posterior/chr1_estimates.txt
cat kum_posterior/kum_chr2_mask.log | grep -P "\t1\t" | cut -f3 > kum_posterior/chr2_estimates.txt
cat kum_posterior/kum_chr3_mask.log | grep -P "\t1\t" | cut -f3 > kum_posterior/chr3_estimates.txt
