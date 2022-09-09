#!/bin/bash
#SBATCH -N1 -c1 -t 2:00:00

cat chr1_mask.vcf | perl vcf2ahmm_ogd.pl > ogd_chr1_mask.inputfile
cat chr2_mask.vcf | perl vcf2ahmm_ogd.pl > ogd_chr2_mask.inputfile
cat chr3_mask.vcf | perl vcf2ahmm_ogd.pl > ogd_chr3_mask.inputfile
cat chr1.vcf | perl vcf2ahmm_ogd.pl > ogd_chr1.inputfile
cat chr2.vcf | perl vcf2ahmm_ogd.pl > ogd_chr2.inputfile
cat chr3.vcf | perl vcf2ahmm_ogd.pl > ogd_chr3.inputfile

Rscript --vanilla fixAhmm.R chr1.gmap.gz ogd_chr1.inputfile ogd_chr1.gmap.inputfile
Rscript --vanilla fixAhmm.R chr2.gmap.gz ogd_chr2.inputfile ogd_chr2.gmap.inputfile
Rscript --vanilla fixAhmm.R chr3.gmap.gz ogd_chr3.inputfile ogd_chr3.gmap.inputfile
Rscript --vanilla fixAhmm.R chr1.gmap.gz ogd_chr1_mask.inputfile ogd_chr1_mask.gmap.inputfile
Rscript --vanilla fixAhmm.R chr2.gmap.gz ogd_chr2_mask.inputfile ogd_chr2_mask.gmap.inputfile
Rscript --vanilla fixAhmm.R chr3.gmap.gz ogd_chr3_mask.inputfile ogd_chr3_mask.gmap.inputfile

ancestry_hmm -v -i ogd_chr1.gmap.inputfile -s ogd.txt -a 2 0.91 0.09 -p 0 10000 0.91 -p 1 -500 0.09 > ogd_posterior/ogd_chr1.log
mv *Ouaga*.viterbi ogd_posterior/chr1
ancestry_hmm -v -i ogd_chr2.gmap.inputfile -s ogd.txt -a 2 0.92 0.08 -p 0 10000 0.92 -p 1 -500 0.08 > ogd_posterior/ogd_chr2.log
mv *Ouaga*.viterbi ogd_posterior/chr2
ancestry_hmm -v -i ogd_chr3.gmap.inputfile -s ogd.txt -a 2 0.90 0.1 -p 0 10000 0.90 -p 1 -500 0.1 > ogd_posterior/ogd_chr3.log
mv *Ouaga*.viterbi ogd_posterior/chr3
ancestry_hmm -v -i ogd_chr1_mask.gmap.inputfile -s ogd.txt -b 80 1000 -a 2 0.94 0.06 -p 0 10000 0.94 -p 1 -500 0.06 > ogd_posterior/ogd_chr1_mask.log
mv *Ouaga*.viterbi ogd_posterior/chr1_mask
ancestry_hmm -v -i ogd_chr2_mask.gmap.inputfile -s ogd.txt -b 80 1000 -a 2 0.95 0.05 -p 0 10000 0.95 -p 1 -500 0.05 > ogd_posterior/ogd_chr2_mask.log
mv *Ouaga*.viterbi ogd_posterior/chr2_mask
ancestry_hmm -v -i ogd_chr3_mask.gmap.inputfile -s ogd.txt -b 80 1000 -a 2 0.91 0.09 -p 0 10000 0.91 -p 1 -500 0.09 > ogd_posterior/ogd_chr3_mask.log
mv *Ouaga*.viterbi ogd_posterior/chr3_mask

cat ogd_posterior/ogd_chr1_mask.log | grep -P "\t1\t" | cut -f3 > ogd_posterior/chr1_estimates.txt
cat ogd_posterior/ogd_chr2_mask.log | grep -P "\t1\t" | cut -f3 > ogd_posterior/chr2_estimates.txt
cat ogd_posterior/ogd_chr3_mask.log | grep -P "\t1\t" | cut -f3 > ogd_posterior/chr3_estimates.txt
