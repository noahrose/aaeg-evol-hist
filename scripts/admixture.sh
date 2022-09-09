#!/bin/bash
#SBATCH -N1 -c4 -t 2:00:00

plink --bfile /tigress/noahr/afrI3.bams/verily_unlinked_norepeat --chr 1 --recode vcf-iid --make-bed --out chr1_adm
plink --bfile /tigress/noahr/afrI3.bams/verily_unlinked_norepeat --chr 2 --recode vcf-iid --make-bed --out chr2_adm
plink --bfile /tigress/noahr/afrI3.bams/verily_unlinked_norepeat --chr 3 --recode vcf-iid --make-bed --out chr3_adm

plink --bfile /tigress/noahr/afrI3.bams/verily_unlinked_norepeat --exclude range outlier_regions.range --chr 1 --recode vcf-iid --make-bed --out chr1_adm_mask
plink --bfile /tigress/noahr/afrI3.bams/verily_unlinked_norepeat --exclude range outlier_regions.range --chr 2 --recode vcf-iid --make-bed --out chr2_adm_mask
plink --bfile /tigress/noahr/afrI3.bams/verily_unlinked_norepeat --exclude range outlier_regions.range --chr 3 --recode vcf-iid --make-bed --out chr3_adm_mask

admixture -j4 chr1_adm.bed 3
admixture -j4 chr2_adm.bed 3
admixture -j4 chr3_adm.bed 3
admixture -j4 chr1_adm_mask.bed 3
admixture -j4 chr2_adm_mask.bed 3
admixture -j4 chr3_adm_mask.bed 3

> ogd_props.txt
> kum_props.txt
> ogd_mask_props.txt
> kum_mask_props.txt
paste unrelated_nomasc.txt chr1_adm_mask.3.Q | grep Kumasi |awk '{sum1+=$2;sum2+=$3;sum3+=$4} END {print sum1/NR,sum2/NR,sum3/NR}' >> kum_mask_props.txt
paste unrelated_nomasc.txt chr2_adm_mask.3.Q | grep Kumasi |awk '{sum1+=$2;sum2+=$3;sum3+=$4} END {print sum1/NR,sum2/NR,sum3/NR}' >> kum_mask_props.txt
paste unrelated_nomasc.txt chr3_adm_mask.3.Q | grep Kumasi |awk '{sum1+=$2;sum2+=$3;sum3+=$4} END {print sum1/NR,sum2/NR,sum3/NR}' >> kum_mask_props.txt
paste unrelated_nomasc.txt chr1_adm_mask.3.Q | grep 'LG_Ouaga\|Tab_Ouaga' |awk '{sum1+=$2;sum2+=$3;sum3+=$4} END {print sum1/NR,sum2/NR,sum3/NR}' >> ogd_mask_props.txt
paste unrelated_nomasc.txt chr2_adm_mask.3.Q | grep 'LG_Ouaga\|Tab_Ouaga' |awk '{sum1+=$2;sum2+=$3;sum3+=$4} END {print sum1/NR,sum2/NR,sum3/NR}' >> ogd_mask_props.txt
paste unrelated_nomasc.txt chr3_adm_mask.3.Q | grep 'LG_Ouaga\|Tab_Ouaga' |awk '{sum1+=$2;sum2+=$3;sum3+=$4} END {print sum1/NR,sum2/NR,sum3/NR}' >> ogd_mask_props.txt
paste unrelated_nomasc.txt chr1_adm.3.Q | grep Kumasi |awk '{sum1+=$2;sum2+=$3;sum3+=$4} END {print sum1/NR,sum2/NR,sum3/NR}' >> kum_props.txt
paste unrelated_nomasc.txt chr2_adm.3.Q | grep Kumasi |awk '{sum1+=$2;sum2+=$3;sum3+=$4} END {print sum1/NR,sum2/NR,sum3/NR}' >> kum_props.txt
paste unrelated_nomasc.txt chr3_adm.3.Q | grep Kumasi |awk '{sum1+=$2;sum2+=$3;sum3+=$4} END {print sum1/NR,sum2/NR,sum3/NR}' >> kum_props.txt
paste unrelated_nomasc.txt chr1_adm.3.Q | grep 'LG_Ouaga\|Tab_Ouaga' |awk '{sum1+=$2;sum2+=$3;sum3+=$4} END {print sum1/NR,sum2/NR,sum3/NR}' >> ogd_props.txt
paste unrelated_nomasc.txt chr2_adm.3.Q | grep 'LG_Ouaga\|Tab_Ouaga' |awk '{sum1+=$2;sum2+=$3;sum3+=$4} END {print sum1/NR,sum2/NR,sum3/NR}' >> ogd_props.txt
paste unrelated_nomasc.txt chr3_adm.3.Q | grep 'LG_Ouaga\|Tab_Ouaga' |awk '{sum1+=$2;sum2+=$3;sum3+=$4} END {print sum1/NR,sum2/NR,sum3/NR}' >> ogd_props.txt

