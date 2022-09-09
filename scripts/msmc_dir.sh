#!/bin/bash
#SBATCH -N1 -c16 -t 4:00:00

#delete any empty input files (potentially generated during bootstrapping of long sparse intervals)
find $1 -name "*multihetsep*.txt" -type f -empty -delete

BASE=$(basename $1)
msmc2 -t 8 -s -I 0-4,0-5,0-6,0-7,1-4,1-5,1-6,1-7,2-4,2-5,2-6,2-7,3-4,3-5,3-6,3-7 -o $1/cross $1/*multihetsep*.txt &
msmc2 -t 2 -s -I 0,1,2,3 -o $1/pop1 $1/*multihetsep*.txt &
msmc2 -t 2 -s -I 4,5,6,7 -o $1/pop2 $1/*multihetsep*.txt &
wait
combineCrossCoal.py $1/cross.final.txt $1/pop1.final.txt $1/pop2.final.txt > $1/$BASE.combined.final.txt
~/MSMC-IM/MSMC_IM.py -mu 1e-8 -N1 1e5 -N2 1e5 -o $1/$BASE.IM --printfittingdetails $1/$BASE.combined.final.txt

