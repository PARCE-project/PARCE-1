#!/bin/bash

#Total Score
complex=`head -1 $1.score | tail -1 | cut -f1`
chainA=`head -1 $2.score | tail -1 | cut -f1`
chainB=`head -1 $3.score | tail -1 | cut -f1`

total=`echo "$complex - $chainA - $chainB" | bc`
echo $total > $1.score_final
rm $1_*.score
rm $1.score

#Total Score BMF
complex=`head -1 $1.score_bmf | tail -1 | cut -f1`
chainA=`head -1 $2.score_bmf | tail -1 | cut -f1`
chainB=`head -1 $3.score_bmf | tail -1 | cut -f1`

totalBMF=`echo "$complex - $chainA - $chainB" | bc`
echo $totalBMF > $1.score_bmf_final
rm $1_*.score_bmf
rm $1.score_bmf

#Total Score BLUUES
complex=`head -1 $1.score_bluues | tail -1 | cut -f1`
chainA=`head -1 $2.score_bluues | tail -1 | cut -f1`
chainB=`head -1 $3.score_bluues | tail -1 | cut -f1`

totalBLUUES=`echo "$complex - $chainA - $chainB" | bc`
echo $totalBLUUES > $1.score_bluues_final
rm $1_*.score_bluues
rm $1.score_bluues

# Delete additional data
rm $1_*.pqr
rm $1.pqr
rm $1_*.gbr
rm $1.gbr
rm $1_*.solv_nrg
rm $1.solv_nrg
rm $1_*.bmf3
rm $1.bmf3
