#!/bin/bash

$2/score_bmf_3 $1.pdb $1.bmf3 >/dev/null
pdb2pqr --ff=CHARMM $1.pdb $1.pqr
$2/bluues_new_2 $1.pqr $1 >/dev/null

# Read BMF
bmf=`awk '{print $2}' $1.bmf3`
bumps=`awk '{print $4}' $1.bmf3`
tors_coars=`awk '{print $8}' $1.bmf3`

# Read BLUUES
tot_energy=`head -n5 $1.solv_nrg | tail -n1 | awk '{print $3}'`

# Final energy
G=`echo "0.17378*$bmf + 0.25789*$bumps + 0.26624*$tors_coars + 0.16446*$tot_energy" | bc`
echo $G > $1.score
echo $bmf > $1.score_bmf
echo $tot_energy > $1.score_bluues
