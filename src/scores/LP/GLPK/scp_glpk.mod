# THIS SOFTWARE IS LICENSED UNDER THE GNU PUBLIC LICENSE VERSION 2.
# See LICENSE file available at http://compbio.cs.princeton.edu/scplp/cks_lp/
# or downloaded with this software. If you use this software please cite the
# paper mentioned at the above website.

#model;

#
# The problem structure: the number of variable residues (num_posn), 
# total number of rotamers (num_nodes), the set of positions (POSN) 
# and the set of rotamers (V).
#
param num_posn integer >= 2;
param num_nodes integer >= num_posn;
set POSN := 1..num_posn;
set V := 1..num_nodes;

#
# The structure of the positions: convert posn_size array into sets of 
# nodes C
#
param posn_size {i in POSN} integer > 0;
check: sum {i in POSN} posn_size[i] = num_nodes;

#?
set C {i in POSN} := {u in V : (sum {j in POSN : j < i} posn_size[j] < u ) and (u <= sum {j in POSN : j <= i} posn_size[j])};
#check: card((union {i in POSN} C[i]) symdiff V) = 0;

#
# The self-energies (costV) and the pairwise energies (costE).
#
set Efull := {(u,v) in {V,V}: u in V and v > u};  # upper tri 
param costV {V} default 0.0;
param costE {Efull} default 0.0;

set negative_pairs := 
    {(i,j) in POSN cross POSN :
        j > i and
        (forall {u in C[i], v in C[j]} costE[u,v] <= 0) and
        (exists {u in C[i], v in C[j]} costE[u,v] < 0)};

set positive_pairs := 
    {(i,j) in POSN cross POSN : 
        j > i and
        exists {u in C[i], v in C[j]} costE[u,v] > 0};

#set NEG_EDGES in Efull := setof {(i,j) in negative_pairs  C[i] cross C[j] } ;
#{(u,v) in C[i] cross C[j] : costE[u,v] != 0.0};

set NEG_EDGES1 in Efull :=  setof{(i,j) in negative_pairs, u in C[i], v in C[j]}(u,v);

set NEG_EDGES in NEG_EDGES1 := 
{(u,v) in NEG_EDGES1 : costE[u,v] != 0.0};

set POS_EDGES  in Efull :=  setof{(i,j) in positive_pairs, u in C[i], v in C[j]}(u,v);


#
# Edge set restricted to edges between allowed pairs.
#
set E in Efull := NEG_EDGES union POS_EDGES;

#
# Ignore the # of interations specified in the file; doesn't make sense for LP
# 

param iterations integer default 1;

#
# The variables
#
var X {V} binary >= 0;  # node variables
var Y {E} binary >= 0;  # edge variables

minimize energy: (sum {v in V} costV[v] * X[v]) +
                 (sum {(u,v) in E} costE[u,v] * Y[u,v]);

subject to column {i in POSN}:
    sum {u in C[i]} X[u] = 1;

#
# Between pairs with only negative edges we have only have
# some flow constraints
#
subject to flow_i_to_j_neg {(i,j) in negative_pairs, u in C[i] :
    exists {v in C[j]} costE[u,v] != 0.0}:
    sum {v in C[j] : (u,v) in E} Y[u,v] <= X[u];

subject to flow_j_to_i_neg {(i,j) in negative_pairs, v in C[j] :
    exists {u in C[i]} costE[u,v] != 0.0}: 
    sum {u in C[i] : (u,v) in E} Y[u,v] <= X[v];

#
# For positive positions, we have all the flow constraints & vars
#
subject to flow_i_to_j_pos {(i,j) in positive_pairs, u in C[i]} :
    sum {v in C[j]} Y[u,v] = X[u];

subject to flow_j_to_i_pos {(i,j) in positive_pairs, v in C[j]} :
    sum {u in C[i]} Y[u,v] = X[v];


