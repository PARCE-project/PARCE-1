# THIS SOFTWARE IS LICENSED UNDER THE GNU PUBLIC LICENSE VERSION 2.
# See LICENSE file available at http://compbio.cs.princeton.edu/scplp/cks_lp/
# or downloaded with this software. If you use this software please cite the
# paper mentioned at the above website.

#
# Set the options
#
option show_stats 1;

#
# 1e-9 is the lowest allowed value for mipgap
# 0 is the default for absmipgap (but listed for sure)
# timing=1 displays how long solving, reading input, and output takes
#
option cplex_options 'primal dualopt absmipgap=0 mipgap=1e-9 timing=1';

#
# Find iteration solutions
#
for {1..iterations} {

  solve;

  #
  # print the solution we got
  #
  print "# energy";
  print energy;
  print "# problem size";
  print card(V);
  print card(POSN);
  print "# primal vertex variables";
  print {v in V : X[v] > 0}: v, X[v];

  #
  # Save solution in model's data
  #
  let M := M + 1;
  for {u in V} let Solns[M,u] := X[u]
}  
