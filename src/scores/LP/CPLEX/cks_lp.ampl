# THIS SOFTWARE IS LICENSED UNDER THE GNU PUBLIC LICENSE VERSION 2.
# See LICENSE file available at http://compbio.cs.princeton.edu/scplp/cks_lp/
# or downloaded with this software. If you use this software please cite the
# paper mentioned at the above website.

option show_stats 1;

#
# timing=1 displays how long solving, reading input, and output takes
#
option cplex_options 'primal dualopt timing=1';

solve;

# print the primal variables (n of them)
print "# energy";
print energy;
print "# problem size";
print card(V);
print card(POSN);
print "# primal vertex variables";
print {v in V : X[v] > 0}: v, X[v];

#
# Uncomment next 3 lines to print out the edge variables too
#
#print "# primal edge variables";
#print card(E);
#print {(u,v) in E}: u, v, Y[u,v];
