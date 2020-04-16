#!/usr/local/bin/perl -w

if ($#ARGV != 0)
{
  print "run_glpk_file.pl <pdb>\n" ;
  exit ;
}

my $pdb = $ARGV[0];
use FindBin;
my $home="$FindBin::Bin";

    print "File: $pdb\n";
    $inputFileName = $pdb.".dat";
    $outputFileName = $pdb.".sol";
    system ("$home/glpsol --cpxlp -m $home/scp_glpk.mod -d $inputFileName -o $outputFileName");




