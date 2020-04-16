#!/usr/local/bin/perl -w
if ($#ARGV != 0)
{
  print "runLPall.pl <pdb>\n" ;
  exit ;
}

my $pdb = $ARGV[0];

use FindBin;
my $home="$FindBin::Bin";



    print "File: $pdb\n";
    $inputFileName = $pdb.".dat";

    $outputFileName = $pdb.".sol";

#system ("$lpscripts/solve_cks_lp  $inputFileName > $outputFileName");
	system ("ampl $home/cks_lp1.mod $inputdir/$filename $home/cks_lp.ampl > $outputdir/$outputFileName");


        open(DATA, "$outputFileName");

        $nonIntegralFound = 0;
        $startCheck = 0;

        while ($line=<DATA> and ($nonIntegralFound == 0)) {
            if ($startCheck == 1) {
                my @tmp=split(/ /,$line);

                if ($#tmp == 1) {

                    my $val=int $tmp[1];
                    if ($val < 1) {
                        $nonIntegralFound = 1;
                    }
                }
            }
            if (($startCheck == 0) and $line eq "# primal vertex variables\n") {
                $startCheck = 1;

            }
        }
	if ($nonIntegralFound == 1){
  print "***************** Run ILP for $pdb ********************\n";
                        $nonIntegralFound = 1;
#			system ("$lpscripts/solve_cks_ilp $inputFileName > $outputFileName");
	system ("ampl $home/cks_ilp.mod $inputdir/$filename $home/cks_ilp.ampl > $outputdir/$outputFileName");
			
}

        close (DATA);







