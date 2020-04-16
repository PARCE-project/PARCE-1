#!/usr/bin/perl
if ($#ARGV != 1)
{
  print "runLPall.pl <inputdir> <outputdir>\n" ;
  exit ;
}

my $inputdir = $ARGV[0];
my $outputdir = $ARGV[1];

opendir (DIR,$inputdir) or die "Couldn't open directory $inputdir, $!";
chdir $inputdir or die "Can't change directory $dir, $!";


use FindBin;
my $home="$FindBin::Bin";

while ($filename = readdir DIR) {

if($filename =~ /.dat$/) {
   
    $var = substr($filename, 0, - 3);

    print "File: $filename\n";

    $outputFileName = $var."sol";

	system ("ampl $home/cks_lp1.mod $inputdir/$filename $home/cks_lp.ampl > $outputdir/$outputFileName");
	
        open(DATA, "$outputdir/$outputFileName");

        $nonIntegralFound = 0;
        $startCheck = 0;

        while ($line=<DATA> and ($nonIntegralFound == 0)) {
            if ($startCheck == 1) {

	chop($line);
                my @tmp=split(/ /,$line);

                if ($#tmp == 1) {


	        my $val = $tmp[1];

                    if (($val ne "1") && ($val ne "0")) {
                       $nonIntegralFound = 1;

                    }
                }
            }
            if (($startCheck == 0) and $line eq "# primal vertex variables\n") {
                $startCheck = 1;

            }
        }
	if ($nonIntegralFound == 1){
  print "***************** Run ILP for $filename ********************\n";
                        $nonIntegralFound = 1;
	system ("ampl $home/cks_ilp.mod $inputdir/$filename $home/cks_ilp.ampl > $outputdir/$outputFileName");

}

        close (DATA);
    }

    }

close DIR;
chdir '..';





