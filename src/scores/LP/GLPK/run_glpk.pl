#!/usr/local/bin/perl -w

if ($#ARGV != 1)
{
  print "run_glpk.pl <inputdir> <outputdir>\n" ;
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

        system ("$home/glpsol --cpxlp -m $home/scp_glpk.mod -d $filename -o $outputFileName");

       
    }

    }

close DIR;
chdir '..';





