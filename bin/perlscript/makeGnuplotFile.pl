#!/usr/bin/perl
#
# By Magnus Isaksson 2011
#

use strict;
use warnings;
use XML::Simple;
use Data::Dumper;

# Subroutine prototypes
sub usage;

my ($next_arg, $inpath, $infile, $outputdir, $plottype);
my $runGnuplot=0;

if(scalar(@ARGV) == 0) { usage(); }
# Parse the command line
while(scalar @ARGV > 0) {
     $next_arg = shift(@ARGV);
     if($next_arg eq "-i") { $infile = shift(@ARGV); }
     elsif($next_arg eq "-o") { $outputdir = shift(@ARGV); }
     elsif($next_arg eq "-p") { $inpath = shift(@ARGV); }
     elsif($next_arg eq "-plot") { $runGnuplot=1; }
     else { print "Invalid argument: $next_arg"; usage(); }
}

# MAIN ----------

# Read XML-file.
my $inputXml = XML::Simple->new();
my $doc = $inputXml->XMLin($infile);

#print Dumper($doc); exit 1;

foreach my $plot (keys (%{$doc->{plot}})){

   # Create gnuplot file.
   open(OUTPUT, "> $outputdir/$doc->{plot}->{$plot}->{output}.gnuplot") or die "Oops, could not open input file: $!";

	# Print general settings	
	print OUTPUT "set terminal postscript portrai enhanced;\n";
	print OUTPUT "set size 1,1;";
	print OUTPUT "set output \"| ps2pdf - '$outputdir/$doc->{plot}->{$plot}->{output}.pdf'\";\n";
	print OUTPUT "set grid;\n";
	print OUTPUT "set yrange[0:100];\n";
	print OUTPUT "set ylabel \"Cumulative base frequency\"\n";	

	# Make multiplot
	print OUTPUT "set multiplot\n";

	# Make Nusbaum plot --------------------------------

	print OUTPUT "set size 1,0.5\n";
	print OUTPUT "set origin 0,0\n";
	print OUTPUT "set xrange[0:1];\n";
	print OUTPUT "set xlabel \"Normalized coverage\";\n";
	print OUTPUT "set notitle;\n";
	print OUTPUT "set nokey";
	
	# Nusbaum data
	my $plotString = "plot";
	foreach my $data ((@{$doc->{plot}->{$plot}->{data}})) {
		
		my $dataFile = $data->{nusbaum};
		if(!($dataFile =~ m/\//g)) { $dataFile = $inpath."/".$dataFile; }

		my $color = $data->{color};
		my $lable = $data->{lable};
	
		$plotString = "$plotString \'$dataFile\' using 1:2 with lines lc rgb \"$color\" lw 4 title \'$lable\',";
	}
	$plotString =~ s/,$/;/;
	print OUTPUT "$plotString\n";

	# Make Stenberg plot --------------------------------

	print OUTPUT "set size 1,0.5\n";
	print OUTPUT "set origin 0,0.5\n";
	print OUTPUT "set auto x;\n";
	print OUTPUT "set logscale x;\n";
	print OUTPUT "set xlabel \"Coverage\";\n";
	print OUTPUT "set title '{/=25 $plot}';\n";
	
	# Stenberg data
	$plotString = "plot";
	foreach my $data ((@{$doc->{plot}->{$plot}->{data}})) {
		
		my $dataFile = $data->{stenberg};
		if(!($dataFile =~ m/\//g)) { $dataFile = $inpath."/".$dataFile; }

		my $color = $data->{color};
		my $lable = $data->{lable};
	
		$plotString = "$plotString \'$dataFile\' using 1:2 with lines lc rgb variable lw 4 title \'$lable\',";
	}
	$plotString =~ s/,$/;/;
	print OUTPUT "$plotString\n";

	# Done plotting	
	print OUTPUT "unset multiplot\n";

   close(OUTPUT);

   if($runGnuplot != 0) { print "Running Gnuplot...\n"; system("gnuplot $outputdir/$doc->{plot}->{$plot}->{output}.gnuplot"); }
}

#----------------

# Sub show how to run this script.
sub usage {
print "\nUsage: $0 -i <input> -o <output_dir>\n 
 -i	Input file format:
 	Plotname {tab} type {tab} datafile1,title,color:datafile2,title,color...
	type = stenberg, nusbaum or ericsson.
	color = Black, Red, Green, Blue, Yellow (as examples)
 .p     Path to find datafiles in if no path is given in input xml-file.
 -o	Output directory.
 -plot  Run Gnuplot on created file.\n\n";   
 exit(0);
}
