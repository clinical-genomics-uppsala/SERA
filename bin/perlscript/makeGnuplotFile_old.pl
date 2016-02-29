#!/usr/bin/perl
#
# By Magnus Isaksson 2010
#

use strict;
use warnings;

# Subroutine prototypes
sub usage;

my ($next_arg, $infile, $outputdir, $plottype);
my $runGnuplot=0;

if(scalar(@ARGV) == 0) { usage(); }
# Parse the command line
while(scalar @ARGV > 0) {
     $next_arg = shift(@ARGV);
     if($next_arg eq "-i") { $infile = shift(@ARGV); }
     elsif($next_arg eq "-o") { $outputdir = shift(@ARGV); }
     elsif($next_arg eq "-p") { $runGnuplot=1; }
     else { print "Invalid argument: $next_arg"; usage(); }
}

# MAIN ----------
open(INPUT, "< $infile") or die "Oops, could not open input file: $!";
while(<INPUT>) {

	chomp;
	if (!(m/^#/)) {
	
		my @columns = split(/\t/, $_);
		my @datafiles = split(/:/, $columns[2]);

		# Create gnuplot file.
		open(OUTPUT, "> $outputdir/$columns[0].gnuplot") or die "Oops, could not open input file: $!";
		
			print OUTPUT "set terminal postscript enhanced;\n";
			print OUTPUT "set output \"| ps2pdf - '$columns[0].pdf'\";\n";
			print OUTPUT "set title '{/=25 $columns[0]}';\n";
			print OUTPUT "set grid;\n";
			print OUTPUT "set yrange[0:100];\n";
			print OUTPUT "set ylabel \"Cumulative base frequency\"\n";
			
			if($columns[1] eq "stenberg") {
				
				print OUTPUT "set auto x;\n";
				print OUTPUT "set logscale x;\n";
				print OUTPUT "set yrange[0:100];\n";
				print OUTPUT "set xlabel \"Coverage\";\n";
#				print OUTPUT "unset key;\n";

			} elsif($columns[1] eq "nusbaum") {

				print OUTPUT "set xrange[0:1];\n";
				print OUTPUT "set xlabel \"Normalized coverage\";\n";
#				print OUTPUT "unset key;\n";

			} elsif($columns[1] eq "eriksson") {

				print OUTPUT "set key inside right bottom vertical Right;\n";
			    	print OUTPUT "set xrange[1:10];\n";
				print OUTPUT "set yrange[0:100];\n";
				print OUTPUT "set xlabel \"Over sequencing factor\";\n";

			} else {
			   print "Error type of plot $columns[1].";
			   exit(1);
			}

			my $plotString = "plot";
			my $c=0;
			foreach my $datafile (@datafiles) {

				$c++;
				print OUTPUT "set style line $c lc $c lt $c lw 6;\n";

				#my ($plotfile, $title, $plotcolor) = $datafile =~ m/(.*?),(.*?),(.*?)/;
				my @plotdata = split(/,/,$datafile);
				$plotString = "$plotString \'$plotdata[0]\' using 1:2 with lines linestyle $c title \'$plotdata[1]\',";
			}

			$plotString =~ s/,$/;/;
			print OUTPUT "$plotString\n";
		
		close(OUTPUT);

		if($runGnuplot != 0) { system("gnuplot $outputdir/$columns[0].gnuplot"); }
	}
}
close(INPUT);
#----------------

# Sub show how to run this script.
sub usage {
print "\nUsage: $0 -i <input> -o <output_dir>\n 
 -i	Input file format:
 	Plotname {tab} type {tab} datafile1,title,color:datafile2,title,color...
	type = stenberg, nusbaum or ericsson.
	color = Black, Red, Green, Blue, Yellow (as examples)
 -o	Output directory.\n\n";   
 exit(0);
}
