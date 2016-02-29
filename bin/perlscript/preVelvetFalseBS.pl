#!/usr/bin/perl
#
# By Magnus Isaksson 2010
#
# Script creates false base-space from color-space.
# This will allow for de novo assembly by Velvet.
#

use strict;
use warnings;

# Subroutine prototypes
sub usage;

my ($next_arg, $infile, $outfile);
my $runGnuplot=0;

if(scalar(@ARGV) == 0) { usage(); }
# Parse the command line
while(scalar @ARGV > 0) {
     $next_arg = shift(@ARGV);
     if($next_arg eq "-i") { $infile = shift(@ARGV); }
     elsif($next_arg eq "-o") { $outfile = shift(@ARGV); }
     else { print "Invalid argument: $next_arg"; usage(); }
}

# MAIN ----------
open(OUTPUT, "> $outfile") or die "Oops, could not create output file: $!";
open(INPUT, "< $infile") or die "Oops, could not open input file: $!";
while(<INPUT>) {

	chomp;
	if (!(m/^#/)) {

		# Looking for read ids
		if (m/@(.*?_F3).*?/) {
			print OUTPUT ">$_\n";
		}

		if (m/^[ACGT][0123]*?$/) {
			
			# Trim zero base and first color call.
			$_=~s/\D\d//; 
			# Make false bases-space.
			$_=~tr/0123/ACGT/;
			print OUTPUT "$_\n";
		}
	
	}
}
close(INPUT);
close(OUTPUT);
#----------------

# Sub show how to run this script.
sub usage {
print "\nUsage: $0 -i <input> -o <output_dir>\n 
 -i	Input file format in color-space fastq (output from MosaikAligner).
 -o	Output file, in false bases-space fasta format.\n\n";   
 exit(0);
}
