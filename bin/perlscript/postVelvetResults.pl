#!/usr/bin/perl
#
# Extracts contigs information into a results file.
# By Magnus Isaksson 2010
#

use strict;
use warnings;
use FileHandle;

sub usage;

my ($next_arg, $infile, $outfile);

if(scalar(@ARGV) == 0){
    &usage();
}
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-i")    { $infile = shift(@ARGV); }
    elsif($next_arg eq "-o") { $outfile = shift(@ARGV); }
    else { print "Invalid argument: $next_arg"; usage(); }
}

if((!$infile)||(!$outfile)) { usage(); }

# Main -----------------------------
open(OUTPUT, "> $outfile") or die "Oops, could not open output file: $!";
print OUTPUT "NodeID\tSize\tCoverage\n";

open(INPUT, "< $infile") or die "Oops, could not open input file: $!";
while(<INPUT>) {
	
	chomp;
	if (/^>/) {
	
		my ($node, $length, $cov) = $_ =~ m/>NODE_(\d*?)_length_(\d*?)_cov_(.*?)$/;
		print OUTPUT "$node \t $length \t $cov \n";
		
	}
}
close(INPUT);
close(OUTPUT);
#-----------------------------------

# Print the usage help for this script.
sub usage {
  print "\nUsage: $0 -i <infile> -o <outfile>\n
  -i	Input file, Velvet contig file.
  -o	Output file.\n\n";
  exit(0);
}