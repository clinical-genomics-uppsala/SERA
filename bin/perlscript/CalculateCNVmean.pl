#!/usr/bin/perl

#  
# Author: Elin Falk 2010

use warnings;
use strict;
#use criticism 'harsh';

use FileHandle;
use Carp;
use Getopt::Long;
#use Readonly;

#Globals
my ($infile, $outfile, $h);

# Subroutine prototypes
sub usage;

# Parse the command line
GetOptions(
	"i=s" => \$infile,
	"o=s" => \$outfile,
	"h=i" => \$h
	
) || usage();

# Error
if (!$infile || !$outfile || !$h) { usage(); }

my $hitColumn = $h - 1;

# Open filehandlers 
my $infh = FileHandle->new($infile) or croak "Couldn't open input file".$infile."!";
my $outfh = FileHandle->new(">".$outfile) or croak "Couldn't open output file".$outfile."!";

my @hits;
my $counter = 0;
my $total = 0;
# Go through input file
while (<$infh>) {
	# Skip comments and empty lines
	if ($_ =~ m/^#/ || $_ eq "" )  { next; }
	chomp;
	# Split on tab
	my @line = split(/\t/, $_);
	
	$hits[$counter] = $line[$hitColumn];
	$counter++;
	$total += $line[$hitColumn];
} # Ending while (<$infh>)

my @sortedHits = sort{$a <=> $b} @hits;
my $median;
if ($counter%2 == 0) {
	print {$outfh} "Median:".(($sortedHits[$counter/2] + $sortedHits[$counter/2-1])/2)."\n"; 
	$median = ($sortedHits[$counter/2] + $sortedHits[$counter/2-1])/2;
}	
else {
	print {$outfh} "Median:".$sortedHits[$counter/2-0.5]."\n";
	$median = $sortedHits[$counter/2-0.5];
}
print {$outfh} "Mean:".$total/$counter."\n";
print "Median - mean diff: ".($median - ($total/$counter))."\n";

# Print the usage help for this script.
sub usage {
  print "\nUsage: $0 \n
 -h Column with hits, first column is 1
 -i Infile, two concatemeriased .map files {Tumor_chrom\tTumor_pos\tTumor_hits\tNormal_chrom\tNormal_pos\tNormal_hits
 -o Output file\n\n";
  exit 1;
}
