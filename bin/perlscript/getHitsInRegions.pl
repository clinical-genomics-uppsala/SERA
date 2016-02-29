#!/usr/bin/perl
#
# Extracts all reads in region from aligned reads
# (sedd file format).
#
# By: Magnus Isaksson 2010
#

use strict;
use warnings;

use File::Spec;
use FileHandle;

# Globals
my($next_arg, $reads_file, $ampregion_file, $hits_file);
my $unique_reads_only = 0;

# Hash for all regions.
my $regions = {};

# Subroutine prototypes
sub checkRead;
sub usage;

if(scalar(@ARGV) == 0) { usage(); }
# Parse the command line
while(scalar @ARGV > 0) {
     $next_arg = shift(@ARGV);
     if($next_arg eq "-i") { $ampregion_file = shift(@ARGV); }
     elsif($next_arg eq "-r") { $reads_file = shift(@ARGV); }
     elsif($next_arg eq "-o") { $hits_file=shift(@ARGV); }
     elsif($next_arg eq "-uniq") { $unique_reads_only=1; }
     else { print "Invalid argument: $next_arg"; usage(); }
}

# Missing options
if ( !($ampregion_file) || !($reads_file)) { usage(); }

# Direct output to stdout if output file not set.
if (!$hits_file) { $hits_file="/dev/stdout"; }

## Main

# Hash ampregion file.
open(REGION, "< $ampregion_file") or die "Oops, could not open region file: $!";
while(<REGION>) {

	chomp;
	if (!/^#/) {
		my @columns = split(/\t/, $_);
		$regions->{ $columns[1] }->{ $columns[0] } = {
														start => $columns[2],
														end   => $columns[3],
													 };
	}
}
close(REGION);


# Filter reads file.
open (READS, "< $reads_file") or die "Oops, could not open reads file: $!";
open (HITS, "> $hits_file") or die "Oops, could not open hits file: $!";
while(<READS>) {

	chomp;
	if (!/^#/) {
		
		my @columns = split(/\t/, $_);
		my ($totNumOfMappedPos) = $columns[5] =~ m/.*?,totNumOfMappedPos=(\d*),/;
		my ($aQual) = $columns[5] =~ m/.*?,aQual=(\d*),/;
		
		if($unique_reads_only == 1) {
			if($totNumOfMappedPos == 1) {
				if(checkRead($columns[1], $columns[2], $columns[3]) == 1) {
					print HITS "$columns[0]\t$columns[1]\t$columns[2]\t$columns[3]\t$aQual\t$totNumOfMappedPos\n";
				}		
			}
		} else {
			
			if(checkRead($columns[1], $columns[2], $columns[3]) == 1) {
					print HITS "$columns[0]\t$columns[1]\t$columns[2]\t$columns[3]\t$aQual\t$totNumOfMappedPos\n";
			}
		}
		
	} 
}	
close(HITS);
close(READS);
## ----

# Sub check i read aligns within region.
sub checkRead {
	
	my($nc, $start, $end) = @_;
	my $hit = 0;
	
	for my $regionId ( keys %{$regions->{ $nc }} ) {
			
		if( ($start <= $regions->{ $nc }->{ $regionId }->{ end }) && ($end >= $regions->{ $nc }->{ $regionId }->{ start })) {
			$hit = 1;
		}
	}
	
	return $hit;	
}

# Sub to tell you how to use this script.
sub usage {
  print "\nUsage: $0 -i <Infile> -o <Output file> -r <reads file> -uniq
 
 -i       Infile, tab-delimited:
          {ID	Chrom	Base	Start	End	Strand
 -o       Output file, tab-delimited (using stdout if not set):
          {ReadID	Chrom	Start End	Score	Hits}
 -r       Aligned reads file in Sedd format.
 -uniq	  Only unique reads are saved, if set.\n\n";
  
  exit(1);
}
