#!/usr/bin/perl -w
#
# By Magnus Isaksson 2010
# Script creates a pileup from a sam file. 
# Id according to SEDD standard.
#

use strict;
use warnings;
use FileHandle;

sub usage;

my $pileup = {}; # Contains all reads.
my ($next_arg, $samFile, $outFile);

if(scalar(@ARGV) == 0){
    &usage();
}
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-i")    { $samFile = shift(@ARGV); }
    elsif($next_arg eq "-o") { $outFile = shift(@ARGV); }
    else { print "Invalid argument: $next_arg"; usage(); }
}

## Error
if(!$samFile) { usage(); }

# If no output file write to stdout.
if (!$outFile) { $outFile="/dev/stdout"; }

## Hashing input reads file.
open(SAM, "< $samFile") or die "Oops, could not open SAM file: $!";
while(<SAM>) {
	if (!/\@/) {
		my @columns = split(/\t/, $_);
		
		my $cigarString = $columns[5];
		my ($matchLength) = $cigarString =~ m/^(\d{2})M$/;
		
		# OBS. this script doesn't support full cigar string yet. 
		if(defined($matchLength)) {
			
			my($ncNumber, $startPos) = $columns[2] =~ m/^.*?#(.*?)#(.*?)#.*?#.*?$/;
			$startPos = $startPos + $columns[3] - 1; # Calculate correct position.
			
			for(my $i = $startPos; $i < ($startPos+$matchLength); $i++) {
				$pileup->{$ncNumber}->{$i}++;
			}
		}	
	}
}
close(SAM);

# Open output
open(OUTPUT, "> $outFile") or die "Oops, could not open output file: $!";

		for my $ncKey ( sort keys %$pileup ) {
        	for my $posKey ( sort keys %{$pileup->{ $ncKey }} ) {
            	print OUTPUT "$ncKey\t$posKey\t$pileup->{ $ncKey }{ $posKey }\n";
        	}
    	}

close(OUTPUT);

# Print the usage help for this script.
sub usage {
  print "\nUsage: $0 -i <in_file> -o <out_file> -e <region_file>\n
  -i	Input SAM file, with reads.
  -o	Output file (stdout is default if not set).\n\n";
  exit(1);
}