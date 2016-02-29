#!/usr/bin/perl -w
#
# By Magnus Isaksson 2010
# Script creates a pileup from a sedd aligned reads file. 
#

use strict;
use warnings;
use FileHandle;

sub usage;

my $pileup = {}; # Contains all reads.
my $roi = {}; # Holds ROI if exist.

my ($next_arg, $seddFile, $roiFile, $outFile);
my $uniqOnly = 0;

if(scalar(@ARGV) == 0){
    &usage();
}
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-i")    { $seddFile = shift(@ARGV); }
    elsif($next_arg eq "-roi") { $roiFile = shift(@ARGV); }
    elsif($next_arg eq "-o") { $outFile = shift(@ARGV); }
    elsif($next_arg eq "-uniq") { $uniqOnly = 1; }
    else { print "Invalid argument: $next_arg"; usage(); }
}

## Error
if((!defined($seddFile))||(!defined($outFile))) { usage(); }

## If ROI file exist hash this regions.
if(defined($roiFile)) {
	open (ROIFILE, "< $roiFile") or die "Oops, could not open roi file: $!";
	while(<ROIFILE>) {
		if(!(/^#/)) {
			chomp;
			my ($ncNumber, $start, $end) = $_ =~ m/.*?\t(NC_.*?)\t(.*?)\t(.*?)\t.*?/;
			for (my $i = $start; $i < ($end+1); $i++) {
				$roi->{$ncNumber}->{$i} = 1;
			}
		}
	}
	close(ROIFILE);
}

## Hashing input reads file.
open(SEDD, "< $seddFile") or die "Oops, could not open Sedd alignment file: $!";
while(<SEDD>) {
	if (!/#/) {
		
		my @columns = split(/\t/, $_);
		
		# Fetch number of mapped positions.
		my $uniq=1;
		if($uniqOnly==1) { ($uniq) = $columns[5] =~ m/totNumOfMappedPos=(\d*)/;	}

		my $matchLength = $columns[3] - $columns[2] + 1;

		if(defined($matchLength)) {
		  if($uniq==1) {
			for(my $i = $columns[2]; $i < $columns[3]; $i++) {
				$pileup->{ $columns[1] }->{ $i }++;
			}
		  }
		} else {
			print "Error: could not calculate match length for ($_).\n";
			exit(1);
		}
	}
}
close(SEDD);

# Open output
open(OUTPUT, "> $outFile") or die "Oops, could not open output file: $!";

	if(defined($roiFile)) {
		
		# Found ROI file output in two separate files.
		for my $ncKey ( sort keys %$pileup ) {
       		for my $posKey ( sort {$a<=>$b} keys %{$pileup->{ $ncKey }} ) {
	        	   	
	        	   	if(defined($roi->{ $ncKey }{ $posKey })) {
	        	   		print OUTPUT "$ncKey\t$posKey\t$pileup->{ $ncKey }{ $posKey }\t1\n";
	        	   	} else {
	        	   		print OUTPUT "$ncKey\t$posKey\t$pileup->{ $ncKey }{ $posKey }\t0\n";
	        	   	}	   	
       		}
    	}
		
	} else {
	
		# No ROI file, just output all hits per base in one file.
		for my $ncKey ( sort keys %$pileup ) {
       		for my $posKey ( sort {$a<=>$b} keys %{$pileup->{ $ncKey }} ) {
	        	   	print OUTPUT "$ncKey\t$posKey\t$pileup->{ $ncKey }{ $posKey }\n";
       		}
    	}
    	
	}

close(OUTPUT);

# Print the usage help for this script.
sub usage {
  print "\nUsage: $0 -i <in_file> -o <out_file> -e <region_file>\n
  -i	Input Sedd alignment file, with reads.
  -roi	ROI file in SEDD format. If defined output file will have
        an extra column with 0=outside or 1=inside roi (Optinal).
  -o	Output file.
  -uniq Only use uniqe mapped reads.\n\n";
  exit(1);
}
