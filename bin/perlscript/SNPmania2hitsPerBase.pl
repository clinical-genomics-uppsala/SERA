#!/usr/bin/perl -w
#
# By Elin Falk Sörqvist, 2011
# Calculate numbers of reads covering each base on- and off target 
# from a SNPmania file
#

use strict;
use warnings;
use FileHandle;

sub usage;

my %bases; # Contains all reads.
my ($next_arg, $inFile, $regionFile, $outFile, $offTarget_hits_file);
my $unique_reads_only = 0;

if(scalar(@ARGV) == 0){
    &usage();
}
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-i")    { $inFile = shift(@ARGV); }
    elsif($next_arg eq "-o") { $outFile = shift(@ARGV); }
    elsif($next_arg eq "-r") { $regionFile = shift(@ARGV); }
	elsif($next_arg eq "-off") { $offTarget_hits_file=shift(@ARGV); }  
    else { print "Invalid argument: $next_arg"; usage(); }
}

## Error
if((!$inFile)||(!$regionFile)) { usage(); }

# If no output file write to stdout.
if (!$outFile) { $outFile="/dev/stdout"; }

my %regions;
open(REGION, "< $regionFile") or die "Oops, could not open input file: $!";
while(<REGION>) {
	if ($_ =~ m/^#/ ||$_ eq "") { next; }
	
	chomp;
	my @r_columns = split(/\t/, $_);
	my $r_ncNumber = $r_columns[1];
	my $r_start = $r_columns[2];
	my $r_end = $r_columns[3]; 
	
	$regions{$r_ncNumber}{$r_start} = $r_end;
}
close(REGION);

## Hashing input reads file.
open(OUTPUT, "> $outFile") or die "Oops, could not open output file: $!";
if ($offTarget_hits_file) {
	open (OFF_HITS, "> $offTarget_hits_file") or die "Oops, could not open off target hits file $offTarget_hits_file: $!";
}
open(INPUT, "< $inFile") or die "Oops, could not open input file: $!";
while(<INPUT>) {
	if ($_ =~ m/^#/ ||$_ eq "") { next; }
	
	chomp;
	my @columns = split(/\t/, $_);
	my $ncNumber = $columns[3];
	my $pos = $columns[4];
	my $hits = $columns[0]; 

	
	foreach my $start (keys %{$regions{$ncNumber}}) {
		if($ncNumber eq $ncNumber && $pos>=$start && $pos<=$regions{$ncNumber}{$start}) {
			print OUTPUT $ncNumber."\t".$pos."\t".$hits."\n";
		}
		else {
			if($offTarget_hits_file) {
				print OFF_HITS $ncNumber."\t".$pos."\t".$bases{$ncNumber}{$pos}."\n";
			}
		}
	}
}
close(INPUT);
close(OUTPUT);
close(OFF_HITS);

# Print the usage help for this script.
sub usage {
  print "\nUsage: $0 -i <in_file> -o <out_file> -e <region_file>\n
  -i    Input file, with reads.
  -o    Output file (stdout is default if not set).
  -r    Region file e.g. amproi.
  -off  If you want the off target reads saved add an off target output file\n\n";
  exit(1);
}