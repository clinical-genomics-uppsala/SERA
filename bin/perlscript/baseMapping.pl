#!/usr/bin/perl -w
#
# By Magnus Isaksson 2011
# Calculate numbers of reads covering each base 
# in a defined region.
#

use strict;
use warnings;
use FileHandle;

sub usage;

my $bases = {}; # Contains all pileups.
my ($next_arg, $inFile, $regionFile, $outFile);

if(scalar(@ARGV) == 0){
    &usage();
}
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-i")    { $inFile = shift(@ARGV); }
    elsif($next_arg eq "-o") { $outFile = shift(@ARGV); }
    elsif($next_arg eq "-r") { $regionFile = shift(@ARGV); }
    else { print "Invalid argument: $next_arg"; usage(); }
}

## Error
if((!$inFile)||(!$regionFile)) { usage(); }

# If no output file write to stdout.
if (!$outFile) { $outFile="/dev/stdout"; }

## Hashing input pileup file.
open(INPUT, "< $inFile") or die "Oops, could not open input file: $!";
while(<INPUT>) {
	if (/^chr/) {	
		chomp;
#		my ($chr, $pos, $depthSeq) = $_ =~ m/^(.*?)\t(\d*?)\t.*?\t(.*?)\t.*?/;

		# adjust whether sequence is chr1 (after corrsel) or chr1#nc#etc if cutadapt, /Linus
		my($chr,$pos,undef,$depthSeq) = split(/\t/,$_);
		if($chr =~ m/#/) {
			my @arr=split(/#/,$chr);
			$chr = $arr[0];
		}

		$bases->{$chr}->{$pos} = $depthSeq;
	}
}
close(INPUT);

## Create output
open(REGION, "< $regionFile") or die "Oops, could not open input file: $!";
open(OUTPUT, "> $outFile") or die "Oops, could not open output file: $!";
while(<REGION>) {
	if (!/#/) {
		
		chomp;
		my ($qChr, $start, $end) = $_ =~ m/.*?\t(NC_.*?)\t(\d+)\t(\d+)/;
		
		# Convert sedd NC-number into chr.
		if($qChr eq "NC_012920.1") {
			$qChr = "chrM";
		} else {		
			$qChr =~ s/NC_0{4,5}(\d{1,2})\..*/chr$1/; 
			if($qChr eq "chr23") { $qChr = "chrX"; }
			if($qChr eq "chr24") { $qChr = "chrY"; }
		}
	
		# Check each base in region file.
		for (my $i = $start; $i < ($end+1); $i++) {
			
			# Do we have coverage?
			if ($bases->{$qChr}->{$i}) {
				print OUTPUT $qChr."\t".$i."\t".$bases->{$qChr}->{$i}."\n";
			} else {
				print OUTPUT $qChr."\t".$i."\t0\n";
			}
			
		}
	}
}
close(OUTPUT);
close(REGION);

# Print the usage help for this script.
sub usage {
  print "\nUsage: $0 -i <in_file> -o <out_file> -e <region_file>\n
  -i	Input file, samtools pileup.
  -o	Output file (stdout is default if not set).
  -r	Region file e.g. amproi.\n\n";
  exit(1);
}
