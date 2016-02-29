#!/usr/bin/perl -w
#
# By Elin Falk Sörqvist, 2011
# Calculate numbers of reads covering each base on- and off target 
# from a pile up file
#

use strict;
use warnings;
use FileHandle;

sub usage;

my %bases; # Contains all reads.
my ($next_arg, $inFile, $regionFile, $outFile, $offTarget_hits_file, $chr2nc);

if(scalar(@ARGV) == 0){
    &usage();
}
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-i")    { $inFile = shift(@ARGV); }
    elsif($next_arg eq "-o") { $outFile = shift(@ARGV); }
    elsif($next_arg eq "-r") { $regionFile = shift(@ARGV); }
    elsif ($next_arg eq "-chr2nc") {$chr2nc = shift(@ARGV); }
    elsif($next_arg eq "-off") { $offTarget_hits_file=shift(@ARGV); }  
    else { print "Invalid argument: $next_arg"; usage(); }
}

## Error
if( !$inFile || !$regionFile || !$chr2nc) { usage(); }

# If no output file write to stdout.
if (!$outFile) { $outFile="/dev/stdout"; }

open(CHR2NC, "< $chr2nc") or die "Oops, could not open chr2nc file: $!";
my %chrhash;
my %nchash;
while (<CHR2NC>) {
	if ($_ =~ m/^#/ ||$_ eq "") { next; }
	chomp;
	my @c_columns = split(/\t/, $_);
	$chrhash{$c_columns[0]} = $c_columns[1];
	$nchash{$c_columns[1]} = $c_columns[0];
}

my $c=0;
## Hashing input reads file.
open(INPUT, "< $inFile") or die "Oops, could not open input file: $!";
while(<INPUT>) {
	if ($_ =~ m/^#/ ||$_ eq "") { next; }
	
	chomp;
	my @columns = split(/\t/, $_);
	my $chr = "";
	if ($columns[0] =~ m/#/) {
                if ($chrhash{$columns[0]}) {
                        $chr = $chrhash{$columns[0]};
           }
        }
	elsif ($columns[0] =~ m/^NC_0{4,5}[1-9]{1,2}\.[0-9]{1,2}$/) {
		$chr = $nchash{$columns[0]};
	}
	else {
		$chr = $columns[0];
	}
	my $pos = $columns[1];
	my $hits = $columns[3]; 

	$c += $hits;
	if ($bases{$chr}{$pos}) {
		$bases{$chr}{$pos} += $hits;
	} 
	else {
		$bases{$chr}{$pos} = $hits;
	}
#	print $chr."\t".$pos."\t".$bases{$chr}{$pos}."\n";
}
print "Total number of bases: ".$c."\n";
close(INPUT);

## Create output
open(REGION, "< $regionFile") or die "Oops, could not open region file: $!";
open(OUTPUT, "> $outFile") or die "Oops, could not open output file: $!";

while(<REGION>) {
	if ($_ =~ m/^#/ ||$_ eq "") { next; }
	
	chomp;
	my @r_columns = split(/\t/, $_);
	my $r_ncNumber = $r_columns[1];
	my $r_start = $r_columns[2];
	my $r_end = $r_columns[3]; 
	
	my $r_chr = "";
	my $ok = 0;

	if ($r_ncNumber =~ m/#/) {
			$r_chr = $chrhash{$r_ncNumber};
			$ok=1;
	}
	elsif ($r_ncNumber =~ m/^NC_0{4,5}[1-9]{1,2}\.[0-9]{1,2}$/) {
		if ($nchash{$r_ncNumber}) {
                        $r_chr = $nchash{$r_ncNumber};
                        $ok=1;
                }
	}
	else {
		$r_chr = $r_ncNumber;
		$ok = 1;
	}

	if($ok == 1) {			
		# Check each base in region file.
		for (my $k = $r_start; $k < ($r_end+1); $k++) {
#			print $r_chr."\t".$k."\t".$bases{$r_chr}{$k}."\n";	
			# Do we have coverage?
			if ($bases{$r_chr}{$k}) {
				print OUTPUT $r_chr."\t".$k."\t".$bases{$r_chr}{$k}."\n";
				delete $bases{$r_chr}{$k};
			}
			else {
				print OUTPUT $r_chr."\t".$k."\t0\n";
			}
		}
	}
	else {
		print "# Unknown chromosome: ".$r_ncNumber."\n";
	}
}

close(OUTPUT);
close(REGION);

if ($offTarget_hits_file) {
	open (OFF_HITS, "> $offTarget_hits_file") or die "Oops, could not open off target hits file $offTarget_hits_file: $!";

	foreach my $chr (keys %bases) {
		foreach my $pos (keys %{$bases{$chr}}) {
			print OFF_HITS $chr."\t".$pos."\t".$bases{$chr}{$pos}."\n";
			delete $bases{$chr}{$pos};
		}
	}
close(OFF_HITS);
}

# Print the usage help for this script.
sub usage {
  print "\nUsage: $0 -i <in_file> -o <out_file> -e <region_file>\n
  -i      Input file, in pile up format.
  -o      Output file (stdout is default if not set).
  -r      Region file e.g. amproi.
  -chr2nc File with conversion chr to nc-number (eg. chr1	NC_000001.10)
  -off    If you want the off target reads saved add an off target output file\n\n";
  exit(1);
}
