#!/usr/bin/perl -w
#

use strict;
use warnings;

# Subroutine prototypes
sub usage;

my ($next_arg, $normalFile, $tumorFile, $nMinRD, $outputFile, $selectionFile);

if(scalar(@ARGV) == 0){ usage(); }
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-n")    { $normalFile = shift(@ARGV); }
    elsif($next_arg eq "-t") { $tumorFile = shift(@ARGV); }
    elsif($next_arg eq "-nMinRD") { $nMinRD = shift(@ARGV); }
    elsif($next_arg eq "-o") { $outputFile = shift(@ARGV); }
    elsif($next_arg eq "-s") { $selectionFile = shift(@ARGV); }
    else { print "Invalid argument: $next_arg"; usage(); }
}
# Defult and errors
if ( !$tumorFile || !$normalFile || !$selectionFile) { &usage(); }
if (!$outputFile) { $outputFile = "/dev/stdout"; }
if (!$nMinRD) { $nMinRD = 0; }

# MAIN ----------

# Hashing data
my %fragmentHash;
my $normalTotalRD = 0;
my $tumorTotalRD = 0;
my %selectors;

open (SELECTION, "< $selectionFile") or die "Oops, could not open selection file: $selectionFile !";
while (<SELECTION>) {
	if($_=~m/^#/ || $_ eq "" || $_ =~m/^FragmentID/) { next; }
	
	# Remove new line
	chomp; 
	# Split line on tab
	my @sLine = split(/\t/, $_);
	$selectors{$sLine[0]} = $sLine[1]."\t".$sLine[2]."\t".$sLine[3]."\t".$sLine[4];
}
close(SELECTION);

open (NORMAL, "< $normalFile") or die "Oops, could not open normal file: $normalFile !";
while (<NORMAL>) {
	if($_=~m/^#/ || $_ eq "" || $_ =~m/^FragmentID/) { next; }
	
	# Remove new line
	chomp; 
	# Split line on tab
	my @nLine = split(/\t/, $_);
	# If normal min RD is given only consider fragments with read depth >= normal min RD
	if($nLine[1]>=$nMinRD) {
		$fragmentHash{$nLine[0]}{'normal'} = $nLine[1];
		$normalTotalRD += $nLine[1];
	}
}
close(NORMAL);

open (TUMOR, "< $tumorFile") or die "Oops, could not open tumor file: $tumorFile !";
while (<TUMOR>) {
	if($_=~m/^#/ || $_ eq "") { next; }
	
	# Remove new line
	chomp; 
	# Split line on tab
	my @tLine = split(/\t/, $_);
	# If the fragment exists in the hash, add info for tumor
	if ($fragmentHash{$tLine[0]}) {
		$fragmentHash{$tLine[0]}{'tumor'} = $tLine[1];
		$tumorTotalRD += $tLine[1];
	}
}
close(TUMOR);

open (OUTPUT, "> $outputFile") or die "Oops, could not open output file: $outputFile !";
for my $frag (sort keys %fragmentHash) {
	if ($fragmentHash{$frag}{'tumor'} == 0) {
		print OUTPUT $frag."\t".$selectors{$frag}."\t-10\n";
	}
	else {
		my $tumorNormRD = $fragmentHash{$frag}{'tumor'}/$tumorTotalRD;
		my $normalNormRD = $fragmentHash{$frag}{'normal'}/$normalTotalRD;
		print OUTPUT $frag."\t".$selectors{$frag}."\t".log($tumorNormRD/$normalNormRD)/log(2)."\n";
	}
	
}

close(OUTPUT);

# ---------------

# Sub show how to run this script.
sub usage {
  print "
************************************************************************************************
  This script calculates the log2ratio(T/N) and adds the information about the selector 
  it comes from. The read depth for a selector is normalized against the total read depth
  of the sample.
************************************************************************************************ 
  \nUsage: $0\n 
 -n      Normal file {fragID [tab] Read depth}
 -t      Tumor file {fragID [tab] Read depth}
 -nMinRD Normal min read depth [optional, default 0]
 -s      Selection file {fragID [tab] chr [tab] start [tab] end [tab] strand}
 -o      Output file [optional, default /dev/stdout.\n\n";   
  exit(1);
}
