#!/usr/bin/perl

use strict;
use warnings;
use FileHandle;

sub revComp;
sub usage;

my ($next_arg, $mate1file, $mate2file, $selectorSeq, $prefix, $readLength);

if(scalar(@ARGV) == 0){
    &usage();
}
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-m1")    { $mate1file = shift(@ARGV); }
    elsif($next_arg eq "-m2") { $mate2file = shift(@ARGV); }
    elsif($next_arg eq "-rl") { $readLength = shift(@ARGV); }
    elsif($next_arg eq "-s") { $selectorSeq = shift(@ARGV); }
    elsif($next_arg eq "-p") { $prefix = shift(@ARGV); }
    else { print "Invalid argument: $next_arg"; usage(); }
}

## Error
if((!$mate1file)||(!$mate2file)||(!$selectorSeq)||(!$prefix)||(!$readLength)) { usage(); }

## Main ---

# Handel selection file.
open(FRAGS, "< $selectorSeq") or die "Oops, could not open fragments file: $!";
 my @frags = <FRAGS>;
close(FRAGS);

# Find minimum fragment size.
my $minimalFragSize = 3000000000;
foreach my $frag (@frags) {
	my ($start, $end) = $frag =~ m/.*?\t.*?\t(\d*?)\t(\d*?)\t/;
	#print $start."\t".$end."\t".($end-$start+1)."\n";
	if(($end-$start+1) < $minimalFragSize) { $minimalFragSize = ($end-$start+1); }
}

if ($readLength < $minimalFragSize) { $minimalFragSize = $readLength; }

# Create fragment fasta file.
open(FRAGOUT, "> $prefix.fragments.fasta") or die "Oops, could not open fragments output file: $!";
foreach my $frag (@frags) {
	
	chomp($frag);
	
	my @col = split(/\t/, $frag);
	my $fragId = $col[0];
	my ($pol) = $fragId =~ m/.*?\|.*?\|([+-])\|.*?/;
	my $fragSeq = "";	

	#if(($col[3]-$col[2]+1) > ($minimalFragSize*2)) { 
		$fragSeq = substr($col[4], 0, $minimalFragSize).substr($col[4], (length($col[4])-$minimalFragSize), $minimalFragSize);
	#} else {
	#	$fragSeq = $col[4];	
	#}

	if($pol eq "+") {	
		print FRAGOUT ">".$fragId."\n".$fragSeq."\n";
	} else {
		print FRAGOUT ">".$fragId."\n".revComp($fragSeq)."\n";
	}
}
close(FRAGOUT);

# Save RAM
@frags = ();

# Hash mate 2 reads.
my $mate2hash = {};
open(M2, "zcat $mate2file |") or die "Oops, could not open mate2 input file: $!";
while(<M2>) {
	
	chomp;
	if(m/^@/) {
	  	
		#my ($readId) = $_ =~ m/^(@.*?)\/2/;
		my ($readId) = $_ =~ m/^(@.*):3:([AGCT]{9}):([AGCT]{10})/;
		my $readSeq = <M2>;
		chomp($readSeq);
		#print "$readId\t$readSeq\t$minimalFragSize\n";		
		$mate2hash->{$readId}->{'part'} = revComp(substr($readSeq,0,$minimalFragSize));
		$mate2hash->{$readId}->{'whole'} = $readSeq;
	}
}
close(M2);

# Create mates fasta file.
open(READSOUT, "> $prefix.reads.fasta") or die "Oops, could not open reads output file: $!";
open(M1, "zcat $mate1file |") or die "Oops, could not open mate2 input file: $!";
while(<M1>) {

	chomp;
	if(m/^@/) {
	  	
		#my ($readId) = $_ =~ m/^(@.*?)\/1/;
		my ($readId) = $_ =~ m/^(@.*):1:([AGCT]{9}):([AGCT]{10})/;
   		my $readSeq = <M1>;
		chomp($readSeq);
#		print $readId."\n";
		if($mate2hash->{$readId}) {
			print READSOUT ">".$readId."_bothMates#".$readSeq."#".$mate2hash->{$readId}->{'whole'}."\n".substr($readSeq,0,$minimalFragSize).$mate2hash->{$readId}->{'part'}."\n";
		} else {
			print "No mate found for $readId...\n";
		}
	}
}
close(M1);
close(READSOUT);

#----------

# Get reverse comp.
sub revComp {
	
	my @inSeq = split('', $_[0]);
	my $outSeq="";
	
	foreach my $base (@inSeq) {
		
		$base =~ tr/ACGTacgt/TGCAtgca/;
		$outSeq = $outSeq.$base;
	}
	
	return reverse($outSeq);
	#return $outSeq;
}

# Print the usage help for this script.
sub usage {
  print "\nUsage: $0 -m1 <mate1> -m2 <mate2> -rl <read_length> -s <fragments> -p <prefix>\n
  -m1	Mate 1 input file in fastq format.
  -m2	Mate 2 input file in fastq format.
  -rl	Read length (bp).
  -s	Selected fragments with sequence.
  -p	Prefix for output files.\n\n";
  exit(1);
}
