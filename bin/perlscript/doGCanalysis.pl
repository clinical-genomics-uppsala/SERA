#!/usr/bin/perl -w
#
# This program takes a tabulated file, and
# calculates %GC and GC-clustering. 
#
# By Magnus Isaksson 2009
#

use strict;
use warnings;

# Subroutine prototypes
sub usage;
sub calculateGC;
sub calculateGCclust;

my ($next_arg, $infile, $outfile, $seqColumn, $clusterThSize, $clusterThGC, $probearms, @probearmlength);

if(scalar(@ARGV) == 0){ usage(); }
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-i")    { $infile = shift(@ARGV); }
    elsif($next_arg eq "-o") { $outfile = shift(@ARGV); }
    elsif($next_arg eq "-c") { $seqColumn = shift(@ARGV); }
    elsif($next_arg eq "-s") { $clusterThSize = shift(@ARGV); }
    elsif($next_arg eq "-n") { $clusterThGC = shift(@ARGV); }
    elsif($next_arg eq "-p") { $probearms = shift(@ARGV); }
    else { print "Invalid argument: $next_arg"; usage(); }
}
# Defult and errors
if ((!$infile)||(!$seqColumn)) { &usage(); }
if (!$outfile) {
	$outfile = $infile;
	$outfile =~ s/(.*?)\..*/$1\.output.txt/;
}
if (!$clusterThSize) { $clusterThSize=2; }
if (!$clusterThGC) { $clusterThGC=3; }
if (!$probearms) {
	$probearmlength[0] = 25;
	$probearmlength[1] = 25;
}
else {
	@probearmlength = split(/,/, $probearms);
}

# MAIN ----------

# Open files.
open (INFILE, "< $infile") or die "Oops, could not open input file: $!";
open (OUTFILE, "> $outfile") or die "Oops, could not open output file: $!";


	while (<INFILE>) {
		
		chomp; # Remove new line.
		my @columns = split(/\t/, $_);
		my $seq = $columns[$seqColumn-1];
		
		# Calculate gc-values.
		my $gc = calculateGC($seq);
		my $gc_clust = calculateGCclust($seq);
		my $gc_probe = calculateProbeGC($seq, $probearmlength[0], $probearmlength[1]);
		
		print OUTFILE $_."\t".$gc."\t".$gc_clust."\t".$gc_probe."\n";
		
	}

close(OUTFILE);
close(INFILE);
# END MAIN ------

# Sub calculates GC-cluster
sub calculateGCclust {
	
	my $seq=$_[0];
	my $length = length($seq);
	my $clusted = 0;
	my $th = $clusterThSize;
	
	for (my $i=0; $i < $length; $i++) {
		
		my $base = substr($seq,$i,1);
		
		# Is base G or C?
		if ($base =~ m/[gc]/ig) {
			
			my $gc = 0;
			
			# Fix subsequence
			my $upStream = $i-$th;
			if ($upStream < 0) { $upStream=0; }
			my $downStream = $i+$th;
			if ($downStream > (length($seq)-1)) { $downStream=length($seq)-1; }
			
			my $seqToAn = substr($seq,$upStream,$downStream-$upStream+1);
			
			# Count G and C in subsequence
			++$gc while $seqToAn =~ m/[gc]/ig; 
			if ( ($gc-1) >= $clusterThGC ) { $clusted++; }
	
			#print $upStream."-".$downStream." ".$seqToAn." GC:".$gc."\n";
			
		}
		
	}
	
	return (($clusted/$length)*100);
}

# Sub calculate %GC
sub calculateGC {
   
   my $gc=0;
   my $seq=$_[0];
   
   my $length = length($seq);
   ++$gc while $seq =~ m/[gc]/ig;

   return (($gc/$length)*100);	
}

# Sub calculate %GC in probearms
sub calculateProbeGC {
   my $seq=$_[0];
   my $probe5 = $_[1];
   my $probe3 = $_[2];
   
   my $gc=0;
   
   my $probeseq = substr($seq, 0, $probe5);
   
   $probeseq .= substr($seq, -$probe3);
   
      
   my $length = length($probeseq);
   ++$gc while $probeseq =~ m/[gc]/ig;

   return (($gc/$length)*100);
}

# Sub show how to run this script.
sub usage {
  print "\nUsage: $0 -i <in-file> -r <out-file> -c <seq_column> -s <search_size> -n <GC_th>\n 
 -i  	Input file in tabular format.
 -o	    Output file with three new columns GC %, GC-clustering, probearm GC.
 -c  	Column with sequence.
 -s     Search window on each side of a found G or C.
 -n     Number of G or C found in the whole search window before count.
 -p     Length of the probearms to analyse, 5'arm,3'arm eg. 25,25\n\n";   
  exit(1);
}