#!/usr/bin/perl

# This script takes a hit per base file and an ampregion file and generates a file with values
# for a Nusbaum plot (with readDepth normalized against ampregion mean read depth)
# and/or cumulative (Stenberg) plot and/or Ericsson plot, depending on the flags.
# The hits in the Nusbaum and Ericsson plot are normalized against the mean read 
# depth in the whole ampregion.
#

use strict;
use FileHandle;

# Subroutine prototypes
sub usage;

my ($next_arg, $infile, $ampregionFile, $cumulativeFile, $ericssonFile, $frequencyFile, $nusbaumFile);

if(scalar(@ARGV) == 0){
    &usage();
}
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-i")    { $infile = shift(@ARGV); }
    elsif($next_arg eq "-a") { $ampregionFile = shift(@ARGV); }
    elsif($next_arg eq "-c") { $cumulativeFile = shift(@ARGV); }
    elsif($next_arg eq "-e") { $ericssonFile = shift(@ARGV); }
    elsif($next_arg eq "-f") { $frequencyFile = shift(@ARGV); }
    elsif($next_arg eq "-n") { $nusbaumFile = shift(@ARGV); }
    else { print "Invalid argument: $next_arg"; usage(); }
}


if ((!$infile)||(!$ampregionFile)||!($nusbaumFile||$ericssonFile||$cumulativeFile||$frequencyFile)) { &usage(); }

# Create output filehandle
my ($outNus, $outCum, $outEric, $outFreq); 


#Reading in the whole ampregion file
my %ampBases = ();
my $hits = 0;
my $counter = 0;

# Put all ampregion bases in a hash
open(ampFile, $ampregionFile);
while (<ampFile>) {
	if (!($_ =~m/^#/)) {
		chomp;
		my @ampLine = split (/\t/, $_);
		$ampBases{$ampLine[0]}{$ampLine[1]} = $ampLine[2];
		$hits += $ampLine[2];
		$counter++;
	}	
}
close(ampFile);


# Calculate mean hit in ampregion
my $meanHit;
if ($counter == 0) {
	$meanHit = 0;
}
else {
	$meanHit = $hits/$counter;
}


my @hitsInRegion;
my $totalBases = 0;
my $inFH_region = new FileHandle($infile);
my $i = 0;
print "Hashing hits\n";
while (<$inFH_region>) {
	if (!($_ =~m/^#/)) {
		chomp;
		my @regionLine = split (/\t/, $_);

		# Convert sedd NC-number into chr.
#		if($regionLine[1] eq "NC_012920.1") {
#			$regionLine[1] = "chrM";
#		} else {		
#			$regionLine[1] =~ s/NC_0{4,5}(\d{1,2})\..*/chr$1/; 
#			if($regionLine[1] eq "chr23") { $regionLine[1] = "chrX"; }
#			if($regionLine[1] eq "chr24") { $regionLine[1] = "chrY"; }
#		}

		# Hash it up.
		$totalBases += ($regionLine[3]-$regionLine[2]+1);
		for my $ampPos (keys %{$ampBases{$regionLine[1]}}) {
			if ($ampPos >= $regionLine[2] && $ampPos <= $regionLine[3]) {
				$hitsInRegion[$i] = $ampBases{$regionLine[1]}{$ampPos};
				$i++;
			}
		}
	}
}

@hitsInRegion = sort(@hitsInRegion);

#If a nusbaum plot is wanted, run the code below
if(defined($nusbaumFile) || defined($frequencyFile)) {
	#Normalize the hits for the Nusbaum plot
	my @normBases;
	for (my $k=0; $k<$#hitsInRegion; $k++) {
		if (!($hitsInRegion[$k] =~m/^#/)) {
			my $normHits = $hitsInRegion[$k]/$meanHit;
			$normBases[$k] = $normHits;
		}
	}
	
	@normBases = sort{$b<=>$a} @normBases;
	
	if(defined($nusbaumFile)) {
		# Nothing underneath 1 hit should be reported since it's uninformative
		my $threshold;
		if ($meanHit == 0) {
			$threshold = 1;
		}
		else {
			$threshold = 1/$meanHit;
		}
		$outNus = new FileHandle(">$nusbaumFile");
		print "Creating nusbaum plot\n";
		#Put the highest hit no. as the last read hit. Calculate the occurences of each hit
		#and print it to file.
				
		my $last = $normBases[0];
		my $count = 1;
		my $freq;
		for ($i=1; $i<=$#normBases; $i++) {
			
			#If the no. of hits on this line is the same as the last line increase
			#the counter with 1.
			if ($normBases[$i] == $last) {
				$count++;			
			}
			#If the hits on this line is lower than the hits on the last line
			#print the number of hits and its frequency. Then set the hits on this
			#line to $last
			elsif ($normBases[$i] < $last) {
				$freq = ($count/$totalBases)*100;
				if ($last >= $threshold) {
					print $outNus "$last\t$freq\n";
				}
				$last = $normBases[$i];
				$count++;
			}
			else {
				print "Bad sorted $normBases[$i] > $last\n";
			}
		}
		
		$freq = ($count/$totalBases)*100;
		if ($last >= $threshold) {
			print $outNus "$last\t$freq\n";
		}
		
		$outNus->close;
	}
	
	if (defined($frequencyFile)) {

		$outFreq = new FileHandle(">$frequencyFile");
		print "Creating frequency plot\n";
		#Put the highest hit no. as the last read hit. Calculate the occurences of each hit
		#and print it to file.
		
		my $last = $normBases[0];
		my $counter = 1;
		my $freq;
		for ($i=0; $i<=$#normBases; $i++) {
					
			#If the no. of hits on this line is the same as the last line increase
			#the counter with 1.
			if ($normBases[$i] == $last) {
				$counter++;			
			}
			#If the hits on this line is lower than the hits on the last line
			#print the number of hits and its frequency. Then set the hits on this
			#line to $last
			elsif ($normBases[$i] < $last) {
				$freq = ($counter/$totalBases)*100;
				print $outFreq "$last\t$freq\n";
				$last = $normBases[$i];
				$counter = 1;
			}
			else {
				print "Bad sorted $normBases[$i] > $last\n";
			}
		}
		$freq = ($counter/$totalBases)*100;
		print $outFreq "$last\t$freq\n";
		$outFreq->close;
	}
	
	
}

# If an ericsson plot is wanted, run the code below
if(defined($ericssonFile)) {
	#print $meanHit."\n";
	$outEric = new FileHandle(">$ericssonFile");
	print "Creating ericsson plot\n";
	# Normalize the hits for the Ericsson plot
	my @normBases;
	for (my $k=0; $k<=$#hitsInRegion; $k++) {
		if (!($hitsInRegion[$k] =~m/^#/)) {
			
			my $normHits;
			if ($hitsInRegion[$k] != 0) {
				$normHits = $meanHit/$hitsInRegion[$k];
			}
			else {
				$normHits = 0;	
			}
			$normBases[$k] = $normHits;
		}
	}
	
	@normBases = sort {$a<=>$b} @normBases;
	
	#Put the highest hit no. as the last read hit. Calculate the occurences of each hit
	#and print it to file.
	my $last = $normBases[0];
	chomp($last);
	my $counter = 0;
	my $freq;
	for ($i=0; $i<=$#normBases; $i++) {
		#If the no. of hits on this line is the same as the last line increase
		#the counter with 1.
		if ($normBases[$i] == $last) {
			$counter++;			
		}
		#If the hits on this line is lower than the hits on the last line
		#print the number of hits and its frequency. Then set the hits on this
		#line to $last
		elsif ($normBases[$i] > $last) {
			$freq = ($counter/ $totalBases)*100;
			if ($last <= $meanHit) {
				print $outEric "$last\t$freq\n";
			}
			$last = $normBases[$i];
			$counter++;
		}
		else {
			print "Bad sorted ".$normBases[$i]." < ".$last."\n";
		}
	}
#	print $last."\n";
	$freq = ($counter/$totalBases)*100;
	if ($last <= $meanHit) {
		print $outEric "$last\t$freq\n";
	}
		
$outEric->close;
}

#If a cumulative plot is wanted, run the code below
if (defined($cumulativeFile)) {
	$outCum = new FileHandle(">$cumulativeFile");
	print "Creating stenberg plot\n";
	my @hitsInRegionCum = sort {$b<=>$a} @hitsInRegion;
	
	#Set the first hit as the last read
	my $last = $hitsInRegionCum[0];
	my $counter = 0;
	my $freq;
	#From second base and ahead count how many times a particular no of hits occur.
	#Make it in to percentage and print to file
	for ($i=1; $i<=$totalBases ; $i++) {
		#If the no. of hits on this line is the same as the last line increase
		#the counter with 1.
		if ($hitsInRegionCum[$i] == $last) {
			$counter++;
		}
		#If the hits on this line is lower than the hits on the last line
		#print the number of hits and its frequency. Then set the hits on this
		#line to $last
		elsif ($hitsInRegionCum[$i] < $last) {
			$freq = ($counter/$totalBases)*100;
			print $outCum "$last\t$freq\n";
			$last = $hitsInRegionCum[$i];
			$counter++;
			
		}
		else {
			print "Bad sorted $hitsInRegionCum[$i] > $last\n";
		}
	}
	
	$freq = ($counter/$totalBases)*100;
	print $outCum "$last\t$freq\n";
	$outCum->close;
}

# Print the usage help for this script.
sub usage {
  print "\nhitsPerBase2ampregionfreqValue.pl Usage: $0 \n 
 -a Ampregion map file:
 	{Chrom	Base pos	Hits}
 -i Region infile:
 	{Id	Chrom	Start	End	Strand} 
 -c If a cumulative output file is wanted:
 	{Hits	Freq}
 -e If an ericsson output file is wanted:
 	{Hits	Freq} 	 
 -f If a frequnecy output file is wanted:
 	{Hits	Freq} 	 
 -n If a nusbaum output file is wanted with read depth normalization against ampregion and \"normal\" region on y-axi region on y-axis:
 	{Hits	Freq}\n\n";
  exit(1);
}
