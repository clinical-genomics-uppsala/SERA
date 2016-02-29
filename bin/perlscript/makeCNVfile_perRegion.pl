#!/usr/bin/perl

#  
# Author: Elin Falk 2010

use warnings;
use strict;
#use criticism 'harsh';

use FileHandle;
use Carp;
use Getopt::Long;
#use Readonly;

#Globals
my ($meanfile_tumor, $meanfile_normal, $infile, $outfile, $pdffile, $tumorid, $normalid, $min, $average);

# Subroutine prototypes
sub usage;
sub calculateRatio;

# Parse the command line
GetOptions(
	"i=s" => \$infile,
	"mt=s" => \$meanfile_tumor,
	"mn=s" => \$meanfile_normal,
	"o=s" => \$outfile,
	"pdf=s" => \$pdffile,
	"t=s"=> \$tumorid,
	"n=s"=> \$normalid,
	"a=s"=> \$average,
	"min=s" => \$min
	
) || usage();

# Error
if (!$infile || !$meanfile_tumor || !$meanfile_normal || !$outfile || !$min || !$average) { usage(); }

# Constants
my $CHROM1 = 0;
my $POS1 = 1;
my $RATIO1 = 2;
my $CHROM2 = 3;
my $POS2 = 4;
my $RATIO2 = 5;

# Open file handler
my $mean_tumorfh = FileHandle->new($meanfile_tumor) or croak "Couldn't open tumor mean file".$meanfile_tumor."!";
my $mean_normalfh = FileHandle->new($meanfile_normal) or croak "Couldn't open tumor mean file".$meanfile_normal."!";
my $infh = FileHandle->new($infile) or croak "Couldn't open input file".$infile."!";
my $outfh = FileHandle->new(">".$outfile) or croak "Couldn't open output file".$outfile."!";


my ($mean_tumor, $mean_normal);
while(<$mean_tumorfh>) {
	if ($average=~m/mean/i && $_ =~ m/Mean/i)  {
		my @tumorLine = split(/:/, $_);
		$mean_tumor = $tumorLine[1];
	}
	elsif ($average=~m/median/i && $_ =~ m/Median/i)  {
		 my @tumorLine = split(/:/, $_);
		 $mean_tumor = $tumorLine[1];
	}
}
while(<$mean_normalfh>) {
	if ($average=~m/mean/i && $_ =~ m/Mean/i)  {
		my @normalLine = split(/:/, $_);
		$mean_normal = $normalLine[1];
	}
	elsif ($average=~m/median/i && $_ =~ m/Median/i)  {
		my @normalLine = split(/:/, $_);
		$mean_normal = $normalLine[1];
	}
}
#my $kvot;
## Go through the input file
#while (<$infh>) {
#	# Skip comments and empty lines
#	if ($_ =~ m/^#/ || $_ eq "" )  { next; }
#	chomp;
#	# Split on tab
#	my @line = split(/\t/, $_);
#	
#	if ($line[5]>=$min) {
#		if ($mean_normal == 0) {
#			$kvot = (($line[2]/$mean_tumor) - 0)/(($line[2]/$mean_tumor) + 0);
#		}
#		elsif($mean_tumor == 0) {
#			$kvot = (0 - ($line[5]/$mean_normal))/(0 + ($line[5]/$mean_normal));
#		}
#		elsif($line[2] + $line[5] == 0) {
#			$kvot = 0;
#		}
#		else {
#			$kvot = (($line[2]/$mean_tumor) - ($line[5]/$mean_normal))/(($line[2]/$mean_tumor) + ($line[5]/$mean_normal));
#		}
#		
#		print $outfh $line[0]."\t".$line[1]."\t".$kvot."\n"; 
#		}
#	
#	}
#close $outfh;


my $firstPos = 1;
my ($pos, $lastpos, $chrom, $lastchrom, $tot_ratio1, $tot_ratio2, $antalpos);
my %pos2print = ();
# Per region
while (<$infh>) {
	# Skip comments and empty lines
	if ($_ =~ m/^#/ || $_ eq "" )  { next; }
	chomp;
	# Split on tab
	my @line = split(/\t/, $_);
	

	$chrom = $line[$CHROM1];
	$pos = $line[$POS1];
	
	if ($firstPos) {
		if ($line[$RATIO2] >= $min){
			$pos = $lastpos = $line[$POS1];
			$chrom = $lastchrom = $line[$CHROM1];
			$tot_ratio1 = $line[$RATIO1];
			$tot_ratio2 = $line[$RATIO2];
			$antalpos++;
			
			$pos2print{$chrom}{$pos}{'ratio1'} = $line[$RATIO1];
			$pos2print{$chrom}{$pos}{'ratio2'} = $line[$RATIO2];
		} # Ending if ($line[$RATIO2] >= $min)
		else {
			$pos = $lastpos = $line[$POS1];
			$chrom = $lastchrom = $line[$CHROM1];
		}
		$firstPos = 0;
	} # Ending if ($firstPos)
	else {
		if ($line[$RATIO2] >= $min) {
			if ($chrom eq $lastchrom) {
				if ($pos == $lastpos+1) {
					$tot_ratio1 += $line[$RATIO1];
					$tot_ratio2 += $line[$RATIO2];
					$antalpos++;
				
					$pos2print{$chrom}{$pos}{'ratio1'} = $line[$RATIO1];
					$pos2print{$chrom}{$pos}{'ratio2'} = $line[$RATIO2];
				
					$lastpos = $pos;
					$lastchrom = $chrom;
				} # Ending if ($line[$RATIO2] >= $min)
				else {
					if (keys(%pos2print) > 0) {
						my $regionkvot1 = calculateRatio("tumor", $tot_ratio1)/$antalpos;
						my $regionkvot2 = calculateRatio("normal", $tot_ratio2)/$antalpos;
				
						my $kvot;
						if (($regionkvot1+$regionkvot2) == 0) {
							$kvot = 0;
						}
						else {
							$kvot = ($regionkvot1-$regionkvot2)/($regionkvot1+$regionkvot2);
						}
					
						for my $printchrom (keys %pos2print) {
							for my $printpos (sort keys %{$pos2print{$printchrom}}) {
#								print $printchrom."\t".$printpos."\t".$regionkvot1."\t".$printchrom."\t".$printpos."\t".$regionkvot2."\n";
								print {$outfh} $printchrom."\t".$printpos."\t".$kvot."\n";
							}
						}
					}
					$tot_ratio1 = $line[$RATIO1];
					$tot_ratio2 = $line[$RATIO2];
					$antalpos = 1;
					%pos2print = ();
				
					$pos2print{$chrom}{$pos}{'ratio1'} = $line[$RATIO1];
					$pos2print{$chrom}{$pos}{'ratio2'} = $line[$RATIO2];
					
					$lastpos = $pos;
					$lastchrom = $chrom;
				} # Ending else ($line[$RATIO2] >= $min)
			} # Ending if ($chrom eq $lastchrom)
			else {
				if (keys(%pos2print) > 0) {
					my $regionkvot1 = calculateRatio("tumor", $tot_ratio1)/$antalpos;
					my $regionkvot2 = calculateRatio("normal", $tot_ratio2)/$antalpos;
					
					my $kvot;
					if (($regionkvot1+$regionkvot2) == 0) {
						$kvot = 0;
					}
					else {
						$kvot = ($regionkvot1-$regionkvot2)/($regionkvot1+$regionkvot2);
					}	
				
					for my $printchrom (keys %pos2print) {
						for my $printpos (sort keys %{$pos2print{$printchrom}}) {
#							print $printchrom."\t".$printpos."\t".$regionkvot1."\t".$printchrom."\t".$printpos."\t".$regionkvot2."\n";
							print {$outfh} $printchrom."\t".$printpos."\t".$kvot."\n";
						}
					}
				}
			
				$tot_ratio1 = $line[$RATIO1];
				$tot_ratio2 = $line[$RATIO2];
				$antalpos = 1;
				%pos2print = ();
				
				$pos2print{$chrom}{$pos}{'ratio1'} = $line[$RATIO1];
				$pos2print{$chrom}{$pos}{'ratio2'} = $line[$RATIO2];
				
				$lastpos = $pos;
				$lastchrom = $chrom;
			} # Ending else ($chrom eq $lastchrom)	
		} # Ending if ($line[$RATIO2] >= $min)
		else {
			if ($chrom eq $lastchrom) {
				if ($pos == $lastpos+1) {
					$lastpos = $pos;
				}
				else {
					if (keys(%pos2print) > 0) {
						my $regionkvot1 = calculateRatio("tumor", $tot_ratio1)/$antalpos;
						my $regionkvot2 = calculateRatio("normal", $tot_ratio2)/$antalpos;
						
						my $kvot;
						if (($regionkvot1+$regionkvot2) == 0) {
							$kvot = 0;
						}
						else {
							$kvot = ($regionkvot1-$regionkvot2)/($regionkvot1+$regionkvot2);
						}
					
						for my $printchrom (keys %pos2print) {
							for my $printpos (sort keys %{$pos2print{$printchrom}}) {
#								print $printchrom."\t".$printpos."\t".$regionkvot1."\t".$printchrom."\t".$printpos."\t".$regionkvot2."\n";
								print {$outfh} $printchrom."\t".$printpos."\t".$kvot."\n";
							}
						}
					}
					$tot_ratio1 = 0;
					$tot_ratio2 = 0;
					$antalpos = 0;
					%pos2print = ();
					
					$lastpos = $pos;
					$lastchrom = $chrom;
				} # Ending else ($chrom eq $lastchrom)
			} # Ending if ($chrom eq $lastchrom)
			else {
				if (keys(%pos2print) > 0) {
					my $regionkvot1 = calculateRatio("tumor", $tot_ratio1)/$antalpos;
					my $regionkvot2 = calculateRatio("normal", $tot_ratio2)/$antalpos;
					
					my $kvot;
					if (($regionkvot1+$regionkvot2) == 0) {
						$kvot = 0;
					}
					else {
						$kvot = ($regionkvot1-$regionkvot2)/($regionkvot1+$regionkvot2);
					}	
			
					for my $printchrom (keys %pos2print) {
						for my $printpos (sort keys %{$pos2print{$printchrom}}) {
#							print $printchrom."\t".$printpos."\t".$regionkvot1."\t".$printchrom."\t".$printpos."\t".$regionkvot2."\n";
							print {$outfh} $printchrom."\t".$printpos."\t".$kvot."\n";
						}
					}
				}
			
				$tot_ratio1 = 0;
				$tot_ratio2 = 0;
				$antalpos = 0;
				%pos2print = ();
				
				$lastpos = $pos;
				$lastchrom = $chrom;
			} # Ending else ($chrom eq $lastchrom)
		} # Ending else ($line[$RATIO2] >= $min)
	} # Ending else ($firstPos)
} # Ending while (<$infh>)

my $regionkvot1 = calculateRatio("tumor", $tot_ratio1)/$antalpos;
my $regionkvot2 = calculateRatio("normal", $tot_ratio2)/$antalpos;
my $kvot;
if (($regionkvot1+$regionkvot2) == 0) {
	$kvot = 0;
}
else {
	$kvot = ($regionkvot1-$regionkvot2)/($regionkvot1+$regionkvot2);
}		

for my $printchrom (keys %pos2print) {
	for my $printpos (sort keys %{$pos2print{$printchrom}}) {
#		print $printchrom."\t".$printpos."\t".$regionkvot1."\t".$printchrom."\t".$printpos."\t".$regionkvot2."\n";
		print {$outfh} $printchrom."\t".$printpos."\t".$kvot."\n";
	}
}

if ($pdffile) {
	
	if (!$tumorid || !$normalid) { usage(); }
	my $code = "<<EOF
		reset
		set terminal postscript  enh color
		set output '| ps2pdf - ".$pdffile."';
		set grid;
		set title \"".$tumorid." ".$normalid."\";
		set yrange [-1:1];
		set xlabel \"Ampregion Position\";
		set ylabel \"(Sick-Carrier)/(Sick+Carrier)\"
		plot '".$outfile."' using 2:3 notitle with points lc 3 pt 6 ps 0.3;
		EOF";
		
		$code =~s/\<\<EOF//;
		$code =~s/EOF//;
		
		my $gnuplotFH = new FileHandle(">gnuplotCNVScript.gp");
		print $gnuplotFH $code."\n";
		`gnuplot gnuplotCNVScript.gp`;
		`rm gnuplotCNVScript.gp`;
	
}

# Print the usage help for this script.
sub usage {
  print "\nUsage: $0 -i <infile> -o <output file>\n
 -i Infile, two concatemeriased .map files {Tumor_chrom\tTumor_pos\tTumor_hits\tNormal_chrom\tNormal_pos\tNormal_hits
 -min Minimum coverage in Normal sample
 -mn File with mean hits for normal {Mean:??}
 -mt File with mean hits for tumor {Mean:??}
 -o Output file
 -pdf Pdf output file
 -t Tumor id
 -n Normal id\n\n";
  exit 1;
}

sub calculateRatio {
	my ($type, $value) = @_;
	if ($type =~ m/normal/) {
		if ($mean_normal == 0) {
			return 0;
		}
		else {
			return $value/$mean_normal;
		}
	}
	else {
		if($mean_tumor == 0) {
			return 0;
		}
		else {
			return $value/$mean_tumor;
		}	
	}
}
