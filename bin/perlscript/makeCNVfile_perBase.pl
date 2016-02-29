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
	"min=i" => \$min
	
) || usage();

# Error
if (!$infile || !$meanfile_tumor || !$meanfile_normal || !$outfile || !$min || !$average) { usage(); }

# Constants
my $SELECTOR => 0;
my $ATTRIBUTE => 1;

# Open file handler
my $mean_tumorfh = FileHandle->new($meanfile_tumor) or croak "Couldn't open tumor mean file".$meanfile_tumor."!";
my $mean_normalfh = FileHandle->new($meanfile_normal) or croak "Couldn't open tumor mean file".$meanfile_normal."!";
my $infh = FileHandle->new($infile) or croak "Couldn't open input file".$infile."!";
my $outfh = FileHandle->new(">".$outfile) or croak "Couldn't open output file".$outfile."!";


print {$outfh} "# Selector_id\tNuber_of_hits\tSelector_length\n";
my %selectors;
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
my $kvot;
# Go through the input file
while (<$infh>) {
	# Skip comments and empty lines
	if ($_ =~ m/^#/ || $_ eq "" )  { next; }
	chomp;
	# Split on tab
	my @line = split(/\t/, $_);
	
	if ($line[5]>=$min) {
		if ($mean_normal == 0) {
			$kvot = (($line[2]/$mean_tumor) - 0)/(($line[2]/$mean_tumor) + 0);
		}
		elsif($mean_tumor == 0) {
			$kvot = (0 - ($line[5]/$mean_normal))/(0 + ($line[5]/$mean_normal));
		}
		elsif($line[2] + $line[5] == 0) {
			$kvot = 0;
		}
		else {
			$kvot = (($line[2]/$mean_tumor) - ($line[5]/$mean_normal))/(($line[2]/$mean_tumor) + ($line[5]/$mean_normal));
		}
		
		print $outfh $line[0]."\t".$line[1]."\t".$kvot."\n"; 
		}
	
	}
close $outfh;

if ($pdffile) {
	
	if (!$tumorid || !$normalid) { usage(); }
	my $code = "<<EOF
		reset
		set terminal postscript  enh color
		set output '| ps2pdf - ".$pdffile."';
		set grid;
		set title \"".$tumorid." ".$normalid."\";
		set yrange [-1:1];
		set xlabel \"Pseudo Position\";
		set ylabel \"(Tumor-Normal)/(Tumor+Normal);\"
		plot '".$outfile."' using 3 notitle with points lc 3 pt 6 ps 0.3;
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
  print "\ngetStartEndPos.pl Usage: $0 -i <infile> -o <output file>\n
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

