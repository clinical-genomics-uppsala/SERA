#!/usr/bin/perl -w
#
# This program takes results from
# jSNPmania and adds SNP info.
#
# By Magnus Isaksson 2010
#

use strict;
use warnings;

# Hold snps from ref. SNP file.
my $snps = {};

# Subroutine prototypes
sub usage;
sub loadSNPs;

my ($next_arg, $infile, $outfile);
my %snpRef;
my $ncToChr = "false";

if(scalar(@ARGV) == 0){ usage(); }
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-i")    { $infile = shift(@ARGV); }
    elsif($next_arg eq "-o") { $outfile = shift(@ARGV); }
    elsif($next_arg eq "-snp_ref") { 
    	($snpRef{'file'}, $snpRef{'chr_col'}, $snpRef{'pos_col'}, $snpRef{'genotype'}) = split(/,/, shift(@ARGV));
    	if( (!$snpRef{'file'}) || (!$snpRef{'chr_col'}) || (!$snpRef{'pos_col'}) || (!$snpRef{'genotype'}) ) 
    	{ usage(); }
    }
    elsif($next_arg eq "-nc2chr") { $ncToChr = "true"; }
    else { print "Invalid argument: $next_arg"; usage(); }
}

if (!$infile) { &usage(); }
if (!$outfile) { $outfile="/dev/sdtout"; }

# MAIN ---------
	
	# Load snps
	print "Hashing reference SNPs ...\n";
	loadSNPs();
	
	# Open files
	print "Processing jSNPmania result file ...\n";
	open(INPUT, "< $infile") or die "Oops, could not open input file: $!";
	open(OUTPUT, "> $outfile") or die "Oops, could not open output file: $!";
	
	while(<INPUT>) {
		
		if(!(m/^#/)) {
		
			chomp;
			my @columns = split(/\t/, $_);
						
			my $chr = $columns[3];
			my $pos = $columns[4];
			
			if($ncToChr eq "true") { 
				$chr =~ s/NC_0{4,5}(\d{1,2})\..*/chr$1/; 
			    if($chr eq "chr23") { $chr = "chrX"; }
			    if($chr eq "chr24") { $chr = "chrY"; }
			}
			
			# Do we have a SNP in this position?		
			if (defined($snps->{$chr}->{$pos})) {
			
				# Is it a hetrozygot or a homozygot
				my $hetHom = 1;
				if(substr($snps->{$chr}->{$pos}->{'genotype'},0,1) eq substr($snps->{$chr}->{$pos}->{'genotype'},1,1)) {
					$hetHom = 2;
				}
			
				print OUTPUT $_."\t".$hetHom."\t".$snps->{$chr}->{$pos}->{'rawRow'}."\n";
			
			} else {
				print OUTPUT $_."\t0\tNo ref. SNP\n";
			}
		}
	}
	
	print "DONE!\n";
	
	# Closing files
	close(OUTPUT);
	close(INPUT);
	
# END MAIN -----

# Sub for loading snps to memory
sub loadSNPs {
	
	# Open dbSNP file.
	open(SNP, "< $snpRef{'file'}") or die "Oops, could not open ref. SNPs file: $!";
			
		while(<SNP>) {
			
			chomp;
			if (!(m/^#/)) {
			
				my @columns = split(/\t/,$_);
				my $chr 	 = $columns[($snpRef{'chr_col'} - 1)];
				my $pos 	 = $columns[($snpRef{'pos_col'} - 1)];
				my $genotype = $columns[($snpRef{'genotype'} - 1)];
							
				$snps->{$chr}->{$pos} = {
										  genotype => $genotype,
										  rawRow   => $_,
										};	
			}
		}
	
	close(SNP);
}

# Sub show how to run this script.
sub usage {
  print "\nUsage: $0 -i <jSNPmania-file> -snp_ref <file,chr_col,pos_col,genotype_col> -o <outputfile>\n 
 -i  	   Input file jSNPmaina format.
 -snp_ref  SNP ref. file and columns with chromosome, position and genotype (separated by \",\".
           e.g. \"file.txt,2,3,4\".
 -o  	   Output file (default /dev/stdout).
 -nc2chr   Converts jSNPmania NC-numbers to chr (e.g. NC_000019.8 -> chr19).
           chr23 and chr24 will be chrX and chrY. !! WARNING !!, this will disregard
           version of the NC-number.\n\n";   
  exit(1);
}