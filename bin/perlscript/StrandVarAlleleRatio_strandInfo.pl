#!/usr/bin/perl

use warnings;
use strict;

use FileHandle;

# Subroutine prototypes
sub usage;
sub getMajorVariantAllele;
sub getInsertions;
sub getDeletions;
sub extractSampleName;
sub countAmplicons;
sub createReferenceHash;

my ( $next_arg, $variantFile, $outputFile, $minRD, $varMinRD, $minVariantAlleleRatio );

if ( scalar(@ARGV) == 0 ) {
	usage();
}

# Parse the command line
while ( scalar @ARGV > 0 ) {
	$next_arg = shift(@ARGV);
	if    ( $next_arg eq "-v" )      { $variantFile = shift(@ARGV); }
	elsif ( $next_arg eq "-o" )      { $outputFile  = shift(@ARGV); }
	elsif ( $next_arg eq "-minRD" )  { $minRD       = shift(@ARGV); }
	elsif ( $next_arg eq "-vMinRD" ) { $varMinRD    = shift(@ARGV); }
	else { print "Invalid argument: $next_arg"; usage(); }
}

# Error
if ( !$variantFile || !$outputFile || !$minRD ) {
	usage();
}

open( OUTPUT, ">", $outputFile )
  or die "Couldn't open output file " . $outputFile . "!";
print OUTPUT
  "#Chr\tPos\tReference_base\tVariant_base\tVariant_allele_ratio_plus\tVariant_allele_ratio_minus\t#_reference_reads_plus\t#_reference_reads_minus\t#_variant_reads_plus\t#_variant_reads_minus\n" ;

open( VARIANT, "<", $variantFile )
  or die "Couldn't open variant file " . $variantFile . "!";
while (<VARIANT>) {
	if ( $_ =~ m/^#/ || $_ =~ m/^$/ ) { next; }
	chomp;

	my @line = split( /\t/, $_ );

	# Check that the total amount of variant is above $minVariantAlleleRatio and that the position has enough read depth
	if ( $line[0] >= $minRD ) {
		getMajorVariantAllele( \@line );
	}
}
close(VARIANT);
close(OUTPUT);

# Print the usage help for this script.
sub usage {
	print "
************************************************************************************************

 This script takes a variation file from SNPmania and outputs the variant allele ratio per
 strand for each position fulfilling the given criteria.

************************************************************************************************
\nUsage: $0\n 
 
 -v             Variation file from SNPmania
 -o             Output file in annovar input format
 -minRD         Minimum total (ref+var) read depth PER strand
 -vMinRD        Minimum read depth for the variant allele
  \n";
	exit 1;
}

sub getMajorVariantAllele {
	my ($lineKey) = @_;
	my @line = @$lineKey;

	# Split column with reads per base on |
	my @alleles = split( /\|/, $line[6] );

	my @alternatives = ( "A", "G", "C", "T" );
	my $ref_plus = my $ref_minus = my $var_plus = my $var_minus = 0;
	my $count    = 0;

	foreach my $ref (@alternatives) {
		if ( $ref eq $line[5] ) {
			( $ref_plus, $ref_minus ) = countStrandReads( $line[ ( $count + 10 ) ], "sub" );
		}
		$count++;
	}

	$count = 0;

	# Go through the base alternatives again
	foreach my $variant (@alternatives) {

		# Check that base is NOT the reference base
		if ( $variant ne $line[5] ) {

			# Count the number of plus and minus amplicons for the variant
			my ( $var_plus, $var_minus ) = countStrandReads( $line[ ( $count + 10 ) ], "sub" );
			my $plusRatio  = my $minusRatio  = 0;
			my $plusStatus = my $minusStatus = "false";
			if ( $var_plus + $var_minus >= $varMinRD ) {

				# Check that there are plus reads
				if ( $ref_plus + $var_plus >= $minRD ) {
					$plusRatio = $var_plus / ( $ref_plus + $var_plus );
					$plusStatus = "true";
				}

				# Check that there are minus reads
				if ( $ref_minus + $var_minus >= $minRD ) {
					$minusRatio = $var_minus / ( $ref_minus + $var_minus );
					$minusStatus = "true";
				}

				if ( $plusStatus eq "true" && $minusStatus eq "true" ) {
					print OUTPUT $line[3] . "\t"
					  . $line[4] . "\t"
					  . $line[5] . "\t"
					  . $variant . "\t"
					  . $plusRatio . "\t"
					  . $minusRatio . "\t"
					  . $ref_plus . "\t"
					  . $ref_minus . "\t"
					  . $var_plus . "\t"
					  . $var_minus . "\n";
				}
			}
		}
		$count++;
	}
}

sub countStrandReads {
	my ( $strandInfo, $type ) = @_;

	# Variables to count plus and minus amplicons
	my $plus  = 0;
	my $minus = 0;

	# Split the strand info in an array
	my @strands = split( /\|/, $strandInfo );
	$plus  = $strands[0] + $strands[3];
	$minus = $strands[1] + $strands[2];

	return ( $plus, $minus );
}
