#!/usr/bin/perl

use warnings;
use strict;

use FileHandle;

# Subroutine prototypes
sub usage;

my ( $next_arg, $variantFile, $outputFile, $regionFile, $mRD, $sampleID, $header );

if ( scalar(@ARGV) == 0 ) {
	usage();
}

# Parse the command line
while ( scalar @ARGV > 0 ) {
	$next_arg = shift(@ARGV);
	if    ( $next_arg eq "-v" ) { $variantFile = shift(@ARGV); }
	elsif ( $next_arg eq "-o" ) { $outputFile  = shift(@ARGV); }
	elsif ( $next_arg eq "-r" ) { $regionFile  = shift(@ARGV); }
	elsif ( $next_arg eq "-m" ) { $mRD         = shift(@ARGV); }
	elsif ( $next_arg eq "-s" ) { $sampleID    = shift(@ARGV); }
	elsif ( $next_arg eq "-h" ) { $header      = 1; }
	else { print "Invalid argument: $next_arg"; usage(); }
}

# Error
if (   !$variantFile
	|| !$regionFile
	|| !$outputFile )
{
	usage();
}

if ( !$mRD ) {
	$mRD = 1;
}

if ( !$sampleID ) {
	$sampleID = $variantFile;
}

my %region    = ();
my %bpCovered = ();

# Split the given min read depths
my @minRDs = split( /\,/, $mRD );

# Create a hash with each rd as key and the number of covered bases as value
my $count = 1;
$bpCovered{$count}{'0'} = 0;
foreach my $rd (@minRDs) {
	$count++;
	$bpCovered{$count}{$rd} = 0;
}

print "Hashing region positions\n";
open( REGION, "<", $regionFile )
  or die "Couldn't open region file " . $regionFile . "!";
while (<REGION>) {
	if ( $_ =~ m/^#/ || $_ =~ m/^$/ ) { next; }
	chomp;

	my @line = split( /\t/, $_ );
	$region{ $line[1] }{ $line[2] } = $line[3];
}
close(REGION);

print "Go through positions\n";
open( VARIANT, "<", $variantFile )
  or die "Couldn't open variant file " . $variantFile . "!";
while (<VARIANT>) {
	if ( $_ =~ m/^#/ || $_ =~ m/^$/ ) { next; }
	chomp;

	my @line = split( /\t/, $_ );

	# Check if pos is within a region
	for my $start ( keys %{ $region{ $line[3] } } ) {
		if ( $start <= $line[4] && $region{ $line[3] }{$start} >= $line[4] ) {
			foreach my $index ( keys %bpCovered ) {
				foreach my $rd ( keys %{ $bpCovered{$index} } ) {
					if ( $line[0] >= $rd ) {
						$bpCovered{$index}{$rd}++;
					}
				}
			}
		}
	}
}
close(VARIANT);

# Print fraction covered
print "Calculate fraction covered\n";
my $print_frac      = $sampleID;
my $print_bp        = "";
my $print_header    = "#Sample";
my $print_header_bp = "";
foreach my $index ( sort keys %bpCovered ) {
	foreach my $rd ( keys %{ $bpCovered{$index} } ) {
		if ( $rd > 0 ) {
			my $frac = $bpCovered{$index}{$rd} / $bpCovered{'1'}{'0'} * 100;
			$print_frac      .= "\t" . $frac;
			$print_bp        .= "\t" . $bpCovered{$index}{$rd};
			$print_header    .= "\tF" . $rd . "x";
			$print_header_bp .= "\tbp_minRD_" . $rd;

		}
	}
}
my $print_str = $print_frac . $print_bp . "\t" . $bpCovered{'1'}{'0'};
$print_header .= $print_header_bp . "\tTot_#_bp";

print "Printing\n";
open( OUTPUT, ">>", $outputFile )
  or die "Couldn't open output file " . $outputFile . "!";
if ($header) {
	print OUTPUT $print_header . "\n";
}
print OUTPUT $print_str . "\n";

# Print the usage help for this script.

sub usage {
	print "
************************************************************************************************

 This script takes the variation from SNPmania and a file with regions and outputs the fraction
 covered at the given read depths.
 
************************************************************************************************
\nUsage: $0\n 
 
 -v            Variation file from SNPmania
 -o            Output file
 -r            Region file, tab-delimited {ID \t chr \t start \t end}
 -m            Read depths to calculate fraction covere for, several RD can be given
                   separated by comma (-m 30,50,100) [optional, default is 1]
 -s            Sample id [optional, default variation file name]
 -h            Set to print header in the file
 \n";
	exit(1);
}
