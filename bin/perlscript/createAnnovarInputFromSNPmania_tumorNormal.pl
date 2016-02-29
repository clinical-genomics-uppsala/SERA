#!/usr/bin/perl

use warnings;
use strict;

use FileHandle;

# Subroutine prototypes
sub usage;
sub getMajorVariantAllele;
sub getInsertions;
sub getDeletions;
sub hashingChrInfo;
sub normalHomozygous;

my ( $next_arg, $tumorVariantFile, $tumorInsertionFile, $tumorDeletionFile, $outputFile, $nc2chr, $nc_col, $chr_col, $tumorMinRD, $minVariantAlleleRatio, $ampliconMapped );

my ( $normalVariantFile, $normalInsertionFile, $normalDeletionFile, $normalMinRD, $normalHomoRatio );

if ( scalar(@ARGV) == 0 ) {
	usage();
}

# Parse the command line
while ( scalar @ARGV > 0 ) {
	$next_arg = shift(@ARGV);
	if    ( $next_arg eq "-tv" )            { $tumorVariantFile      = shift(@ARGV); }
	elsif ( $next_arg eq "-ti" )            { $tumorInsertionFile    = shift(@ARGV); }
	elsif ( $next_arg eq "-td" )            { $tumorDeletionFile     = shift(@ARGV); }
	elsif ( $next_arg eq "-nv" )            { $normalVariantFile     = shift(@ARGV); }
	elsif ( $next_arg eq "-ni" )            { $normalInsertionFile   = shift(@ARGV); }
	elsif ( $next_arg eq "-nd" )            { $normalDeletionFile    = shift(@ARGV); }
	elsif ( $next_arg eq "-o" )             { $outputFile            = shift(@ARGV); }
	elsif ( $next_arg eq "-nc2chr" )        { $nc2chr                = shift(@ARGV); }
	elsif ( $next_arg eq "-nc_col" )        { $nc_col                = shift(@ARGV); }
	elsif ( $next_arg eq "-chr_col" )       { $chr_col               = shift(@ARGV); }
	elsif ( $next_arg eq "-tminRD" )        { $tumorMinRD            = shift(@ARGV); }
	elsif ( $next_arg eq "-nminRD" )        { $normalMinRD           = shift(@ARGV); }
	elsif ( $next_arg eq "-am" )            { $ampliconMapped        = shift(@ARGV); }
	elsif ( $next_arg eq "-tminVarRatio" )  { $minVariantAlleleRatio = shift(@ARGV); }
	elsif ( $next_arg eq "-nHomoRefRatio" ) { $normalHomoRatio       = shift(@ARGV); }
	else { print "Invalid argument: $next_arg"; usage(); }
}

# Error
if (   !$tumorVariantFile
	|| !$normalVariantFile
	|| !$outputFile
	|| !$tumorMinRD
	|| !$normalMinRD
	|| !$normalHomoRatio
	|| !$minVariantAlleleRatio )
{
	usage();
}

# Create hash with NC-number as key and chr number as value
my %ncToChr   = ();
my %positions = ();
my %tumorRef  = ();
my %normalRef = ();

hashingChrInfo();

print "Hashing normal variants\n";
open( VARIANT, "<", $normalVariantFile )
  or die "Couldn't open normal variant file " . $normalVariantFile . "!";
while (<VARIANT>) {
	if ( $_ =~ m/^#/ || $_ =~ m/^$/ ) { next; }
	chomp;

	my @variantLine = split( /\t/, $_ );

	#	$normalRef{$variantLine[3]}{$variantLine[4]}=$variantLine[5];
	if ( $variantLine[0] >= $normalMinRD ) {
		if ( $variantLine[2] >= $normalHomoRatio ) {
			createNormalHash( \@variantLine, $variantLine[5] );
		}
		elsif ( $variantLine[2] <= ( 1 - $normalHomoRatio )
			&& normalHomozygous( $variantLine[6], $normalHomoRatio, $variantLine[0] ) )
		{
			createNormalHash( \@variantLine, $variantLine[7] );
		}
	}
}
close(VARIANT);

print "Hashing tumor variants\n";
open( VARIANT, "<", $tumorVariantFile )
  or die "Couldn't open tumor variant file " . $tumorVariantFile . "!";
while (<VARIANT>) {
	if ( $_ =~ m/^#/ || $_ =~ m/^$/ ) { next; }
	chomp;

	# Split line on tab
	my @variantLine = split( /\t/, $_ );

	# Hashing chr, pos, referenceBase, referenceBase_RD, ampInfo to be able to check later on
	if ( $normalRef{ $variantLine[3] }{ $variantLine[4] } ) {
		createReferenceHash( \@variantLine );
	}

	# If the read depth is higher than tumorMinRD and the reference allele ratio is less than 1-minVariantAlleleRatio
	if (   $variantLine[0] >= $tumorMinRD
		&& $variantLine[2] <= ( 1 - $minVariantAlleleRatio )
		&& $normalRef{ $variantLine[3] }{ $variantLine[4] }{'allele'} )
	{
		getMajorVariantAllele( \@variantLine );
	}
}
close(VARIANT);

print "\tPrinting variants!\n";
open( OUTPUT, ">", $outputFile )
  or die "Couldn't open output file " . $outputFile . "!";

# Print parameter info
print OUTPUT "#Min_tumor_read_depth=" . $tumorMinRD . "\n";
print OUTPUT "#Min_normal_read_depth=" . $normalMinRD . "\n";
print OUTPUT "#Min_tumor_variant_allele_ratio=" . $minVariantAlleleRatio . "\n";
print OUTPUT "#Min_reference_allele_ratio_to call_normal_position_homozygot=" . $normalHomoRatio . "\n";
print OUTPUT "#Min_amplicon_read_depth=" . $ampliconMapped . "\n";

#Print all variants found
for my $key ( keys %positions ) {
	print OUTPUT $positions{$key} . "\n";
}

%positions = ();

# If insertion file is given
my %normalInsRef = ();
if ($tumorInsertionFile) {

	open( INSERTION, "<", $normalInsertionFile )
	  or die "Couldn't open insertion file " . $normalInsertionFile . "!";
	while (<INSERTION>) {
		if ( $_ =~ m/^#/ || $_ =~ m/^$/ ) { next; }
		chomp;

		my @insertionLine = split( /\t/, $_ );

		# Check that the position is present in the normalRef hash and that it doesn't have an insertion and save ok positions
		if (   $insertionLine[4] <= ( 1 - $normalHomoRatio )
			&& $normalRef{ $insertionLine[1] }{ $insertionLine[2] }{'allele'} )
		{
			$normalInsRef{ $insertionLine[1] }{ $insertionLine[2] } = 1;
		}
	}

	close(INSERTION);

	print "Hashing tumor insertions\n";

	# Open tumor insertion file
	open( INSERTION, "<", $tumorInsertionFile )
	  or die "Couldn't open insertion file " . $tumorInsertionFile . "!";
	while (<INSERTION>) {
		if ( $_ =~ m/^#/ || $_ =~ m/^$/ ) { next; }
		chomp;

		# Split the line by tab
		my @insertionLine = split( /\t/, $_ );

		# If the read depth is higher than tumorMinRD and the insertion ratio is higher than minVariantAlleleRatio
		if (   $insertionLine[0] >= $tumorMinRD
			&& $insertionLine[4] >= $minVariantAlleleRatio
			&& $normalInsRef{ $insertionLine[1] }{ $insertionLine[2] } )
		{
			getInsertions( \@insertionLine );
		}
	}
	close(INSERTION);

	# Print all insertions found
	print "\tPrinting insertions!\n";
	for my $key ( keys %positions ) {
		print OUTPUT $positions{$key} . "\n";
	}
	%positions    = ();
	%normalInsRef = ();
}

# If deletion file is given
my %normalDelRef = ();
if ($tumorDeletionFile) {
	print "Hasing normal deletions\n";
	open( DELETION, "<", $normalDeletionFile )
	  or die "Couldn't open deletion file " . $normalDeletionFile . "!";
	while (<DELETION>) {
		if ( $_ =~ m/^#/ || $_ =~ m/^$/ ) { next; }
		chomp;

		my @deletionLine = split( /\t/, $_ );

		# Check that the position is present in the normalRef hash and that it doesn't have a deletion and save ok positions
		if (   $normalRef{ $deletionLine[1] }{ $deletionLine[2] }{'allele'}
			&& $deletionLine[4] <= ( 1 - $normalHomoRatio ) )
		{
			$normalDelRef{ $deletionLine[1] }{ $deletionLine[2] } = 1;
		}
	}
	close(DELETION);

	print "Hasing tumor deletions\n";
	open( DELETION, "<", $tumorDeletionFile )
	  or die "Couldn't open deletion file " . $tumorDeletionFile . "!";
	while (<DELETION>) {
		if ( $_ =~ m/^#/ || $_ =~ m/^$/ ) { next; }
		chomp;

		my @deletionLine = split( /\t/, $_ );

		# If the read depth is higher than tumorMinRD, the deletion ratio is higher than minVariantAlleleRatio and the positions is present in the normalDelRef extract deletion
		if (   $deletionLine[0] >= $tumorMinRD
			&& $deletionLine[4] >= $minVariantAlleleRatio
			&& $normalDelRef{ $deletionLine[1] }{ $deletionLine[2] } )
		{
			getDeletions( \@deletionLine );
		}
	}
	close(DELETION);

	# Print all deletions
	print "\tPrinting deletions!\n";
	for my $key ( keys %positions ) {
		print OUTPUT $positions{$key} . "\n";
	}
}

print "DONE!\n";
close(OUTPUT);

sub normalHomozygous {
	my ( $normalVariants, $minHomo, $readDepth ) = @_;

	my @alleles = split( /\|/, $normalVariants );

	foreach my $var (@alleles) {
		if ( ( $var / $readDepth ) >= $minHomo ) {
			return 1;
		}
	}
	return 0;
}

sub getMajorVariantAllele {
	my ($lineKey) = @_;
	my @line = @$lineKey;

	my @alleles = split( /\|/, $line[6] );
	my @alternatives = ( "A", "G", "C", "T" );
	my @tumorAmp = split( /\|/, $line[14] );

	my $count = 0;

	my $tumorRefRD = $tumorRef{ $line[3] }{ $line[4] }{'refRD'};

	$count = 0;
	foreach my $variant (@alternatives) {

		# Check that the variant is not the same as the reference in the normal
		if ( $variant ne $normalRef{ $line[3] }{ $line[4] }{'allele'} ) {

			if ( ( $alleles[$count] / $line[0] ) >= $minVariantAlleleRatio ) {
				my $key = $ncToChr{ $line[3] } . "#" . $line[4] . "#" . $line[5] . "#" . $variant;

				my $sample = extractSampleName($tumorVariantFile);

				# Check if the SNPmania file has ampliconmapping info
				if ($ampliconMapped) {

					# Count the number of plus and minus amplicons for the variant
					my ( $amp_plus, $amp_minus ) = countAmplicons( $tumorAmp[$count], "sub" );

					$positions{$key} =
					    $ncToChr{ $line[3] } . "\t"
					  . $line[4] . "\t"
					  . $line[4] . "\t"
					  . $normalRef{ $line[3] }{ $line[4] }{'allele'} . "\t"
					  . $variant
					  . "\tcomments: sample="
					  . $sample
					  . " variantAlleleRatio="
					  . ( $alleles[$count] / $line[0] )
					  . " alleleFreq="
					  . $tumorRefRD . ","
					  . $alleles[$count]
					  . " readDepth="
					  . $line[0]
					  . " Tumor_A="
					  . $line[10]
					  . " Tumor_G="
					  . $line[11]
					  . " Tumor_C="
					  . $line[12]
					  . " Tumor_T="
					  . $line[13] . " "
					  . $normalRef{ $line[3] }{ $line[4] }{'strandCov'}
					  . " Tumor_var_plusAmplicons="
					  . $amp_plus
					  . " Tumor_var_minusAmplicons="
					  . $amp_minus
					  . " Tumor_ref_plusAmplicons="
					  . $tumorRef{ $line[3] }{ $line[4] }{'amp+'}
					  . " Tumor_ref_minusAmplicons="
					  . $tumorRef{ $line[3] }{ $line[4] }{'amp-'}
					  . " Normal_ref_plusAmplicons="
					  . $normalRef{ $line[3] }{ $line[4] }{'amp+'}
					  . " Normal_ref_minusAmplicons="
					  . $normalRef{ $line[3] }{ $line[4] }{'amp-'}
					  . " Tumor_var_ampliconInfo="
					  . $tumorAmp[$count]
					  . " Tumor_ref_ampliconInfo="
					  . $tumorRef{ $line[3] }{ $line[4] }{'ampInfo'}
					  . " Normal_ref_ampliconInfo="
					  . $normalRef{ $line[3] }{ $line[4] }{'ampInfo'};
				}
				else {
					$positions{$key} =
					    $ncToChr{ $line[3] } . "\t"
					  . $line[4] . "\t"
					  . $line[4] . "\t"
					  . $normalRef{ $line[3] }{ $line[4] }{'allele'} . "\t"
					  . $variant
					  . "\tcomments: sample="
					  . $sample
					  . " variantAlleleRatio="
					  . ( $alleles[$count] / $line[0] )
					  . " alleleFreq="
					  . $tumorRefRD . ","
					  . $alleles[$count]
					  . " readDepth="
					  . $line[0]
					  . " Tumor_A="
					  . $line[10]
					  . " Tumor_G="
					  . $line[11]
					  . " Tumor_C="
					  . $line[12]
					  . " Tumor_T="
					  . $line[13] . " "
					  . $normalRef{ $line[3] }{ $line[4] }{'strandCov'};
				}
			}
		}
		$count++;
	}
}

sub getInsertions {
	my ($lineKey) = @_;
	my @line = @$lineKey;

	# Check that the line has an insert
	if ( $line[5] =~ m/[ACGTN]/ ) {

		# Split the insertion on |
		my @insertions = split( /\|/, $line[5] );

		# Retrieve the read depth for the reference
		my $tumorRefRD = $tumorRef{ $line[1] }{ $line[2] }{'refRD'};

		# Go through the array of insertions
		for ( my $i = 0 ; $i < scalar(@insertions) ; $i++ ) {

			# Find the read depth of this particular insertion by substitute all letters with nothing
			my $insertRD = $insertions[$i];
			$insertRD =~ s/[ACGTN]+//;

			# Find what is inserted by substituting all numbers with nothing
			my $insert = $insertions[$i];
			$insert =~ s/[0-9]+//;

			# Check if the insertion ratio is above $minVariantAlleleRatio
			if ( ( $insertRD / $line[0] ) >= $minVariantAlleleRatio ) {

				# If so create key with chr {tab} pos {tab} pos {tab} insert
				my $key = $ncToChr{ $line[1] } . "#" . $line[2] . "#" . $line[2] . "#" . $insert;

				# Extract sample name from file name
				my $sample = extractSampleName($tumorVariantFile);

				if ($ampliconMapped) {

					my @ampliconOptions = split( /\|/, $line[7] );
					my $ampInfo = my $amp_plus = my $amp_minus = my $amp = "-";
					foreach my $am (@ampliconOptions) {
						my $match = "^" . $insert . ":";

						# Check if the amplicon info matches the insert
						if ( $am =~ m/$match/ ) {

							# Count the number of plus and minus amplicons for the variant
							my ( $amp_plus, $amp_minus ) = countAmplicons( $am, "ins" );
							$amp = $am;
						}
					}

					$positions{$key} =
					    $ncToChr{ $line[1] } . "\t"
					  . $line[2] . "\t"
					  . $line[2] . "\t-\t"
					  . $insert
					  . "\tcomments: sample="
					  . $sample
					  . " variantAlleleRatio="
					  . ( $insertRD / $line[0] )
					  . " alleleFreq="
					  . $tumorRefRD . ","
					  . $insertRD
					  . " readDepth="
					  . $line[0]
					  . " Tumor_Ins="
					  . $line[6] . " "
					  . $normalRef{ $line[1] }{ $line[2] }{'strandCov'}
					  . " Tumor_var_plusAmplicons="
					  . $amp_plus
					  . " Tumor_var_minusAmplicons="
					  . $amp_minus
					  . " Tumor_ref_plusAmplicons="
					  . $tumorRef{ $line[1] }{ $line[2] }{'amp+'}
					  . " Tumor_ref_minusAmplicons="
					  . $tumorRef{ $line[1] }{ $line[2] }{'amp-'}
					  . " Normal_ref_plusAmplicons="
					  . $normalRef{ $line[1] }{ $line[2] }{'amp+'}
					  . " Normal_ref_minusAmplicons="
					  . $normalRef{ $line[1] }{ $line[2] }{'amp-'}
					  . " Tumor_var_ampliconInfo="
					  . $amp
					  . " Tumor_ref_ampliconInfo="
					  . $tumorRef{ $line[1] }{ $line[2] }{'ampInfo'}
					  . " Normal_ref_ampliconInfo="
					  . $normalRef{ $line[1] }{ $line[2] }{'ampInfo'};
				}
				else {
					$positions{$key} =
					    $ncToChr{ $line[1] } . "\t"
					  . $line[2] . "\t"
					  . $line[2] . "\t-\t"
					  . $insert
					  . "\tcomments: sample="
					  . $sample
					  . " variantAlleleRatio="
					  . ( $insertRD / $line[0] )
					  . " alleleFreq="
					  . $tumorRefRD . ","
					  . $insertRD
					  . " readDepth="
					  . $line[0]
					  . " Tumor_Ins="
					  . $line[6] . " "
					  . $normalRef{ $line[1] }{ $line[2] }{'strandCov'};
				}
			}
		}
	}
}

sub getDeletions {
	my ($lineKey) = @_;
	my @line = @$lineKey;

	# Check that the position contains deletions
	if ( $line[5] =~ m/\(/ ) {

		# Split the deletion column on |
		my @deletions = split( /\|/, $line[5] );

		# Retrieve the read depth for the reference
		my $tumorRefRD = $tumorRef{ $line[1] }{ $line[2] }{'refRD'};

		# Go through all deletions
		for ( my $i = 0 ; $i < scalar(@deletions) ; $i++ ) {

			# Split the deletion on ( which gives the read depth for the deletion and the rest
			my ( $deleteRD, $tmp ) = split( /\(/, $deletions[$i] );

			# Substitute ) with nothing
			$tmp =~ s/\)//;

			# By splitting on , the from and to position for the deletion are extracted
			my ( $fromBP, $toBP ) = split( /,/, $tmp );

			# Extract the bases which are deleted from the ref-hash
			my $deleted    = "";
			my $deletionOK = "true";
			for ( my $pos = $fromBP ; $pos <= $toBP ; $pos++ ) {

				# Check that all positions are ok in the normal
				if ( $normalRef{ $line[1] }{ ( $line[2] + $pos ) }{'allele'} ) {
					$deleted .= $normalRef{ $line[1] }{ ( $line[2] + $pos ) }{'allele'};
				}

				# If not the deletion is not OK and should be excluded from the analysis
				else {
					$deletionOK = "false";
				}
			}

			# Only continue with deletions that are OK!
			if ( $deletionOK eq "true" ) {

				# Create a key string with chr#posFrom {tab} posTo#-
				my $key = $ncToChr{ $line[1] } . "#" . ( $line[2] + $fromBP ) . "\t" . ( $line[2] + $toBP ) . "#-";

				# Check if the deletion has a variant allele ratio that is above $minVariantAlleleRatio
				if (   ( $deleteRD / $line[0] ) >= $minVariantAlleleRatio
					&& ( $deleteRD + $tumorRefRD ) >= $tumorMinRD )
				{

					# Check if the deletion already is found, if not add it
					if ( !$positions{$key} ) {
						my $sample = extractSampleName($tumorVariantFile);

						# Check if there are ampliconmapping info in the file
						if ($ampliconMapped) {

							my @ampliconOptions = split( /\|/, $line[7] );
							my $ampInfo = my $amp_plus = my $amp_minus = my $amp = "-";
							foreach my $am (@ampliconOptions) {

								# Create pattern to match ex. ^\\(0,0\\)
								my $match = "^\\(" . $fromBP . "," . $toBP . "\\):";

								# Check if the amplicon info matches the deletion
								if ( $am =~ m/$match/ ) {

									# Count the number of plus and minus amplicons for the variant
									my ( $amp_plus, $amp_minus ) = countAmplicons( $am, "del" );
									$amp = $am;
								}
							}

							$positions{$key} =
							    $ncToChr{ $line[1] } . "\t"
							  . ( $line[2] + $fromBP ) . "\t"
							  . ( $line[2] + $toBP ) . "\t"
							  . $deleted
							  . "\t-\tcomments: sample="
							  . $sample
							  . " variantAlleleRatio="
							  . ( $deleteRD / $line[0] )
							  . " alleleFreq="
							  . $tumorRefRD . ","
							  . $deleteRD
							  . " readDepth="
							  . $line[0]
							  . " Tumor_Del="
							  . $line[6] . " "
							  . $normalRef{ $line[1] }{ $line[2] }{'strandCov'}
							  . " Tumor_var_plusAmplicons="
							  . $amp_plus
							  . " Tumor_var_minusAmplicons="
							  . $amp_minus
							  . " Tumor_ref_plusAmplicons="
							  . $tumorRef{ $line[1] }{ $line[2] }{'amp+'}
							  . " Tumor_ref_minusAmplicons="
							  . $tumorRef{ $line[1] }{ $line[2] }{'amp-'}
							  . " Normal_ref_plusAmplicons="
							  . $normalRef{ $line[1] }{ $line[2] }{'amp+'}
							  . " Normal_ref_minusAmplicons="
							  . $normalRef{ $line[1] }{ $line[2] }{'amp-'}
							  . " Tumor_var_ampliconInfo="
							  . $amp
							  . " Tumor_ref_ampliconInfo="
							  . $tumorRef{ $line[1] }{ $line[2] }{'ampInfo'}
							  . " Normal_ref_ampliconInfo="
							  . $normalRef{ $line[1] }{ $line[2] }{'ampInfo'};
						}
						else {
							$positions{$key} =
							    $ncToChr{ $line[1] } . "\t"
							  . ( $line[2] + $fromBP ) . "\t"
							  . ( $line[2] + $toBP ) . "\t"
							  . $deleted
							  . "\t-\tcomments: sample="
							  . $sample
							  . " variantAlleleRatio="
							  . ( $deleteRD / $line[0] )
							  . " alleleFreq="
							  . $tumorRefRD . ","
							  . $deleteRD
							  . " readDepth="
							  . $line[0]
							  . " Tumor_Del="
							  . $line[6] . " "
							  . $normalRef{ $line[1] }{ $line[2] }{'strandCov'};
						}
					}
					else {

						# Split the line on tab
						my @hashLine = split( /\t/, $positions{$key} );

						# Extract the comments by splitting on white space
						my @comments =
						  split( /\s/, $hashLine[ ( scalar(@hashLine) - 1 ) ] );

						# Go through all comments
						for ( my $i = 1 ; $i < scalar(@comments) ; $i++ ) {

							# and split on =
							my @commentLine = split( /=/, $comments[$i] );

							# If the name of the comment is variantAlleleRatio and the variant allele ratio is higher than the one saved update info
							if ( $commentLine[0] =~ m/variantAlleleRatio/i
								&& ( $deleteRD / $line[0] ) > $commentLine[1] )
							{
								my $sample = extractSampleName($tumorVariantFile);

								# Check if there are ampliconmapping info in the file
								if ($ampliconMapped) {

									my @ampliconOptions =
									  split( /\|/, $line[7] );
									my $ampInfo = my $amp_plus = my $amp_minus = my $amp = "-";
									foreach my $am (@ampliconOptions) {

										# Create pattern to match ex. ^\\(0,0\\)
										my $match = "^\\(" . $fromBP . "," . $toBP . "\\):";

										# Check if the amplicon info matches the deletion
										if ( $am =~ m/$match/ ) {

											# Count the number of plus and minus amplicons for the variant
											my ( $amp_plus, $amp_minus ) = countAmplicons( $am, "del" );
											$amp = $am;
										}
									}
									$positions{$key} =
									    $ncToChr{ $line[1] } . "\t"
									  . ( $line[2] + $fromBP ) . "\t"
									  . ( $line[2] + $toBP ) . "\t"
									  . $deleted
									  . "\t-\tcomments: sample="
									  . $sample
									  . " variantAlleleRatio="
									  . ( $deleteRD / $line[0] )
									  . " alleleFreq="
									  . $tumorRefRD . ","
									  . $deleteRD
									  . " readDepth="
									  . $line[0]
									  . " Tumor_Del="
									  . $line[6] . " "
									  . $normalRef{ $line[1] }{ $line[2] }{'strandCov'}
									  . " Tumor_var_plusAmplicons="
									  . $amp_plus
									  . " Tumor_var_minusAmplicons="
									  . $amp_minus
									  . " Tumor_ref_plusAmplicons="
									  . $tumorRef{ $line[1] }{ $line[2] }{'amp+'}
									  . " Tumor_ref_minusAmplicons="
									  . $tumorRef{ $line[1] }{ $line[2] }{'amp-'}
									  . " Normal_ref_plusAmplicons="
									  . $normalRef{ $line[1] }{ $line[2] }{'amp+'}
									  . " Normal_ref_minusAmplicons="
									  . $normalRef{ $line[1] }{ $line[2] }{'amp-'}
									  . " Tumor_var_ampliconInfo="
									  . $amp
									  . " Tumor_ref_ampliconInfo="
									  . $tumorRef{ $line[1] }{ $line[2] }{'ampInfo'}
									  . " Normal_ref_ampliconInfo="
									  . $normalRef{ $line[1] }{ $line[2] }{'ampInfo'};
								}
								else {
									$positions{$key} =
									    $ncToChr{ $line[1] } . "\t"
									  . ( $line[2] + $fromBP ) . "\t"
									  . ( $line[2] + $toBP ) . "\t"
									  . $deleted
									  . "\t-\tcomments: sample="
									  . $sample
									  . " variantAlleleRatio="
									  . ( $deleteRD / $line[0] )
									  . " alleleFreq="
									  . $tumorRefRD . ","
									  . $deleteRD
									  . " readDepth="
									  . $line[0]
									  . " Tumor_Del="
									  . $line[6] . " "
									  . $normalRef{ $line[1] }{ $line[2] }{'strandCov'};
								}
							}
						}
					}
				}
			}
		}
	}
}

sub hashingChrInfo {
	print "Hashing accession and chromosome conversion\n";
	if ($nc2chr) {
		open( NC, "<", $nc2chr )
		  or die "Couldn't open nc2chr file " . $nc2chr . "!";
		while (<NC>) {
			if ( $_ =~ m/^#/ || $_ =~ m/^$/ ) { next; }

			my @nc_line = split( /\t/, $_ );
			$nc_line[ ( $chr_col - 1 ) ] =~ s/chr//;
			$ncToChr{ $nc_line[ ( $nc_col - 1 ) ] } = $nc_line[ ( $chr_col - 1 ) ];
		}
	}
}

# Print the usage help for this script.
sub usage {
	print "
************************************************************************************************

 This script takes the variation file (and optionally insertion and deletion file) from SNPmania
 for a tumor and and a normal sample and outputs the result as an annovar input file.
 It can also filter on read depth and min variant allele ratio, minimum reference allele ratio,
 in the normal sample.
 
************************************************************************************************
\nUsage: $0\n 
 
 -tv             Tumor variation file from SNPmania
 -ti             Tumor insertion file from SNPmania (optional)
 -td             Tumor deletion file from SNPmania (optional)
 -nv             Normal variation file from SNPmania
 -ni             Normal insertion file from SNPmania (optional)
 -nd             Normal deletion file from SNPmania (optional)
 -o              Output file in annovar input format
 -nc2chr         File with NC and corresponding chr number
 -nc_col         Column with NC number in nc2chr file (first column is 1)
 -chr_col        Column with chr number in nc2chr file (first column is 1)
 -tminRD         Minimum tumor read depth for a position to be considered (optional, default: 1)
 -nminRD         Minimum normal read depth for a position to be considered (optional, default: 1)
 -nHomoRefRatio  Minimum REFERENCE allele ratio for a position to be called homozygouz
 -tminVarRatio   Minimum tumor variant allele ratio for a position to be considered (optional, default: 0)
 -am             Minimum read depth per amplicon to count as covered
 \n";
	exit 1;
}

sub createNormalHash {
	my ( $lineKey, $allele ) = @_;
	my @line = @$lineKey;

	if ($ampliconMapped) {
		my @alleles = split( /\|/, $line[6] );
		my @alternatives = ( "A", "G", "C", "T" );

		my $count = 0;

		# Go through the alternatives array base by base
		foreach my $var (@alternatives) {

			# If the base is the same as $allele
			if ( $var eq $allele ) {

				# Add info to normal hash
				$normalRef{ $line[3] }{ $line[4] } = {
					'allele'    => $allele,
					'strandCov' => "Normal_A=" . $line[10] . " Normal_G=" . $line[11] . " Normal_C=" . $line[12] . " Normal_T=" . $line[13]
				};

				# If amplicon mapping information exists
				if ($ampliconMapped) {

					# Split amplicon information
					my @ampInfo = split( /\|/, $line[14] );

					# Count the number of amplicons per strand
					my ( $amp_plus, $amp_minus ) = countAmplicons( $ampInfo[$count], "sub" );

					# Add amplicon information to normal hash
					$normalRef{ $line[3] }{ $line[4] }{'amp+'}    = $amp_plus;
					$normalRef{ $line[3] }{ $line[4] }{'amp-'}    = $amp_minus;
					$normalRef{ $line[3] }{ $line[4] }{'ampInfo'} = $ampInfo[$count];
				}
			}
			$count++;
		}
	}

	return ();
}

sub createReferenceHash {
	my ( $lineKey, $allele ) = @_;
	my @line = @$lineKey;

	my @alleles = split( /\|/, $line[6] );
	my @alternatives = ( "A", "G", "C", "T" );

	my $count = 0;

	# Go through the alternatives array base by base
	foreach my $var (@alternatives) {

		# Check if the base is the same base as the reference base in the normal, if so save the read depth
		if ( $var eq $normalRef{ $line[3] }{ $line[4] }{'allele'} ) {
			$tumorRef{ $line[3] }{ $line[4] }{'ref'}   = $normalRef{ $line[3] }{ $line[4] }{'allele'};
			$tumorRef{ $line[3] }{ $line[4] }{'refRD'} = $alleles[$count];
			if ($ampliconMapped) {
				my @ampInfo = split( /\|/, $line[14] );
				my ( $amp_plus, $amp_minus ) = countAmplicons( $ampInfo[$count], "sub" );
				$tumorRef{ $line[3] }{ $line[4] }{'amp+'}    = $amp_plus;
				$tumorRef{ $line[3] }{ $line[4] }{'amp-'}    = $amp_minus;
				$tumorRef{ $line[3] }{ $line[4] }{'ampInfo'} = $ampInfo[$count];
			}
		}
		$count++;
	}
}

sub extractSampleName {
	my $variantFile = shift;

	# Extract the sample name from the file name
	my @tmpSample = split( /\//, $variantFile );
	my @sample =
	  split( /\./, $tmpSample[ ( scalar(@tmpSample) - 1 ) ] );

	#return ( \@sample, \@tmpSample );
	return $sample[0];
}

sub countAmplicons {
	my ( $ampliconInfo, $type ) = @_;

	# Variables to count plus and minus amplicons
	my $a_plus  = 0;
	my $a_minus = 0;

	if ( $ampliconInfo =~ m/^0$/ ) {
		return ( 0, 0 );
	}
	else {

		# Split the amplicon string into amplicons in an array
		my @amplicons = split( /\#/, $ampliconInfo );

		# Go through all amplicons
		foreach my $amplicon (@amplicons) {

			# Split the amplicon string on :
			my @info = split( /:/, $amplicon );

			# If it is a substitution
			if ( $type eq "sub" ) {

				# Check that the number of reads for the amplicon is above ampliconMapped
				if ( $info[3] >= $ampliconMapped ) {

					# Check if it is on the plus or minus strand
					if ( $info[2] =~ m/\+/ ) {
						$a_plus++;
					}
					else {
						$a_minus++;
					}
				}
			}

			# If it is an insertion or deletion
			else {

				# Check that the number of reads for the amplicon is above ampliconMapped
				if ( $info[4] >= $ampliconMapped ) {

					# Check if it is on the plus or minus strand
					if ( $info[3] =~ m/\+/ ) {
						$a_plus++;
					}
					else {
						$a_minus++;
					}
				}
			}

		}
		return ( $a_plus, $a_minus );
	}
}
