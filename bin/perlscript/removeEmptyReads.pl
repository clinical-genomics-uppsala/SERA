#!/usr/bin/perl

use warnings;
use strict;

use FileHandle;

# Subroutine prototypes
sub usage;

my (
	$next_arg, $readsFile1, $readsFile2, $indexFile,
	$output1,  $output2,    $outindex
);

if ( scalar(@ARGV) == 0 ) {
	usage();
}

# Parse the command line
while ( scalar @ARGV > 0 ) {
	$next_arg = shift(@ARGV);
	if    ( $next_arg eq "-r1" ) { $readsFile1 = shift(@ARGV); }
	elsif ( $next_arg eq "-r2" ) { $readsFile2 = shift(@ARGV); }
	elsif ( $next_arg eq "-in" ) { $indexFile  = shift(@ARGV); }
	elsif ( $next_arg eq "-o1" ) { $output1    = shift(@ARGV); }
	elsif ( $next_arg eq "-o2" ) { $output2    = shift(@ARGV); }
	elsif ( $next_arg eq "-oi" ) { $outindex   = shift(@ARGV); }
	else { print "Invalid argument: $next_arg"; usage(); }
}

# Error
if ( !$readsFile1 || !$output1) { usage(); }

# Open output file handlers
my $read1out = new FileHandle(">$output1");

# Check if read2 file is given
my $read2out;
if ($readsFile2) {
	if ($output2) {
		$read2out = new FileHandle(">$output2");
	}
	else {
		print "-o2 has to be set when -r2 is used!\n";
		usage();
	}
}

# Check if index file is given
my $indexout;
if ($indexFile) {
	if ($outindex) {
		$indexout = new FileHandle(">$outindex");

	}
	else {
		print "-oi has to be set when -in is used!\n";
		usage();
	}
}

# Open read1 files, possible in gzip format
if ( $readsFile1 =~ m/\.gz$/ ) {
	open( READ1, "gunzip -c $readsFile1 |" )
	  or die "Couldn't open reads file 1 " . $readsFile1 . "!";
}
#else {
#	open( READ1, '<', $readsFile1 )
#	  or die "Couldn't open reads file " . $readsFile1 . "!";
#}

# Check if the index file exists, if so open
if ($indexFile) {
	if ( $indexFile =~ m/\.gz$/ ) {
		open( INDEX, "gunzip -c $indexFile |" )
		  or die "Couldn't open reads file 1 " . $indexFile . "!";
	}
#	else {
#		open( INDEX, '<', $indexFile )
#		  or die "Couldn't open reads file " . $indexFile . "!";
#	}
}

# Check if it's paired end reads, if so open
if ($readsFile2) {
	if ( $readsFile2 =~ m/\.gz$/ ) {
		open( READ2, "gunzip -c $readsFile2 |" )
		  or die "Couldn't open reads file 2 " . $readsFile2 . "!";
	}
#	else {
#		open( READ2, '<', $readsFile2 )
#		  or die "Couldn't open reads file " . $readsFile2 . "!";
#	}
}

# Since all three files opened above are sorted the same way we can read all three
# simultaneously line by line
while (<READ1>) {

	#	if($_ =~m/^#/ || $_ =~ m/^$/) { next; }
	my $read1_id   = $_;
	my $read1_seq  = <READ1>;
	my $read1_plus = <READ1>;
	my $read1_qual = <READ1>;

	# Define variables for index read
	my ( $index_id, $index_seq, $index_plus, $index_qual );

	if ($indexFile) {
		$index_id   = <INDEX>;
		$index_seq  = <INDEX>;
		$index_plus = <INDEX>;
		$index_qual = <INDEX>;
	}
	chomp($read1_seq);

	# Define variables for read file 2
	my ( $read2_id, $read2_seq, $read2_plus, $read2_qual );

	# if read file 2 is defined
	if ($readsFile2) {
		$read2_id   = <READ2>;
		$read2_seq  = <READ2>;
		$read2_plus = <READ2>;
		$read2_qual = <READ2>;
	}
	chomp($read2_seq);

	if ($readsFile2) {
		if ( !$read1_seq =~ m/^$/ && !$read2_seq =~ m/^$/ ) {

			#			print $read2_id.$read2_seq.$read2_plus.$read2_qual;
			print $read1out $read1_id
			  . $read1_seq . "\n"
			  . $read1_plus
			  . $read1_qual;
			print $read2out $read2_id
			  . $read2_seq . "\n"
			  . $read2_plus
			  . $read2_qual;

			if ($indexFile) {
				print $indexout $index_id
				  . $index_seq
				  . $index_plus
				  . $index_qual;
			}
		}
	}

	else {
		if ( !$read1_seq =~ m/^$/ ) {
			print $read1out $read1_id
			  . $read1_seq . "\n"
			  . $read1_plus
			  . $read1_qual;
		}
		if ($indexFile) {
			print $indexout $index_id . $index_seq . $index_plus . $index_qual;
		}
	}
}
close(READ1);
close(READ2);
close(INDEX);
close($read1out);
close($read2out);
close($indexout);

# Print the usage help for this script.
sub usage {
	print "
************************************************************************************************

This script removes reads which are empty after being run through for example cutAdapt. 
 
************************************************************************************************
\nUsage: $0\n 
 
 -r1 Read file for the first read in fastq-format (gzipped is ok)
 -r2 File with paired read [optional]
 -in File with index read [optional]
 -o1 Output file for read 1
 -o2 Output file for read 2, if -r2 is given
 -oi Output file for index, if -in is given\n\n";
	exit 1;
}
