#!/usr/bin/perl
#
# Extract reads of quality score from SOLiD sequencing from its name. 
# Author Elin Falk, 2009 
# Modified Magnus Isaksson, 2010 (Working with Illumina fastq files also.)
#

use strict;
use warnings;
use FindBin;
use FileHandle;

sub usage;

my ($next_arg, $infile, $extractFile, $outputFile, $illumina, $solid);

if(scalar(@ARGV) == 0){
    &usage();
}
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-i")    { $infile = shift(@ARGV); }
    elsif($next_arg eq "-o") { $outputFile = shift(@ARGV); }
    elsif($next_arg eq "-e") { $extractFile = shift(@ARGV); }
    elsif($next_arg eq "-illumina") { $illumina = 1; }
    elsif($next_arg eq "-solid") { $solid = 1; }
    else { print "Invalid argument: $next_arg"; usage(); }
}

## Error
if((!$infile)||(!$extractFile)||(!($solid||$illumina))) { usage(); }
if(($illumina)&&($solid)) { print "Flaggs Illumina and SOLiD can not be both set.\n\n"; usage(); }

# If no output file write to stdout.
if (!$outputFile) { $outputFile="/dev/stdout"; }

## Hashing reads file.
my %extract = ();
my $key;
my $extractFH = new FileHandle($extractFile) or die "Couldn't open input file";
if ($solid) {
	while (<$extractFH>) {
		
		# Skip "#" and empty rows
		if ($_ =~ m/^#/ || $_ eq "") { next;}
		
		chomp;
		
		# Hash up reads with read id as key.
		if ($_ =~ m/^>/) {
			$key = $_;
			$extract{$key} = <$extractFH>;
		}
	}
}
elsif($illumina) {
	while (<$extractFH>) {
		# Skip "#" and empty rows
		if ($_ =~ m/^#/ || $_ eq "") { next;}
		
		chomp;
		
		# Hash up reads with read id as key.
		# (Illumina readid always ends with "/1" or "/2")
		if ($_ =~ m/^([@\+].*?\/[12]).*?$/) {
			$extract{$1} = <$extractFH>;
		}
	}
}
$extractFH->close();

# Open output file for writing and name file for input.
my $outputFH = new FileHandle(">".$outputFile);
my $inFH = new FileHandle($infile);
while (<$inFH>) {
	
	# Skip "#" and empty rows
	if ($_ =~ m/^#/ || $_ eq "") { next; }
	
	chomp;
	
	if ($solid) {
		print $outputFH ">".$_."\n".$extract{">".$_}."\n";
	}
	elsif($illumina) {
		print $outputFH "@".$_."\n".$extract{"@".$_};
		print $outputFH "+".$_."\n".$extract{"+".$_};
	}
	
}
$inFH->close();
$outputFH->close;


# Print the usage help for this script.
sub usage {
  print "\nUsage: $0 -i <infile> -o <output file> -e <extract file>\n
  -i	Input file, list of read names to extract without \">\" or \"@\" in the begining.
  -o	Output file (stdout is default if not set).
  -e	File to extract from (either *.csfasta, *.qual or *.fastq)
  
  -solid	If it's a solid run.
  -illumina	If it's an illumina run.\n\n";
  exit(1);
}
