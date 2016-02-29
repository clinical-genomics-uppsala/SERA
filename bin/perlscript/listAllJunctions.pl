#!/usr/bin/perl

# This script takes a list with ID and number of junction hits and a list of all junctions
# and outputs all junctions in the design 

use warnings;
use strict;
#use criticism 'harsh';

use FileHandle;
use Carp;

use Readonly;

# Subroutine prototypes
sub usage;

my ($next_arg, $junctionFile, $selectionFile, $outputFile);

if(scalar(@ARGV) == 0){
    usage();
}
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-j")    { $junctionFile = shift(@ARGV); }
    elsif($next_arg eq "-s")    { $selectionFile = shift(@ARGV); }
    elsif($next_arg eq "-o")    { $outputFile = shift(@ARGV); }
    else { print "Invalid argument: $next_arg"; usage(); }
}

# Error
if (!$junctionFile || !$selectionFile || !$outputFile) { usage(); }

Readonly my $ID => 0;

# Open junction list file
open(my $junction_fh, '<', $junctionFile) or croak "Couldn't open junction hit file ".$junctionFile."!";
my %hash;

# Go through the file and put ID as key
while (<$junction_fh>) {
	# Check if it's a comment or a an empty line
	if ($_ =~ m/^#/ || $_ =~ m/^\n/ || $_ =~m/^\r/) { next; }
	chomp;
	my @jh_line = split /\t/, $_;
	$hash{$jh_line[$ID]} = $jh_line[1];
}
close $junction_fh or croak "Couldn't close junction hit file ".$junctionFile."!";

# Open target- and output file handler 
my $selection_fh = FileHandle->new($selectionFile) or croak "Couldn't open selection file ".$selectionFile."!";
my $outputFH = FileHandle->new(">".$outputFile) or croak "Couldn't open output file ".$outputFile."!";
# Go through the target file
while (<$selection_fh>) {
	
	# If it's a comment print directly to file
	if ($_ =~ m/^\#/) {
		print $outputFH $_;
	}
	elsif ($_ =~ m/^\n/ || $_ =~m/^\r/) { next; }
	# If it's not a blank line check if the id on the line exists in the hash. If that is the case
	# print the line to file
	else {
		chomp;
		my @line = split /\t/, $_;
		
		if (exists $hash{$line[$ID]}) {
			print $outputFH $line[$ID]."\t".$hash{$line[$ID]}."\n";
		} 
		else {
			print $outputFH $line[$ID]."\t0\n";
		}
	}
}

close $selection_fh or croak "Couldn't close selection file ".$selectionFile."!";
close $outputFH or croak "Couldn't close output file ".$outputFile."!";


# Print the usage help for this script.
sub usage {
  print "\nUsage: $0 -j <junction hits file> -o <output file> -s <selection file>  \n 
 
 -j junction hits file with the same ID:s as in the selection file, tab delimited: {ID\tNumberOfHits}
 -o output file
 -s selection file\n\n";
  exit 1;
}