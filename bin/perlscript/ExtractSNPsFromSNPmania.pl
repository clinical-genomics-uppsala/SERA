#!/usr/bin/perl

use warnings;
use strict;
#use criticism 'harsh';

use FileHandle;
#use Carp;



# Subroutine prototypes
sub usage;

my ($next_arg, $regionFile, $SNPmaniaFile, $outputFile);

if(scalar(@ARGV) == 0){
    usage();
}
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-r")    { $regionFile = shift(@ARGV); }
    elsif($next_arg eq "-o")    { $outputFile = shift(@ARGV); }
    elsif($next_arg eq "-s")    { $SNPmaniaFile = shift(@ARGV); }
    else { print "Invalid argument: $next_arg"; usage(); }
}

# Error
if (!$regionFile || !$SNPmaniaFile) { usage(); }

if ( !$outputFile ) { $outputFile = "/dev/stdout"; }

my %region;
open(REGION,"<", $regionFile) or die "Couldn't open region file ".$regionFile."!";
while (<REGION>) {
	
	if($_ =~m/^#/ || $_ =~ m/^$/) { next; }
	
	chomp;
	my @line = split(/\t/, $_);
	$region{$line[1]}{$line[2]} = $line[3];
}
close(REGION) or die "Couldn't close region file ".$regionFile."!";

# Open region list file
open(SNPmania, "<", $SNPmaniaFile) or die "Couldn't open SNPmania file ".$SNPmaniaFile."!";
open(OUTPUT, "> ".$outputFile) or die "Oops, could not open output file: $!";

while (<SNPmania>) {

	if($_ =~m/^#/ || $_ =~ m/^$/) { next; }

	chomp;
	my @line = split(/\t/, $_);
	for my $start (keys %{$region{$line[3]}}) {
		if ($line[4]>=$start && $line[4]<=$region{$line[3]}{$start}) {
			print OUTPUT $_."\n";
		}
	}
}

# Print the usage help for this script.
sub usage {
  print "\nUsage: $0\n 
 
 -o Output file
 -r Region file
 -s SNPmania output file\n\n";
  exit 1;
}