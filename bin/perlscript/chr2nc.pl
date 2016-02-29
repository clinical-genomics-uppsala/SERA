#!/usr/bin/perl

# This script transforms chr to NC nomenclature or vice versa.
# The transformation file has to be in format chr{tab}NC_.....
# Elin Falk Sorqvist

use warnings;
use strict;

#use criticism 'harsh';
use FileHandle;
use Carp;

# Subroutine prototypes
sub usage;

my ($next_arg, $chr2nc, $infile, $col, $outfile, $nc2chr);


if(scalar(@ARGV) == 0){
    usage();
}
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-chr2nc")    { $chr2nc = shift(@ARGV); }
    elsif($next_arg eq "-nc2chr")    { $nc2chr = shift(@ARGV); }
    elsif($next_arg eq "-i")    { $infile = shift(@ARGV); }
    elsif($next_arg eq "-c")    { $col = shift(@ARGV); }
    elsif($next_arg eq "-o")    { $outfile = shift(@ARGV); }
    else { print "\nInvalid argument: $next_arg\n"; usage(); }
}

# Error
if (!$infile || !$outfile || !$col || !($chr2nc||$nc2chr) ) { usage(); }

my $in_fh = FileHandle->new($infile) or croak "Couldn't open input file ".$infile."!";
my $out_fh = FileHandle->new("> ".$outfile) or croak "Couldn't open output file ".$outfile."!";

my $chr2nc_fh;
if ($chr2nc) {
	$chr2nc_fh = FileHandle->new($chr2nc) or croak "Couldn't open chr2nc file ".$chr2nc."!";	
}
elsif($nc2chr) {
	$chr2nc_fh = FileHandle->new($nc2chr) or croak "Couldn't open chr2nc file ".$nc2chr."!";
}


my %chr=();
while (<$chr2nc_fh>) {
	if ($_ =~ m/^#/ || $_ eq "" )  { next; }
	chomp;
	
	my @chrLine = split(/\t/, $_);
	if ($chr2nc) {
		$chr{$chrLine[0]} = $chrLine[1];
	}
	elsif ($nc2chr) {
		$chr{$chrLine[1]} = $chrLine[0];
	}
	

} 

while(<$in_fh>) {
	if ($_ =~ m/^#/ || $_ eq "" )  { next; }
	chomp;
	
	my @line = split(/\t/, $_);
	
	if ($chr{$line[($col-1)]}) {
		
		my $tmp = $chr{$line[($col-1)]};
		$line[($col-1)] = $tmp;
		
		my $c=0;
		foreach (@line) {
			if($c == 0) {
				print {$out_fh} $_;
				$c++;
			}
			else {
				print {$out_fh} "\t".$_;
				$c++;
			}
		}
		print {$out_fh} "\n";
	}
	else {
		print "Unknown chromosome: ".$_."\n";
	}
}	
	
	
	


# Print the usage help for this script.
sub usage {
  print "
************************************************************************************************

 This script transforms the chromosome annotation from chr to NC nomenclature or vice versa.

************************************************************************************************

  \nUsage: ".$0." \n 

 -i Infile, tab-delimited
 -c Column in input file with chromosome
 -o Output file
 -chr2nc tab-delimited file: {chr1\tNC_000001.10} (can't be used together with -nc2chr)
 -nc2chr tab-delimited file: {chr1\tNC_000001.10} (can't be used together with -chr2nc)\n\n";
 exit 1;
}
1;