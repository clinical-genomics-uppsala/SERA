#!/usr/bin/perl -w

# Uses a SAM file and output the amount of reads aligning to each reference.
# Useful when extracting junction hits.
#
# By: Linus, 2010
#

use strict;
use warnings;
use Getopt::Long;
use FileHandle;

# parse SAM file and fetch files
my($file,$output);
GetOptions("i=s" => \$file,"o=s" => \$output);
if (!$output) { $output="/dev/stdout"; }

my $fh = new FileHandle($file);
my $output_fh = new FileHandle($output,">");
my (@refs, %reads, %hash);

# go through SAM file and add to array/hash
while(<$fh>) {
	chomp;
	if ($_ =~ m/^\@SQ/) {
		my @parse = split(/\t/, $_);
		push(@refs, substr($parse[1],3));
	}
	if ($_ !~ m/^\@/) {
		my @parse = split(/\t/,$_);
		my @parse2 = split(/#/, $parse[2]);
		push(@{$reads{$parse2[1]}},getRef($parse[2]));
	}
}

# calculate how many reads for each reference
for my $ref (@refs) {
	my @parse = split(/#/, $ref);
	my $chr = $parse[1];

	for my $read (@{$reads{$chr}}) {
#		print "Comparing ".getRef($ref)." and $read\n";
		if(getRef($ref) eq $read) {
			$hash{$ref}++;
		}
	}
}

# print output for all fragments
for my $ref (@refs) {
	my $reads = 0;
	if($hash{$ref}) { $reads = $hash{$ref}; }
	print $output_fh getRef($ref)."\t$reads\n";
}

$fh->close();
$output_fh->close();

sub getRef {
	my $inp = shift(@_);
	my @arr = split(/#/,$inp);
	return $arr[0];
}
