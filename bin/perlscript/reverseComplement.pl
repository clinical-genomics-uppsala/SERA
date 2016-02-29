#!/usr/bin/perl
#
# Script to return the reverse complement of a dna sequence
# 

use warnings;
use strict;

print revComp(shift);

sub revComp {
	my $seq = reverse(lc(shift));
	$seq =~ tr/acgt/tgca/;
	return uc($seq);
}
