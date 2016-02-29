#!/usr/bin/perl -w

# Converts NC chromosome ID from reference file with chr id everywhere in 
# ampregion file
#
# By: Linus Forsmark, Jun 2011

use strict;
use warnings;
use FileHandle;

#my ($ampFh,$refFh,$selFh)=(new FileHandle($ARGV[0]),new FileHandle($ARGV[1]),new FileHandle($ARGV[2]));
my ($ampFh,$refFh)=(new FileHandle($ARGV[0]),new FileHandle($ARGV[1]));
my @conv;

# generate array with entries to convert
while(<$refFh>) {
	chomp;

	# parse line
	my(@arr)=split(/\t/,$_);

	# add to convert array
	push(@conv,\@arr);
}
close($refFh);

# loop input file to convert
while(<$ampFh>) {
	chomp;

	my($out,undef)=split(/\t/,$_);

	# loop over all lines and replace all matches
	for my $ref (@conv) {
		my($chrNr,$NCNR,$alignRef)=@{$ref};
		$out =~ s/$NCNR/$chrNr/;

		#replace ref id with correct chr
#		my($chr,undef)=split(/\#/,$alignRef);
#		$out =~ s/$chr/$chrNr/;
	}
	print $out." ";
}
close($ampFh);

#while(<$selFh>) {
#	chomp;
#	my($id,$chr,$start,$end)=split(/\t/,$_);
#	print "$id#$chr#$start#$end#".(abs($end-$start)+1)." ";
#}
#close($selFh);

# ref: chr1 NC_000001.10 chr1#NC_000001.10#1#249250621#-1 249250621

