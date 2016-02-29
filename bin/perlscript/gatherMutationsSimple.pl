#!/usr/bin/perl -w

# Script for gathering variations from several files and generate output for each gene, position and sample
#
# By: Linus Forsmark, Aug 2011

# testline on uppmax:
# perl /bubo/home/h20/elinfalk/program/SERA/bin/perlscript/gatherMutatedRegions.pl -s /bubo/home/h20/elinfalk/tmp/LarryRebeqa.withUce.roi -r /bubo/home/h20/elinfalk/glob/private/b2011080/refFiles/45-001-1.ampregion -i /bubo/home/h20/elinfalk/glob/private/b2011080/effectPrediction/ -o /bubo/home/h20/elinfalk/testGather/ -g /bubo/home/h20/elinfalk/CLL_analysis/samples.txt -h /bubo/home/h20/elinfalk/CLL_analysis/names.txt -e roi,severe,variations,indels

use warnings;
use strict;
use FileHandle;
use Getopt::Long;
use List::Util qw(max min);

# define globals
my($filepath,$outpath,$groupfile,%hash,%covered,%extra,%hgncHash);

# parse command line
GetOptions(
	"i=s" => \$filepath,
	"o=s" => \$outpath,
	"g=s" => \$groupfile,
	) || usage();

sub usage {
  print "\nUsage: $0 \n
 -i Path to variant files
 -o Path for output txt-files
 -g External information (format: sample_id\textra_info) (Optional)\n\n";
  exit 1;
}

# open directory
opendir(DIR,$filepath);
my @file_arr=readdir(DIR);

# prepare file with groups
getGroups($groupfile);

# generate hash with all variants for all files
hashAllVariations();

# print output
printOutput();


sub printOutput {

	# open files for output
#	open(POSITIONS,">$outpath/positions.txt");
	open(SAMPLES,">$outpath/samples.txt");
	open(GENES,">$outpath/genes.txt");
	open(INDELS,">$outpath/indels.txt");

	# print out all results for each ROI
	for my $chr (sort keys %hash) {
	
		# loop over all positions in ROI
		for my $pos (sort keys %{$hash{$chr}}) {

			# loop over all different types
			for my $type (keys %{$hash{$chr}{$pos}}) {

				# loop over all samples for this found at this position
				for my $sample (keys %{$hash{$chr}{$pos}{$type}}) {

					# get info about recurrent samples and indels
					my($re_snvs,$re_indels,$re_total,$different)=(0,0,0,0);

					$re_snvs=$covered{$chr}{$pos}{snvs} if($covered{$chr}{$pos}{snvs});
					$re_indels=$covered{$chr}{$pos}{indels} if($covered{$chr}{$pos}{indels});
					$re_total=$re_snvs+$re_indels;
					$different=keys(%{$hash{$chr}{$pos}});

					# get additional information from external file
					my $extra="";
					$extra=$extra{$sample} if($extra{$sample});


					my $var=$hash{$chr}{$pos}{$type}{$sample};

=pod
					print $$var{'sample'}."\n";
					print $$var{'hgnc'}."\n";
					print $$var{'chr'}."\n";
					print $$var{'pos'}."\n";
					print $$var{'location'}."\n";
					print $$var{'depth'}."\n";
					print $$var{'variantDepth'}."\n";
					print $$var{'refBase'}."\n";
					print $$var{'variantBase'}."\n";
					print $$var{'info'}."\n";
					print $$var{'annotations'}."\n";
					print $$var{'allelefrequency'}."\n";
					print $$var{'cons'}."\n";
					print $$var{'codonChange'}."\n";
					print $$var{'aaChange'}."\n";
					print $$var{'polyphenPrediction'}."\n";
					print $$var{'polyphenScore'}."\n";
					print $$var{'siftPrediction'}."\n";
					print $$var{'siftScore'}."\n";
					print $$var{'condelPrediction'}."\n";
					print $$var{'condelScore'}."\n";
					print $$var{'locationDomain'}."\n";

					print "$re_snvs\t$re_indels\t$re_total\t$different\n";

					print $extra."\n";

					my $out="";
=cut

					my $out = 
						$$var{'sample'}."\t".
						$$var{'hgnc'}."\t".
						$$var{'chr'}."\t".
						$$var{'pos'}."\t".
						$$var{'end'}."\t".
						$$var{'location'}."\t".
						$$var{'depth'}."\t".
						$$var{'variantDepth'}."\t".
						$$var{'refBase'}."\t".
						$$var{'variantBase'}."\t".
						$$var{'info'}."\t".
						$$var{'annotations'}."\t".
						$$var{'allelefrequency'}."\t".
						$$var{'cons'}."\t".

						$$var{'codonChange'}."\t".
						$$var{'aaChange'}."\t".
						$$var{'polyphenPrediction'}."\t".
						$$var{'polyphenScore'}."\t".
						$$var{'siftPrediction'}."\t".
						$$var{'siftScore'}."\t".
						$$var{'condelPrediction'}."\t".
						$$var{'condelScore'}."\t".
						$$var{'locationDomain'}."\t".

						"$re_snvs\t$re_indels\t$re_total\t$different\t".

						"$extra\n"
					;

					if($type =~ m/^[I|D]/) { print INDELS $out; }
					else { print SAMPLES $out; }

				}
			}
		}
	}

	# print gene file
	for my $hgnc (keys(%hgncHash)) {
		print GENES "$hgnc\t$hgncHash{$hgnc}\n";
	}

	close(GENES);
	close(INDELS);
	close(SAMPLES);
#	close(POSITIONS);
}



# returns the sampleids corresponding groupid from hash
sub getGroupId {
	my $sid=shift;
	return $extra{$sid};
}

# for dividing all samples into subgroups (requires input file in format: sample_id\tgroup_id);
sub getGroups {

	my $file=shift;

	my $fh = new FileHandle($file);
	while(<$fh>) {
		chomp;
		my($sampleid,$rest) = $_ =~ m/^(\S+)\s+(.+)$/g;

		$extra{$sampleid}=$rest;
	}
	$fh->close();
}



# hash all variations in all files
sub hashAllVariations {

	for my $file (@file_arr) {

		# loop over relevant files
		if($file =~ m/variants$/) {

			# get samplename id
			my($sampleid,undef) = split(/\./,$file);

			my $fh=new FileHandle("$filepath/$file");

			# parse previous output file
			while(<$fh>) {

				chomp;
				my($chr,$pos,$end,$depth,$variantDepth,$refBase,$variantBase,$allelefrequency,$location,$hgnc,$annotations,$info,$cons,$codonChange,$aaChange,$polPred,$polScore,$siftPred,$siftScore,$condPred,$condScore,$locationDomain)=split(/\t/,$_);

				# adjust to compensate for problems in the main program, some values not defined apparently
				$locationDomain="-" if(!$locationDomain);

				# adjust the hash key if deletion
				my $keyBase=$variantBase;
				$keyBase="D$refBase" if($variantBase eq "-");
				$keyBase="I$variantBase" if($refBase eq "-");

				# add everything to hash, prev marks it as belonging to old hash for adjustments later
				$hash{$chr}{$pos}{$keyBase}{$sampleid} = {
					'sample' => $sampleid,
					'chr' => $chr,
					'pos' => $pos,
					'end' => $end,
					'refBase' => $refBase,
					'variantBase' => $variantBase,
					'variantDepth' => $variantDepth,
					'depth' => $depth,
					'location' => $location,
					'info' => $info,
					'hgnc' => $hgnc,
					'annotations' => $annotations,
					'allelefrequency' => $allelefrequency,
					'cons' => $cons,
					'codonChange' => $codonChange,
					'aaChange' => $aaChange,
					'polyphenPrediction' => $polPred,
					'polyphenScore' => $polScore,
					'siftPrediction' => $siftPred,
					'siftScore' => $siftScore,
					'condelPrediction' => $condPred,
					'condelScore' => $condScore,
					'locationDomain' => $locationDomain
				};

				# fix insertion positions
				my($s,$e)=($pos,$end);
				($s,$e)=($end,$pos) if($pos>$end);

				# add to counters
				if($keyBase =~ m/^[I|D]/) {
					$covered{$chr}{$pos}{indels}++ for($s..$e);
				}
				else {
					$covered{$chr}{$pos}{snvs}++ for($s..$e);
				}

				# increase hgnc count for gene list
				$hgncHash{$hgnc}++;

			}
			$fh->close();

		}
	}
}
