#!/usr/bin/perl
#
# Input SNPmania files, download ensembl information, output filtered variants
# By: Linus Forsmark, Aug 2011
#
# ideas:
# option not to download ensembl information if we just want basic filters
# option to force consitency (e.g. so that roi and 1k allele information is removed if -m or -r is not specified)
# make a * to mark if variant is in exon, in hapmap sequence or in hapmap/normal sequencing file (not removing), good if updating to database later on
# option to overwrite the file (everything in hash when outputting, should be ok, perhaps already working?) and update log with that information
# take a folder as input and take all SNPmania files and perform all of the filters and generate a SAMPLE-list instead
# change to filtering at snpmania ref allele ratio, and also report allele frequency
# accept VCF as input and output

# filters:
# check for a third allele above a threshold, e.g. 10%
# filtering away normal variants based on a threshold span (e.g. ratio 0.25 in sample, remove if 0.15-0.35 in normal, to remove some recurrent positions)

#load uppmax modules
#module load bioinfo-tools
#module load BioPerl/1.6.1_PERL5.10.1
#PERL5LIB=${PERL5LIB}:${HOME}/modules/lib64/perl5
#PERL5LIB=${PERL5LIB}:${HOME}/modules/ensembl/modules
#PERL5LIB=${PERL5LIB}:${HOME}/modules/ensembl-compara/modules
#PERL5LIB=${PERL5LIB}:${HOME}/modules/ensembl-variation/modules
#PERL5LIB=${PERL5LIB}:${HOME}/modules/ensembl-functgenomics/modules
#export PERL5LIB
#
#for testing:
#zcat /bubo/home/h20/elinfalk/data/1000Genomes/alleleFrequencies_interim_phase1.frq.gz | $HOME/program/SERA/bin/perlscript/snpManiaToVariationFile.pl -i $HOME/glob/private/b2011080/SNPmania/C1.h.sapiens.uniq.insertions -d $HOME/glob/private/b2011080/SNPmania/C1.h.sapiens.uniq.deletions -v $HOME/glob/private/b2011080/SNPmania/C1.h.sapiens.uniq.variations -t 0.1 -m 0.1 -n 0.2 -h 20 -r /bubo/home/h20/elinfalk/data/HapMap/NA12802.genotypes_CEU_r27_nr.b36_fwd.lifted2hg19.txt -a 0.05 -b /dev/stdin


use FileHandle;
use strict;
use warnings;
use Getopt::Long;
#use Bio::EnsEMBL::Registry;

# declare variables
my($insFile,$delFile,$varFile,$thrFilter,$depVarFilter,$depIndFilter,$hapmapFilter,$allelicThrFilter,$allelicFile,$normalFile,$roiFile,$prevOutFile,$prevOutNormFile,$locationFilter,$severityFilter,$logFile,$outFile,$coveredVariantsFile,$ignoreFetch,$clusterThr);
my(%ref,%out,%ignore,%hap,%af,%roi,@medians,%roiNames,%hgncNames,%tmpRoiNames);

GetOptions(
	"i=s" => \$insFile,
	"d=s" => \$delFile,
	"v=s" => \$varFile,
	"t=s" => \$thrFilter,
	"ds=s" => \$depVarFilter,
	"di=s" => \$depIndFilter,
	"h=s" => \$hapmapFilter,
	"m=s" => \$allelicThrFilter,
	"n=s" => \$normalFile,
	"a=s" => \$allelicFile,
	"r=s" => \$roiFile,
	"s!" => \$severityFilter,
	"p=s" => \$prevOutFile,
	"q=s" => \$prevOutNormFile,
	"o=s" => \$outFile,
	"f=s" => \$logFile,
	"k=s" => \$coveredVariantsFile,
	"l!" => \$locationFilter,
	"c=s" => \$clusterThr,
	"--ignore_fetch" => \$ignoreFetch
) || usage("Illegal parameters provided");

# specify default values
$logFile="/dev/stdout" if($outFile && !$logFile);
$logFile="filter.log" if(!$outFile && !$logFile);
$outFile="/dev/stdout" if(!$outFile);
my $afExist=0;	# keeping track whether allelefrequencies have been reported in old hash

sub usage {
  my $warn=shift;
  print "\nUsage: $0, Warning: $warn\n
 -v jSNPmania variation file
 -d jSNPmania deletion file [Optional]
 -i jSNPmania insertion file [Optional]
 -p 
 output file [Optional]
 -h HapMap reference sequence (remove if variant found) [Optional]
 -r ROI file (marking locations) [Recommended]
 -l Output only ROI locations
 -s Output only severe predictions
 -n SNPmania files for removing identical SNVs (tumor normal / hashmap sample) [Optional]
 -q Previous output files for eliminating identical SNVs (tumor normal / hashmap sample) [Optional]
 -t Threshold value (ratio variant_depth/total_depth) [Optional]
 -m Allelic frequency threshold [Optional]
 -a Allelic frequency input file (frq format) [Optional]
 -c Window size for SNP clusters
 -ds Minimum read depth for SNVs [Optional]c
 -di Minimum read depth for indels [Optional]
 -o Output file [Default: /dev/stout]
 -f Log file [Default: /dev/stout | filter.log]\n\n";
  exit 1;
}

# error handling
usage("ROI-file required if set to filter ROI") if($locationFilter && !$roiFile);	# roifile must be specified for location
usage("Indel files must be specified if threshold set") if($depIndFilter && !$insFile && !$delFile);	# depth threshold set for indels but missing indel files
usage("Found no input file") if(!$varFile && !$prevOutFile); # found no input file, neither SNPmania nor previous outfile
usage("Specified files not found") if( ($varFile && !-e $varFile) || ($insFile && !-e $insFile) || ($delFile && !-e $delFile) || ($prevOutFile && !-e $prevOutFile) || ($allelicFile && !-e $allelicFile) || ($roiFile && !-e $roiFile) || ($hapmapFilter && !-e $hapmapFilter) || ($prevOutNormFile && !-e $prevOutNormFile) );
usage("Insertion file not found for previous normal output file") if( $insFile && $normalFile && !-e "$normalFile.insertion" );
usage("Deletion file not found for previous normal output file") if( $delFile && $normalFile && !-e "$normalFile.deletion" );
#usage("Program must fetch allelic frequencies from file or ensembl if allelic threshold set") if( $allelicThr && $ignoreFetch && !$allelicFile);

# prepare files for output
open(LOG,">$logFile") or die("Can't output log to $logFile");

#print LOG "Input line: $0";
#print $_ for(@ARGV);
#print LOG "\n\n";

# execute program:
print LOG "Starting program at ".returnDate()."\n";
print LOG "===========\n";

print LOG returnTime()." Marking locations and validating HGNC names with $roiFile\n" if($roiFile);
my $changedHgnc=hashROI($roiFile) if($roiFile);

# check input files and parse hash correspondingly
print LOG returnTime()." Regenerating previous output data from $prevOutFile\n" if($prevOutFile);
regenerateHash($prevOutFile,\%out) if($prevOutFile);

# add data from all snpmania files
print LOG returnTime()." Hashing variations from $varFile\n" if($varFile);
my($inRoiTotal,$roiTotal)=hashVariations($varFile,\%out) if($varFile);

print LOG returnTime()." Hashing insertions from $insFile\n" if($insFile);
hashIndels($insFile,"insertion",\%out) if($insFile);

print LOG returnTime()." Hashing deletions from $delFile\n" if($delFile);
hashIndels($delFile,"deletion",\%out) if($delFile);

# these should only be run if a snpmaniafile is provided unless they don't have values in maf etc
print LOG returnTime()." Removing hapmap variants from $hapmapFilter\n" if($hapmapFilter);
hashHapmap($hapmapFilter) if($hapmapFilter);

print LOG returnTime()." Hashing allelic frequencies from $allelicFile\n" if($allelicFile);
hashAF($allelicFile) if($allelicFile);

# apply basic filters (sequencing filters, dependent on sequencing data)
print LOG returnTime()." Applying filters on sequencing data\n";
print LOG indent()."Removing variations at a sequence depth of $depVarFilter \n" if($depVarFilter);
print LOG indent()."Removing indels at a sequence depth of $depIndFilter\n" if($depIndFilter);
print LOG indent()."Removing variants with $thrFilter ratio threshold\n" if($thrFilter);
print LOG indent()."Removing ampregion variants\n" if($locationFilter && $roiFile);
applySeqFilters(\%out);

# get the coverage of 1k variants
my $covered1k;
$covered1k=reportCovered1kVariants(0,$coveredVariantsFile) if($coveredVariantsFile);

# remove variations from normal/hapmap sample
print LOG returnTime()." Removing variants from normal file\n" if($prevOutNormFile || $normalFile);
purgeIdenticalSNVs($normalFile,$prevOutNormFile) if($normalFile || $prevOutNormFile);

# fetches data from ensembl and updates the output hash, only if either maf should be fetched or
#if(!$ignoreFetch) {
print LOG returnTime()." Fetching information from ensembl\n";
print LOG indent()."Predicting consequences\n";
#fetchEnsemblData();
#}
# fill with empty and set status not to download
#else {

# apply filters based on position data
print LOG returnTime()." Applying position-based filters\n";
print LOG indent()."Removing variants present in hapmap sequence\n" if($hapmapFilter);
print LOG indent()."Removing non-severe predictions\n" if($severityFilter);
print LOG indent()."Removing variants with allelic frequencies above $allelicThrFilter\n" if($allelicThrFilter);
applyPosFilters(\%out);

# print
print LOG returnTime()." Printing output to $outFile\n\n";
my $c=printOutput();

# Extra information
print LOG "Additional results:\n" if($changedHgnc || $roiTotal);
print LOG "===========\n";
print LOG "Total ROI coverage with $depVarFilter minimum depth: ".($inRoiTotal/$roiTotal)."\n" if($roiTotal);
print LOG "Covered 1k Genome variations: ".($covered1k)."\n" if($covered1k);
print LOG "$changedHgnc\n" if($changedHgnc);

print LOG "===========\n";
print LOG "Ending program at ".returnTime().", $c variants reported, log available at $logFile\n";

# close files
close(LOG);






# fetches domain features for this domain
sub getDomainFeatures {

	# define transcript variation
	my($tv)=shift;

#	my @exons = @{ $tran->get_all_Exons() };
	
	# use only coding transcripts and if we have a start and end
	if($tv->transcript()->translation() && $tv->translation_start && $tv->translation_end) {

		for my $domain (grep {$_->analysis()->logic_name() =~ m/^pf/} @{$tv->transcript()->translation()->get_all_DomainFeatures()}) {

			# return the name if found, could be extended to include the location and favor pfam before pfscan
			# returns as long as at least one end is covered by the transcript
			if( ( $tv->translation_start => $domain->start && $tv->translation_start <= $domain->end ) || ( $tv->translation_end >= $domain->start && $tv->translation_end <= $domain->end ) ) {
				return $domain->idesc();
			}
		}
	}
	return;
}


# subroutine to extract sift and polyphen scores from transcript, send in transcriptvariation
sub getSiftPolyphen {

	my($tv,$var)=@_;

	# get all variation alleles for this transcript variation
#	my $tva = $tv->get_TranscriptVariationAllele_for_allele_seq($$var{'variantBase'});
	my @tvas=$tv->get_all_alternate_TranscriptVariationAlleles();

#	die("Found multiple TranscriptVariationAlleles for this TranscriptVariation using a VariationFeature") if(scalar(@tvas)>1);

#	print "Transcript: ".$tva->transcript()->display_id()."\n";

#	print "Polyphen".$tva->polyphen_score().", Sift ".$tva->sift_score().", Condel ".$tva->condel_score()."\n";

	my($codonChange,$aaChange,$polScore,$siftScore,$condScore,$polPred,$condPred,$siftPred)=("-","-","-","-","-","-","-","-");

	for my $tvArr (@tvas) {

#		die("Found multiple TranscriptVariationAlleles in return array even though it should get only the alternate") if(scalar(@$tvArr)>1);

#		print $("Found multiple TranscriptVariationAlleles in return array even though it should get only the alternate") if(scalar(@$tvArr)>1);
		# get the only element in array
#		my $tva=shift(@$tvArr);

		# loop over all
		for my $tva (@{$tvArr}) {

			# safety check that we are not looking at erronous alleles and transcripts
#			die("When checking alleles for SIFT, Condel and PolyPhen, variant ".$$var{'refBase'}."/".$$var{'variantBase'}." do not agree with ".$tva->allele_string()) if($$var{'refBase'}."/".$$var{'variantBase'} ne $tva->allele_string());

			if($tva && $$var{'refBase'}."/".$$var{'variantBase'} eq $tva->allele_string()) {

				$aaChange=$tva->pep_allele_string() if($tva->pep_allele_string());
				$codonChange=$tva->display_codon_allele_string() if($tva->display_codon_allele_string());

				$polScore=$tva->polyphen_score() if($tva->polyphen_score());
				$siftScore=$tva->sift_score() if($tva->sift_score());
				$condScore=$tva->condel_score() if($tva->condel_score());

				$polPred=$tva->polyphen_prediction() if($tva->polyphen_prediction());
				$siftPred=$tva->sift_prediction() if($tva->sift_prediction());
				$condPred=$tva->condel_prediction() if($tva->condel_prediction());

	#			print "amino acid change: ", $aaChange, "\n";
	#			print "resulting codon: ", $codonChange, "\n";
	#			print "PolyPhen prediction: ", $polPred, "\n";
	#			print "PolyPhen score: ", $polScore, "\n";
	#			print "SIFT prediction: ".$siftPred."\n";
	#			print "SIFT score: ".$siftScore."\n";
	#			print "CONDEL prediction: ".$condPred."\n";
	#			print "CONDEL score: ".$condScore."\n";
	#			print "DISPLAY_CODON: ".$tva->display_codon."\n";
	#			print "DISPLAY_CODON_ALLELE_STRING: ".$tva->display_codon_allele_string."\n";
	#			print "CODON ALLELEL STRING: ".$tva->codon_allele_string."\n";
	#			print "allele_string: ".$tva->allele_string."\n";
	#			print "\n";
			}
		}
	}

	return($codonChange,$aaChange,$polScore,$siftScore,$condScore,$polPred,$siftPred,$condPred);
}


# parses through a region file and add information to all variations found
sub hashROI {

	# check the hgnc name for all rois from ensembl to avoid non-consistency
	my $reg = 'Bio::EnsEMBL::Registry';
	$reg->load_registry_from_db(-host => 'ensembldb.ensembl.org', -user => 'anonymous');

	# fetch adaptors
	my $ga = $reg->get_adaptor('human', 'core', 'gene');
	my $sa = Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "slice" );

	# loop over ROI file
	my $fh = new FileHandle(shift);
	while(<$fh>) {

		chomp;

		my($roiId,$ncChr,$start,$end)=split(/\t/,$_);

		# parse roiid to get the roi name (assuming format "GENE_ROI_NR#")
		my ($roiname,undef,$roinumber)=split(/_/,$roiId);

#		print "$roiname $roinumber\n";

		my $chr=convertChrName($ncChr);

		# add roi information for all variations found in roi
		$roi{$chr}{$_}=1 for($start..$end);

#		print "Checking $roiname\n";

		# add to tmp hash for checking afterwards
		$tmpRoiNames{$roiname}=1;

		# check if this hgnc name exists in their name nomenclature
		if(!$roiNames{$roiname}) {
			
			# define the slices and get the gene
			my $slice = $sa->fetch_by_region('chromosome', $chr, $start, $end);
			my @genes = @{$ga->fetch_all_by_Slice($slice)};

#			print "$roiname found ".scalar(@genes)." genes\n";

			if(!scalar(@genes)) { die("Found no genes for ROI"); }

			elsif(scalar(@genes)==1) {

#				print "Entered first\n";

				my $hgnc;
				my $gene=shift(@genes);
				my @names = @{$gene->get_all_DBEntries("HGNC")};
				$hgnc=shift(@names)->display_id if(@names);

				# declare value for tmp array to identify contradicting reported hgnc for this roi
#				$tmp{$roiname}=$hgnc if(!$tmp{$roiname});

				# kill program if contradicting found
#				die("Found contradicting HGNC names in ROI for $roiname") if($tmp{$roiname} ne $hgnc);

				# update hgnc
				$roiNames{$roiname}=$hgnc if($hgnc);
				$hgncNames{$hgnc}=$roiname if($hgnc);

#				print "FOUND $hgnc\n";
			
			}
			# if we found many, check if any of them are specified as hgnc, else continue looking in nextcomming rois
			elsif(scalar(@genes)>=2) {

				# look in all genes for a hgnc that matches
				for my $gene (@genes) {

					my $hgnc;
					my @names = @{$gene->get_all_DBEntries("HGNC")};

					$hgnc=shift(@names)->display_id if(@names);

					# exit if hgnc name found
					if($hgnc && $hgnc eq $roiname) {
						$roiNames{$roiname}=$hgnc;
						$hgncNames{$hgnc}=$roiname if($hgnc);
					}
				}
			}
		}
	}
#	print LOG "\nVerifying ROI with HGNC names:\n";
	# loop over all and report hgnc

	my $changedHgnc;

	for my $id (keys %tmpRoiNames) {
#		print "tmproinames: ".$id."\t".$roiNames{$id}."\n";
		if(!$roiNames{$id}) { $roiNames{$id}="-"; $changedHgnc=$changedHgnc."HGNC for $id not found, all transcripts will be reported"; }
		elsif($id ne $roiNames{$id}) { $changedHgnc=$changedHgnc."$id changed to $roiNames{$id}\n"; }
	}

	return $changedHgnc;
}

# returns ampregion or roi from chr and pos
sub getLocation {
	my($chr,$pos)=@_;
	return "roi" if($roi{$chr}{$pos});
	return "ampregion" if($roiFile);
	return "";
}

# print timestamp
sub returnTime {
	my($sec,$min,$hour)=localtime(time);
	$hour="0$hour" if($hour<10);
	$min="0$min" if($min<10);
	$sec="0$sec" if($sec<10);
	return "$hour:$min:$sec";
}

sub returnDate {
	my(undef,undef,undef,$mday,$mon,$year)=localtime(time);
	$mday="0$mday" if($mday<10);
	$mon="0$mon" if($mon<10);
	return ($year+1900)."-".($mon+1)."-$mday";
}

sub indent {
	return "         - ";
}


# filters based on settings provided
sub applySeqFilters {
	
	# get hash to add output
	my $key=shift;
	my %hash=%$key;

	# remove with threshold, at specified sequence depth, hapmap sequence
	for my $chr (keys %hash) {
		for my $pos (keys %{$hash{$chr}}) {
			for my $b (keys %{$hash{$chr}{$pos}}) {

				# get variation key
				my $var=$hash{$chr}{$pos}{$b};
				my $hapPos=$hap{$chr}{$pos};

				# remove if not enough sequence depth
				if( ($depVarFilter && $$var{'depth'} < $depVarFilter && $b.$$var{'refBase'} !~ m/\-/) || ($depIndFilter && $$var{'depth'} < $depIndFilter && $b.$$var{'refBase'} =~ m/\-/) ) { delete($hash{$chr}{$pos}{$b}); }

				# remove based on read depth ratio, now no upper threshold
#				elsif($thrFilter && $$var{'variantDepth'}/$$var{'depth'} > 1-$thrFilter || $thrFilter && $$var{'variantDepth'}/$$var{'depth'} < $thrFilter) { delete($hash{$chr}{$pos}{$b}); }
				elsif($thrFilter && $$var{'variantDepth'}/$$var{'depth'} < $thrFilter) { delete($hash{$chr}{$pos}{$b}); }

				# remove based on location (dependent on sequence, not sequencing)
				elsif($roiFile && $$var{'location'} && $$var{'location'} ne "roi" && $locationFilter) { delete($hash{$chr}{$pos}{$b}); }

			}
		}
	}
}

sub applyPosFilters {

	# get hash to add output
	my $key=shift;
	my %hash=%$key;

	# remove with threshold, at specified sequence depth, hapmap sequence
	for my $chr (keys %hash) {
		for my $pos (keys %{$hash{$chr}}) {
			for my $b (keys %{$hash{$chr}{$pos}}) {

				my $var = $out{$chr}{$pos}{$b};
				my $hapPos=$hap{$chr}{$pos};

				# remove if in hapmap sequence (dependent on sequence, not sequencing)
				if($hapmapFilter && $hapPos && $hapPos =~ m/$b/) { delete($hash{$chr}{$pos}{$b}); }

				# check for snp clusters at this position
#				snpClusterFilter($chr,$pos,$clusterThr) if($clusterThr);

				# remove based on severity of consequences
#				if($severityFilter && $$var{'cons'} && $$var{'cons'} !~ m/STOP/ && $$var{'cons'} !~ m/FRAMESHIFT/ && $$var{'cons'} !~ m/NON_SYNONYMOUS/ && $$var{'cons'} !~ m/SPLICE/) { delete($hash{$chr}{$pos}{$b}); }

				# remove based on severity of consequences, also removes if hgnc can not be found (thus no relevant transcripts for severity)
				if($severityFilter && ( !$$var{'hgnc'} || ( $$var{'cons'} && $$var{'cons'} !~ m/STOP/ && $$var{'cons'} !~ m/FRAMESHIFT/ && $$var{'cons'} !~ m/NON_SYNONYMOUS/ && $$var{'cons'} !~ m/SPLICE/ ) ) ) { delete($hash{$chr}{$pos}{$b}); }

				# remove based on MAF (dependent on sequence, not sequencing)
				if($allelicThrFilter && $af{$chr}{$pos}{$b} && $af{$chr}{$pos}{$b}>$allelicThrFilter && $af{$chr}{$pos}{$b}>$allelicThrFilter) { delete($hash{$chr}{$pos}{$b}); }

				# check for snp clusters at this position
#				snpClusterFilter($chr,$pos,$clusterThr) if($clusterThr);

			}
		}
	}
}


# covnert chr name from NC to chr format
sub convertChrName {

	my $inp = shift;

	# fix some variables
	my @chr = $inp =~ m/NC\_0*(\d+)\./g;
	my $chr=shift(@chr);

	# adjust chrX chrY
	if($chr==23) { $chr="X"; }
	elsif($chr==24) { $chr="Y"; }

	return $chr;
}

# takes snpmania files from tumor normal or hapmap samples and removes identical snvs and indels in multiple samples
sub purgeIdenticalSNVs {

#	print "Entered purge identical\n";

	# parse input
	my($stub,$prev)=@_;

	# get snpmania files
	my($var,$ins,$del)=("$stub.variations","$stub.insertion","$stub.deletion") if($stub);

	# define hash for normal sample
	my %normal;

	# add everything from snpmaniafile, adds unneccessarily much, should be adjusted to increase performance? accepts both snpmania and filtered format

	print LOG indent()."Regenerating previous output data from $prev \n" if($prev);
	regenerateHash($prev,\%normal) if($prev);

	print LOG indent()."Hashing variations from $var \n" if($var);
	(undef,undef)=hashVariations($var,\%normal) if($var);

	print LOG indent()."Hashing insertions from $insFile \n" if($insFile && $var);
	print LOG indent()."Hashing deletions from $delFile \n" if($delFile && $var);
	hashIndels($del,"deletion",\%normal) if($delFile && $del);
	hashIndels($ins,"insertions",\%normal) if($insFile && $ins);

	# sort according to user criteria

	# this is the problem, applying sequence filter should not be performed if old prev out file?
#	print "NOT APPLYING SEQFILTERS\n" if(!$stub);
	print LOG indent()."Applying filters on sequencing data (as previous)\n" if($stub);
	applySeqFilters(\%normal) if($stub);

	print LOG indent()."Removing identical variants\n";
	# check for identical variations in %normal and compare to %out
	for my $chr (keys %out) {
		for my $pos (keys %{$out{$chr}}) {
			for my $b (keys %{$out{$chr}{$pos}}) {
#				print "Removing $chr $pos $b\n" if($normal{$chr}{$pos}{$b});
				delete($out{$chr}{$pos}{$b}) if($normal{$chr}{$pos}{$b});
			}
		}
	}
}

# hash all allelic frequncies reported from 1k Genomes from local file instead of using ensembl api
# should preferably only hash from regions where variations are found to spare some memory
# add fallback mode to gather from ensembl api
sub hashAF {

	my $fh=new FileHandle(shift);
	while(<$fh>) {

		chomp;

		# remove comment lines
		if(substr($_,0,5) ne "CHROM") {

			# parse line
			my($chr,$pos,$freqs) = $_ =~ m/^(\S+)\s+(\S+)\s+\S+\s+\S+\s+(.+)/g;

			# check that this position should be hashed
			if($ref{$chr}{$pos}) {

				# parse alleles and frequncies from freqs
				my @alleles=split(/\s+/,$freqs);
				for my $alleleStr (@alleles) {
					my($b,$fr)=split(/:/,$alleleStr);

					# add to hash
					$af{$chr}{$pos}{$b}=$fr;
				}
			}
		}
	}

	addAF();
}

# add alleleic frequencis to output
sub addAF {
	for my $chr (keys %out) {
		for my $pos (keys %{$out{$chr}}) {
			for my $b (keys %{$out{$chr}{$pos}}) {

				# update af and overwrites
				$out{$chr}{$pos}{$b}{'allelefrequency'}=$af{$chr}{$pos}{$b} if($af{$chr}{$pos}{$b});

			}
		}
	}
}


# fetch ensembl data for each variation, fetch maf and or consequences
sub fetchEnsemblData {

#	print "Fetching ensembl data\n";

	# connect to api
	my $reg = 'Bio::EnsEMBL::Registry';
	$reg->load_registry_from_db(-host => 'ensembldb.ensembl.org', -user => 'anonymous');

	# fetch classes consequences
	my $vfa = $reg->get_adaptor('human', 'variation', 'variationfeature');
	my $sa = $reg->get_adaptor('human', 'core', 'slice');
	my $va = $reg->get_adaptor('human', 'variation', 'variation');
	my $ga = $reg->get_adaptor('human', 'core', 'gene');
#	my $vaa = $reg->get_adaptor('homo_sapiens', 'variation', 'variationannotation');

	# loop over all positions
	for my $chr (keys %out) {
		for my $pos (keys %{$out{$chr}}) {
			for my $b (keys %{$out{$chr}{$pos}}) {

#				print "Fetching for $chr $pos $b\n";

				# get variation key
				my $v = $out{$chr}{$pos}{$b};

				# donwload for regenerated varations only if alle frequencies not repfetchEnsorted last time
				if(!$$v{'cons'} || ($$v{'cons'} && !$afExist && $allelicThrFilter && !$allelicFile) || $$v{'cons'} eq "-") {

					# define slices
					my $vfs = $sa->fetch_by_region('chromosome', $$v{'chr'});
					my $vfsPos = $sa->fetch_by_region('chromosome', $$v{'chr'}, $$v{'pos'}, $$v{'end'});

					# get variations from slice (variation api)
					my @vfs = @{$vfa->fetch_all_by_Slice($vfsPos)};

					# get cosmic froms lice
					push(@vfs,@{$vfa->fetch_all_somatic_by_Slice($vfsPos)});

					# get allele frequency from 1000 genome project if set, otherwise get from hash
					my $afr;
					$afr = gatherAF($v,\@vfs) if(!$allelicFile && $allelicThrFilter);

#					$af = $af{$chr}{$pos}{$b} if($allelicThr && $allelicFile && $af && $af{$chr}{$pos}{$b});

					# get the list of predicted consequnences (should a sort based on MAF fetch be implemented to skip these?)
					if(!$out{$chr}{$pos}{$b}{'cons'} || $out{$chr}{$pos}{$b}{'cons'} eq "-") {
						my($cons,$id,$hgnc,$trKey,$domain) = predictConsequences($v,\@vfs,$vfs,$vfa,$ga);
						
						# update output array with all of the results
						$out{$chr}{$pos}{$b}{'annotations'}=$id if($id);
						$out{$chr}{$pos}{$b}{'cons'}=$cons;
						$out{$chr}{$pos}{$b}{'hgnc'}=$hgnc;

#						print $_."\n" for(@$trKey);

						# fet all of the predicted tr consequences 
						$out{$chr}{$pos}{$b}{'codonChange'}=$$trKey[0];
						$out{$chr}{$pos}{$b}{'aaChange'}=$$trKey[1];
						$out{$chr}{$pos}{$b}{'polyphenScore'}=$$trKey[2];
						$out{$chr}{$pos}{$b}{'polyphenPrediction'}=$$trKey[5];
						$out{$chr}{$pos}{$b}{'siftScore'}=$$trKey[3];
						$out{$chr}{$pos}{$b}{'siftPrediction'}=$$trKey[6];
						$out{$chr}{$pos}{$b}{'condelScore'}=$$trKey[4];
						$out{$chr}{$pos}{$b}{'condelPrediction'}=$$trKey[7];
						$out{$chr}{$pos}{$b}{'locationDomain'}=$domain;

					}

					# update output array with all of the results
					$out{$chr}{$pos}{$b}{'allelefrequency'}=$afr if($afr);

				}
			}
		}
	}
}


# de novo prediction of all consequences for this variation
# pass variation, list of variationfeatureobjects and a slice and adaptor for creating a new object if novel
sub predictConsequences {

	my($var,$vfs,$slice,$adaptor,$ga)=@_;

	# define to hold current variant
	my $variation;

	# str for merging annotations and consequences
	my $ann="";
	my $conslist="";
	my $hgnc;
	my $domain;
	my @severestTrCons;

	# check for variations containing the correct variation and add to valid variants
	for my $vf (@$vfs) {

		# det kan vara flera på en allele_string, te.x. A/G/G, borde checka båda

		# check if the allele_string have the same reference as our reported, please observe that for some there are weird strings
		 if($vf->allele_string && $vf->allele_string =~ m/\//) {

			my(undef,$varBase)=split(/\//,$vf->allele_string);

			# get the variation name
			$ann=$ann.$vf->variation_name."," if($vf->variation_name);

			# add the first entry (should be the rs-number)
			$variation=$vf if($varBase && $varBase eq $$var{'variantBase'} && !$variation);

		}
		# add the weird one too just bevause it is included in the official variant effect predictor script, but as it do not have a correct string, dont use it for the prediction
		elsif($vf->allele_string && $vf->allele_string !~ m/\// && $vf->variation_name) {
			$ann=$ann.$vf->variation_name.",";
		}
	}

	# if no variant was found, add a new one to array
	if(!$variation) {

		# adjust for indels
		my $str = $$var{'refBase'}."/".$$var{'variantBase'};
		$str = "-/".$$var{'variantBase'} if($$var{'refBase'} eq "-");
		$str = $$var{'refBase'}."/-" if($$var{'variantBase'} eq "-");

#		print "Creating new for ".$$var{'refBase'}."/".$$var{'variantBase'}."\n";
#		print "with $str\n";

		$variation = Bio::EnsEMBL::Variation::VariationFeature->new(
			-start => $$var{'pos'},
			-end => $$var{'end'},
			-slice => $slice,
			-allele_string => $str,
			-strand => 1,
			-map_weight => 1,
			-adaptor => $adaptor,
			-variation_name => "-"
		);
	}

	# get all consequences and transcripts 
#	my $cons="";
#	$cons="$cons$_," for(@{$variation->consequence_type()});

	foreach my $tv (@{$variation->get_all_TranscriptVariations}) {

		# get the associated hgnc names by fetching transcript and gene
		my @names = grep {$_->database eq 'HGNC'} @{$ga->fetch_by_transcript_stable_id($tv->transcript->stable_id)->get_all_DBEntries()};

		# if nothing found, grep all
		@names=@{$ga->fetch_by_transcript_stable_id($tv->transcript->stable_id)->get_all_DBEntries()} if(!@names);

		# update name and consequences if available in our defined roi list OR if locations have not been set
		if($hgncNames{$names[0]->display_id} || !$roiFile) {

			# remove all from conslist if they have been defined by elsif
			$conslist="" if(!$hgnc);

			# get the HGNC name if set
			$hgnc=$names[0]->display_id;

			# update list of consequences
			foreach my $cons (@{$tv->consequence_type}) {
				$conslist="$conslist$cons," if($conslist !~ m/$cons/);
			}

			# get the SIFT and PolyPhen predictions, in array for later parsing
			my @trCons=getSiftPolyphen($tv,$var);

			# check if the condel score is higher and change the transcript if it is
			@severestTrCons=@trCons if( (!$severestTrCons[0] || $severestTrCons[4] eq "-") || ($trCons[4] ne "-" && $trCons[4] > $severestTrCons[4]) || ( !$severestTrCons[0] && $severestTrCons[4] eq "-" && ( $severestTrCons[2] ne "-" || $severestTrCons[3] ne "-") ) );

			# check the protein domain for this transcript, should preferentially be done when defining genes and hashed instead
			$domain=getDomainFeatures($tv) if(!$domain);

		}
		# if there are no hgnc name for this roi, then just add all transcripts, as long as NONE of the following transcripts have no
		elsif(!$hgnc) {
			foreach my $cons (@{$tv->consequence_type}) {
				$conslist="$conslist$cons," if($conslist !~ m/$cons/);
			}
		}
	}

	# adjust domain if nothing found
	$domain="-" if(!$domain);

	# adjust hgnc to read - if not set
	$hgnc="-" if(!$hgnc);

	return(substr($conslist,0,length($conslist)-1),substr($ann,0,length($ann)-1),$hgnc,\@severestTrCons,$domain);

}

# pass the variant and a list of variationfeature objects and extract MAFs from the 1000 genome project
# should look after a more suitable adaptor to fetch only the 1000 genomes and not having to sort manually for non-relevant ones
sub gatherAF() {

	print LOG indent()."Gathering allelic frequencies\n";

	my($varKey,$vfs)=@_;
	my $af="-";

	# loop over all variations found at this position
	foreach my $vf (@$vfs){
		foreach my $allele (@{$vf->variation->get_all_Alleles()}) {

			$af=$allele->frequency if($$varKey{'variantBase'} eq $allele->allele && $vf->strand()==1 && $allele->population && $allele->population->name =~ /1000GENOMES.+CEU/);

			# check that this is the correct variation at this position and return if it is
			return $af;
		}
	}
}


# get slices for contigous regions for ensembl API (ampregions)
sub getSlices {

	# store slices
	my @slices;

	# loop over reference region
	for my $chr (keys %ref) {

		my($start,$end);
		for my $pos (sort {$a<=>$b} keys %{$ref{$chr}}) {

			# set start position
			$start=$pos if(!$start);

			# store contig in slices and reset
			if(!$ref{$chr}{$pos+1}) {

				$end=$pos;
				my @contig=($chr,$start,$end);

				push(@slices,\@contig);
				$start=0;

			}
		}
	}
	return @slices;
}

# hash hapmap variations
sub hashHapmap {

	my $fh = new FileHandle(shift);

	while(<$fh>) {

		if(substr($_,0,3) ne "rs#") { 

			# parse hapmap file
			chomp;
			my($chr,$pos,$a1,$a2)=($_ =~ m/chr(\w+)\s(\d+)\s.+\s.+\s.+\s(\D)(\D)/g);

			# add to reference hash
			$hap{$chr}{$pos}="$a1$a2" if($ref{$chr}{$pos});
		}
	}
}

# add all variations to hash
sub hashVariations {

	my($file,$key)=@_;
	my $fh=new FileHandle($file);

	# array to keep track of ROI coverage
	my(@inRoi)=(0,0);

	while(<$fh>) {

		chomp;
		if(substr($_,0,1) ne "#") { 

			#split line
			my($depth,undef,$refAlleleRatio,$ncChr,$pos,$refBase,$variants,undef,$corrAlleleRatio)=split(/\t/,$_);

			# fix some variables
			my $chr=convertChrName($ncChr);

			# check all variants
			my @base=("A","G","C","T");
			my @hits=(split(/\|/,$variants));

			for my $i (0..scalar(@base)-1) {

				my $b=$base[$i];
				my $h=$hits[$i];

				# do not add but skip if this base is annotated in reference or hapmap samples
				next if($b eq $refBase);

				# add to hash if there are enough depth and are not already hashed from regeneratedhash
				$$key{$chr}{$pos}{$b} = {
					'chr' => $chr,
					'pos' => $pos,
					'end' => $pos,
					'refBase' => $refBase,
					'variantBase' => $b,
					'variantDepth' => $h,
					'depth' => $depth,
					'location' => getLocation($chr,$pos),
					'info' => $variants,
					'alleleRatio' => $corrAlleleRatio,
					'refAlleleRatio' => $refAlleleRatio
				} if($h && !$$key{$chr}{$pos}{$b});

			}

			# keep track of ROI coverage with given thr
			if($depVarFilter) {
				$inRoi[1]++ if($roi{$chr}{$pos});
				$inRoi[0]++ if($roi{$chr}{$pos} && $depth>=$depVarFilter);
			}

			# save all reference positions
			$ref{$chr}{$pos}=$refBase;
		}

	}
	return($inRoi[0],$inRoi[1]) if($inRoi[1]);
	return;
}


# hash all indels from files
sub hashIndels {

	my($file,$type,$key)=@_;

	# get hash to add output
	my %hash=%$key;

	my $fh = new FileHandle($file);
	while(<$fh>) {
		
		chomp;
		if(substr($_,0,1) ne "#") {	

			# split line
			my($depth,$ncChr,$pos,$amount,undef,$indel) = split(/\t/,$_);

			# convert chr name
			my $chr=convertChrName($ncChr);

			# if indel is insertion and insertion file is specified in input
			if($type eq "insertion" && $indel && $depth) {

				# get all insertions found
				my @arr=split(/\|/,$indel);
				for my $ins (@arr) {

					# get depth and sequence
					my($indDepth,$indSeq) = $ins =~ m/(\d+)(\D+)/;

#					print "$indDepth $indSeq\n" if($indSeq && $indDepth && !$$key{$chr}{$pos}{"I$indSeq"});
#					print "Added $chr $pos ".($pos+1)."-/$indSeq with I$indSeq from file\n" if($indSeq && $indDepth && !$$key{$chr}{$pos}{"I$indSeq"});

					# add to hash
					$$key{$chr}{$pos+1}{"I$indSeq"} = {
						'chr' => $chr,
						'pos' => $pos+1,
						'end' => $pos,
						'refBase' => "-",
						'variantBase' => $indSeq,
						'variantDepth' => $indDepth,
						'depth' => $depth,
						'location' => getLocation($chr,$pos),
						'info' => $indel
					} if($indSeq && $indDepth && !$$key{$chr}{$pos}{"I$indSeq"});

				}
			}

			# if indel is deletion and deletion file provided
			elsif($type eq "deletion" && $indel && $depth) {
					
				my @arr=split(/\|/,$indel);
				for my $del (@arr) {

					my($indDepth,$offsetStart,$offsetEnd) = $del =~ m/(\d+)\((\S+),(\S+)\)/;

					# get the reference sequence of this deletion
					my $seq = getRefSeq($chr,$pos+$offsetStart,$offsetEnd-$offsetStart);

					# check that this should be added
					if($seq && $indDepth && !$ignore{$chr}{$pos} && !$$key{$chr}{$pos}{"D$seq"}) {

						$$key{$chr}{$pos+$offsetStart}{"D$seq"} = {
							'chr' => $chr,
							'pos' => $pos+$offsetStart,
							'end' => $pos+$offsetEnd,
							'refBase' => $seq,
							'variantBase' => "-",
							'variantDepth' => $indDepth,
							'depth' => $depth,
							'location' => getLocation($chr,$pos),
							'info' => $indel
						};

#						print "Added $chr ".($pos+$offsetStart)." ".($pos+$offsetEnd)." $seq/- with D$seq from file\n";

						# add the whole deletion to ignore list to avoid duplicates in this region
						$ignore{$chr}{$_}=1 for($pos+$offsetStart..$pos+$offsetEnd);
					}
				}
			}
		}
	}	
}


# get reference positions from hashed variation file
sub getRefSeq {
	my($chr,$pos,$length)=@_;
	my $seq="";

	for($pos..$pos+$length) {

		# return empty to exclude this deletion if it stretches beyond ampregion
		return if(!$ref{$chr}{$_});

		# concatenate sequence 
		$seq=$seq.$ref{$chr}{$_};
	}

	return $seq;
}

# very much to be done!
# print lines with given thresholds etc
sub printOutput {

	# print information about sorting, what was included etc, time and etc...
	# to be included...

	open(OUT,">$outFile") or die("Can't output data to $outFile");

	my $c=0;

	for my $chr (keys %out) {
		for my $pos (keys %{$out{$chr}}) {
			for my $b (keys %{$out{$chr}{$pos}}) {

				$c++;

				# get variant key
				my $var=$out{$chr}{$pos}{$b};

				print OUT $$var{'chr'}."\t" if($$var{'chr'});
				print OUT $$var{'pos'}."\t" if($$var{'pos'});
				print OUT $$var{'end'}."\t" if($$var{'pos'});
#				print OUT $$var{'end'}."\t" if($$var{'end'});
				print OUT $$var{'depth'}."\t" if($$var{'depth'});
				print OUT $$var{'variantDepth'}."\t" if($$var{'variantDepth'});
				print OUT $$var{'refBase'}."\t" if($$var{'refBase'});
				print OUT $$var{'variantBase'}."\t" if($$var{'variantBase'});
				if($$var{'allelefrequency'}) { print OUT $$var{'allelefrequency'}."\t"; } else { print OUT "-\t" }
				if($$var{'location'}) { print OUT $$var{'location'}."\t"; } else { print OUT "-\t"; }
				if($$var{'hgnc'}) { print OUT $$var{'hgnc'}."\t"; } else { print OUT "-\t"; }
				if($$var{'annotations'}) { print OUT $$var{'annotations'}."\t"; } else { print OUT "-\t"; }
				print OUT $$var{'info'}."\t" if($$var{'info'});
				print OUT $$var{'cons'}."\t" if($$var{'cons'});

				print OUT $$var{'codonChange'}."\t";;
				print OUT $$var{'aaChange'}."\t";
				print OUT $$var{'polyphenPrediction'}."\t";
				print OUT $$var{'polyphenScore'}."\t";
				print OUT $$var{'siftPrediction'}."\t";
				print OUT $$var{'siftScore'}."\t";
				print OUT $$var{'condelPrediction'}."\t";
				print OUT $$var{'condelScore'}."\t";
				print OUT $$var{'locationDomain'}."";

				print OUT "\n";
			}
		}
	}

	close(OUT);
	return $c;
}

# parses old input file previously generated as output from program, parse the new file and extract a subhash with values to update
sub regenerateHash {
	
	my($file,$key)=@_;

#	print "Entered regeneratehash with $file...\n";

	my $fh=new FileHandle($file);

	# parse previous output file
	while(<$fh>) {

		chomp;
		my($chr,$pos,$end,$depth,$variantDepth,$refBase,$variantBase,$refAlleleRatio,$corrAlleleRatio,$allelefrequency,$location,$hgnc,$annotations,$info,$cons,$codonChange,$aaChange,$polPred,$polScore,$siftPred,$siftScore,$condPred,$condScore,$locationDomain)=split(/\t/,$_);

		# adjust the hash key if deletion
		my $keyBase=$variantBase;
		$keyBase="D$refBase" if($variantBase eq "-");
		$keyBase="I$variantBase" if($refBase eq "-");

		# remove old values to remain consistency with new sort
#		$allelefrequency="-" if(!$allelicThrFilter);
#		$location="-" if(!$roiFile);

		# update new values for all old if additional optional changes not taken care of in sort
		$location=getLocation($chr,$pos) if($roiFile);

		# add everything to hash, prev marks it as belonging to old hash for adjustments later
		$$key{$chr}{$pos}{$keyBase} = {
			'chr' => $chr,
			'pos' => $pos,
			'end' => $end,
			'refBase' => $refBase,
			'variantBase' => $variantBase,
			'variantDepth' => $variantDepth,
			'depth' => $depth,
			'refAlleleRatio' => $refAlleleRatio,
			'alleleRatio' => $corrAlleleRatio,
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
			'alleleRatio' => $corrAlleleRatio,
			'locationDomain' => $locationDomain
		};

#		print "Regenerated $chr $pos $end $refBase/$variantBase with with $keyBase\n" if($refBase eq "-" || $variantBase eq "-");

		# adjust "globals"
		$afExist=1 if($allelefrequency ne "-");
	
		# save all reference positions
		$ref{$chr}{$pos}=$refBase;
		
#		print "Adding $chr $pos $variantBase\n";

#		print "In regenerating: value ".$$key{$chr}{$pos}{$variantBase}{'info'}." is set";

	}

	$fh->close();
}

# check a specified window for multiple SNPs at this position (takes HapMap variants into account)
# could be even more preferable to remove snp clusters PRIOR to all filters (in case of bad region, one SNP could easily escape if all other are filtered away by chance when filtering like this below)
sub snpClusterFilter {

	my($chr,$pos,$window)=@_;

	my(@variants);

	# check for multiple snps in this window
	my $skip;
	for my $i ($pos-$window..$pos+$window) {

		for my $b (keys %{$out{$chr}{$i}}) {

			# skip insertions and deletions that are more then 1bp 
			 if($b !~ m/[I|D]/ || ($b =~ m/[I|D]/ && length($out{$chr}{$i}{$b}{'variantBase'}.$out{$chr}{$i}{$b}{'refBase'}) <= 2)) {

				# to remove later
				push(@variants,$i) if($out{$chr}{$i});

				# check for known snps
				if($out{$chr}{$i}{$b}{'allelefrequency'} ne "-") { 
					$skip=1;
					last;
				}
			}
		}
	}

	# now, evaluate which to remove, we must make sure that none of them are SNPs if we did not remove them
	if(!$skip && scalar(@variants)>1) {
		for my $i (@variants) {
			delete($out{$chr}{$i});
		}
	}
}
	

# special subroutine to report a file with the amount of 1000k variants are detected
sub reportCovered1kVariants {

	my($thr,$outputfile)=@_;

	open(AF,">$outputfile") or die("Can't output 1k genome project data to $outputfile");
	open(AFPOS,">$outputfile.positions") or die("Can't output 1k genome project data to $outputfile.positions");

	# get percentage
	my @c=(0,0);

	# loop over all of the hashed 1k genome files
	for my $chr (keys %af) {
		for my $pos (keys %{$af{$chr}}) {
			for my $b (keys(%{$af{$chr}{$pos}})) {

				# only consider locations in ROI, enough frequent and not of the reference
				if($roi{$chr}{$pos} && $af{$chr}{$pos}>$thr && $b ne $ref{$chr}{$pos}) {

					$c[1]++;

					print AF "$chr\t$pos\t".$ref{$chr}{$pos}."\t$b\t".$af{$chr}{$pos}{$b};

					my $v=$out{$chr}{$pos}{$b};

					# get the variant if found in our analysis
					if($out{$chr}{$pos}{$b}) {

						$c[0]++;

						print AF "\t";
						print AF $$v{'variantDepth'}/$$v{'depth'}."\t";
						print AF $$v{'variantDepth'}."\t";
						print AF $$v{'depth'}."\t";

					}

					# just print the threshold at this position
					print AFPOS $$v{'variantDepth'}/$$v{'depth'}."\n";

					print AF "\n";
				}
			}
		}
	}

	close(AF);

	return($c[0]/$c[1]) if($c[1]);
	return;
}
