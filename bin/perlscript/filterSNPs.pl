#!/usr/bin/perl
#
# Input SNPmania files, download ensembl information, output filtered variants
# By: Linus Forsmark, Aug 2011
#
# ideas:
# option not to download ensembl information if we just want basic filters
# option to force consistency (e.g. so that roi and 1k allele information is removed if -m or -r is not specified)
# flags to mark if variants are in exons, hapmap sequence or in hapmap/normal sequencing file (not removing), good if updating to database later on
# take a folder as input and take all SNPmania files and perform all of the filters and generate a SAMPLE-list instead?
# filtering and reporting at snpmania ref allele and allele ratio rather then variant depth ratio
# accept VCF as input and output

# filters:
# check for a third allele above a threshold, e.g. 10%

# run in shell for testing without SLURM
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


#######################################
#
#	 PROGRAM SETUP
#
######################################

use FileHandle;
use strict;
use warnings;
use Getopt::Long;
#use Bio::EnsEMBL::Registry;

# declare variables for inputfile
my($insFile,$delFile,$varFile,$thrFilter,$depVarFilter,$depIndFilter,$hapmapFilter,$allelicThrFilter,$allelicFile,$normalFile,$roiFile,$prevOutFile,$prevOutNormFile,$locationFilter,$severityFilter,$logFile,$outFile,$coveredVariantsFile,$ignoreFetch,$clusterThr,$ratioDiff,$geneListFile);

# declare various lists
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
	"g=s" => \$geneListFile,
	"o=s" => \$outFile,
	"f=s" => \$logFile,
	"k=s" => \$coveredVariantsFile,
	"l!" => \$locationFilter,
	"c=s" => \$clusterThr,
	"rd=s" => \$ratioDiff,
	"--ignore_fetch" => \$ignoreFetch	# not implemented yet
) || usage("Illegal parameters provided");

# specify default values
$logFile="/dev/stdout" if($outFile && !$logFile);
$logFile="filter.log" if(!$outFile && !$logFile);
$outFile="/dev/stdout" if(!$outFile);
$ratioDiff=0 if(!$ratioDiff);
# print user information
sub usage {
  my $warn=shift;
  print "\nUsage: $0, Warning: $warn\n
 -v jSNPmania variation file
 -d jSNPmania deletion file [Optional]
 -i jSNPmania insertion file [Optional]
 -p Previous output file [Optional]
 -h HapMap reference sequence (remove if variant found) [Optional]
 -r ROI file (marking locations) [Recommended]
 -g Gene file (validating transcripts) [Recommended]
 -l Output only ROI locations
 -s Output only severe predictions
 -n SNPmania files for removing identical SNVs (tumor normal / hapmap sample) [Optional]
 -q Previous output files for eliminating identical SNVs (tumor normal / hapmap sample) [Optional]
 -t Threshold value (ratio variant_depth/total_depth) [Optional]
 -m Allelic frequency threshold [Optional]
 -a Allelic frequency input file (frq format) [Optional]
 -c Window size for SNP clusters [Optional]
 -ds Minimum read depth for SNVs [Optional]
 -di Minimum read depth for indels [Optional]
 -rd Filtering Tumor Normal on min difference in reference allele ratio [optional]
 -k Output file for reporting covered 1000 Genome Project variants [Optional]
 -o Output file [Default: /dev/stout]
 -f Log file [Default: /dev/stout | filter.log]\n\n";
  exit 1;
}

# error handling
usage("ROI-file required if set to filter ROI") if($locationFilter && !$roiFile);	# roifile must be specified for location
usage("Indel files must be specified if threshold set") if($depIndFilter && !$insFile && !$delFile);	# depth threshold set for indels but missing indel files
usage("Found no input file") if(!$varFile && !$prevOutFile); # found no input file, neither SNPmania nor previous outfile
usage("Specified files not found") if( ($varFile && !-e $varFile) || ($insFile && !-e $insFile) || ($delFile && !-e $delFile) || ($prevOutFile && !-e $prevOutFile) || ($allelicFile && !-e $allelicFile) || ($roiFile && !-e $roiFile) || ($hapmapFilter && !-e $hapmapFilter) || ($prevOutNormFile && !-e $prevOutNormFile) );
usage("Insertion file not found for previous normal output file") if( $insFile && $normalFile && !-e "$normalFile.insertions" );
usage("Deletion file not found for previous normal output file") if( $delFile && $normalFile && !-e "$normalFile.deletions" );
#usage("Program must fetch allelic frequencies from file or ensembl if allelic threshold set") if( $allelicThr && $ignoreFetch && !$allelicFile);

# prepare files for output
open(LOG,">$logFile") or die("Can't output log to $logFile");



#######################################
#
#	 PROGRAM EXECUTION
#
######################################
print LOG "Starting program at ".returnDate()."\n";
print LOG "===========\n";


# checking for HGNC names in ROI file
my $changedHgnc;
if($roiFile) {
	print LOG returnTime()." Marking locations and validating HGNC names with $roiFile\n";

	# variable for printing to log
	if($roiFile && $geneListFile) {
		$changedHgnc=hashROIWithGeneFile($roiFile,$geneListFile);
	}
	elsif($roiFile) {
		$changedHgnc=hashROI($roiFile);
	}
}


# check if previous output data was provided and hash it
if($prevOutFile) {
	print LOG returnTime()." Regenerating previous output data from $prevOutFile\n";
	regenerateHash($prevOutFile,\%out);
}


# add data from all input snpmania files for tumor
my($inRoiTotal,$roiTotal);
if($varFile) {
	print LOG returnTime()." Hashing variations from $varFile\n";
	($inRoiTotal,$roiTotal)=hashVariations($varFile,\%out) if($varFile);
}

if($insFile) {
	print LOG returnTime()." Hashing insertions from $insFile\n";
	hashIndels($insFile,"insertions",\%out);
}

if($delFile) {
	print LOG returnTime()." Hashing deletions from $delFile\n";
	hashIndels($delFile,"deletions",\%out) if($delFile);
}

# remove hapmap variants from known sequence file
if($hapmapFilter) {
	print LOG returnTime()." Removing hapmap variants from $hapmapFilter\n";
	hashHapmap($hapmapFilter) if($hapmapFilter);
}

# hash allele frequencies from input file
if($allelicFile) {
	print LOG returnTime()." Hashing allelic frequencies from $allelicFile\n";
	hashAF($allelicFile);
}


# apply basic filters (sequencing filters, dependent on sequencing data)
print LOG returnTime()." Applying filters on sequencing data\n";
print LOG indent()."Removing variations at a sequence depth of $depVarFilter \n" if($depVarFilter);
print LOG indent()."Removing indels at a sequence depth of $depIndFilter\n" if($depIndFilter);
print LOG indent()."Removing variants with $thrFilter ratio threshold\n" if($thrFilter);
print LOG indent()."Removing ampregion variants\n" if($locationFilter && $roiFile);
applySeqFilters(\%out);

# save data for covered variants in 1000 Genome Project to specified file
reportCovered1kVariants(0,$coveredVariantsFile) if($coveredVariantsFile);

# remove variations from normal sample (could also be used with hapmap sample)
if($prevOutNormFile || $normalFile) {
	print LOG returnTime()." Removing variants from normal file\n";
	purgeIdenticalSNVs($normalFile,$prevOutNormFile);
}

# fetches data from ensembl and updates the output hash
if(!$ignoreFetch) {
	print LOG returnTime()." Fetching information from ensembl\n";
	print LOG indent()."Predicting consequences\n";
	fetchEnsemblData();
}


# apply filters based on position data
print LOG returnTime()." Applying position-based filters\n";
print LOG indent()."Removing variants present in hapmap sequence\n" if($hapmapFilter);
print LOG indent()."Removing non-severe predictions\n" if($severityFilter);
print LOG indent()."Removing variants with allelic frequencies above $allelicThrFilter\n" if($allelicThrFilter);
applyPosFilters(\%out);


# print data to file and check the amount of variants
print LOG returnTime()." Printing output to $outFile\n\n";
my $c=printOutput();


# Append extra information to log
print LOG "Additional results:\n" if($changedHgnc || $roiTotal);
print LOG "===========\n";
print LOG "Total ROI coverage with $depVarFilter minimum depth: ".($inRoiTotal/$roiTotal)."\n" if($roiTotal);
print LOG "$changedHgnc\n" if($changedHgnc);

print LOG "===========\n";
print LOG "Ending program at ".returnTime().", $c variants reported, log available at $logFile\n";

# close files
close(LOG);




########################################
#
#	  SUBROUTINES
#
#########################################


#
# Fetches domain features for a given TranscriptVariation and return the domain name if found
#
sub getDomainFeatures {

	# define transcript variation
	my($tv)=shift;

	# use only coding transcripts and if we have a translation start and end
	if($tv->transcript()->translation() && $tv->translation_start && $tv->translation_end) {

		# loop over all domains for this translated protein
		for my $domain (grep {$_->analysis()->logic_name() =~ m/^pf/} @{$tv->transcript()->translation()->get_all_DomainFeatures()}) {

			# returns the domain name if found, could be extended to include the location and favor pfam before pfscan
			if( ( $tv->translation_start => $domain->start && $tv->translation_start <= $domain->end ) || ( $tv->translation_end >= $domain->start && $tv->translation_end <= $domain->end ) ) {
				return $domain->idesc();
			}
		}
	}
	return;
}



#
# Subroutine to extract and return sift and polyphen scores from a TranscriptVariation object
#
sub getSiftPolyphen {

	my($tv,$var)=@_;

	# get all variation alleles for this transcript variation
	my @tvas=$tv->get_all_alternate_TranscriptVariationAlleles();

	# define variables to overwrite if adjusted
	my($codonChange,$aaChange,$polScore,$siftScore,$condScore,$polPred,$condPred,$siftPred)=("-","-","-","-","-","-","-","-");

	# loop over all TranscriptVariationAlleles
	for my $tvArr (@tvas) {

		# loop over all TranscriptVariations
		for my $tva (@{$tvArr}) {

			# check that the object corresponds to our variation
			if($tva && $$var{'refBase'}."/".$$var{'variantBase'} eq $tva->allele_string()) {

				# get and update all of the sift, condel, polyphen values, could be extended to include a lot more information
				$aaChange=$tva->pep_allele_string() if($tva->pep_allele_string());
				$codonChange=$tva->display_codon_allele_string() if($tva->display_codon_allele_string());

				$polScore=$tva->polyphen_score() if($tva->polyphen_score());
				$siftScore=$tva->sift_score() if($tva->sift_score());
				$condScore=$tva->condel_score() if($tva->condel_score());

				$polPred=$tva->polyphen_prediction() if($tva->polyphen_prediction());
				$siftPred=$tva->sift_prediction() if($tva->sift_prediction());
				$condPred=$tva->condel_prediction() if($tva->condel_prediction());

			}
		}
	}

	# return the values
	return($codonChange,$aaChange,$polScore,$siftScore,$condScore,$polPred,$siftPred,$condPred);
}



#
# Parses through a roi file, validates (via ensembl) and adjustes the HGNC names, and marks location as ROI or ampregion by updating the %roi hash
# Returns a text string with information about the adjustments.
#
sub hashROI {

	# make ensembl connection
	my $reg = 'Bio::EnsEMBL::Registry';
	$reg->load_registry_from_db(-host => 'ensembldb.ensembl.org', -user => 'anonymous');

	# fetch adaptors
	my $ga = $reg->get_adaptor('human', 'core', 'gene');
	my $sa = Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "slice" );

	# loop over ROI file
	my $fh = new FileHandle(shift);
	while(<$fh>) {

		if ($_ =~ m/^#/ || $_ eq "") { next; }
		chomp;

		# parse line
		my($roiId,$ncChr,$start,$end)=split(/\t/,$_);

		# parse roiid to get the roi name (assuming format "GENE_ROI_NR#")
		my ($roiname,undef,$roinumber)=split(/_/,$roiId);

		# convert the chr name from nc to chr
		my $chr=convertChrName($ncChr);

		# mark position as belonging to roi
		$roi{$chr}{$_}=1 for($start..$end);

		# add to tmp hash for post-checks in subroutine
		$tmpRoiNames{$roiname}=1;

		# check if this hgnc name exists in their name nomenclature
		if(!$roiNames{$roiname}) {
			
			# define the ensembl slices and get the genes for this ROI
			my $slice = $sa->fetch_by_region('chromosome', $chr, $start, $end);
			my @genes = @{$ga->fetch_all_by_Slice($slice)};

			# basic error handling, every ROI must belong to at least one gene
			if(!scalar(@genes)) { die("Found no genes for ROI"); }

			# loop over the found genes
			elsif(scalar(@genes)==1) {

				# get the hgnc name
				my $hgnc;
				my $gene=shift(@genes);
				my @names = @{$gene->get_all_DBEntries("HGNC")};
				$hgnc=shift(@names)->display_id if(@names);

				# update hash with hgnc for easy access later in program
				$roiNames{$roiname}=$hgnc if($hgnc);
				$hgncNames{$hgnc}=$roiname if($hgnc);

			}

			# if we found many genes, check if ensembl hgnc name agrees
			elsif(scalar(@genes)>=2) {

				# look in all genes for a hgnc that matches
				for my $gene (@genes) {

					# get hgnc name
					my $hgnc;
					my @names = @{$gene->get_all_DBEntries("HGNC")};
					$hgnc=shift(@names)->display_id if(@names);

					# if hgnc agree with name in roi, update hashes 
					if($hgnc && $hgnc eq $roiname) {
						$roiNames{$roiname}=$hgnc;
						$hgncNames{$hgnc}=$roiname if($hgnc);
					}
				}
			}
		}
	}

	# declare object for storing return string with information about changes
	my $changedHgnc;

	# check the status of each roiname
	for my $id (keys %tmpRoiNames) {

		# check whether we found the roi name and update information
		if(!$roiNames{$id}) {
			$roiNames{$id}="-";
			$changedHgnc=$changedHgnc."HGNC for $id not found, all transcripts will be reported";
		}
		
		# update if the name was changed
		elsif($id ne $roiNames{$id}) {
			$changedHgnc=$changedHgnc."$id changed to $roiNames{$id}\n";
		}
	}

	return $changedHgnc;
}


# Parses through a ROI without provided ROI-names, also takes a list of genenames and tries to adjust the gene names. Also marks positions as belonging to ROI.
sub hashROIWithGeneFile {

	# make ensembl connection
	my $reg = 'Bio::EnsEMBL::Registry';
	$reg->load_registry_from_db(-host => 'ensembldb.ensembl.org', -user => 'anonymous');

	# fetch adaptors
	my $ga = $reg->get_adaptor('human', 'core', 'gene');
	my $sa = Bio::EnsEMBL::Registry->get_adaptor( "human", "core", "slice" );

	# input
	my($roifile,$genefile)=@_;

	# parsed gene names
	my %geneFileNames;
	
	# string for storing return information
	$changedHgnc="";

	# loop over gene file
	my $fh = new FileHandle($genefile);
	while(<$fh>) {

		chomp;

		# extract only the first field
		my ($genename,undef,undef,undef)=split(//,$_);

		# store for later adjustments
		$geneFileNames{$genename}=1;

		# connect and see whether this name agree with hgnc name in ensembl
		my $gene;
		$gene = @{$ga->fetch_by_display_label($genename)};

		# print information that this gene was not found
		$changedHgnc=$changedHgnc."HGNC for $gene not found, all transcripts will be reported" if(!$gene);

		# add to hash if it was found
		$hgncNames{$genename}=1 if($gene);

	}
	$fh->close();

	# loop over ROI file and try to find roi numbers with unique genes
	$fh = new FileHandle(shift);
	while(<$fh>) {

		chomp;

		# parse line
		my($roiId,$ncChr,$start,$end)=split(/\t/,$_);

		# convert the chr name from nc to chr
		my $chr=convertChrName($ncChr);

		# mark position as belonging to roi
		$roi{$chr}{$_}=1 for($start..$end);

		# define the ensembl slices and get the genes for this ROI
		my $slice = $sa->fetch_by_region('chromosome', $chr, $start, $end);
		my @genes = @{$ga->fetch_all_by_Slice($slice)};

		# basic error handling, every ROI must belong to at least one gene
		if(!scalar(@genes)) { die("Found no genes for ROI"); }

		# add to hgnc list if we found unique genes at one ROI
		elsif(scalar(@genes)==1) {

			# get the hgnc name
			my $hgnc;
			my $gene=shift(@genes);
			my @names = @{$gene->get_all_DBEntries("HGNC")};
			$hgnc=shift(@names)->display_id if(@names);

			# add to hash list if gene was found and not present in gene file
			if($hgnc && !$hgncNames{$hgnc}) {
				$hgncNames{$hgnc}=1;
				$changedHgnc=$changedHgnc."$hgnc appended, unique for at least one of the ROIs\n";
			}
		}
	}

	return $changedHgnc;
}



#
# Returns "ampregion" or "roi" from roi hash given chr and pos
#
sub getLocation {
	my($chr,$pos)=@_;
	return "roi" if($roi{$chr}{$pos});
	return "ampregion" if($roiFile);
	return "";
}

#
# Returns the timestamp as HH:MM:SS
#
sub returnTime {
	my($sec,$min,$hour)=localtime(time);
	$hour="0$hour" if($hour<10);
	$min="0$min" if($min<10);
	$sec="0$sec" if($sec<10);
	return "$hour:$min:$sec";
}

#
# Retuns the current date as YYYY-MM-DD
#
sub returnDate {
	my(undef,undef,undef,$mday,$mon,$year)=localtime(time);
	$mday="0$mday" if($mday<10);
	$mon="0$mon" if($mon<10);
	return ($year+1900)."-".($mon+1)."-$mday";
}

#
# Makes a simple indent for the log file
#
sub indent {
	return "         - ";
}

#
# Covnert chr name from NC to chr format
#
sub convertChrName {

	my $inp = shift;

	# fix variable
	my @chr = $inp =~ m/NC\_0*(\d+)\./g;
	my $chr=shift(@chr);

	# adjust chrX chrY
	if($chr==23) { $chr="X"; }
	elsif($chr==24) { $chr="Y"; }

	return $chr;
}



#
# Subroutine for performing filters based on settings provided by input file, requires a reference key to a hash to perform the operations (e.g. hash for storing tumor or normal)
#
sub applySeqFilters {
	
	# get hash to add output
	my $key=shift;
	my %hash=%$key;

	# loop over all current variants in hash
	for my $chr (keys %hash) {
		for my $pos (keys %{$hash{$chr}}) {
			for my $b (keys %{$hash{$chr}{$pos}}) {

				# get variation key
				my $var=$hash{$chr}{$pos}{$b};
				my $hapPos=$hap{$chr}{$pos};

				# remove variant if not enough sequence depth
				if( ($depVarFilter && $$var{'depth'} < $depVarFilter && $b.$$var{'refBase'} !~ m/\-/) || ($depIndFilter && $$var{'depth'} < $depIndFilter && $b.$$var{'refBase'} =~ m/\-/) ) { delete($hash{$chr}{$pos}{$b}); }

				# remove variant if not enough read depth ratio (no upper thr)
				elsif($thrFilter && $$var{'variantDepth'}/$$var{'depth'} < $thrFilter) { delete($hash{$chr}{$pos}{$b}); }

				# remove based on location (in roi or ampregion)
				elsif($roiFile && $$var{'location'} && $$var{'location'} ne "roi" && $locationFilter) { delete($hash{$chr}{$pos}{$b}); }

			}
		}
	}
}



#
# Subroutine for filters based on static positions in the targeted region (e.g. hapmap sequence, allele frequency in 1000 genome project etc). Requires ta reference key to a hash to perform the operations on.
#
sub applyPosFilters {

	# get hash to add output
	my $key=shift;
	my %hash=%$key;

	# loop over all variants in hash
	for my $chr (keys %hash) {
		for my $pos (keys %{$hash{$chr}}) {
			for my $b (keys %{$hash{$chr}{$pos}}) {

				# stores the key to the variant
				my $var = $out{$chr}{$pos}{$b};
				my $hapPos=$hap{$chr}{$pos};

				# remove if variant found in hapmap sequence
				if($hapmapFilter && $hapPos && $hapPos =~ m/$b/) { delete($hash{$chr}{$pos}{$b}); }

				# snp clusters should more approprtiately be removed after all variants have been iterated, not to affect the situation for adjacent variants

				# check for snp clusters, more appropriate to remove before additional filters?
#				snpClusterFilter($chr,$pos,$clusterThr) if($clusterThr);

				# remove based on severity of consequences, also removes if hgnc can not be found (thus no relevant transcripts for severity)
				if($severityFilter && ( !$$var{'hgnc'} || ( $$var{'cons'} && $$var{'cons'} !~ m/STOP/ && $$var{'cons'} !~ m/FRAMESHIFT/ && $$var{'cons'} !~ m/NON_SYNONYMOUS/ && $$var{'cons'} !~ m/SPLICE/ ) ) ) { delete($hash{$chr}{$pos}{$b}); }

				# remove based on allelic frequency in 1000 Genome Project
				if($allelicThrFilter && $af{$chr}{$pos}{$b} && $af{$chr}{$pos}{$b}>$allelicThrFilter && $af{$chr}{$pos}{$b}>$allelicThrFilter) { delete($hash{$chr}{$pos}{$b}); }

				# check for snp clusters at this position, based on variants after all filters have been performed
#				snpClusterFilter($chr,$pos,$clusterThr) if($clusterThr);
			}
		}
	}
}



#
# Subroutine to remove variants found in normal sample from tumor sample. Requires a file stub to SNPmania input files and an optional file of previous output data for normal sample. Performs all of the filters and compares to the tumor, removes all variants that coincide with the normal within a given threshold span for the variant depth ratio.
#
sub purgeIdenticalSNVs {

	# parse input
	my($stub,$prev)=@_;

	# get snpmania files
	my($var,$ins,$del)=("$stub.variations","$stub.insertions","$stub.deletions") if($stub);

	# define temp hash for normal sample
	my %normal;

	# add everything from snpmaniafile (could be adjusted to only add variants present in tumor sample)
	if($prev) {
		print LOG indent()."Regenerating previous output data from $prev \n";
		regenerateHash($prev,\%normal);
	}

	if($var) {
		print LOG indent()."Hashing variations from $var\n";
		(undef,undef)=hashVariations($var,\%normal);
	}

	if($insFile && $var) {
		print LOG indent()."Hashing insertions from $insFile\n" ;
		hashIndels($ins,"insertions",\%normal);
	}

	if($delFile && $var) {
		print LOG indent()."Hashing deletions from $delFile\n";
		hashIndels($del,"deletions",\%normal);
	}

	# apply filters only if previous outputfile was not specified as input
	if($stub) {
		print LOG indent()."Applying filters on sequencing data (as previous)\n";
		applySeqFilters(\%normal) if($stub);
	}

	# loop over all variants in the tumor hash
	print LOG indent()."Removing identical variants\n";
	for my $chr (keys %out) {
		for my $pos (keys %{$out{$chr}}) {
			for my $b (keys %{$out{$chr}{$pos}}) {

				# remove variant from tumor file if depth ratio for the normal is within the defined ratio threshold
				if ($normal{$chr}{$pos}{$b}
					&& ($out{$chr}{$pos}{$b}{'variantDepth'}/$out{$chr}{$pos}{$b}{'depth'}) <= ($normal{$chr}{$pos}{$b}{'variantDepth'}/$normal{$chr}{$pos}{$b}{'depth'}+$ratioDiff) 
					&& ($out{$chr}{$pos}{$b}{'variantDepth'}/$out{$chr}{$pos}{$b}{'depth'}) >= ($normal{$chr}{$pos}{$b}{'variantDepth'}/$normal{$chr}{$pos}{$b}{'depth'}-$ratioDiff)) {
					delete($out{$chr}{$pos}{$b}) if($normal{$chr}{$pos}{$b});
				}
			}
		}
	}
}


#
# Hash all allelic frequncies for roi reported from 1k Genomes from local file (instead of ensembl api)
#
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

#
# Add alleleic frequencies to output
#
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


#
# Fetch ensembl data for each variation (protein family, sift, condel, polyphen, cosmic, dbsnp, predicted consequence, 1000 Genome Project allele frequency)
#
sub fetchEnsemblData {

	# connect to api
	my $reg = 'Bio::EnsEMBL::Registry';
	$reg->load_registry_from_db(-host => 'ensembldb.ensembl.org', -user => 'anonymous');

	# fetch classes consequences
	my $vfa = $reg->get_adaptor('human', 'variation', 'variationfeature');
	my $sa = $reg->get_adaptor('human', 'core', 'slice');
	my $va = $reg->get_adaptor('human', 'variation', 'variation');
	my $ga = $reg->get_adaptor('human', 'core', 'gene');

	# loop over all positions
	for my $chr (keys %out) {
		for my $pos (keys %{$out{$chr}}) {
			for my $b (keys %{$out{$chr}{$pos}}) {

				# get variation key
				my $v = $out{$chr}{$pos}{$b};

				# download for regenerated varations only if alle frequencies not fetched last time
				if(!$$v{'cons'} || ($$v{'cons'} && $allelicThrFilter && !$allelicFile) || $$v{'cons'} eq "-") {

					# define slices
					my $vfs = $sa->fetch_by_region('chromosome', $$v{'chr'});
					my $vfsPos = $sa->fetch_by_region('chromosome', $$v{'chr'}, $$v{'pos'}, $$v{'end'});

					# get variations from slice (variation api)
					my @vfs = @{$vfa->fetch_all_by_Slice($vfsPos)};

					# get cosmic from slice
					push(@vfs,@{$vfa->fetch_all_somatic_by_Slice($vfsPos)});

					# get allele frequency from 1000 genome project if set, otherwise get from hash
					my $afr;
					$afr = gatherAF($v,\@vfs) if(!$allelicFile && $allelicThrFilter);

					# update variant with allelefrequency
					$out{$chr}{$pos}{$b}{'allelefrequency'}=$afr if($afr);

					# update variants if values are not present
					if(!$out{$chr}{$pos}{$b}{'cons'} || $out{$chr}{$pos}{$b}{'cons'} eq "-") {

						# get predictions and its corresponding hgnc gene and protein domain
						my($cons,$id,$hgnc,$trKey,$domain) = predictConsequences($v,\@vfs,$vfs,$vfa,$ga);
						
						# update output array with all of the results
						$out{$chr}{$pos}{$b}{'annotations'}=$id if($id);
						$out{$chr}{$pos}{$b}{'cons'}=$cons;
						$out{$chr}{$pos}{$b}{'hgnc'}=$hgnc;

						# fetch all of the predicted transcript consequences 
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
				}
			}
		}
	}
}


#
# Returns ensembl information for a variation. Requires the variation key, list of VariationFeature objects, slice and adaptor
#
sub predictConsequences {

	my($var,$vfs,$slice,$adaptor,$ga)=@_;

	# define to hold current variant
	my $variation;

	# str for merging annotations and consequences for all transcripts
	my $ann="";
	my $conslist="";
	my $hgnc;
	my $domain;
	my @severestTrCons;

	# check for variations containing the correct variation and add to valid variants
	for my $vf (@$vfs) {

		# check if the variation exists in ensembl, (notice the regexp, present to avoid unexpected format for strings that are sometimes returned)
		 if($vf->allele_string && $vf->allele_string =~ m/\//) {

			my(undef,$varBase)=split(/\//,$vf->allele_string);

			# get the variation name
			$ann=$ann.$vf->variation_name."," if($vf->variation_name);

			# add the first entry (should be the rs-number)
			$variation=$vf if($varBase && $varBase eq $$var{'variantBase'} && !$variation);

		}

		# add only the annotation for the unexpected string if found (it is present in the official variantEffectPrediction)
		elsif($vf->allele_string && $vf->allele_string !~ m/\// && $vf->variation_name) {
			$ann=$ann.$vf->variation_name.",";
		}
	}

	# if no variant was found in ensembl, create a new local object
	if(!$variation) {

		# adjust for indels
		my $str = $$var{'refBase'}."/".$$var{'variantBase'};
		$str = "-/".$$var{'variantBase'} if($$var{'refBase'} eq "-");
		$str = $$var{'refBase'}."/-" if($$var{'variantBase'} eq "-");

		# create object and store in variation variable
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

	# get all TranscriptVariations in order to extract consequence and additional information
	foreach my $tv (@{$variation->get_all_TranscriptVariations}) {

		# get the associated hgnc names by fetching transcript and gene
		my @names = grep {$_->database eq 'HGNC'} @{$ga->fetch_by_transcript_stable_id($tv->transcript->stable_id)->get_all_DBEntries()};

		# if nothing found, grep all
		@names=@{$ga->fetch_by_transcript_stable_id($tv->transcript->stable_id)->get_all_DBEntries()} if(!@names);

		# update name and consequences if the HGNC name is defined in ROI (considers automatically corrected ROI-names)
		if($hgncNames{$names[0]->display_id} || !$roiFile) {

			# remove all consequences from conslist if they have been defined by fallback elsif
			$conslist="" if(!$hgnc);

			# get the HGNC name if set
			$hgnc=$names[0]->display_id;

			# update list of consequences
			foreach my $cons (@{$tv->consequence_type}) {
				$conslist="$conslist$cons," if($conslist !~ m/$cons/);
			}

			# get SIFT and PolyPhen predictions and store in array for later processing
			my @trCons=getSiftPolyphen($tv,$var);

			# check if the condel score is higher and switch transcripts if it is (program only report values for ONE transcript)
			@severestTrCons=@trCons if( (!$severestTrCons[0] || $severestTrCons[4] eq "-") || ($trCons[4] ne "-" && $trCons[4] > $severestTrCons[4]) || ( !$severestTrCons[0] && $severestTrCons[4] eq "-" && ( $severestTrCons[2] ne "-" || $severestTrCons[3] ne "-") ) );

			# check the protein domain for this transcript, should preferentially be done when defining genes and hashed instead
			$domain=getDomainFeatures($tv) if(!$domain);

		}

		# if there are no hgnc name for this roi, use all transcripts as long as NONE of the following transcripts have defined hgnc names
		elsif(!$hgnc) {
			foreach my $cons (@{$tv->consequence_type}) {
				$conslist="$conslist$cons," if($conslist !~ m/$cons/);
			}
		}
	}

	# adjust domain if nothing was found
	$domain="-" if(!$domain);

	# adjust hgnc to read - if not set
	$hgnc="-" if(!$hgnc);

	return(substr($conslist,0,length($conslist)-1),substr($ann,0,length($ann)-1),$hgnc,\@severestTrCons,$domain);

}

#
# Returns the allelic frequency from the 1000 Genome Project from ensembl given a VariationFeature object and variation key. Could easilly be manually changed to include populations from the pilot studies etc.
# (a more appropriate adaptor may be available to avoid fetching additional non-relevant data)
#
sub gatherAF() {

	print LOG indent()."Gathering allelic frequencies\n";

	my($varKey,$vfs)=@_;
	my $af="-";

	# loop over all variations found at this position
	foreach my $vf (@$vfs){
		foreach my $allele (@{$vf->variation->get_all_Alleles()}) {

			# extract the allelic frequency
			$af=$allele->frequency if($$varKey{'variantBase'} eq $allele->allele && $vf->strand()==1 && $allele->population && $allele->population->name =~ /1000GENOMES.+CEU/);

			return $af;
		}
	}
}

#
# Hash hapmap variations from input file and store in hash
#
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

#
# Add all variations from input SNPmania file to hash. Requires a path to input file and a key to a hash to store data. Returns the amount of bases in file (ampregion) and amount of bases in ROI.
#
sub hashVariations {

	my($file,$key)=@_;
	my $fh=new FileHandle($file);

	# array to keep track of ROI coverage
	my(@inRoi)=(0,0);

	# loop over snpmania file
	while(<$fh>) {

		chomp;
		if(substr($_,0,1) ne "#") { 

			#split line
			my($depth,undef,undef,$ncChr,$pos,$refBase,$variants,undef,undef)=split(/\t/,$_);

			# skip if this position have no sequencing depth
			next if($depth==0);

			# fix chr name
			my $chr=convertChrName($ncChr);

			# check all bases at this position
			my @base=("A","G","C","T");
			my @hits=(split(/\|/,$variants));

			for my $i (0..scalar(@base)-1) {

				my $b=$base[$i];
				my $h=$hits[$i];

				# skip if the base is annotated in the reference genome
				next if($b eq $refBase);

				# add to hash if variant is not already hashed from the regenerated file
				$$key{$chr}{$pos}{$b} = {
					'chr' => $chr,
					'pos' => $pos,
					'end' => $pos,
					'refBase' => $refBase,
					'variantBase' => $b,
					'variantDepth' => $h,
					'depth' => $depth,
					'location' => getLocation($chr,$pos),
					'info' => $variants
				} if($h && !$$key{$chr}{$pos}{$b});

			}

			# keep track of ROI coverage at specified depth threshold
			if($depVarFilter && $roiFile) {
				$inRoi[1]++ if($roi{$chr}{$pos});
				$inRoi[0]++ if($roi{$chr}{$pos} && $depth>=$depVarFilter);
			}

			# save all reference positions (ampregion)
			$ref{$chr}{$pos}=$refBase;
		}

	}
	return($inRoi[0],$inRoi[1]) if($inRoi[1]);
	return;
}


#
# Hash all indels from SNP-mania file, requires path to file, type of indel (insertion, deletion) and key to hash for storing output
#
sub hashIndels {

	my($file,$type,$key)=@_;

	# get hash to add output
	my %hash=%$key;

	# loop over snpmania-file
	my $fh = new FileHandle($file);
	while(<$fh>) {
		
		chomp;
		if(substr($_,0,1) ne "#") {	

			# split line
			my($depth,$ncChr,$pos,$amount,undef,$indel) = split(/\t/,$_);

			# skip if no sequencing depth
			next if($depth==0);

			# convert chr name
			my $chr=convertChrName($ncChr);

			# if indel is insertion and insertion file is specified in input
			if($type eq "insertions" && $indel) {

				# get all insertions found
				my @arr=split(/\|/,$indel);
				for my $ins (@arr) {

					# get depth and sequence
					my($indDepth,$indSeq) = $ins =~ m/(\d+)(\D+)/;

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

			# if indel is deletion and deletion file is provided
			elsif($type eq "deletions" && $indel && $depth) {
					
				my @arr=split(/\|/,$indel);
				for my $del (@arr) {

					# parse get all types of deletions at this position
					my($indDepth,$offsetStart,$offsetEnd) = $del =~ m/(\d+)\((\S+),(\S+)\)/;

					# get the reference sequence of this deletion
					my $seq = getRefSeq($chr,$pos+$offsetStart,$offsetEnd-$offsetStart);

					# check that this deletion should be added (ignore list to avoid duplicates when checking subsequent positions)
					if($seq && $indDepth && !$ignore{$chr}{$pos} && !$$key{$chr}{$pos}{"D$seq"}) {

						# add to hash
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

						# add the whole deletion to ignore list to avoid duplicates in this region
						$ignore{$chr}{$_}=1 for($pos+$offsetStart..$pos+$offsetEnd);
					}
				}
			}
		}
	}	
}


#
# Returns the reference sequence (from input snpmania-files) given chr, pos and length of the sequence to extract
#
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


#
# Prints the final and processed hash of variations
#
sub printOutput {

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

#
# Parses a previous output file and regenerates the hash, given a input file and key to hash for storing output
#
sub regenerateHash {
	
	my($file,$key)=@_;

	my $fh=new FileHandle($file);

	# parse previous output file
	while(<$fh>) {

		chomp;
		my($chr,$pos,$end,$depth,$variantDepth,$refBase,$variantBase,$allelefrequency,$location,$hgnc,$annotations,$info,$cons,$codonChange,$aaChange,$polPred,$polScore,$siftPred,$siftScore,$condPred,$condScore,$locationDomain)=split(/\t/,$_);

		# adjust the hash key if deletion
		my $keyBase=$variantBase;
		$keyBase="D$refBase" if($variantBase eq "-");
		$keyBase="I$variantBase" if($refBase eq "-");

		# remove old values to remain consistency with new sort
#		$allelefrequency="-" if(!$allelicThrFilter);
#		$location="-" if(!$roiFile);

		# update new values for all old positions if the optional choises was not specified last time
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


		# save all reference positions in case no snpmania file is provided as input
		$ref{$chr}{$pos}=$refBase;

	}

	$fh->close();
}

#
# Removes variants if multiple SNPs are found in a specified window spanning the variant (ignores windows with known SNPs from 1000 Genome Project)
#
sub snpClusterFilter {

	# get input
	my($chr,$pos,$window)=@_;

	# store variants encountered
	my(@variants);

	# check for snps in window
	my $skip;
	for my $i ($pos-$window..$pos+$window) {

		# check if variant is found in our hash of variants
		for my $b (keys %{$out{$chr}{$i}}) {

			# skip insertions and deletions that are longer then 1bp 
			 if($b !~ m/[I|D]/ || ($b =~ m/[I|D]/ && length($out{$chr}{$i}{$b}{'variantBase'}.$out{$chr}{$i}{$b}{'refBase'}) <= 2)) {

				# mark variant as present in window
				push(@variants,$i) if($out{$chr}{$i});

				# check if this position is a known SNP 
				if($out{$chr}{$i}{$b}{'allelefrequency'} ne "-") { 
					$skip=1;
					last;
				}
			}
		}
	}

	# if multiple SNPs were found, we must make sure that none of them are SNPs if we did not remove them
	if(!$skip && scalar(@variants)>1) {
		for my $i (@variants) {
			delete($out{$chr}{$i});
		}
	}
}
	
#
# Prints the read depth ratio for all detected variants that coincide with SNPs in the 1000 Genome Project. Output to separate file specified in command line.
#
sub reportCovered1kVariants {

	my($thr,$outputfile)=@_;

	open(AF,">$outputfile") or die("Can't output 1000 Genome Project data to $outputfile");

	# loop over all snps in the 1000 Genome Project hash
	for my $chr (keys %af) {
		for my $pos (keys %{$af{$chr}}) {
			for my $b (keys(%{$af{$chr}{$pos}})) {

				# only consider locations in ROI, enough frequent and not of the reference
				if($roi{$chr}{$pos} && $af{$chr}{$pos}>$thr && $b ne $ref{$chr}{$pos}) {

					# extract the variant from our variation hash
					my $v=$out{$chr}{$pos}{$b};

					# prints information for the SNP
					print AF "$chr\t$pos\t".$ref{$chr}{$pos}."\t$b\t".$af{$chr}{$pos}{$b};
					
					# print sequencing data information if variant found in our variation hash
					if($out{$chr}{$pos}{$b}) {
						print AF "\t";
						print AF $$v{'variantDepth'}/$$v{'depth'}."\t";
						print AF $$v{'variantDepth'}."\t";
						print AF $$v{'depth'}."\t";
					}

					print AF "\n";
				}
			}
		}
	}

	close(AF);
}

#
# Get slices for all contigous regions for ensembl API (ampregions)
# ABANDONED!
#
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
