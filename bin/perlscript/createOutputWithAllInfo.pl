#!/usr/bin/perl

use warnings;
use strict;

use FileHandle;

# Subroutine prototypes
sub usage;

my ($next_arg, $exonicVariantFile, $variantFile, $outputDir, $outputFile, $File1000g, $snpFile, $snpFlagFile );

if(scalar(@ARGV) == 0){
    usage();
}
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-e")    { $exonicVariantFile = shift(@ARGV); }
    elsif($next_arg eq "-v")    { $variantFile = shift(@ARGV); }
    elsif($next_arg eq "-fo")    { $outputDir = shift(@ARGV); }
    elsif($next_arg eq "-o")    { $outputFile = shift(@ARGV); }
    elsif($next_arg eq "-g")    { $File1000g = shift(@ARGV); }
    elsif($next_arg eq "-s")    { $snpFile = shift(@ARGV); }
    elsif($next_arg eq "-sf")    { $snpFlagFile = shift(@ARGV); }
    else { print "Invalid argument: $next_arg"; usage(); }
}

# Error
if ( !$variantFile ) { usage(); }

if ($outputDir && $outputFile) {
	print "ONLY one of -o and -fo can be set!\n";
	usage();
}
elsif (!$outputDir && !$outputFile) {
	print "One of -o and -fo has to be set!\n";
	usage();
}

my %infoHash;


if ($variantFile) {
	open(VARIANT, "<", $variantFile) or die "Couldn't open variant file ".$variantFile."!";
	
	while (<VARIANT>) {
		if($_ =~m/^#/ || $_ =~ m/^$/) { next; }
		chomp;
		my @variantLine = split(/\t/, $_);
		my @comments = split(/\s/, $variantLine[(scalar(@variantLine)-1)]);
		my $sample;
	
		my $transcripts;
		if ($variantLine[1] =~ m/\(/) {
			my ($tmp, $tran) = split(/\(/, $variantLine[1]);
			$variantLine[1] = $tmp;
			$tran =~ s/\)//;
			my @trans = split(/,/, $tran);
			$transcripts = join("\t", @trans); 
			
		}
		else {
			$transcripts = "-";
		}
		for (my $i=1; $i<scalar(@comments); $i++) {
			my @commentLine = split(/=/, $comments[$i]);
			if ($commentLine[0] =~ m/sample/ && !$infoHash{$commentLine[1]}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}) {
				$sample = $commentLine[1];
				$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]} = {
																					'end' => $variantLine[4],
																					'variantType' => $variantLine[0],
																					'gene' => $variantLine[1],
																					'refAllele' => $variantLine[5],
																					'transcripts' => $transcripts };
			}
			elsif ($commentLine[0] =~ m/variantAlleleRatio/i && !$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'variantAlleleRatio'}) {
				$commentLine[1] =~ s/\./,/;
				$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'variantAlleleRatio'} = $commentLine[1];
			}
			elsif ($commentLine[0] =~ m/alleleFreq/i && !$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'alleleFreq'}) {
				my @alleleFreqs = split(/,/, $commentLine[1]);
				$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'alleleFreq'} = join("\t", @alleleFreqs);
			}
			elsif ($commentLine[0] =~ m/filter/i && !$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'filter'}) {
				$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'filter'} = $commentLine[1];
			}
			elsif ($commentLine[0] =~ m/readDepth/i && !$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'readDepth'}) {
				$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'readDepth'} = $commentLine[1];
			}
			elsif ($commentLine[0] =~ m/StrandsByDesign/i && !$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'StrandsByDesign'}) {
				$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'StrandsByDesign'} = $commentLine[1];
			}
			elsif ($commentLine[0] =~ m/normalAfterSequencing/i && !$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'normalAfterSequencing'}) {
				$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'normalAfterSequencing'} = $commentLine[1];
			}
			elsif ($commentLine[0] =~ m/tumorAfterSequencing/i && !$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'tumorAfterSequencing'}) {
				$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'tumorAfterSequencing'} = $commentLine[1];
			}
			elsif ($commentLine[0] =~ m/Tumor_A/ && !$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Tumor_A'}) {
                                $infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Tumor_A'} = $commentLine[1];
                        }
			elsif ($commentLine[0] =~ m/Tumor_G/ && !$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Tumor_G'}) {
                                $infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Tumor_G'} = $commentLine[1];
                        }
			elsif ($commentLine[0] =~ m/Tumor_C/ && !$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Tumor_C'}) {
                                $infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Tumor_C'} = $commentLine[1];
                        }
			elsif ($commentLine[0] =~ m/Tumor_T/ && !$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Tumor_T'}) {
                                $infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Tumor_T'} = $commentLine[1];
                        }
			elsif ($commentLine[0] =~ m/Normal_A/ && !$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Normal_A'}) {
				$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Normal_A'} = $commentLine[1];
			}
			elsif ($commentLine[0] =~ m/Normal_G/ && !$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Normal_G'}) {
				$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Normal_G'} = $commentLine[1];
			}
			elsif ($commentLine[0] =~ m/Normal_C/ && !$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Normal_C'}) {
				$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Normal_C'} = $commentLine[1];
			}
			elsif ($commentLine[0] =~ m/Normal_T/ && !$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Normal_T'}) {
				$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Normal_T'} = $commentLine[1];
			}
			elsif ($commentLine[0] =~ m/Normal_Ins/ && !$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Normal_Ins'}) {
				$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Normal_Ins'} = $commentLine[1];
			}
			elsif ($commentLine[0] =~ m/Tumor_Ins/ && !$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Tumor_Ins'}) {
				$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Tumor_Ins'} = $commentLine[1];
			}
			elsif ($commentLine[0] =~ m/Normal_Del/ && !$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Normal_Del'}) {
				$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Normal_Del'} = $commentLine[1];
			}
			elsif ($commentLine[0] =~ m/Tumor_Del/ && !$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Tumor_Del'}) {
				$infoHash{$sample}{$variantLine[2]}{$variantLine[3]}{$variantLine[5]}{$variantLine[6]}{'Tumor_Del'} = $commentLine[1];
			}
		}
	}
}
close(VARIANT);

# If an exonic_variant_function file is given
if ($exonicVariantFile) {
	
	# Open filehandler and read through file
	open(EXONIC_VARIANT, "<", $exonicVariantFile) or die "Couldn't open variant file ".$exonicVariantFile."!";
	
	while (<EXONIC_VARIANT>) {
		
		if($_ =~m/^#/ || $_ =~ m/^$/) { next; }
		chomp;
		# Split every line on tab and split the last column on comma
		my @exonicVariantLine = split(/\t/, $_);
		my @exonicComments = split(/\s/, $exonicVariantLine[(scalar(@exonicVariantLine)-1)]);
		
		# Variable to store the exonicSample name in
		my $exonicSample;
		
		# Go through all fields in the comment column
		for (my $i=1; $i<scalar(@exonicComments); $i++) {
			my @exonicCommentLine = split(/=/, $exonicComments[$i]);
			
			# If the name of the field is exonicSample save the name in $exonicSample			
			if ($exonicCommentLine[0] =~ m/sample/) {
				$exonicSample = $exonicCommentLine[1];
				# Split the given transcripts on comma
				my @transcript = split(",", $exonicVariantLine[2]);

				# Check if the comment info already is stored if not - store it
#				print $exonicCommentLine[1]."\t".$exonicVariantLine[3]."\t".$exonicVariantLine[4]."\t".$exonicVariantLine[6]."\t".$exonicVariantLine[7]."\n";
				if (!$infoHash{$exonicCommentLine[1]}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}) {
					$infoHash{$exonicCommentLine[1]}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]} = {
																					'end' => $exonicVariantLine[5],
																					'refAllele' => $exonicVariantLine[6],
																					'exonicType' => $exonicVariantLine[1],
																					'transcripts' => join("\t", @transcript)
																					 };												
				}
				# If the comment info already exists
				else {
					# Save info about the exonic type and the transcripts in the hash
					$infoHash{$exonicCommentLine[1]}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'exonicType'} = $exonicVariantLine[1];
					$infoHash{$exonicCommentLine[1]}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'transcripts'} = join("\t", @transcript);
				}
			}
			elsif ($exonicCommentLine[0] =~ m/variantAlleleRatio/i && !$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'variantAlleleRatio'}) {
				$exonicCommentLine[1] =~ s/\./,/;
				$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'variantAlleleRatio'} = $exonicCommentLine[1];
			}
			elsif ($exonicCommentLine[0] =~ m/alleleFreq/i && !$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'alleleFreq'}) {
				my @alleleFreqs = split(/,/, $exonicCommentLine[1]);
				$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'alleleFreq'} = join("\t", @alleleFreqs);
			}
			elsif ($exonicCommentLine[0] =~ m/filter/i && !$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'filter'}) {
				$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'filter'} = $exonicCommentLine[1];
			}
			elsif ($exonicCommentLine[0] =~ m/readDepth/i && !$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'readDepth'}) {
				$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'readDepth'} = $exonicCommentLine[1];
			}
			elsif ($exonicCommentLine[0] =~ m/StrandsByDesign/i && !$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'StrandsByDesign'}) {
				$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'StrandsByDesign'} = $exonicCommentLine[1];
			}
			elsif ($exonicCommentLine[0] =~ m/normalAfterSequencing/i && !$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'normalAfterSequencing'}) {
				$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'normalAfterSequencing'} = $exonicCommentLine[1];
			}
			elsif ($exonicCommentLine[0] =~ m/tumorAfterSequencing/i && !$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'tumorAfterSequencing'}) {
				$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'tumorAfterSequencing'} = $exonicCommentLine[1];
			}
			elsif ($exonicCommentLine[0] =~ m/Tumor_A/ && !$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'Tumor_A'}) {
                                $infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'Tumor_A'} = $exonicCommentLine[1];
                        }
                        elsif ($exonicCommentLine[0] =~ m/Tumor_G/ && !$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'Tumor_G'}) {
                                $infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'Tumor_G'} = $exonicCommentLine[1];
                        }
                        elsif ($exonicCommentLine[0] =~ m/Tumor_C/ && !$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'Tumor_C'}) {
                                $infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'Tumor_C'} = $exonicCommentLine[1];
                        }
                        elsif ($exonicCommentLine[0] =~ m/Tumor_T/ && !$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'Tumor_T'}) {
                                $infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'Tumor_T'} = $exonicCommentLine[1];
                        }
                        elsif ($exonicCommentLine[0] =~ m/Normal_A/ && !$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'Normal_A'}) {
                                $infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'Normal_A'} = $exonicCommentLine[1];
                        }
                        elsif ($exonicCommentLine[0] =~ m/Normal_G/ && !$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'Normal_G'}) {
                                $infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'Normal_G'} = $exonicCommentLine[1];
                        }
                        elsif ($exonicCommentLine[0] =~ m/Normal_C/ && !$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'Normal_C'}) {
                                $infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'Normal_C'} = $exonicCommentLine[1];
                        }
                        elsif ($exonicCommentLine[0] =~ m/Normal_T/ && !$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'Normal_T'}) {
                                $infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'Normal_T'} = $exonicCommentLine[1];
                        }
                        elsif ($exonicCommentLine[0] =~ m/Normal_Ins/ && !$infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'Normal_Ins'}) {
                                $infoHash{$exonicSample}{$exonicVariantLine[3]}{$exonicVariantLine[4]}{$exonicVariantLine[6]}{$exonicVariantLine[7]}{'Normal_Ins'} = $exonicCommentLine[1];
                        }	
		}
	}
}

my $G1000_log = "";
# If an exonic_variant_function file is given
if ($File1000g) {
	
	# Open filehandler and read through file
	open(FILE1000G, "<", $File1000g) or die "Couldn't open 1000 genome file ".$File1000g."!";
	
	while (<FILE1000G>) {
		
		if($_ =~m/^#/ || $_ =~ m/^$/) { next; }
		chomp;
		# Split every line on tab and split the last column on comma
		my @G1000Line = split(/\t/, $_);
		my @G1000Comments = split(/\s/, $G1000Line[(scalar(@G1000Line)-1)]);
		
		$G1000_log = $G1000Line[0];
		
		# Variable to store the G1000Sample name in
		my $G1000Sample;
		
		# Go through all fields in the comment column
		for (my $i=1; $i<scalar(@G1000Comments); $i++) {
			my @G1000CommentLine = split(/=/, $G1000Comments[$i]);
			
			# If the name of the field is G1000Sample save the name in $G1000Sample			
			if ($G1000CommentLine[0] =~ m/sample/) {
				$G1000Sample = $G1000CommentLine[1];
				# Split the given transcripts on comma
				my @transcript = split(",", $G1000Line[2]);

				# Check if the comment info already is stored if not - store it
#				print $G1000CommentLine[1]."\t".$G1000Line[2]."\t".$G1000Line[3]."\t".$G1000Line[6]."\n";
				if (!$infoHash{$G1000CommentLine[1]}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}) {
					$G1000Line[1] =~ s/\./,/;
					$infoHash{$G1000CommentLine[1]}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]} = {
																					'end' => $G1000Line[4],
																					'refAllele' => $G1000Line[5],
																					'1000g' => $G1000Line[1]
																					};												
				}
				# If the comment info already exists
				else {
					# Save info about the exonic type and the transcripts in the hash
					$G1000Line[1] =~ s/\./,/;
					$infoHash{$G1000CommentLine[1]}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'1000g'} = $G1000Line[1];
				}
			}
			elsif ($G1000CommentLine[0] =~ m/variantAlleleRatio/i && !$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'variantAlleleRatio'}) {
				$G1000CommentLine[1] =~ s/\./,/;
				$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'variantAlleleRatio'} = $G1000CommentLine[1];
			}
			elsif ($G1000CommentLine[0] =~ m/alleleFreq/i && !$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'alleleFreq'}) {
				my @alleleFreqs = split(/,/, $G1000CommentLine[1]);
				$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'alleleFreq'} = join("\t", @alleleFreqs);;
			}
			elsif ($G1000CommentLine[0] =~ m/filter/i && !$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'filter'}) {
				$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'filter'} = $G1000CommentLine[1];
			}
			elsif ($G1000CommentLine[0] =~ m/readDepth/i && !$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'readDepth'}) {
				$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'readDepth'} = $G1000CommentLine[1];
			}
			elsif ($G1000CommentLine[0] =~ m/StrandsByDesign/i && !$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'StrandsByDesign'}) {
				$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'StrandsByDesign'} = $G1000CommentLine[1];
			}
			elsif ($G1000CommentLine[0] =~ m/normalAfterSequencing/i && !$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'normalAfterSequencing'}) {
				$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'normalAfterSequencing'} = $G1000CommentLine[1];
			}
			elsif ($G1000CommentLine[0] =~ m/tumorAfterSequencing/i && !$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'tumorAfterSequencing'}) {
				$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'tumorAfterSequencing'} = $G1000CommentLine[1];
			}
			elsif ($G1000CommentLine[0] =~ m/Tumor_A/ && !$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Tumor_A'}) {
                                $infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Tumor_A'} = $G1000CommentLine[1];
                        }
                        elsif ($G1000CommentLine[0] =~ m/Tumor_G/ && !$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Tumor_G'}) {
                                $infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Tumor_G'} = $G1000CommentLine[1];
                        }
			elsif ($G1000CommentLine[0] =~ m/Tumor_C/ && !$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Tumor_C'}) {
                                $infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Tumor_C'} = $G1000CommentLine[1];
                        }
                        elsif ($G1000CommentLine[0] =~ m/Tumor_T/ && !$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Tumor_T'}) {
                                $infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Tumor_T'} = $G1000CommentLine[1];
                        }
                        elsif ($G1000CommentLine[0] =~ m/Normal_A/ && !$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Normal_A'}) {
                                $infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Normal_A'} = $G1000CommentLine[1];
                        }
                        elsif ($G1000CommentLine[0] =~ m/Normal_G/ && !$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Normal_G'}) {
                                $infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Normal_G'} = $G1000CommentLine[1];
                        }
                        elsif ($G1000CommentLine[0] =~ m/Normal_C/ && !$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Normal_C'}) {
                                $infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Normal_C'} = $G1000CommentLine[1];
                        }
                        elsif ($G1000CommentLine[0] =~ m/Normal_T/ && !$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Normal_T'}) {
                                $infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Normal_T'} = $G1000CommentLine[1];
                        }
                        elsif ($G1000CommentLine[0] =~ m/Normal_Ins/ && !$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Normal_Ins'}) {
                                $infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Normal_Ins'} = $G1000CommentLine[1];
                        }
                        elsif ($G1000CommentLine[0] =~ m/Tumor_Ins/ && !$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Tumor_Ins'}) {
                                $infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Tumor_Ins'} = $G1000CommentLine[1];
                        }
                        elsif ($G1000CommentLine[0] =~ m/Normal_Del/ && !$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Normal_Del'}) {
                                $infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Normal_Del'} = $G1000CommentLine[1];
                        }
                        elsif ($G1000CommentLine[0] =~ m/Tumor_Del/ && !$infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Tumor_Del'}) {
                                $infoHash{$G1000Sample}{$G1000Line[2]}{$G1000Line[3]}{$G1000Line[5]}{$G1000Line[6]}{'Tumor_Del'} = $G1000CommentLine[1];
			}
		}
	}
}

my $dbSNP_log = "";
# If an exonic_variant_function file is given
if ($snpFile) {
	
	# Open filehandler and read through file
	open(SNP, "<", $snpFile) or die "Couldn't open dbSNP file ".$snpFile."!";
	
	while (<SNP>) {
		
		if($_ =~m/^#/ || $_ eq "") { next; }
		chomp;
		# Split every line on tab and split the last column on comma
		my @snpLine = split(/\t/, $_);
		my @snpComments = split(/\s/, $snpLine[(scalar(@snpLine)-1)]);
		
		$dbSNP_log = $snpLine[0];
		
		# Variable to store the snpSample name in
		my $snpSample;
		
		# Go through all fields in the comment column
		for (my $i=1; $i<scalar(@snpComments); $i++) {
			my @snpCommentLine = split(/=/, $snpComments[$i]);
			
			# If the name of the field is snpSample save the name in $snpSample			
			if ($snpCommentLine[0] =~ m/sample/) {
				$snpSample = $snpCommentLine[1];
				# Split the given transcripts on comma
				my @transcript = split(",", $snpLine[2]);

				# Check if the comment info already is stored if not - store it
#				print $snpCommentLine[1]."\t".$snpLine[2]."\t".$snpLine[3]."\t".$snpLine[6]."\n";
				if (!$infoHash{$snpCommentLine[1]}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}) {
					$infoHash{$snpCommentLine[1]}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]} = {
																					'end' => $snpLine[4],
																					'refAllele' => $snpLine[5],
																					'dbsnp' => $snpLine[1]
																					};												
				}
				# If the comment info already exists
				else {
					# Save info about the exonic type and the transcripts in the hash
					$infoHash{$snpCommentLine[1]}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'dbsnp'} = $snpLine[1];
				}
			}
			elsif ($snpCommentLine[0] =~ m/variantAlleleRatio/i && !$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'variantAlleleRatio'}) {
				$snpCommentLine[1] =~ s/\./,/;
				$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'variantAlleleRatio'} = $snpCommentLine[1];
			}
			elsif ($snpCommentLine[0] =~ m/alleleFreq/i && !$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'alleleFreq'}) {
				my @alleleFreqs = split(/,/, $snpCommentLine[1]);
				$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'alleleFreq'} = join("\t", @alleleFreqs);;
			}
			elsif ($snpCommentLine[0] =~ m/filter/i && !$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'filter'}) {
				$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'filter'} = $snpCommentLine[1];
			}
			elsif ($snpCommentLine[0] =~ m/readDepth/i && !$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'readDepth'}) {
				$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'readDepth'} = $snpCommentLine[1];
			}
			elsif ($snpCommentLine[0] =~ m/StrandsByDesign/i && !$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'StrandsByDesign'}) {
				$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'StrandsByDesign'} = $snpCommentLine[1];
			}
			elsif ($snpCommentLine[0] =~ m/normalAfterSequencing/i && !$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'normalAfterSequencing'}) {
				$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'normalAfterSequencing'} = $snpCommentLine[1];
			}
			elsif ($snpCommentLine[0] =~ m/tumorAfterSequencing/i && !$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'tumorAfterSequencing'}) {
				$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'tumorAfterSequencing'} = $snpCommentLine[1];
			}
			elsif ($snpCommentLine[0] =~ m/Tumor_A/ && !$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Tumor_A'}) {
                                $infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Tumor_A'} = $snpCommentLine[1];
                        }
                        elsif ($snpCommentLine[0] =~ m/Tumor_G/ && !$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Tumor_G'}) {
                                $infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Tumor_G'} = $snpCommentLine[1];
                        }
			elsif ($snpCommentLine[0] =~ m/Tumor_C/ && !$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Tumor_C'}) {
                                $infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Tumor_C'} = $snpCommentLine[1];
                        }
                        elsif ($snpCommentLine[0] =~ m/Tumor_T/ && !$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Tumor_T'}) {
                                $infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Tumor_T'} = $snpCommentLine[1];
                        }
                        elsif ($snpCommentLine[0] =~ m/Normal_A/ && !$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Normal_A'}) {
                                $infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Normal_A'} = $snpCommentLine[1];
                        }
                        elsif ($snpCommentLine[0] =~ m/Normal_G/ && !$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Normal_G'}) {
                                $infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Normal_G'} = $snpCommentLine[1];
                        }
                        elsif ($snpCommentLine[0] =~ m/Normal_C/ && !$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Normal_C'}) {
                                $infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Normal_C'} = $snpCommentLine[1];
                        }
                        elsif ($snpCommentLine[0] =~ m/Normal_T/ && !$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Normal_T'}) {
                                $infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Normal_T'} = $snpCommentLine[1];
                        }
                        elsif ($snpCommentLine[0] =~ m/Normal_Ins/ && !$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Normal_Ins'}) {
                                $infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Normal_Ins'} = $snpCommentLine[1];
                        }
                        elsif ($snpCommentLine[0] =~ m/Tumor_Ins/ && !$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Tumor_Ins'}) {
                                $infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Tumor_Ins'} = $snpCommentLine[1];
                        }
                        elsif ($snpCommentLine[0] =~ m/Normal_Del/ && !$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Normal_Del'}) {
                                $infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Normal_Del'} = $snpCommentLine[1];
                        }
                        elsif ($snpCommentLine[0] =~ m/Tumor_Del/ && !$infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Tumor_Del'}) {
                                $infoHash{$snpSample}{$snpLine[2]}{$snpLine[3]}{$snpLine[5]}{$snpLine[6]}{'Tumor_Del'} = $snpCommentLine[1];
			}
		}
	}
}
close(SNP);


my $dbSNP_Flag_log = "";
# If an exonic_variant_function file is given
if ($snpFlagFile) {
	
	# Open filehandler and read through file
	open(SNP_FLAG, "<", $snpFlagFile) or die "Couldn't open dbSNP_FLAG file ".$snpFlagFile."!";
	
	while (<SNP_FLAG>) {
		
		if($_ =~m/^#/ || $_ eq "") { next; }
		chomp;
		# Split every line on tab and split the last column on comma
		my @snpFlagLine = split(/\t/, $_);
		my @snpFlagComments = split(/\s/, $snpFlagLine[(scalar(@snpFlagLine)-1)]);
		
		$dbSNP_Flag_log = $snpFlagLine[0];
		
		# Variable to store the snpFlagSample name in
		my $snpFlagSample;
		
		# Go through all fields in the comment column
		for (my $i=1; $i<scalar(@snpFlagComments); $i++) {
			my @snpFlagCommentLine = split(/=/, $snpFlagComments[$i]);
			
			# If the name of the field is snpFlagSample save the name in $snpFlagSample			
			if ($snpFlagCommentLine[0] =~ m/sample/) {
				$snpFlagSample = $snpFlagCommentLine[1];
				# Split the given transcripts on comma
				my @transcript = split(",", $snpFlagLine[2]);

				# Check if the comment info already is stored if not - store it
#				print $snpFlagCommentLine[1]."\t".$snpFlagLine[2]."\t".$snpFlagLine[3]."\t".$snpFlagLine[6]."\n";
				if (!$infoHash{$snpFlagCommentLine[1]}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}) {
					$infoHash{$snpFlagCommentLine[1]}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]} = {
																					'end' => $snpFlagLine[4],
																					'refAllele' => $snpFlagLine[5],
																					'dbsnpFlag' => "Yes"
																					};												
				}
				# If the comment info already exists
				else {
					# Save info about the exonic type and the transcripts in the hash
					$infoHash{$snpFlagCommentLine[1]}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'dbsnpFlag'} = "Yes";
				}
			}
			elsif ($snpFlagCommentLine[0] =~ m/variantAlleleRatio/i && !$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'variantAlleleRatio'}) {
				$snpFlagCommentLine[1] =~ s/\./,/;
				$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'variantAlleleRatio'} = $snpFlagCommentLine[1];
			}
			elsif ($snpFlagCommentLine[0] =~ m/alleleFreq/i && !$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'alleleFreq'}) {
				my @alleleFreqs = split(/,/, $snpFlagCommentLine[1]);
				$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'alleleFreq'} = join("\t", @alleleFreqs);;
			}
			elsif ($snpFlagCommentLine[0] =~ m/filter/i && !$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'filter'}) {
				$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'filter'} = $snpFlagCommentLine[1];
			}
			elsif ($snpFlagCommentLine[0] =~ m/readDepth/i && !$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'readDepth'}) {
				$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'readDepth'} = $snpFlagCommentLine[1];
			}
			elsif ($snpFlagCommentLine[0] =~ m/StrandsByDesign/i && !$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'StrandsByDesign'}) {
				$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'StrandsByDesign'} = $snpFlagCommentLine[1];
			}
			elsif ($snpFlagCommentLine[0] =~ m/normalAfterSequencing/i && !$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'normalAfterSequencing'}) {
				$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'normalAfterSequencing'} = $snpFlagCommentLine[1];
			}
			elsif ($snpFlagCommentLine[0] =~ m/tumorAfterSequencing/i && !$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'tumorAfterSequencing'}) {
				$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'tumorAfterSequencing'} = $snpFlagCommentLine[1];
			}
			elsif ($snpFlagCommentLine[0] =~ m/Tumor_A/ && !$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Tumor_A'}) {
                                $infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Tumor_A'} = $snpFlagCommentLine[1];
                        }
                        elsif ($snpFlagCommentLine[0] =~ m/Tumor_G/ && !$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Tumor_G'}) {
                                $infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Tumor_G'} = $snpFlagCommentLine[1];
                        }
			elsif ($snpFlagCommentLine[0] =~ m/Tumor_C/ && !$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Tumor_C'}) {
                                $infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Tumor_C'} = $snpFlagCommentLine[1];
                        }
                        elsif ($snpFlagCommentLine[0] =~ m/Tumor_T/ && !$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Tumor_T'}) {
                                $infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Tumor_T'} = $snpFlagCommentLine[1];
                        }
                        elsif ($snpFlagCommentLine[0] =~ m/Normal_A/ && !$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Normal_A'}) {
                                $infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Normal_A'} = $snpFlagCommentLine[1];
                        }
                        elsif ($snpFlagCommentLine[0] =~ m/Normal_G/ && !$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Normal_G'}) {
                                $infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Normal_G'} = $snpFlagCommentLine[1];
                        }
                        elsif ($snpFlagCommentLine[0] =~ m/Normal_C/ && !$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Normal_C'}) {
                                $infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Normal_C'} = $snpFlagCommentLine[1];
                        }
                        elsif ($snpFlagCommentLine[0] =~ m/Normal_T/ && !$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Normal_T'}) {
                                $infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Normal_T'} = $snpFlagCommentLine[1];
                        }
                        elsif ($snpFlagCommentLine[0] =~ m/Normal_Ins/ && !$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Normal_Ins'}) {
                                $infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Normal_Ins'} = $snpFlagCommentLine[1];
                        }
                        elsif ($snpFlagCommentLine[0] =~ m/Tumor_Ins/ && !$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Tumor_Ins'}) {
                                $infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Tumor_Ins'} = $snpFlagCommentLine[1];
                        }
                        elsif ($snpFlagCommentLine[0] =~ m/Normal_Del/ && !$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Normal_Del'}) {
                                $infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Normal_Del'} = $snpFlagCommentLine[1];
                        }
                        elsif ($snpFlagCommentLine[0] =~ m/Tumor_Del/ && !$infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Tumor_Del'}) {
                                $infoHash{$snpFlagSample}{$snpFlagLine[2]}{$snpFlagLine[3]}{$snpFlagLine[5]}{$snpFlagLine[6]}{'Tumor_Del'} = $snpFlagCommentLine[1];
			}
		}
	}
}
close(SNP_FLAG);


foreach my $sample (keys %infoHash) {
	my $out = "";
	if ($outputFile) {
		$out = new FileHandle(">$outputFile");
	}
	else {
		$out = new FileHandle(">$outputDir/$sample.txt");	
	}
	print {$out} "#Chromosome\tStart\tEnd\tReference_base\tVariant_base\tFilter\tGene\tType\tExonic_type\tVariant_allele_ratio\t\#_ref_allele\t\#_var_allele\tRead_depth\tdbSNP_id\tClinically_flagged\tRatio_in_1000Genome\tPossible_strands_by_design\tStrands_for_homozygous_in_normal\tStrands_for_variant_in_tumor\tTumor_A(F+|F-|S+|S-)\tTumor_G(F+|F-|S+|S-)\tTumor_C(F+|F-|S+|S-)\tTumor_T(F+|F-|S+|S-)\tNormal_A(F+|F-|S+|S-)\tNormal_G(F+|F-|S+|S-)\tNormal_C(F+|F-|S+|S-)\tNormal_T(F+|F-|S+|S-)\tTumor_Del(F+|F-|S+|S-)\tNormal_Del(F+|F-|S+|S-)\tTumor_Ins(F+|F-|S+|S-)\tNormal_Ins\tTranscripts\n";
	if ($G1000_log !~ m/^$/) { print {$out} "##1000_genome_project_file_used=".$G1000_log."\n"; }
	if ($dbSNP_log !~ m/^$/) { print {$out} "##dbSNP_file_used=".$dbSNP_log."\n"; }
	 
	foreach my $chr (keys %{$infoHash{$sample}}) {
		foreach my $start (keys %{$infoHash{$sample}{$chr}}) {
			foreach my $refAllele (keys %{$infoHash{$sample}{$chr}{$start}}) {
				foreach my $variantAllele (keys %{$infoHash{$sample}{$chr}{$start}{$refAllele}}) {
					print {$out} $chr."\t".$start."\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'end'}."\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'refAllele'}."\t".$variantAllele;
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'filter'}) {print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'filter'}; } else { print $out "\t-"; }
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'gene'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'gene'}; } else { print $out "\t-"; }
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'variantType'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'variantType'}; } else { print $out "\t-"; }
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'exonicType'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'exonicType'}; } else { print $out "\t-"; }
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'variantAlleleRatio'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'variantAlleleRatio'}; } else { print {$out} "\t-";}
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'alleleFreq'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'alleleFreq'}; } else { print {$out} "\t-";}
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'readDepth'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'readDepth'}; } else { print {$out} "\t-"; }
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'dbsnp'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'dbsnp'}; } else { print {$out} "\t-";}
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'dbsnpFlag'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'dbsnpFlag'}; } else { print {$out} "\t-";}
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'1000g'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'1000g'}; } else { print {$out} "\t-";}
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'StrandsByDesign'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'StrandsByDesign'}; } else { print {$out} "\t-";}
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'normalAfterSequencing'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'normalAfterSequencing'}; } else { print {$out} "\t-";}
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'tumorAfterSequencing'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'tumorAfterSequencing'}; } else { print {$out} "\t-";}
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Tumor_A'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Tumor_A'}; } else { print {$out} "\t-";}
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Tumor_G'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Tumor_G'}; } else { print {$out} "\t-";}
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Tumor_C'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Tumor_C'}; } else { print {$out} "\t-";}
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Tumor_T'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Tumor_T'}; } else { print {$out} "\t-";}
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Normal_A'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Normal_A'}; } else { print {$out} "\t-";}
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Normal_G'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Normal_G'}; } else { print {$out} "\t-";}
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Normal_C'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Normal_C'}; } else { print {$out} "\t-";}
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Normal_T'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Normal_T'}; } else { print {$out} "\t-";}
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Tumor_Del'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Tumor_Del'}; } else { print {$out} "\t-";}
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Normal_Del'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Normal_Del'}; } else { print {$out} "\t-";}
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Tumor_Ins'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Tumor_Ins'}; } else { print {$out} "\t-";}
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Normal_Ins'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'Normal_Ins'}; } else { print {$out} "\t-";}
					if ($infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'transcripts'}) { print {$out} "\t".$infoHash{$sample}{$chr}{$start}{$refAllele}{$variantAllele}{'transcripts'}; } else { print $out "\t-"; }
				
					print {$out} "\n";
				} 
			}
		}
	}
	close($out);
}
## Close all FileHandles
#foreach my $key (keys %indexHash) {
#	close ($indexHash{$key});
#	close ($read1Hash{$key});
#	close ($read2Hash{$key});
#}

# Print the usage help for this script.
sub usage {
  print "
******************************************************************************************************************

 This script 
 
******************************************************************************************************************
\nUsage: $0\n 
 
 -fo  Output directory [use when input file consists of more than one sample, can NOT be used together with -o]
 -o   Output file [use when input file consits of ONLY one sample, can NOT be used together with -fo]
 -v   Variant function output from annovar (file ending .variant_function)
 -e   Exonic variant function output from annovar (file ending .exonic_variant_function) [optional]
 -g   Dropped 1000 genome project output from annovar (file ending _dropped [optional]
 -s   Dropped dbSNP output from annovar (file ending _dropped [optional]
 -sf  Dropped dbSNP flagged output from annovar (file ending Flagged_dropped [optional]
 \n";

}
