#!/usr/bin/perl

#This script gives a file with GC bins of the given size with
#corresponding no. mean no. of hits per base within the bin, 
# the std and the number of bases within the bin

use FileHandle;
use strict;
use Carp;

# Subroutine prototypes
sub usage;


my ($next_arg, $infile, $bin, $column, $errorbars, $hitsC, $outfile, $pdf, $title, $xlabel, $ymax, $y2max, $type);
if(scalar(@ARGV) == 0){
    &usage();
}
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-i")    { $infile = shift(@ARGV); }
    elsif($next_arg eq "-b")    { $bin = shift(@ARGV); }
    elsif($next_arg eq "-c")    { $column = shift(@ARGV); }
    elsif($next_arg eq "-errorbars")    { $errorbars = 1; }
    elsif($next_arg eq "-h")    { $hitsC = shift(@ARGV); }
    elsif($next_arg eq "-o")    { $outfile = shift(@ARGV); }
    elsif($next_arg eq "-pdf")    { $pdf = shift(@ARGV); }
    elsif($next_arg eq "-t")    { $title = shift(@ARGV); }
    elsif($next_arg eq "-xlabel")    { $xlabel = shift(@ARGV); }
    elsif($next_arg eq "-ymax")    { $ymax = shift(@ARGV); }
    elsif($next_arg eq "-y2max")    { $y2max = shift(@ARGV); }
    elsif($next_arg eq "-type")    { $type = shift(@ARGV); }
    else { print "Invalid argument: $next_arg"; usage(); }
}

# Error
if (!$infile || !$outfile || !$column || !$hitsC) { &usage(); }

#If no bin size is given it's set to default, 10.
if (!$bin) {
	$bin = 10;
}

# Create input & output filehandle
my $in_fh = FileHandle->new($infile) or croak "Couldn't open input file ".$infile."!";
my $out_fh = FileHandle->new(">".$outfile) or croak "Couldn't open output file ".$outfile."!";

print $out_fh "#Bin_number	GC_bins	Mean_hits	Std	No_of_junction\n";

# Add important info in hash
my %hash = ();
my $hitsColumn = $hitsC - 1;
my $gcColumn = $column -1;
#As long as there are lines to read, continue 
while (<$in_fh>) {
	
	if ($_ =~ m/^#/ || $_ =~ m/^\n/)  { next; }
	chomp;
	# Split on tab
	my @line = split(/\t/, $_);
	# If the fragment gc are found before increase antal with one and add the hits	
	if ($hash{$line[$gcColumn]}) {
		$hash{$line[$gcColumn]}{'antal'}++;
		$hash{$line[$gcColumn]}{'hits'} += $line[$hitsColumn];
		$hash{$line[$gcColumn]}{'squaredHits'} += ($line[$hitsColumn] * $line[$hitsColumn]);
	}
	# Else add the gc to the hash
	else {
		$hash{$line[$gcColumn]}{'antal'} = 1;
		$hash{$line[$gcColumn]}{'hits'} = $line[$hitsColumn];
		$hash{$line[$gcColumn]}{'squaredHits'} = ($line[$hitsColumn] * $line[$hitsColumn]);
	}	
}

# Print the bins to file
my ($counter, $totHits, $totSquaredHits);
my $c = 0;
for (my $i=0; $i<100; $i+=$bin) {
	$counter = 0;
	$totHits = 0;
	$totSquaredHits = 0;
	
	foreach my $gc (sort keys %hash) {
		if ($gc < $i+$bin && $gc >= $i) {
			$counter += $hash{$gc}{'antal'};
			$totHits += $hash{$gc}{'hits'};
			$totSquaredHits += $hash{$gc}{'squaredHits'};	
		}
	}
	# If there aren't any fragments/junctions in the bin print zeros
	if ($counter == 0) {
		print {$out_fh} $c."\t".$i."-".($i+$bin)."\t0\t0\t0\n";
	}
	# Else print mean hits, std and number of fragments per bin
	else {
		my $std = sqrt(($totSquaredHits - (($totHits * $totHits)/$counter))/$counter);
		print {$out_fh} $c."\t".$i."-".($i+$bin)."\t".($totHits/$counter)."\t".$std."\t".$counter."\n";
	}
	# Counting number of bins
	$c++;
}

# If a pdf plot is wanted continue run the code below
if ($pdf) {
	
	# Set xtics
	my $xtic = "set xtics nomirror scale 0.3 rotate by -45 ('";
	my $a = 0;
	for (my $i=0; $i<100; $i+=$bin) {
		# If it's the first bin don't add a comma before
		if ($i == 0) {
			$xtic .= $i." - ".($i+$bin)."' ".$a;
		}
		else{
			# If bin end is larger than 100 set xtic to 100
			if (($i+$bin) > 100) {
				$xtic .= ", '".$i." - 100' ".$a;
			}
			else {
				$xtic .= ", '".$i." - ".($i+$bin)."' ".$a;
			}
		}
		$a++;
	}
	$xtic .= ")";
	
	my $plot;
	if(defined($errorbars)) {
		$plot = "plot '".$outfile."' using 3 lt 2 title 'Mean hits/".$type."', '' using (\$1-0.175):3:4 lc rgb 'black' notitle with errorbars, '' using 5 axes x1y2 lt 3 title 'No. of ".$type."'";
	}
	else {
		$plot = "plot '".$outfile."' using 3 lt 2 title 'Mean hits/".$type."', '' using 5 axes x1y2 lt 3 title 'No. of ".$type."'";
	}
	
	my $gc = "<<EOF
		reset
		set terminal postscript enh color
		set output '| ps2pdf - ".$pdf."

		set bar 1.000000
		set boxwidth 0.9 absolute
		set style fill solid 0.5  border -1 
		set style histogram gap 1
		set style data histograms
		".$xtic."
		set title '".$title."'
		set ylabel  'Mean hits/".$type."'
		set ytics nomirror
		set y2tics nomirror
		set y2label 'Number of ".$type."'
		set xlabel '".$xlabel."'
		set yrange[0:".$ymax."]
		set y2range[0:".$y2max."]
		set xrange[-0.4:".$a."]
		".$plot."
		EOF";
		
		$gc =~s/\<\<EOF//;
		$gc =~s/EOF//;
		
		my $gnuplotFH = new FileHandle(">gnuplotGcScript.gp");
		print $gnuplotFH $gc."\n";
		`gnuplot gnuplotGcScript.gp`;
		system ("rm gnuplotGcScript.gp");
}
		
# Print the usage help for this script.
sub usage {
  print "\ngcBins.pl Usage: $0 -b <bin size> -c <mean GC column> -h <hits column> -i <infile> -o <outfile> -pdf <pdf file> -t <title> -xlabel <x label> -ymax\n 
 
 -b Bin size (Default size is 10)
 -c Column with mean GC, first column is number one
 -errorbars To plot errorbars (Optional)
 -h Column with hits, first column is number one
 -i Infile
 -o Outfile:
 	{Bin	mean hits per base	std}
 -pdf If a pdf file is wanted, give the path to the file
 -t If a title on the plot is wanted, give it here
 -type Either base or junction, depending on the input file type
 -xlabel The xlabel wanted
 -ymax Y max (Default is auto scale)
 -y2max Y2 max (Default is auto scale)\n\n";
  exit(1);
}
