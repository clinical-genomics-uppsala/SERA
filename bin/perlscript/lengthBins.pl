#!/usr/bin/perl

#This script gives a file with length bins of given size, the mean hits per base
#and the std for each of the bins
 
use strict;
use FileHandle;
use Carp;

# Subroutine prototypes
sub usage;

my ($next_arg, $infile, $bin, $column, $maxBin, $errorbars, $hitsC, $outfile, $title, $xlabel, $pdf, $ymax, $y2max, $type);

if(scalar(@ARGV) == 0){
    &usage();
}
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-i")    { $infile = shift(@ARGV); }
    elsif($next_arg eq "-b")    { $bin = shift(@ARGV); }
	elsif($next_arg eq "-c")    { $column = shift(@ARGV); } 
    elsif($next_arg eq "-e")    { $maxBin = shift(@ARGV); }
    elsif($next_arg eq "-errorbars")    { $errorbars = 1; }
    elsif($next_arg eq "-h")    { $hitsC = shift(@ARGV); }
    elsif($next_arg eq "-o")    { $outfile = shift(@ARGV); }
    elsif($next_arg eq "-t")    { $title = shift(@ARGV); }
    elsif($next_arg eq "-xlabel")    { $xlabel = shift(@ARGV); }
    elsif($next_arg eq "-pdf")    { $pdf = shift(@ARGV); }
    elsif($next_arg eq "-ymax")    { $ymax = shift(@ARGV); }
    elsif($next_arg eq "-y2max")    { $y2max = shift(@ARGV); }
    elsif($next_arg eq "-type")    { $type = shift(@ARGV); }
    else { print "Invalid argument: $next_arg"; usage(); }
}

# Error
if ((!$infile)||!($outfile)||!($column)||!($hitsC)||!($type)) { &usage(); }

#If no bin size is given it's set to default, 50.
if (!$bin) {
	$bin = 50;
}
if (!$maxBin) {
	$maxBin = 1000;
}

# Create input & output filehandle
my $in_fh = FileHandle->new($infile) or croak "Couldn't open input file ".$infile."!";
my $out_fh = FileHandle->new(">".$outfile) or croak "Couldn't open output file ".$outfile."!";

#Print header in ouput file
print $out_fh "#Bin_number Length_bin	Mean_hits	Std	No_in_bin\n";

my %hash = ();
my $hitsColumn = $hitsC - 1;
my $lengthColumn = $column -1;

#As long as there are lines to read, continue 
while (<$in_fh>) {
	
	if ($_ =~ m/^#/ || $_ =~ m/^\n/)  { next; }
	chomp;
	# Split on tab
	my @line = split(/\t/, $_);
	
	if ($hash{$line[$lengthColumn]}) {
		$hash{$line[$lengthColumn]}{'antal'}++;
		$hash{$line[$lengthColumn]}{'hits'} += $line[$hitsColumn];
		$hash{$line[$lengthColumn]}{'squaredHits'} += ($line[$hitsColumn] * $line[$hitsColumn]);
	}
	else {
		$hash{$line[$lengthColumn]}{'antal'} = 1;
		$hash{$line[$lengthColumn]}{'hits'} = $line[$hitsColumn];
		$hash{$line[$lengthColumn]}{'squaredHits'} = ($line[$hitsColumn] * $line[$hitsColumn]);
	}	
}

my ($counter, $totHits, $totSquaredHits);
my $c = 0;
for (my $i=0; $i<$maxBin; $i+=$bin) {
	$counter = 0;
	$totHits = 0;
	$totSquaredHits = 0;
	
	foreach my $length (sort keys %hash) {
		if ($length < $i+$bin && $length >= $i) {
			$counter += $hash{$length}{'antal'};
			$totHits += $hash{$length}{'hits'};
			$totSquaredHits += $hash{$length}{'squaredHits'};	
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

# If a pdf-plot is wanted run the code below
my $y;
if ($pdf) {
	
	# Set xtics
	my $xtic = "set xtics nomirror scale 0.3 rotate by -45 ('";
	my $a = 0;
	for (my $i=0; $i<$maxBin; $i+=$bin) {
		if ($i == 0) {
			$xtic .= $i." - ".($i+$bin)."' ".$a;
		}
		else{
			$xtic .= ", '".$i." - ".($i+$bin)."' ".$a;
		}
		$a++;
	}
	$xtic .= ")";
	
	# Xlabel aren't set, set it
	if (! $xlabel) {
		$xlabel = "Fragment length bins"
	}
		
	my $plot;
	if(defined($errorbars)) {
		$plot = "plot '".$outfile."' using 3 lt 2 title 'Mean hits/".$type."', '' using (\$1-0.175):3:4 lc rgb 'black' notitle with errorbars, '' using 5 axes x1y2 lt 3 title 'No. of ".$type."'";
	}
	else {
		$plot = "plot '".$outfile."' using 3 lt 2 title 'Mean hits/".$type."', '' using 5 axes x1y2 lt 3 title 'No. of ".$type."'";
	}
	
#	print $xtic."\n";
	my $gpScript = "<<EOF
		reset
		set terminal postscript  enh color
		set output '| ps2pdf - ".$pdf."
		set bar 1.000000
		set boxwidth 0.9 absolute
		set style fill solid 0.5  border -1
		set style histogram gap 1
		set style data histograms
		".$xtic."
		set title '".$title."'
		set ylabel  'Mean hits/".$type."'
		set y2label 'Number of ".$type."'
		set ytics nomirror
		set y2tics nomirror
		set xlabel '".$xlabel."'
		set yrange[0:".$ymax."]
		set y2range[0:".$y2max."]
		set xrange[-0.4:".$a."]
		".$plot."
		EOF";
		
		$gpScript =~s/\<\<EOF//;
		$gpScript =~s/EOF//;
		
		my $gnuplotFH = new FileHandle(">gnuplotLengthScript.gp");
		print $gnuplotFH $gpScript."\n";
		`gnuplot gnuplotLengthScript.gp`;
		
		system ("rm gnuplotLengthScript.gp");
		
	}	

# Print the usage help for this script.
sub usage {
  print "\nlengthBins.pl Usage: $0 -c <mean length column> -i <infile> -o <outfile> -pdf <pdf outfile> -t <plot title> -ymax <\n 
 
 -i Infile
 -o Outfile:
 	{Bin	mean hits per base	std}
 -b Bin size (Default size is 50)
 -c Column with mean length, first column is number one
 -h Column with hits, first column is number one 
 -e End point for bins (Default is 1000)
 -errorbars To plot errorbars (Optional)
 -pdf If a pdf file is wanted, give the path to the file
 -t If a title on the plot is wanted, give it here
 -type Either base or junction, depending on the input file type 
 -xlabel Set xlabel, Default is 'Fragment length bins'
 -ymax Y max (Default is autoscale)
 -y2max Y2 max (Default is autoscale)\n\n";
  exit(1);
}
