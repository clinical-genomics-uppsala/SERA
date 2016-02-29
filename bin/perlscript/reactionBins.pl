 #!/usr/bin/perl
 
 #This script gives a file with length bins of given size, the mean hits per base
 #and the std for each of the bins
 
use strict;
use FileHandle;


# Subroutine prototypes
sub usage;

my ($next_arg, $file, $hitsC, $outfile, $reactionFile, $title, $xlabel, $pdf, $ymax, $type, $y2max, $errorbars, $perc);

if(scalar(@ARGV) == 0){
    &usage();
}
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-i")    { $file = shift(@ARGV); }
    elsif($next_arg eq "-errorbars")    { $errorbars = 1; }
    elsif($next_arg eq "-h")    { $hitsC = shift(@ARGV); }
    elsif($next_arg eq "-o")    { $outfile = shift(@ARGV); }
    elsif($next_arg eq "-r")    { $reactionFile = shift(@ARGV); }
    elsif($next_arg eq "-t")    { $title = shift(@ARGV); }
    elsif($next_arg eq "-xlabel")    { $xlabel = shift(@ARGV); }
    elsif($next_arg eq "-pdf")    { $pdf = shift(@ARGV); }
    elsif($next_arg eq "-ymax")    { $ymax = shift(@ARGV); }
    elsif($next_arg eq "-y2max")    { $y2max = shift(@ARGV); }
    elsif($next_arg eq "-type")    { $type = shift(@ARGV); }
    elsif($next_arg eq "-p")    { $perc = 1; }
    else { print "Invalid argument: $next_arg"; usage(); }
}

# Error
#if ((!$file)||!($outfile)||!($hitsC)||!($type)||!($reactionFile)) { &usage(); }
if ((!$file)||!($outfile)||!($hitsC)||!($type)) { &usage(); }


# Hits column in array
my $hitsColumn = $hitsC-1;

=pod
# Put all reactions in a hash with value 0;
my $reaction = {};
my $reactionFH = new FileHandle($reactionFile) or die "Couldn't open reaction file ".$reactionFile."!\n";
my $boolean = 0;
my $counter = 0;
while (<$reactionFH>) {
	chomp;
	if ($boolean == 1) {
		$counter ++;
		my @reactionLine = split(/\t/, $_);
		# If only one enzyme is used in the reaction use only that one as reaction id...
		if (@reactionLine == 1) {
			$reaction->{$reactionLine[0]} = {	'no' => 0,
												'hits' => 0,
												'squaredHits' => 0};
		}
		# ...otherwise add a / between the two enzymes
		else {
			$reaction->{$reactionLine[0]."/".$reactionLine[1]} = {	'no' => 0,
																	'hits' => 0,
																	'squaredHits' => 0};
		}
		
	}
	if ($_ =~m/^\[Reactions\]/) {
		$boolean = 1;
	}
		
}
=cut

# Put all reactions from a selection file in a hash with value 0, ignoring reaction file
# Modified: Linus, May 2011
my $reaction={};
my $reactionFH = new FileHandle($reactionFile) or die "Couldn't open selectionfile!";
while (<$reactionFH>) {
	chomp;
	my @parse = split(/\|/,$_);
	my ($enz1,$enz2) = split(/\//,$parse[1]);
	if(!$enz2) {
		$reaction->{$enz1} = {'no' => 0,'hits' => 0,'squaredHits' => 0};
	}
	else {
		$reaction->{$enz1."/".$enz2} = {'no' => 0,'hits' => 0,'squaredHits' => 0};
	}
}
my $counter = keys %$reaction;

my $infileFH = new FileHandle($file) or die "Couldn't open file ".$file."!\n";
my $outfileFH = new FileHandle(">".$outfile) or die "Couldn't open output file ".$outfile."!\n";
print $outfileFH "#Reaction	Mean hits	Std	No in bin\n";

my $totalHits = 0;
while ( <$infileFH> ) {
	if (!($_=~m/#/)) {
		chomp;
		my @line = split(/\t/, $_);
		# Split ID to reach the reactions used to produce the fragment
		my @idParts = split (/\|/, $line[0]);
		my $enzymes = $idParts[1];
				
		$reaction->{$enzymes}->{'no'}++;
		$reaction->{$enzymes}->{'hits'} += $line[$hitsColumn];
		$reaction->{$enzymes}->{'squaredHits'} += ($line[$hitsColumn]*$line[$hitsColumn]);
		
		$totalHits += $line[$hitsColumn];
	} # Slut på if (!($_=~m/#/))
	
} # Slut på while ( <$infileFH> )

my $c = 0;
my @outputID;
my $abs = 0;
for my $id (keys %{$reaction}) {
	my $meanHit=0;
	my $std=0;
	my $percOfTot = 0;
	
	if ($perc) {
		if ( $reaction->{$id}->{'no'} == 0) {
			$meanHit = 0;
			$std = 0;
			$percOfTot = 0;
		}
		else {
			$meanHit = $reaction->{$id}->{'hits'}/$reaction->{$id}->{'no'};
			$std = sqrt(( $reaction->{$id}->{'squaredHits'} - (( $reaction->{$id}->{'hits'} *  $reaction->{$id}->{'hits'})/ $reaction->{$id}->{'no'}))/ $reaction->{$id}->{'no'});
			$percOfTot = $reaction->{$id}->{'hits'}/$totalHits*100;
		}
		print $outfileFH $c."\t".$id."\t".$meanHit."\t".$std."\t". $reaction->{$id}->{'no'}."\t".$percOfTot."\n";
		$abs += abs($percOfTot-12.5);
	}
	else {
		if ( $reaction->{$id}->{'no'} == 0) {
			$meanHit = 0;
			$std = 0;
		}
		else {
			$meanHit = $reaction->{$id}->{'hits'}/$reaction->{$id}->{'no'};
			$std = sqrt(( $reaction->{$id}->{'squaredHits'} - (( $reaction->{$id}->{'hits'} *  $reaction->{$id}->{'hits'})/ $reaction->{$id}->{'no'}))/ $reaction->{$id}->{'no'});
		}
		
		print $outfileFH $c."\t".$id."\t".$meanHit."\t".$std."\t". $reaction->{$id}->{'no'}."\n";	
	}
	
	$outputID[$c] = $id;
	$c++;
}
#print $title."\t".$abs."\n";
$infileFH->close;
$outfileFH->close;

  

my ($y);
if ($pdf) {
		
	# Set xtic
	my $xtic = "set xtics nomirror scale 0.3 rotate by -45 ('".$outputID[0]."' 0";
	for (my $i=1; $i < $counter; $i++) {
		$xtic .=", '".$outputID[$i]."' ".$i;
	}
	$xtic = $xtic.")";
	
	if (!(defined($xlabel))) {
		$xlabel = "Restriction enzymes";
	}
	
	# Set upper point of xrange
	my $xmax = $counter - 0.5;
	
	my $plot;
	if(defined($errorbars)) {
		$plot = "plot '".$outfile."' using 3 lt 2 title 'Mean hits/".$type."', '' using (\$1-0.175):3:4 lc rgb 'black' notitle with errorbars, '' using 5 axes x1y2 lt 3 title 'No. of ".$type."'";
	}
	else {
		if($perc) {
			$plot = "plot '".$outfile."' using 6 lt 2 title 'Percentage of total hits/".$type."', '' using 5 axes x1y2 lt 3 title 'No. of ".$type."'";
		}
		else {
			$plot = "plot '".$outfile."' using 3 lt 2 title 'Mean hits/".$type."', '' using 5 axes x1y2 lt 3 title 'No. of ".$type."'";
		}
		
	}
	
	# Put the gnuplot script together
	my $length = "<<EOF
		reset;

		set terminal postscript  enh color;
		set output '| ps2pdf - ".$pdf.";
		set bar 1.0;
		set boxwidth 0.9 absolute;
		set style fill solid 0.5  border -1;
		set style histogram gap 1;
		set style data histograms;
		".$xtic.";
		set title '".$title."';
		set ylabel  'Mean hits/".$type."';
		set y2label 'Number of ".$type."';
		set ytics nomirror;
		set y2tics nomirror;
		set xlabel '".$xlabel."';
		set yrange[0:".$ymax."];
		set y2range[0:".$y2max."];
		set xrange[-0.5:".$xmax."];
		".$plot.";
		EOF";
		
		$length =~s/\<\<EOF//;
		$length =~s/EOF//;
		
		my $gnuplotFH = new FileHandle(">gnuplotReactionScript.gp");
		print $gnuplotFH $length."\n";
		`gnuplot gnuplotReactionScript.gp`;
		system ("rm gnuplotReactionScript.gp");
}	




# Print the usage help for this script.
sub usage {
  print "\nreactionBins.pl Usage: $0 \n 
 
 -i Infile
 -o Outfile:
 	{Bin	mean hits per base	std}
 -h Column with hits, first column is number one 
 -errorbars To plot errorbars (Optional)
 -pdf If a pdf file is wanted, give the path to the file
 -r Selection file
 -t If a title on the plot is wanted, give it here
 -type Either base or junction, depending on the input file type 
 -ymax Y max (Default is auto scale)
 -y2max Y2 max (Default is auto scale)\n\n";
 
  exit(1);
}
