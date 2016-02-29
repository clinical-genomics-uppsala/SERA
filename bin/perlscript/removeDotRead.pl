#!/usr/bin/perl
#
# This script removes reads that contains one or more a dots (.)
# and it can also remove reads from the file with base quality scores
# that contains one or more -1 (the read has a dot in colorspace).
#
use FileHandle;

# Subroutine prototypes
sub usage;

if(scalar(@ARGV) == 0){
    &usage();
}
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-i")    { $infile = shift(@ARGV); }
    elsif($next_arg eq "-o") { $outfile = shift(@ARGV); }
    else { print "Invalid argument: $next_arg"; usage(); }
}

# Error
if (!$infile) { &usage(); }

# If no output write to stdout
if (!$outfile) { $outfile="/dev/stdout"; }

# Create input & output filehandle
my $in_fh = new FileHandle("$infile") or die "Oops, couldn't open input file ".$infile."!\n";
my $out_fh = new FileHandle(">$outfile") or die "Oops, couldn't open output file ".$outfile."!\n"; 

#read the first line
$line = <$in_fh>;

#as long as $line isn't null continue 
while ((!$line eq '')) {
	
	#if there is a comment write direct to file
	if ($line =~m/^#/)
	{
		print $out_fh "$line";
		$line = <$in_fh>;
	}
	
	else
	{
		#if $line starts with >, set $id to $line 
		#and read a new line (contains the read in colorspace) 
		if ($line =~m/^>/)
		{
			$id = $line;
			$line = <$in_fh>;
			
			#if the line contains one or more dots or -1 remove
			#from the read list 
			if (!($line =~m/\.+|-1/)) {
				print $out_fh $id;
				print $out_fh $line;
				$id = "null";
				$line = <$in_fh>;
			}
			else {
				$id = "null";
				$line = <$in_fh>;
			}
		}
		else {
			$line = <$in_fh>;
		}
	}
}
			
$in_fh->close;
$out_fh->close;

# Print the usage help for this script.
sub usage {
  print "\nremoveDotRead.pl Usage: $0 -i <infile> -o <outfile>\n 
 
 -i Infile:
 	{>ID
 	seq}
 	
 -o Outfile (writes to stdout if not set):
 	>ID 
	seq\n\n";
  exit(1);
}
