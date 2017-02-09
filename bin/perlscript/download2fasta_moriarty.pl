#!/usr/bin/perl
#
# This script takes a indexfile with format:
#
# Id {tab} Chr {tab} Start {tab} Stop
# Seq1  NC_000017.9 10000   10010
#
# And download each sequence to a fasta/table file by using fastacmd.
# By: Magnus Isaksson 2008
#

use strict;
use warnings;

use threads;
use threads::shared;

use Data::Dumper;

my($next_arg, $index_file, $outputfile, $blastdb, $tableformat,$numberOfThreads);

# Subroutine prototypes
sub usage;
sub save_to_fasta_file;
sub save_to_table_file;
sub indexfile2fastacmdOpt;
sub extractFasta;

if(scalar(@ARGV) == 0){
    &usage();
}
# Parse the command line
while(scalar @ARGV > 0){
    $next_arg = shift(@ARGV);
    if($next_arg eq "-s")    { $index_file = shift(@ARGV); }
    elsif($next_arg eq "-o") { $outputfile = shift(@ARGV); }
    elsif($next_arg eq "-d") { $blastdb = shift(@ARGV); }
    elsif($next_arg eq "-t") { $tableformat = shift(@ARGV); }
    elsif($next_arg eq "-th") { $numberOfThreads = shift(@ARGV); }
    else { print "Invalid argument: $next_arg"; usage(); }
}

# Error
if ((!$blastdb)||(!$index_file)) { &usage(); }
# Defaults
if (!$outputfile) { $outputfile="/dev/stdout"; }
if(!$numberOfThreads) { $numberOfThreads = 5; }
# Make fastacmd options for all posts in the index file.
our %fastaCmd:shared;
our %fastaOpt:shared;
our %fastaSeq:shared;

# Open output file.
open(OUTPUT, "> $outputfile" ) or die "Can't save to output $outputfile : $!";
# Open input file.
open( FILE, "< $index_file" ) or die "Can't open $index_file : $!";
while (<FILE>) { 
        
    chomp; # Remove all \n in line.
    if (!(/^#/)) {
        
        my @columns = split(/\t/, $_);
        my $id = $columns[0];
        my $chr = $columns[1];
        my $start = $columns[2];
        my $stop = $columns[3];

        #fastacmd -d human_genomic -s NC_000017.9 -L 7520186,7520668
        $fastaOpt{$id} = $id."\t".$chr."\t".$start."\t".$stop;
        $fastaCmd{$id} = "-d ".$blastdb." -s ".$chr." -L ".$start.",".$stop;
        #my $cmdLine = "-d ".$blastdb." -s ".$chr." -L ".$start.",".$stop;
        #my %content: shared; = {'table' => $table, 'cmdline' => $cmdLine};
        
    }
}
close(FILE);

# Create an array contaning all keys
my @keys;
foreach my $key ( keys %fastaCmd ) {
    push(@keys,$key);
}

# counter used to keep track of the next available key
our $indexCounter:shared = 0;

# start threads
for(my $i = 0; $i < $numberOfThreads; $i++){
    threads->new(\&extractFasta, \@keys);
}

foreach my $thr (threads->list) { 
    # Don't join the main thread or ourselves 
    if ($thr->tid && !threads::equal($thr, threads->self)) { 
        $thr->join; 
    } 
}

# Print all fasta data to a file.
for( my $i = 0; $i <= $#keys; $i++) {  
   if (defined($tableformat)) {
        save_to_table_file($keys[$i], $fastaOpt{$keys[$i]}, $fastaSeq{$keys[$i]});
   } else {
        save_to_fasta_file($keys[$i], $fastaSeq{$keys[$i]});
   }
}
close(OUTPUT);

# Function used to iterate through all input values.
sub extractFasta{
    my @keys = @{(shift)};
    my $index = 0;
    
    {
         
        lock($indexCounter);
        $index = $indexCounter;
        $indexCounter++;
    }

    while($index <= $#keys){
        my $key = $keys[$index];
        my $ans = `blastdbcmd $fastaCmd{$key}`;
        #my $ans = `fastacmd $fastaCmd{$key}`;
        $ans =~ s/^>.*?\n//;  # Remove fasta header from ans.
        $ans =~ s/[\n|\s]//g; # Remove all new lines and speces from ans.
        $fastaSeq{$key} = $ans;
        {
            lock($indexCounter);
            $index = $indexCounter;
            $indexCounter++;
        }
    }
    return 1;
}

# Saves to fasta file append
sub save_to_fasta_file {

    my($id, $seq) = @_;
    my $seq_length = length($seq);
    my $count = 0;
    
    # Print header
    print OUTPUT ">".$id."\n";
    # Print sequence
    for(my $i = 0; $i < ($seq_length); $i++){
      print OUTPUT substr($seq, $i, 1);
      $count++; 
      if($count == 70){ print OUTPUT "\n"; $count=0; }
    }
    print OUTPUT "\n";
}

# Saves to table file append
sub save_to_table_file {

    my($smallid, $fullid, $seq) = @_;

    # Print to file
    if($tableformat eq "full") {
        print OUTPUT $fullid."\t".$seq."\n";
    } elsif($tableformat eq "small") {
        print OUTPUT $smallid."\t".$seq."\n";
    } else {
        print "ERROR: Table format should be small of full.\n";
        usage(); exit(1);
    }
}

# Print the usage help for this script.
sub usage {
  print "\nUsage: $0 -s <index_file> -o <fasta_outputfile> -d <blastdb> -t <full/small>\n 
 You need to have <fastacmd> from the BLAST package in your path for this
 script to work.
 
 -s Index file in format:
    
    # Id {tab} Chr {tab} Start {tab} Stop
    Seq1    NC_000017   10001   10010
    Seq2    NC_000017   10011   10020

 -o Outputfile saved in fasta or table format.

 -d Blast database
  
 -th number of threads that should be created.

 -t Table format output full or small.
     full  = id {tab} nc-number {tab} start [tab} end {tab} sequence
     small = id {tab} seq\n\n";
    
  exit(1);
}