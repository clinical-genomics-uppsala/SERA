#!/bin/awk -f
# Takes an ampregion file (merged selection) and output it
# as a download2fasta file.

#BEGIN { File=1; }
{
  if (substr($1,0,1)!="#") {
  	print $1"#"$2"#"$3"#"$4"#-1\t"$2"\t"$3"\t"$4"\t"$5
  } 
} 
