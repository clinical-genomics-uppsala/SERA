#!/bin/awk -f
# 
# Script concatanets each selected fragment and creates
# an artificial sequence.
#
# By: Magnus Isaksson 2011
#

BEGIN { FS="\t"; }
{
  if (substr($1,0,1)!="#") {
  	print $1"#"$2"#"$3"#"$4"#"($4-$3+1)"\n"$6$6;
  } 
} 
