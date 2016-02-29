#!/bin/awk -f
# 
# Script concatanets each selected fragment and creates
# an artificial sequence.
#
# By: Magnus Isaksson 2011
#

BEGIN { FS="\t"; oldNc="Non"; oldStart=0; oldEnd=0;}
{
  if (substr($1,0,1)!="#") {

	# Did we add this fragment from an other reaction already?
	if(($2 != oldNc) || ($3 != oldStart) || ($4 != oldEnd)) {
	  	print ">"$1"#"$2"#"$3"#"$4"#"($4-$3+1)"\n"$5$5;
	}
	
	oldNc = $2;
	oldStart = $3;
	oldEnd = $4;
  } 
} 
