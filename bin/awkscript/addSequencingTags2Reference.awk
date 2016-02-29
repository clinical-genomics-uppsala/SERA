#!/bin/awk -f
# 
# Script adds sequencing tags to each selected fragments.
# By: Magnus Isaksson 2011
#
# -v fivePtag="seq"
# -v threePtag="seq"
# -v r=readLength

BEGIN { FS="\t"; oldNc="Non"; oldStart=0; oldEnd=0;}
{
  if (substr($1,1,1)!="#") {
	
	# Check if fragment is shorter than read length.
	if (($4-$3+1) < r) {
		
		# Did we add this fragment from an other reaction already?
		if(($2 != oldNc) || ($3 != oldStart) || ($4 != oldEnd)) {
			print ">"$1"#"$2"#"$3"#"$4"#"length(fivePtag)"\n"fivePtag $5 threePtag;
		}
	}

	oldNc = $2;
	oldStart = $3;
	oldEnd = $4;
  } 
} 
