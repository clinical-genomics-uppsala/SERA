#!/bin/awk -f
# 
# Script calculates the number used for -mfl flagg
# in MosaikBuild, later used for mate pair mapping.
#
# By: Magnus Isaksson 2011
#
# -v mfl=max_frag_length

BEGIN { FS="\t"; min=3000000; max=0; }
{
  if (substr($1,0,1)!="#") {

	l=$4-$3+1;
	if((l < mfl)||(mfl==-1)) {
		if(l < min) { min = l; }
		if(l > max) { max = l; }
	}
  } 
}
END {
	printf("%d\n", (((max-min)/2) + min));
}
