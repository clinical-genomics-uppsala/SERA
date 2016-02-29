#!/usr/bin/awk -f
BEGIN {FS="\t"}
{
	refAllel=0;
	for(i=1; i<=length($5); i++) {
		if((substr($5,i,1)==".")||(substr($5,i,1)==",")) {
			refAllel++;
		}
	}

	print $0 "\t" (refAllel/$4);
}
