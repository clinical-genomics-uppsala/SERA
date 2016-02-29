#!/bin/awk -f
#
# -v n="Name"
# -v d="Description"
#
BEGIN { 
	FS="\t"; 
	print "## Hits per base BED-file format.";
	print "## Olink Genomics, http://www.olinkgenomics.com";
	print "browser full altGraph";
	print "track type=bedGraph name=\"" n "\" description=\"" d "\" visibility=full color=205,95,152 altColor=95,152,205 priority=20 windowingFunction=maximum smoothingWindow=2 alwaysZero=on autoScale=on";
}
{
	split($1, s, ".");
	
	if(s[i] ~ /NC_00000/) {
		chr = substr(s[1],9,1);
	} else {
		chr = substr(s[1],8,2);
	}

	if (chr == 23) { chr = "X"; }
	if (chr == 24) { chr = "Y"; }

	print "chr" chr "\t" $2 "\t" ($2+1) "\t" $3;
}
