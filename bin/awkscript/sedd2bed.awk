#!/usr/bin/awk -f
#
# Convert SEDD to BED, By Magnus Isaksson 2012
#
# Example:
# ./sedd2bed.awk -v name="test" -v desc="This is a test track" hg19.chrdata testregion.selection;
#

BEGIN {

	# Define fields as separated by tabs
	FS  = "\x09";
	OFS = "\x09";

	if(length(name) == 0) { name=ARGV[2]; }
	if(length(desc) == 0) { desc="Produced using input file "ARGV[2]; }

	print "track name=\""name"\" description=\""desc"\"";
}

# Make lookup array.
FNR==NR{chr[$2]=$1;next}

{
	if($0!~/#/) {
		pol="+"; if($5=="-1") { pol="-" };
		start = ($3-1);
		end   = $4;

		print chr[$2], start, end, $1, "0", pol;
	}
}
