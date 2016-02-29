#!/bin/awk -f
# This script extracts the read names from a sam fil
# Input is an alignment file in SAM format.

BEGIN { File=1; }
{
	if (substr($1,0,1)!="@") {
		print $1;
	}
}
