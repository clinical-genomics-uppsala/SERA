#!/bin/awk -f
# It takes a .map file with hits as third column to calculate mean number of hits
# Input file: {Chrom    Pos     Hits}

BEGIN { File=1; c=0; sum=0;}
{
        c++;
	sum+=$k;
} 
END { print "Mean:"sum/c; }
