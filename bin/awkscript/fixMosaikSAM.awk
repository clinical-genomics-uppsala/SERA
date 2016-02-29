#!/usr/bin/awk -f
{
  if(substr($0,0,1)!="@") { 
	print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t0\t0\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14;
  } else {
	print $0;
  }
}
