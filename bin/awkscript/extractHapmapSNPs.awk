#!/usr/bin/awk -f

# This script takes a goSNPmania2.pl output file and keep the bases which are in Hapmap
# and aren't annotated as NN

{
	if(substr($0,0,1)!="#") {
  		if($15>0 && $23!="NN") {
			print $0;
		}
	}	
}
