#!/bin/awk -f
# This script only saves the positions with higher or equal coverage
# to min

BEGIN { File=1; }
{
  if(substr($1,0,1)!="#") {
	if($11>=min) {
		print $0;
	}
  }
}
