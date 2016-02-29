BEGIN {
	FILE=1;
	bp=0;
	min30=0;
	min100=0;
	min300=0;
	min500=0;
	print "#Sample\tF30x\tF100x\tF300x\tF500x\tminRD30\tminRD100\tminRD300\tminRD500\ttot#bp"
}
{
if($1!~/^#/) {
	bp++;
	if($1>29) {
		min30++;
		if ($1>99) {
			min100++;
			if($1>299) {
				min300++;
				if($1>499) {
					min500++;
					}
				}
			}
		}
	}
} 
END {
	print FILENAME"\t"min30/bp"\t"min100/bp"\t"min300/bp"\t"min500/bp"\t"min30"\t"min100"\t"min300"\t"min500"\t"bp
} 
