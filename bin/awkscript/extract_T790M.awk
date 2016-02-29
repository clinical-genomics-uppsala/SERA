#!/bin/awk -f

# -v fpath="file path to the sample running"
# -v tumor="the tumour type of the sample"
# -v cds="cds change"
# -v aa="amino acid change"
# -v r="reference base"
# -v v="variant base"

BEGIN{
	FS="\t";
#	print "Run\tSample\tTumour\tVaf\tRef_RD\tVar_RD\tTot_RD\t#Ref_amp\t#Var_amp\tChr\tPos\tRef\tVar\tCDS_change\tAA_change\tRef_amp\tVar_amp";
	base["A"] = 1;
	base["G"] = 2;
	base["C"] = 3;
	base["T"] = 4;
}

{
	rbase = base[r];
	vbase = base[v];
	split($7,a,"|");
	split(fpath,b,"/");
	split(b[length(b)],c,".");
	split($15,d,"|");
	refAmp=0;
	varAmp=0; 
	split(d[rbase],ref,"#");
	split(d[vbase],var,"#");

	for (i in ref) {
		split(ref[i],s,":");
		if (s[4]>=5) {
			refAmp++;
		}
	}

	for (k in var) {
		split(var[k],x,":");
		if(x[4]>=5) {
			varAmp++;
		}
	}

	print b[length(b)-2]"\t"c[1]"\t"tumor"\t"a[vbase]/(a[vbase]+a[rbase])"\t"a[rbase]"\t"a[vbase]"\t"a[rbase]+a[vbase]"\t"refAmp"\t"varAmp"\t"$4"\t"$5"\t"r"\t"v"\t"cds"\t"aa"\t"d[rbase]"\t"d[vbase]
}

