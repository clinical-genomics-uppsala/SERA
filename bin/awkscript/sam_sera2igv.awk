#!/usr/bin/awk -f
BEGIN { FS="\t"; }
{
	# Header of not?
	if(substr($0,1,1)=="@") {
	
		# Reference post?
		if(substr($0,1,3)=="@SQ") {
		
			split($2,s,"#");
			if (s[2]=="NC_012920.1") {
				chr="M"; # Special NC-number.
			} else {
			
				chr=substr(s[2],8,2);
				
				if(substr(chr,1,1)=="0") { chr = substr(chr,2,1); }
				if(chr=="23") { chr="X"; }
				if(chr=="24") { chr="Y"; }
			}

			# Print new @SQ row.
			print $1 "\t" "SN:chr" chr "\t" $3 "\t" $4;

		} else {
			print $0;
		}

	} else {
	
		# Fix a read row.
		split($3,s,"#");
		if (s[2]=="NC_012920.1") {
			chr="M"; # Special NC-number.
		} else {
			
			chr=substr(s[2],8,2);

			if(substr(chr,1,1)=="0") { chr = substr(chr,2,1); }
			if(chr=="23") { chr="X"; }
			if(chr=="24") { chr="Y"; }
		}

		if($14=="") { $14="PG:Z:Mosaik"; }

		# Print new read row.
		print $1 "\t" $2 "\t" "chr" chr  "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14;
		#printf("%s\t%s\tchr%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, chr, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14);

	}
}
