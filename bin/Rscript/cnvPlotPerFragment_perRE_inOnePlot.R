#!/usr/bin/env Rscript
#arguments: log2ratio_file	Sorted_ampregion_file	Pdf_output_file

args <- commandArgs(TRUE);
cnvT <-( read.table(args[1], header=F));
ampr <-( read.table(args[2], header=F));
enzymes <-(read.table(args[3], header=F));
pdf(args[4])
#pdf(args[4], width=4, height=6);


cnvTable <- as.data.frame(cnvT);
ampregion <- as.data.frame(ampr);
RE_array <- as.data.frame(enzymes)

# Add names to the columns
names(cnvTable) <- c("fragID", "chr", "start", "end", "strand", "cnvRatio");
names(ampregion) <- c("ampID", "chr", "start", "end");
names(RE_array) <- c("name");

noRE <- length(RE_array$name)

# Set ymax for the plot
ymax <- 5;

#Set xmax for the plot
xmax <- 0;
for (a in 1:nrow(ampregion)) {
	xmax <- xmax+ampregion[a,4]-ampregion[a,3]+1+10;
}

#par(mfrow=c(4,1), mar=(c(3,4,3,2)+0.1))
mtitle <-strsplit(as.character(args[1]), "/")[[1]];
#plot(1, type="n", axes=F, ylab="", xlab="", xlim=c(1,xmax), ylim=c(-ymax,ymax), main=mtitle[length(mtitle)]);
#par(new=TRUE)

par(fig=c(0, 1, 0.8, 1), mar=c(0,4,2,2), oma=c(0,0,0,0))
plot(1,type="n", axes=F,ylab="", xlab="");
farg <- rainbow(8)
legend_name <- (c(as.character(RE_array$name)));
legend("topleft", legend_name, text.col=farg, cex=0.4, horiz=TRUE)
title(main=mtitle[length(mtitle)])

par(fig=c(0, 1, 0, 0.9), mar=c(4,4,0,2), oma=c(0,0,0,0), new=TRUE)
plot(1, type="n", axes=T, ylab="log2(T/N)", xlab="Pseudo genome position", xlim=c(1,xmax), ylim=c(-ymax,ymax));
#legend("bottomleft", title="Restriction combinations", legend_name, text.col=farg, horiz_TRUE)
#RE_array <- c("DdeI/BfaI", "Hpy188I/NlaIII", "HpyCH4III/AluI", "HpyCH4III/Hpy188I", "HpyCH4V/HaeIII", "MseI/BfaI", "MseI/Bsp1286I", "NlaIII/DraI")
#r <-1
#for (RE in RE_array$name) {
#	print (RE)
#	RE_set <- cnvTable[grep(RE, cnvTable[,1]),]
	#print(RE_set)


	#frame[grep("string",frame[,1]),]

	#Create an empty plot with correct axes
#	mtitle <-strsplit(as.character(args[1]), "/")[[1]];
	
#	plot(1, type="n", axes=T, ylab="log2(T/N)", xlab="Pseudo genome position", xlim=c(1,xmax), ylim=c(-ymax,ymax));
#	mtext(RE, 3)

	#Set the pseudo start position for the first ampregion 
	start <- 1;

	# Set parameters for creating grey backgrounds for every second new chromosome
	lastChr <- "a";
	rectStart <- 1;
	col = 0;

	# Go through all lines in the ampregion file
	for (a in 1:nrow(ampregion)) {
#		par(new=TRUE)

		# If it is the first line in the ampregion file set lastChr to the chr of the first line
		if (a==1) {
			lastChr <- ampregion[a,2]
		}
 		# If the chr for the last line is not the same as for this line change background color
		if (lastChr != ampregion[a,2]) {
    			if (col == 0) {
      				col <- 1;      
    			}
    			# If col was set to 1 make a grey background
    			else {
				rect(start, -ymax, start+ampregion[a,4]-ampregion[a,3]+1+10, ymax, border="grey", col="grey");
      				col <- 0;
    			}
			# Set lastChr to chr of this line
			lastChr <- ampregion[a,2];
		}
		# If lastChr and chr of this line is the same and the col is set to 1 continue draw a grey background
		else {
			if (col==1) {
				rect(start, -ymax, start+ampregion[a,4]-ampregion[a,3]+1+10, ymax, border="grey", col="grey");
			}
		}
 
		r <-1
		for (RE in RE_array$name) {
			#print (RE)
			RE_set <- cnvTable[grep(RE, cnvTable[,1]),]

 
				# Extract all fragments within the ampregion
			selInAmp <- subset(RE_set, (as.character(ampregion[a,2])==as.character(RE_set[,2]) & as.numeric(ampregion[a,3])<=as.numeric(RE_set[,3]) & as.numeric(ampregion[a,4])>=as.numeric(RE_set[,4])), select=c("fragID", "chr", "start", "end", "strand", "cnvRatio"));
			# Go through all the selected fragments
			for (f in 1:nrow(selInAmp)) {
				#print (paste("A:", ampregion[a,2], ampregion[a,3],ampregion[a,4], "S:", selInAmp[f,2], selInAmp[f,3], selInAmp[f,4], sep=" "));
				par(fig=c(0, 1, 0, 0.9),  mar=c(4,4,0,2), oma=c(0,0,0,0), new=TRUE);
				# Calculate the start and end position in pseudo positions
				selStart <- start+selInAmp[f,3]-ampregion[a,3];
				selEnd <- start+selInAmp[f,4]-ampregion[a,3];
   	 
				# Draw the fragment
				lines(c(selStart,selEnd), c(selInAmp[f,6],selInAmp[f,6]), col=farg[r])
				r<-r+1
			}
		}
  		# Calculate the pseudo start position for the next ampregion
		start <- start+ampregion[a,4]-ampregion[a,3]+1+10;}
	}
	dev.off();
