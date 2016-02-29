#!/usr/bin/env Rscript

args <- commandArgs(TRUE);
file <- file(args[1], "r");
pdffile <- args[2];
readDepth <-args[3];
title <- args[4];

hapmap <- read.table(file, header=F);

pdf(pdffile, width=9, height=6);
#jpeg(pdffile, width=1280, height=1024);

	# Fix layout
	zones = matrix(c(1,2), ncol=2, byrow=TRUE);
	layout( zones, widths=c(9/10,1/10));

	# Main plot
	par(mar=c(5,5,3,1));
	plot(hapmap$V1, hapmap$V3, col=(hapmap$V11+1), pch=20, cex=0.6, ylim=c(0,1), xlim=c(20,max(hapmap$V1)), log="x", xlab="Sequencing Depth", ylab="Reference Allele Ratio", main=title);
	#grid(col="black", lty="dotted");
	
#	# Calculate density for heterozygotes
#	heteroz <- hapmap[1,]; hetCount=0;
#	homoz <- hapmap[1,]; homCount=0;
#	noRef <- hapmap[1,]; noRefCount=0;
#
#	for (i in 1:(NROW(hapmap))) {
#	  if(hapmap$V1[i] >= readDepth) {
#	     if(hapmap$V11[i] == 1) { heteroz[(hetCount=hetCount+1),] <- hapmap[i,]; }
#	     if(hapmap$V11[i] == 2) { homoz[(homCount=homCount+1),] <- hapmap[i,]; }
#	     if(hapmap$V11[i] == 0) { noRef[(noRefCount=noRefCount+1),] <- hapmap[i,]; }
#	  }			
#	}
#	
#	dHet <- density(heteroz$V3);
#	dHom <- density(homoz$V3);
#	dnoRef <- density(noRef$V3);
	
	# Density plots
#	par(mar=c(5,0,3,1));
#	plot(dHet$y, dHet$x, pch="", ylim=c(0,1), yaxt='n', xaxt='n', ylab="", xlab="Density"); lines(dHet$y, dHet$x, col=(heteroz$V11+1));
#	par(new=TRUE);
#	plot(dnoRef$y, dnoRef$x, pch="", ylim=c(0,1), yaxt='n', xaxt='n', ylab="", xlab=""); lines(dnoRef$y, dnoRef$x, col=(noRef$V11+1));
#	par(new=TRUE);
#	plot(dHom$y, dHom$x, pch="", ylim=c(0,1), yaxt='n', xaxt='n', ylab="", xlab=""); lines(dHom$y, dHom$x, col=(homoz$V11+1));
#	grid(col="black", lty="dotted");
#	legend("top", box.lty=0, text.col="black", cex=0.6, c(paste("AA: ",homCount), paste("Ab: ",hetCount))); 
		
dev.off();


