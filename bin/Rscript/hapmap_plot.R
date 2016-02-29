hapmap <- read.table("Sample5.only.ROI.hapmap.noNN", header=F);

pdf("Design1_roi.pdf", width=9, height=6);

	# Fix layout
	zones = matrix(c(1,2), ncol=2, byrow=TRUE);
	layout( zones, widths=c(9/10,1/10));

	# Main plot
	par(mar=c(5,5,1,1));
	plot(hapmap$V1, hapmap$V3, col=(hapmap$V11+1), pch=20, cex=0.6, ylim=c(0,1), xlim=c(10,2000), log="x", xlab="Sequencing Depth", ylab="Reference Allele Ratio");
	#grid(col="black", lty="dotted");
	
	# Calculate density for heterozygotes
	heteroz <- hapmap[1,]; hetCount=0;
	homoz <- hapmap[1,]; homCount=0;

	for (i in 1:(NROW(hapmap))) {
	  if(hapmap$V1[i] >= 10) {
	     if(hapmap$V11[i] == 1) { heteroz[(hetCount=hetCount+1),] <- hapmap[i,]; }
	     if(hapmap$V11[i] == 2) { homoz[(homCount=homCount+1),] <- hapmap[i,]; }
	  }			
	}
	
	dHet <- density(heteroz$V3);
	dHom <- density(homoz$V3);
	
	# Density plots
	par(mar=c(5,0,1,1));
	plot(dHet$y, dHet$x, pch="", ylim=c(0,1), yaxt='n', xaxt='n', ylab="", xlab="Density"); lines(dHet$y, dHet$x, col=(heteroz$V11+1));
	par(new=TRUE);
	plot(dHom$y, dHom$x, pch="", ylim=c(0,1), yaxt='n', xaxt='n', ylab="", xlab=""); lines(dHom$y, dHom$x, col=(homoz$V11+1));
	grid(col="black", lty="dotted");
	
dev.off();


