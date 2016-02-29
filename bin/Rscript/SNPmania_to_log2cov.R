#!/usr/bin/env Rscript
#arguments: log2ratio_file	Sorted_ampregion_file	Pdf_output_file

args <- commandArgs(TRUE);
SNPmania <- read.table(args[1], header=F)
pdf(args[2]);
tit<-(args[3])
end<-(args[4])
step<-(args[5])

SNPframe <- as.data.frame(SNPmania)

names(SNPframe) <- c("rd", "a", "b", "chr", "pos", "e", "f", "g", "h", "i", "j", "k", "l", "m");

log2val <- subset(SNPframe, SNPframe$rd>0, select=c("chr","pos","rd"))
log2val$log2 <- log(log2val$rd)/log(2)

nonlog <- subset(SNPframe, SNPframe$rd==0, select=c("chr","pos","rd"))
log2all <- 0

if (nrow(nonlog)>0) {
	nonlog$log2 <- rep(0,nrow(nonlog))
	log2all <- rbind(log2val, nonlog)
} else {
	log2all <- logval
}

bins <- seq(0,as.numeric(end),by=as.numeric(step))

hist(log2all$log2, breaks=bins, xlab="log2 of read depth", ylab="Percentage of positions", ylim=c(0,0.4), freq=FALSE, main=tit)
dev.off();
