#!/usr/bin/env Rscript
#arguments: SNPmania_tumor	SNPmania_normal	Gene_file	Min_read_depth	Avgerage_method(mean or median)	Lower_ratio_limit	Upper_ratio_limit	Pdf_plot_file	Output_file

args <- commandArgs(TRUE);
h1 <- read.table(args[1], skip=3);
h2 <- read.table(args[2], skip=3);
minReadDepth <- (args[3]);
pdf(args[4]);
title <- (args[5]);

# Add names to the two SNPmania input files
names(h1) <-c("readHits","a", "b", "chr", "pos", "c", "d", "e", "f", "g","h","i","j","k");
names(h2) <-c("moleculeHits", "a", "b", "chr", "pos", "c", "d", "e", "f", "g","h","i","j","k");

# Combine the two SNPmania input files on chromosome and position
hd <- merge(h1, h2, by=c("chr", "pos"));
# Only use the positions with a normal read depth above minReadDepth
hits <- subset(hd, hd$readHits>=as.numeric(minReadDepth) & hd$moleculeHits>=as.numeric(minReadDepth), select=c(chr, pos, readHits, moleculeHits));
hits$rdPerMol <- (hits$readHits/hits$moleculeHits)

xrange<-c(1,(max(hits$readHits)+100));
yrange<-c(1,(max(hits$rdPerMol)+3));
plot(hits$moleculeHits, hits$rdPerMol, xlim=xrange, ylim=yrange, pch = 20, cex = 0.5, col = "darkgreen" , log="x", ylab = "Reads/molecule", xlab = "Molecule depth");
title(title);
noPos <- NROW(hits$readHits);
dev.off();
