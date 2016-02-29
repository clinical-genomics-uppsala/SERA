#!/usr/bin/env Rscript
#arguments: SNPmania_tumor	SNPmania_normal	Gene_file	Min_read_depth	Avgerage_method(mean or median)	Lower_ratio_limit	Upper_ratio_limit	Pdf_plot_file	Output_file

args <- commandArgs(TRUE);
h1 <- read.table(args[1]);
minReadDepth <- (args[2]);
pdf(args[3]);
title <- (args[4]);

# Add names to the two SNPmania input files
names(h1) <-c("chr","start", "end", "amplicon", "a", "strand", "read", "molecule");

# Only use the positions with a normal read depth above minReadDepth
hits <- subset(h1, h1$read>=as.numeric(minReadDepth) & h1$molecule>=as.numeric(minReadDepth), select=c(chr, start, end, amplicon, strand, read, molecule));
hits$rdPerMol <- (hits$read/hits$molecule)

xrange<-c(1,(max(hits$read)+100));
#yrange<-c(1,(max(hits$rdPerMol)+1));
yrange<-c(1,3);
plot(hits$read, hits$rdPerMol, xlim=xrange, ylim=yrange, pch = 20, cex = 0.5, col = "darkgreen" , log="x", ylab = "Reads/molecule", xlab = "Read depth");
title(title);
noPos <- NROW(hits$readHits);
dev.off();
