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

xrange<-c(1,50000);
yrange<-c(1,50000);
plot(hits$readHits, hits$moleculeHits, xlim=xrange, ylim=yrange, pch = 20, cex = 0.5, col = "darkgreen" , log="xy", xlab = "Reads", ylab = "Molecules");
title(title);
noPos <- NROW(hits$readHits);
reg <- lm(hits$moleculeHits ~ readHits, hits);
adjustedR <- paste("Adjusted R-squared: ", round(summary(reg)$adj.r.squared, digits = 4), "\nNo of positions: ", noPos);
print(adjustedR);
legend("topleft",  box.lty = 0, text.col = "darkgreen", c(adjustedR));
dev.off();
