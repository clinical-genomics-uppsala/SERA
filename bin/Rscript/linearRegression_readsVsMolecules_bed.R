#!/usr/bin/env Rscript
#arguments: SNPmania_tumor	SNPmania_normal	Gene_file	Min_read_depth	Avgerage_method(mean or median)	Lower_ratio_limit	Upper_ratio_limit	Pdf_plot_file	Output_file

args <- commandArgs(TRUE);
h1 <- read.table(args[1]);
minReadDepth <- (args[2]);
pdf(args[3]);
title <- (args[4]);

# Add names to the two SNPmania input files
names(h1) <-c("chr", "start", "end", "id", "a", "strand","reads","molecules");

hits <- subset(h1, h1$readHits>=as.numeric(minReadDepth) & h1$moleculeHits>=as.numeric(minReadDepth), select=c(chr, start, end, readHits, moleculeHits));

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
