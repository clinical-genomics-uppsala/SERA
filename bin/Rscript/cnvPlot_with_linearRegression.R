#!/usr/bin/env Rscript
#arguments: SNPmania_tumor	SNPmania_normal	Gene_file	Min_read_depth	Avgerage_method(mean or median)	Lower_ratio_limit	Upper_ratio_limit	Pdf_plot_file	Output_file

args <- commandArgs(TRUE);
h1 <- read.table(args[1], skip=3);
h2 <- read.table(args[2], skip=3);
g <- read.table(args[3]);
minReadDepth <- (args[4]);
avg <- (args[5]);
lowerRatio <- as.numeric(args[6]);
upperRatio <- as.numeric(args[7]);
pdf(args[8]);
out <- (args[9]);
title <- (args[10]);
# Add names to the two SNPmania input files
#names(h1)<-c("tumorHits","a", "b", "chr", "pos", "c", "d", "e", "f", "g");
#names(h2) <-c("normalHits", "a", "b", "chr", "pos", "c", "d", "e", "f", "g");

names(h1) <-c("tumorHits","a", "b", "chr", "pos", "c", "d", "e", "f", "g","h","i","j","k");
names(h2) <-c("normalHits", "a", "b", "chr", "pos", "c", "d", "e", "f", "g","h","i","j","k");


# Open channel for output writing
sink(out);
cat(paste("# Min normal read depth: ", minReadDepth, "\n", sep=" "));
cat(paste("# Normalization method: ", avg, "\n", sep=" "));

# Add names to the columns in the gene file
genes <- as.data.frame(g);
names(genes) <- c("gene", "chr", "start", "end");
genes <- genes[order(genes$chr, genes$start),];

# Combine the two SNPmania input files on chromosome and position
hd <- merge(h1, h2, by=c("chr", "pos"));
# Only use the positions with a normal read depth above minReadDepth
hits <- subset(hd, hd$normalHits>=as.numeric(minReadDepth), select=c(chr, pos, tumorHits, normalHits));

# Sort all positions in the "hits" file on chromosome and position
hits <- (hits[order(hits$chr, hits$pos),]);

# Put all unique gene names in a list and create a list to save the cnv for all positions in a gene
geneList <- unique(genes$gene);
cnvList <- list();

# Create variables to save tumor and normal avg in
tumorMean <- numeric(0);
normalMean <- numeric(0);

# Calculate mean or median only including positions where the normal read depth is higher than minReadDepth
if (tolower(avg) == "mean") {
	tumorMean <- mean(hits$tumorHits);
	normalMean <- mean(hits$normalHits);
	cat(paste("# Tumor mean:", tumorMean, "\n", sep=" "));
	cat(paste("# Normal mean:", normalMean, "\n", sep=" "));
	 
}
if (tolower(avg) == "median") {
	tumorMean <- median(hits$tumorHits);
	normalMean <- median(hits$normalHits);
	cat(paste("# Tumor median:", tumorMean, "\n", sep=" "));
	cat(paste("# Normal median:", normalMean, "\n", sep=" "));
}	

# Calculate cnv only for positions with a normal read depth higher than minReadDepth and combine with its chromosome and position
cnv <- ((as.numeric(hits$tumorHits)/as.numeric(tumorMean)) - (as.numeric(hits$normalHits)/as.numeric(normalMean))) / ((as.numeric(hits$tumorHits)/as.numeric(tumorMean)) + (as.numeric(hits$normalHits)/as.numeric(normalMean)));

# Combine chromosome, position and the calculated cnv ratio
kvot <- cbind(as.character(hits$chr), as.numeric(hits$pos), as.numeric(cnv)); 

# Add names to the columns
names(kvot) <- c("chr", "pos", "cnv");

# Create a variable to save the bases within a gene
cnvRoi <- numeric(0);

# Go through all genes in the gene file and select cnv ratios only for the bases which are within the gene
for (k in 1:nrow(genes)) {
	cnvRow <- (kvot[which(kvot[,1]==as.character(genes[k,2]) & kvot[,2]>=as.numeric(genes[k,3]) & kvot[,2]<=as.numeric(genes[k,4])), ]);
	cnvRoi <- rbind(cnvRoi, cnvRow);
	cnvList[[genes[k,1]]] <- cnvRow[,3];
}


# Plot cnv ratio for each position
plot(seq(1:length(cnvRoi[,3])), cnvRoi[,3], ylim=c(-1,1), xlab="Pseudo position", ylab="((Tumor-Normal)/(Tumor+Normal))", main = title, pch = 20, cex = 0.5, col="blue");

# s is the start pos for the linear regression line
s<-1;

cat(paste("# Lower ratio: ", lowerRatio, "\n", sep=" "));
cat(paste("# Upper ratio: ", upperRatio, "\n", sep=" "));
# Print header to the output file
cat("#Gene_name	Gene_chr	Gene_start	Gene_end	CNV_ratio_start	CNV_ratio_end\n");

# Go through all genes in the gene file
for (k in 1:nrow(genes)) {
	# If there are some bases with cnv ratios in the gene continue
	if (length((cnvList[[genes[k,1]]]))) {
		# Calculating end of linear regression area
		end <- s+length((cnvList[[genes[k,1]]]))-1;
		# Creating vector with x-values
		x <- seq(as.numeric(s), as.numeric(end), by=1 );
		# Doing linear regression for the gene
		fit <- lm(as.numeric(cnvList[[genes[k,1]]]) ~ x);
		# Draw a regression line only covering the gene
		z <- predict(fit, as.data.frame(x));
		lines(x, z, col="red");
		if((z[1]<=lowerRatio && z[length(z)]<=lowerRatio) || (z[1]>=upperRatio && z[length(z)]>=upperRatio)) {
			cat(paste (genes[k,1], genes[k,2], genes[k,3], genes[k,4], z[1], z[length(z)], sep="    "));
			cat("\n");
		}

		# Calculating the starting point for the next gene
		s<- end+1;
	}	
}
sink();
dev.off();
