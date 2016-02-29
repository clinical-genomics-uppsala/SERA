#!/usr/bin/env Rscript
#arguments: log2ratio_file	Sorted_ampregion_file	Pdf_output_file

args <- commandArgs(TRUE);
cnvT <-( read.table(args[1], header=F));
ampr <-( read.table(args[2], header=F));
pdf(args[3]);

cnvTable <- as.data.frame(cnvT);
ampregion <- as.data.frame(ampr);

# Add names to the columns
names(cnvTable) <- c("fragID", "chr", "start", "end", "strand", "cnvRatio");
names(ampregion) <- c("ampID", "chr", "start", "end");

# Set ymax for the plot
ymax <- 5;
#if (max(cnvTable$cnvRatio) >= ymax ) {
#  ymax <-max(cnvTable$cnvRatio) + 1;
#}

#Set xmax for the plot
xmax <- 0;
for (a in 1:nrow(ampregion)) {
  xmax <- xmax+ampregion[a,4]-ampregion[a,3]+1+10;
}

#Create an empty plot with correct axes
mtitle <-strsplit(as.character(args[1]), "/")[[1]];

plot(1, type="n", axes=T, ylab="log2(T/N)", xlab="Pseudo genome position", main=mtitle[length(mtitle)], xlim=c(1,xmax), ylim=c(-ymax,ymax));

#Set the pseudo start position for the first ampregion 
start <- 1;

# Set parameters for creating grey backgrounds for every second new chromosome
lastChr <- "a";
rectStart <- 1;
col = 0;

# Go through all lines in the ampregion file
for (a in 1:nrow(ampregion)) {
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
      rect(start, -11, start+ampregion[a,4]-ampregion[a,3]+1+10, ymax, border="grey", col="grey");
      col <- 0;
    }
    # Set lastChr to chr of this line
    lastChr <- ampregion[a,2];
  }
  # If lastChr and chr of this line is the same and the col is set to 1 continue draw a grey background
  else {
    if (col==1) {
      rect(start, -11, start+ampregion[a,4]-ampregion[a,3]+1+10, ymax, border="grey", col="grey");
    }
  }
  
  # Extract all fragments within the ampregion
  selInAmp <- subset(cnvTable, (ampregion[a,2]==cnvTable[,2] & as.numeric(ampregion[a,3])<=as.numeric(cnvTable[,3]) & as.numeric(ampregion[a,4])>=as.numeric(cnvTable[,4])), select=c("fragID", "chr", "start", "end", "strand", "cnvRatio"));
  # Go through all the selected fragments
  for (f in 1:nrow(selInAmp)) {
    #     print (paste("A:", ampregion[a,2], ampregion[a,3],ampregion[a,4], "S:", selInAmp[f,2], selInAmp[f,3], selInAmp[f,4], sep=" "));
    par(new=TRUE);
    # Calculate the start and end position in pseudo positions
    selStart <- start+selInAmp[f,3]-ampregion[a,3];
    selEnd <- start+selInAmp[f,4]-ampregion[a,3];
    
    # Draw the fragment
    lines(c(selStart,selEnd), c(selInAmp[f,6],selInAmp[f,6]));
  }
  # Calculate the pseudo start position for the next ampregion
  start <- start+ampregion[a,4]-ampregion[a,3]+1+10;
}
dev.off();
