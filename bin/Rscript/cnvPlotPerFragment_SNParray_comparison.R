#!/usr/bin/env Rscript
#arguments: log2ratio_file	Sorted_ampregion_file	Pdf_output_file

args <- commandArgs(TRUE);
cnvT <-( read.table(args[1], header=F));
ampr <-( read.table(args[2], header=F));
pdf(args[3]);
snpArrayTable <- (read.table(args[4], header=F));
selection <- ( read.table(args[5], header=F));


ampregion <- as.data.frame(ampr);
cnvTable <- as.data.frame(cnvT);
snpArray <- as.data.frame(snpArrayTable);
selTable <- as.data.frame(selection);

# Add names to the columns

names(cnvTable) <- c("fragID", "chr", "start", "end", "strand", "cnvRatio");
names(ampregion) <- c("ampID", "chr", "start", "end");
names(snpArray) <- c("fragID", "chr", "start", "end", "cnvRatio", "Cn","mCn");
names(selTable) <- c("fragID", "chr", "start", "end", "strand", "cnvRatio");

# Set ymax for the plot
ymax <- 10;
#if (max(cnvTable$cnvRatio) >= ymax ) {
#  ymax <-max(cnvTable$cnvRatio) + 1;
#}

#Set xmax for the plot
xmax <- 0;
for (a in 1:nrow(ampregion)) {
  xmax <- xmax+ampregion[a,4]-ampregion[a,3]+1+10;
}

ytics <- c(-9, -5, 0, 5, 9);

#Create an empty plot witih correct axes
mtitle <-strsplit(as.character(args[1]), "/")[[1]];
layout(matrix(c(1,2),2,1,byrow=TRUE), height=c(0.6,0.4), width=1);
par(mar=c(0,4,2,2)); 
plot(1, type="n", ylab="log2(T/N)", xlab="", axes=TRUE, main=mtitle[length(mtitle)], xlim=c(1,xmax), ylim=c(-ymax, ymax), xaxt="n", yaxt="n");
axis(1, tck=0, labels=FALSE);
axis(2, at=ytics, labels=NULL);
#plot(1, type="n", axes=T, ylab="log2(T/N)", xlab="Pseudo genome position", main=mtitle[length(mtitle)], xlim=c(1,xmax), ylim=c(-ymax, ymax));

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
      par(mar=c(0,4,2,2), new=TRUE);
      rect(start, -ymax, start+ampregion[a,4]-ampregion[a,3]+1+10, ymax, border="grey", col="grey");
      col <- 0;
    }
    # Set lastChr to chr of this line
    lastChr <- ampregion[a,2];
  }
  # If lastChr and chr of this line is the same and the col is set to 1 continue draw a grey background
  else {
    if (col==1) {
      par(mar=c(0,4,2,2), new=TRUE);
      rect(start, -ymax, start+ampregion[a,4]-ampregion[a,3]+1+10, ymax, border="grey", col="grey");
    }
  }
 par(mar=c(0,4,2,2), new=TRUE);
 lines(c(start,start+ampregion[a,4]), c(0,0), lty=2, col="black");
# plot(c(start,start+ampregion[a,4]), c(0,0), type="l",lty=2, col="blue")

  # Extract all fragments within the ampregion
  cnvInAmp <- subset(cnvTable, (as.character(ampregion[a,2])==as.character(cnvTable[,2]) & as.numeric(ampregion[a,3])<=as.numeric(cnvTable[,3]) & as.numeric(ampregion[a,4])>=as.numeric(cnvTable[,4])), select=c("fragID", "chr", "start", "end", "cnvRatio"));
  selInAmp <- subset(selTable, (as.character(ampregion[a,2])==as.character(selTable[,2]) & as.numeric(ampregion[a,3])<=as.numeric(selTable[,3]) & as.numeric(ampregion[a,4])>=as.numeric(selTable[,4])), select=c("fragID", "chr", "start", "end", "cnvRatio"));

 for (f in 1:nrow(selInAmp)) {
    #     print (paste("A:", ampregion[a,2], ampregion[a,3],ampregion[a,4], "S:", cnvInAmp[f,2], cnvInAmp[f,3], cnvInAmp[f,4], sep=" "));
    par(mar=c(0,4,2,2), new=TRUE);
    # Calculate the start and end position in pseudo positions
    selStart <- start+selInAmp[f,3]-ampregion[a,3];
    selEnd <- start+selInAmp[f,4]-ampregion[a,3];

    # Draw the fragment
    lines(c(selStart,selEnd), c(selInAmp[f,5],selInAmp[f,5]));
#    plot(c(selStart,selEnd), c(cnvInAmp[f,5],cnvInAmp[f,5]), type="l")
  }
  # Go through all the selected fragments 
 for (f in 1:nrow(cnvInAmp)) {
    #     print (paste("A:", ampregion[a,2], ampregion[a,3],ampregion[a,4], "S:", cnvInAmp[f,2], cnvInAmp[f,3], cnvInAmp[f,4], sep=" "));
    par(mar=c(0,4,2,2), new=TRUE);
    # Calculate the start and end position in pseudo positions
    selStart <- start+cnvInAmp[f,3]-ampregion[a,3];
    selEnd <- start+cnvInAmp[f,4]-ampregion[a,3];
    
    # Draw the fragment
    lines(c(selStart,selEnd), c(cnvInAmp[f,5],cnvInAmp[f,5]), col="red");
#    plot(c(selStart,selEnd), c(cnvInAmp[f,5],cnvInAmp[f,5]), type="l")
  }
# Calculate the pseudo start position for the next ampregion
  start <- start+ampregion[a,4]-ampregion[a,3]+1+10;
}
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# Plot lower graph
ymax <- 1
ytics <- c(-1,0,1)

par(mar=c(4,4,0,2));
plot(1, type="n", ylab="log2(T/N)", xlab="Pseudo genome position", axes=TRUE, xlim=c(1,xmax), ylim=c(-ymax, ymax), yaxt="n");
axis(2, at=ytics, labels=NULL)
  
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
     par(mar=c(4,4,0,2), new=TRUE);
      rect(start, -ymax, start+ampregion[a,4]-ampregion[a,3]+1+10, ymax, border="grey", col="grey");
      col <- 0;
    }

    # Set lastChr to chr of this line
    lastChr <- ampregion[a,2];
  }
  # If lastChr and chr of this line is the same and the col is set to 1 continue draw a grey background
  else {
    if (col==1) {
     par(mar=c(4,4,0,2), new=TRUE);
      rect(start, -ymax, start+ampregion[a,4]-ampregion[a,3]+1+10, ymax, border="grey", col="grey");
    }
  }
 par(mar=c(4,4,0,2), new=TRUE);
 lines(c(start,start+ampregion[a,4]), c(0,0), lty=2, col="black", axes=FALSE);

# Extract all fragments within the ampregion
 arrayInAmp <- subset(snpArray, (as.character(ampregion[a,2])==as.character(snpArray[,2]) & as.numeric(ampregion[a,3])<=as.numeric(snpArray[,3]) & as.numeric(ampregion[a,4])>=as.numeric(snpArray[,4])), select=c("fragID", "chr", "start", "end", "cnvRatio", "Cn"))

 for (g in 1:nrow(arrayInAmp)) {
     par(mar=c(4,4,0,2), new=TRUE);
     
     #Calculate the start and end position in pseudo positions
    arrayStart <- start+arrayInAmp[g,3]-ampregion[a,3];
    arrayEnd <- start+arrayInAmp[g,4]-ampregion[a,3];
    
    # Draw the fragment
  lines(c(arrayStart,arrayEnd), c(arrayInAmp[g,5],arrayInAmp[g,5]), col="blue", ylab="", xlab="") ;

     #lines( c(arrayStart,arrayEnd), c((log(arrayInAmp[g,6])/log(2)-1),(log(arrayInAmp[g,6])/log(2)-1)), col="red", axes=FALSE, ylab="", xlab="")
  }

  # Calculate the pseudo start position for the next ampregion
  start <- start+ampregion[a,4]-ampregion[a,3]+1+10;
}
dev.off();
