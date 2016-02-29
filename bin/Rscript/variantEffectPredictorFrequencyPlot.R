#!/usr/bin/env Rscript
library("ggplot2")

#args <- commandArgs(TRUE)
#data <- read.table(args[1]);
data <- read.table("../testGather/variants_rfile.txt");
names(data)<-c("roi","pos","count","type","group")
#names(data)<-c("roi","roinr","pos","count","type")

# roi vector
rois<-unique(data$roi);

# bin size
bin <- 100;

# prepare plots
#pdf("out.pdf")
jpeg("out.jpg", width=5000, height=300);

pseudoData <- data.frame();

# loop over each
for (r in rois) {

	sub <- subset(data,roi==r);

	start<-min(sub$pos)
	max<-min(sub$pos)

	if(dim(sub)[1]>0) {
		pseudoVec <- c(1:dim(sub)[1])
		sub <- cbind(sub,pseudoVec);
	}

	pseudoData <- rbind(pseudoData,sub)

	# overwrite exclude all samples
	pseudoData <- subset(pseudoData,!(type=="Frequency" & group=="All") & group!="HM")

#	psudoD <- subset(psuedoData,!(type==all && 

	plot<-ggplot(data=pseudoData,aes(pseudoData$pseudo))+geom_point(aes(y=count,colour=roi),size=1)+facet_grid(group~roi, space="free",scales="free")+opts(legend.position="none",title="",axis.text.x=theme_blank())

	# normalize each position

#	plot<-ggplot(data=subvar,aes(pseudo))+geom_point(aes(y=factor(count)))+opts(legend.position="none",title=r,axis.text.x=theme_text(angle=90))+facet_grid(type~.,scales="free",space="free")
#	print(plot)

}

dev.off()
