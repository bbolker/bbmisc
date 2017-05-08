x <- read.table("timings.txt",header=TRUE,
                dec=",",
                stringsAsFactors=FALSE,
                na.strings="-")
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
xg <- gather(x,taxon,time,-Model) %>%
  mutate(Model=reorder(factor(Model),time),
         taxon=reorder(factor(taxon),time))
ggplot(xg,aes(time,Model,colour=taxon,group=taxon))+
  geom_point()+
  geom_smooth(se=FALSE)+
  scale_colour_brewer(palette="Dark2")+
  scale_x_log10()+
  labs(x="time (seconds)",y="")
ggsave("pix/timings.png")


pdf("compTime.pdf",width=8, height=5)
compTime <- x[,-1]
rownames(compTime) <- x[,1]
for(i in 1:6){
	par(mar=c(10,5,0.5,0.5),oma=c(0,4,0,0))
	plot(compTime[,i],ylim=c(0,max(compTime[,i],na.rm=TRUE)),
	     type="l",axes=FALSE,frame.plot=TRUE,ylab="",xlab="",lwd=3)

	mtext("Time in minutes",side=2, outer=TRUE,cex=2.75,line=0)
	axis(1,1:nrow(compTime),labels=rownames(compTime),las=2,cex.axis=1.5)
	axis(2,las=2,cex.axis=1.5)
}

dev.off()
