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
