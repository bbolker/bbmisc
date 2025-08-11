library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)

if (interactive() && !file.exists("tab.csv")) {

  ## select only response-values rectangle
  ## n.b. I haven't tried this since I cropped the PDF
  tt4 <- extract_areas("tab_crop.pdf")[[1]]
  ## convert all vals to numeric, reshape to matrix
  m <- matrix(readr::parse_number(tt4),ncol=ncol(tt4))
  ## add left-hand indicator columns, and column names, manually
  nms <- c("dataset","r","type","MGHD.ERR","MGHD.ARI",
           "MST.ERR","MST.ARI","MI/MGHD.ERR","MI/MGHD.ARI","MI/MST.ERR","MI/MST.ARI")
  xcols <- expand.grid(type=c("est","sd"),r=c(0.05,0.1,0.2,0.3),
                       dataset=paste0("sim",1:6))[3:1]
  dd <- setNames(data.frame(xcols,m),nms)
  write.csv(dd, file="tab.csv", row.names=FALSE)
} else {
  dd <- read_csv("tab.csv")
}

## rearrange
dd2 <- (dd 
    |> pivot_longer(-c(dataset,r,type))
    |> separate(name,into=c("model","stat"),sep="\\.") 
    |> pivot_wider(names_from=type,values_from=value)
    |> mutate(across(dataset, factor))
)

## independent information on simulation codes
simtab <- read.table(header=TRUE,text="
dataset distribution covstruc separation
sim1 MGHD VEE well-separated
sim2 MGHD VEE overlapping
sim3 MST VEI well-separated
sim4 MST VEI overlapping
sim5 GMM VEE well-separated
sim6 GMM VEE overlapping
")

labs <- with(simtab,sprintf("%s/%s\n%s",distribution,covstruc,separation))

## slightly dangerous, assumes matching order ...
levels(dd2$dataset) <- labs

## make the picture
ggplot(dd2,aes(factor(r),est,colour=model)) + 
  ## points and lines
    geom_point()+geom_line(aes(group=model)) +
  ## transparent ribbons, +/- 1 SD
    geom_ribbon(aes(ymin=est-sd,ymax=est+sd,group=model,fill=model),
                colour=NA,
     alpha=0.3)+
  facet_grid(stat~dataset)+
  ## the rest is cosmetic
  theme_bw()+labs(x="r (proportion missing)",y="")+
  theme(panel.spacing=grid::unit(0,"lines"))+
  scale_colour_brewer(palette="Set1")+
    scale_fill_brewer(palette="Set1")
    
## TO DO: order datasets (MGHD,MST,MI/MGHD,MI/MST)
## include characteristics of sims? order differently?
## scale="free_y" ?

## could adjust panel backgrounds, e.g. to distinguish between
## well-separated and overlapping cases

## direct labels?

ggsave("fig.png", width = 10, height = 6)

## ANOVA table
## factor plot
## interaction.plot()
