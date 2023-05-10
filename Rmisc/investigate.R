#https://dynamicecology.wordpress.com/2021/03/15/how-long-do-institutional-investigations-into-accusations-of-serious-scientific-misconduct-typically-take-heres-some-data/
## https://secretariat.mcmaster.ca/app/uploads/Research-Integrity-Policy.pdf
## https://www-nature-com.libaccess.lib.mcmaster.ca/articles/d41586-020-00287-y

library(tidyverse)
dd <- read.table(header=TRUE, text="
subject duration year
Fuji 24 2010
Fuji 3 2012
Boldt 21 2010
Sata 12 2017
Sata 12 2017
Sata 12 2017
Nazari 24 2019
Stapel 11 2011
Chen 4 2009
Reuben 10 2008
Schön 4 2002
Lönnstedt 12 2016
Lönnstedt 9 2019
Møller 19 2001
")

labs <- c("", letters)
library(ggplot2); theme_set(theme_bw())
dd2 <- (dd
  %>% group_by(subject, year)
  %>% mutate(fyear=paste0(year,labs[1:n()]),
             label = factor(sprintf("%s (%s)", subject, fyear)))
  %>% mutate(across(label, ~forcats::fct_reorder(., duration)))
)
## why is this necessary?

dd2$label <- factor(dd2$label, levels=dd2$label[order(dd2$duration)])

pdur <- as.numeric(Sys.Date() - as.Date("2020-02-05"))/30.5
ggplot(dd2, aes(x=duration, y=label)) + geom_point() +
    geom_vline(xintercept=pdur,col=2,lty=2) +
    expand_limits(x=c(0,50), y = 15) +
    geom_text(y=Inf, x=pdur, hjust=-0.2,  vjust=2, col=2, label="Pruitt (2020)") +
    labs(y="", x="duration (months)", title = "Duration of investigations of scientific misconduct",
         subtitle = "Data from https://tinyurl.com/dynecology-pruittpost")
ggsave("pruitt.png")

ggsave("pruitt.pdf")

## alt text
"Duration of investigations of scientific misconduct, data from https://tinyurl.com/dynecology-pruittpost. Fuji (2021) is minimum at 4 months or so, Pruitt is max at 39 months (next longest are Nazari 2019 and Fuji 2010 at 24 months)
"

