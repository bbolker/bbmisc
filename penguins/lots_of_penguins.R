## Claude prompt:

## I am looking for a clean, freely available data set that includes morphometric measurements of at least three traits of a large number of species (at least 10), with at least 10 individuals measured for each species.

## Claude response:

## > The AVONET dataset would be my top recommendation as it's the most comprehensive, recent, and globally representative, with excellent sample sizes across species that far exceed your requirements.

## retrieve data from
## https://figshare.com/s/b990722d72a26b5bfeadu

## doi: https://doi.org/10.111/ele.13898
## Tobias, Joseph A., Catherine Sheard, Alex L. Pigot, et al. 2022. “AVONET: Morphological, Ecological and Geographical Data for All Birds.” Ecology Letters 25 (3): 581–97. https://doi.org/10.1111/ele.13898.

## filter to penguin species only

library(tidyverse)
library(tinyplot) 
library(ggalt)
library(directlabels)

if (!file.exists("penguin_traits.rds")) {
  spp <- read_csv("../ELEData/TraitData/AVONET_Extant_Species_List.csv") |>
    filter(grepl("Spheniscidae", Family.name)) |>
    select(c(Species.name, Family.name)) |>
    distinct()


  dd <- (read_csv("../ELEData/TraitData/AVONET_Raw_Data.csv")
    |> right_join(spp, by = c("Species1_BirdLife" = "Species.name"))
    |> select(Species = Species1_BirdLife,
              Sex, Age, Locality,
              Country, Beak.Length_Culmen, Beak.Length_Nares,  
              Beak.Width, Beak.Depth, Tarsus.Length,
              Wing.Length, Kipps.Distance, Secondary1,
              "Hand-wing.Index", Tail.Length)
    |> mutate(genus = stringr::str_extract(Species, "^[[:alpha:]]+"), .after = 1,
              shortname = stringr::str_replace(Species, "^([[:alpha:]])[[:alpha:]]+", "\\1\\. "))
  )

  saveRDS(dd, "penguin_traits.rds")
} else {
  dd <- readRDS("penguin_traits.rds")
}

with(dd, table(Species))
plt(Beak.Depth ~ Beak.Width | Species, data = dd)

## with ggplot so we can do ellipses
theme_set(theme_bw())
gg0 <- ggplot(dd, aes(Beak.Width, Beak.Depth, color = Species, shape = genus)) +
  geom_point(size=4)

## gg0 + stat_ellipse(aes(fill=Species), geom = "polygon", alpha = 0.1)

gg0 + geom_encircle(aes(fill=Species), alpha = 0.1, expand=0.01) + geom_dl(method="smart.grid", aes(label=shortname)) +
  guides(color="none", fill="none")

ggsave("penguin_traits.pdf", width=8, height=6)
