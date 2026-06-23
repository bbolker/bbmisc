## https://opentraits.org/datasets/avonet.html
## Download from https://figshare.com/ndownloader/files/34480865?private_link=b990722d72a26b5bfead
library(tidyverse)

f0 <- c("AVONET_Raw_Data.csv", "AVONET_Extant_Species_List.csv")
files <- paste0("ELEData/TraitData/", f0) |>
    setNames(c("trait_data", "species_list"))
unzip("ELEData.zip", files)

spp <- read_csv(files[["species_list"]], col_types = "c") |>
    filter(grepl("Spheniscidae", Family.name)) |>
    select(c(Species.name, Family.name)) |>
    distinct()

dd <- (read_csv(files[["trait_data"]])
    |> right_join(spp, by = c("Species1_BirdLife" = "Species.name"))
    |> select(Species = Species1_BirdLife,
              Sex, Age, Locality,
              Country, Beak.Length_Culmen, Beak.Length_Nares,  
              Beak.Width, Beak.Depth, Tarsus.Length,
              Wing.Length, Kipps.Distance, Secondary1,
              "Hand-wing.Index", Tail.Length)
    |> mutate(genus = stringr::str_extract(Species, "^[[:alpha:]]+"),
              .after = 1,
              shortname =
                  stringr::str_replace(Species, "^([[:alpha:]])[[:alpha:]]+", "\\1\\. "))
  )
