start_date <- "2023-06-21"
end_date <- "2023-07-01"
refdate  <- "2023-06-21"

start_date <- "2024-12-01"
end_date <- "2024-12-12"
refdate <- "2024-12-06"
season <- "winter"

library(ggplot2); theme_set(theme_bw())
library(colorspace)
library(purrr)
library(dplyr)
## need archived package for sunrise/sunset calcs
##  (no compiled code)
while (!require("RAtmosphere")) {
    remotes::install_version("RAtmosphere", "1.1")
}
datevec <- seq.Date(as.Date(start_date), as.Date(end_date), by = "1 day")
## julian days (day-of-year)
dvec <- as.numeric(format(datevec, "%j"))

lonlat <- list(
    Hamilton_ON = c(-79.844, 43.255),
    Somerville_MA = c(-71.1022, 42.388),
    Dover_NH = c(-71.0284,43.3267),
    Newton_MA = c(-71.2053,  42.3612))

## compute sunrise/sunset times for all location
dd <- (purrr::map_dfr(lonlat,
                     ~ data.frame(
    date = datevec,
    suncalc(dvec, Lat = .[2], Lon = .[1])),
    .id = "place")
)

sdiff_fun <- function(x, time = "sunset", season = "summer") {
    if ((season == "summer" && time == "sunset")  ||
        (season == "winter" && time == "sunrise"))
        return(max(x) - x)
    ## else
    return(x - min(x))
}

## scale to 0 at latest sunset
dd2 <- (dd
    |> group_by(place)
    |> mutate(sunset_diff = sdiff_fun(sunset, "sunset", season)*60,
              sunrise_diff = sdiff_fun(sunrise, "sunrise", season)*60)
    |> ungroup()
)
    
                                                 
ggplot(dd2, aes(date, sunset_diff, colour = place)) +
    geom_line(aes(linetype = place)) + 
    scale_color_discrete_qualitative() +
    geom_vline(xintercept = as.Date(refdate), lty = 2) +
    labs(y = "sunset diff (minutes)")
ggsave("analemma.png", width = 5, height = 5)
