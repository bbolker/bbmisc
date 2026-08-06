library(tidyverse)
library(mgcv)
library(Matrix)

theme_set(theme_bw())
zmargin <- theme(panel.spacing = grid::unit(0, "lines"))

dd <- read_csv("metfin.csv") |>
  select(STA_PROV:WYear,TMINfin) |>
  mutate(extremecold = as.numeric(TMINfin < -18))

length(unique(dd$Site_ID)) ## 37 sites

gg0 <- ggplot(dd, aes(doy2, TMINfin)) +
  geom_line(aes(colour = WYear)) +
  facet_wrap(~Site_ID) +
  zmargin +
  scale_colour_viridis_c()

## glance at a subset of data
gg0 %+% filter(dd, substr(Site_ID, 1, 1) == "3", WYear > 1990, WYear < 1995)

## not sure what doy22 is ... ?
dd_sum <- dd |>
  group_by(WYear, doy2) |>
  summarise(
    extremecold = sum(extremecold, na.rm = TRUE),
    .groups = "drop")

tt <- with (dd, table(WYear, doy2))
image(Matrix(as.matrix(tt)), useRaster = TRUE, useAbs = FALSE)
tt2 <- tt
tt2[tt2>37] <- NA
image(Matrix(as.matrix(tt2)), useRaster = TRUE, useAbs = FALSE)
## one day in the middle (which one?) with n>37

image(Matrix(as.matrix(tt>37)))

colnames(tt)[colSums(tt>37)>0] ## day 61 is double-counted? what day is that?

with(dd, table(table(WYear, doy2))) ## shouldn't have values > 37 ?

gg1 <- ggplot(dd_sum, aes(doy2, extremecold)) +
  geom_line() + 
  facet_wrap(~WYear) +
  zmargin

gg1 + aes(y = qlogis((extremecold+0.5)/38)) +
  geom_smooth(method = "lm", formula = y ~ poly(x,2))

gg1 + aes(y = qlogis((extremecold+0.5)/38)) +
  geom_smooth(method = "lm", formula = y ~ sin(2*pi*x/365) + cos(2*pi*x/365))



## double-check WYear, doy2 ...

## mgcv? glmmTMB?
## cyclic spline x year spline?
## binomial? beta-binomial?
## go back and fit a spatiotemporal model (changes in seasonality in space and time?)

glm(extremecold/37 ~ (sin(2*pi*doy2/365) + cos(2*pi*doy2/365))*WYear,
    family = binomial, data = dd_sum)
