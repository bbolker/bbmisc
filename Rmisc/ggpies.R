## floating pies

## https://stackoverflow.com/questions/9684371/create-floating-pie-charts-with-ggplot
## https://stackoverflow.com/questions/47718439/multiple-pies-with-scatterpie-or-ggforce
## https://stackoverflow.com/questions/35532683/r-inserting-a-picture-into-faceted-ggplot-using-annotation-custom
## https://stackoverflow.com/questions/10673481/add-a-geom-layer-for-a-single-panel-in-a-faceted-plot
## https://clarewest.github.io/blog/post/2019-08-27-combining-inset-plots-with-facets-using-ggplot/

library(tidyverse)
library(gtools)
library(ggplot2)
library(ggimage)
source("annotation.R")   ## hacked geom_subview/annotation_custom

rd <- function(n, alpha) {
    rdirichlet(n, alpha) |> as.data.frame() |> setNames(paste0("c", 1:length(alpha)))
}

## fake data
set.seed(101)
df <- (expand_grid(x = LETTERS[1:4],
                  f = letters[1:5])
    |> mutate(y = rgamma(length(x), scale = 3, shape = 5),
              s = rgamma(length(x), scale = 3, shape = 5))
    ## would like to use dynamic nrow(.data) but ?? (need %>% ?)
    |> bind_cols(rd(20, alpha = c(3,2,1)))
    |> mutate(across(where(is.character), factor))
)

## construct a single pie chart from a one-row tibble
## FIXME: add colour specification?
pie_fun <- function(x) {
    (x
        |> pivot_longer(everything())
        |> ggplot(aes(x=1, y=value, fill = name)) +
           geom_bar(stat="identity", width=1) +
           coord_polar(theta="y") +
           theme_void() +
           theme(legend.position="none") +
           theme_transparent()
    )
}

## construct pie chart
dfpie <- (df
    |> nest(data = c(c1, c2, c3))
    |> mutate(pie = map(data, pie_fun))
    ## all mappings, factor/numeric conversion etc. must be done up front!
    ##  (geom_subview() does arithmetic directly on x, extracts mapped columns
    ##   from data **by name**)
    |> mutate(s = s/10, x = as.numeric(x))
)

## base plot
p <- p0 <- (ggplot(dfpie, aes(x=x, y=y))
    + geom_point()
    + facet_wrap(~f))

## add subviews
p <- p + geom_subview(data = dfpie,
                      aes(x = x,
                          y = y,
                          width = s,
                          height = s,
                          subview = pie))
                          
