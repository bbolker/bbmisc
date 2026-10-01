## Monthly international visitor arrivals to the US by country of residence,
## 2000-present (NTTO I-94 Arrivals Program)
## https://www.trade.gov/i-94-arrivals-program
library(tidyverse)
library(readxl)
library(janitor)

url <- paste0("https://www.trade.gov/sites/default/files/2024-06/",
              "Monthly%20Arrivals%202000%20to%20Present%20%E2%80%93%20",
              "Country%20of%20Residence%20%28COR%29_1.xlsx")
fn <- tempfile(fileext = ".xlsx")
download.file(url, fn, mode = "wb")

## row 1 holds the dates; read everything as text because the date header
## mixes Excel serial numbers (through 2016-12) with "YYYY-MM" strings
## (single-digit months in 2022 are written "YYYY-M"), some tagged
## "\r\nPreliminary"
raw <- read_excel(fn, sheet = "Monthly", col_names = FALSE,
                  col_types = "text", na = c("", "-"),
                  .name_repair = "minimal")
hdr <- unlist(raw[1, ])

ym_rx <- "^\\d{4}-\\d{1,2}"
dates <- tibble(col = seq_along(hdr), hdr = hdr) |>
  slice(-(1:3)) |>
  filter(str_detect(hdr, paste0("^\\d{5}$|", ym_rx))) |>
  ## the Excel serials fall at arbitrary days within each month
  mutate(date = convert_to_date(hdr, character_fun = \(x) ym(str_extract(x, ym_rx))) |>
           floor_date("month"),
         preliminary = str_detect(hdr, "Preliminary"))
## every month from the first to the last must be present exactly once
stopifnot(identical(dates$date,
                    seq(min(dates$date), max(dates$date), by = "month")))

arrivals <- raw[-1, c(2, 3, dates$col)] |>
  set_names(c("country", "region", dates$col)) |>
  filter(!is.na(country)) |>
  pivot_longer(-c(country, region), names_to = "col", values_to = "arrivals",
               names_transform = as.integer) |>
  left_join(select(dates, -hdr), by = "col") |>
  transmute(country = str_trim(country), region, date, preliminary,
            arrivals = as.numeric(arrivals))

## 12-month trailing mean; the pivot gives every country a complete monthly
## series, so a 12-row window is 12 months. Missing months are dropped from
## the window; an all-missing window gives NaN, recoded to NA.
arrivals <- arrivals |>
  arrange(country, date) |>
  mutate(arrivals_12mo = zoo::rollmeanr(arrivals, k = 12, fill = NA,
                                        na.rm = TRUE) |>
           na_if(NaN),
         n_months = zoo::rollsumr(!is.na(arrivals), k = 12, fill = NA),
         .by = country)

## averages over a full 12 months only (a missing peak or trough month
## biases the average); `arrivals` keeps every monthly observation
arrivals_smooth <- filter(arrivals, n_months == 12)

## plot
theme_set(theme_bw())
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
opts <- options(ggplot2.discrete.colour = okabe_ito)

## log-axis ticks placed as in base R's axisTicks(..., log = TRUE)
base_breaks <- function(n = 10) {
  function(x) axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
}

countries <- c("Canada", "Mexico", "United Kingdom", "Japan", "France")
p <- arrivals_smooth |>
  filter(country %in% countries) |>
  mutate(country = factor(country, levels = countries)) |>
  ggplot(aes(date, arrivals_12mo / 1000, colour = country)) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(trans = "log10", breaks = base_breaks(),
                     labels = scales::label_comma()) +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y",
               expand = expansion(mult = 0.02)) +
  labs(x = NULL,
       y = "Monthly visits (×1000, 12-month trailing mean)",
       caption = "Source: NTTO I-94 Arrivals Program, US Department of Commerce") +
  ## direct labels are drawn outside the panel, in a widened right margin
  coord_cartesian(clip = "off") +
  ## L-shaped frame (left and bottom axes only), so the labels aren't
  ## separated from their lines by the panel border
  theme(panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "grey20"),
        plot.margin = margin(5.5, 90, 5.5, 5.5, "pt"))
## labels start 0.1 cm to the right of each line's end. "last.qp" fails
## because the lines end on different dates (Canada and Mexico have no
## preliminary months); "bumpup" separates overlapping labels instead
p <- directlabels::direct.label(p, list("last.points", cex = 0.8,
                                        directlabels::dl.trans(x = x + 0.1),
                                        "bumpup"))
ggsave("i94_arrivals.png", p, width = 8, height = 5)
options(opts)
