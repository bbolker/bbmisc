## from Wolfgang Huber
## https://bsky.app/profile/wkhuber.bsky.social/post/3m7pzim5gbk2z
## https://github.com/wolfganghuber/demos/blob/master/sunrise%2Bset.R
library("suncalc")
library("ggplot2")
library("dplyr")

# Heidelberg coordinates
lat = 49.38
lon = 8.71

## Hamilton
lat = 43.26
lon = -79.87
tz = "America/Toronto"

## Generate dates for this year
datelims <- as.Date(c("2026-01-03", "2026-01-10"))
dates = seq(datelims[1], datelims[2], by = "day")

# Calculate sunrise and sunset times
suntimes = suncalc::getSunlightTimes(
  date = dates,
  lat = lat, lon = lon,
  tz = tz
)

suntimes |> select(c(sunrise, sunset))

## compare with
## https://www.timeanddate.com/sun/canada/hamilton
## for Jan 1-6
read.table(header=TRUE, text = "
1 7:51 am 4:54 pm 9:03:05 +0:46
2 7:51 am 4:55 pm 9:03:56 +0:51
3 7:51 am 4:56 pm 9:04:52 +0:55
4 7:51 am 4:57 pm 9:05:51 +0:59
5 7:51 am 4:58 pm 9:06:55 +1:03
6 7:51 am 4:59 pm 9:08:02 +1:07
"

# Extract hour of day as decimal (e.g., 6.5 = 6:30 AM)
hour = function(x) as.numeric(format(x, "%H")) + as.numeric(format(x, "%M")) / 60 + as.numeric(format(x, "%S")) / 3600

suntimes = mutate(suntimes, 
   sunriseHour = hour(sunrise),
   sunsetHour  = hour(sunset))

EarliestNight = slice_min(suntimes, sunsetHour)
LatestDay     = slice_max(suntimes, sunriseHour)

ggplot(suntimes, aes(x = sunriseHour, y = sunsetHour, label = format(date))) + 
  geom_point(col = "#b0b0b0") +
  geom_path(col = "lightblue") +
  geom_point(data = EarliestNight, col = "red") +
  geom_point(data = LatestDay, col = "red") +
  geom_text(data = EarliestNight, vjust = "top") +
  geom_text(data = LatestDay, hjust = "right") +
  xlab("Sunrise (hour)") + ylab("Sunset (hour)") +
  ggtitle("Sunrise vs Sunset Times in Heidelberg") 

ggsave("sunrise+set.png", width = 160, height = 160, units = "mm")
