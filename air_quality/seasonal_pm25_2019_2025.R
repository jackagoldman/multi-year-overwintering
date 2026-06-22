# Air quality analysis: seasonal PM2.5 trends in Fort Nelson region, 2019–2025
# Author: Jack A. Goldman
# Date: 2024-07-17

library(tidyverse)
library(lubridate)
library(here)
library(sf)
library(exactextractr)
library(ggplot2)

# get sdd average in Fort Nelson for 2023
sdd_23 <- raster::raster("/Volumes/top-secret/PostDoc/multi-year-overwintering/data/snow_dd/SDD_2023.tif")
# for nelson
nelson_sf <- st_read("/Volumes/top-secret/PostDoc/multi-year-overwintering/data/study_aoi/study_firezones.shp") |> st_transform(st_crs(sdd_23)) |> dplyr::filter(MFFRZNNM == "Fort Nelson Fire Zone")
# calculate average mean SDD in nelson
nelson_sdd_23 <- exact_extract(sdd_23, nelson_sf, 'mean')

AQ_raw <- read.csv(here('data', 'air_quality', 'Fort_Nelson_Region_AQ_2019-2026.csv'))

# Parse datetime, filter 2019–2025 
AQ_raw$Date_Time <- ymd_hms(paste(AQ_raw$Date_stamp, AQ_raw$Time_stamp))
AQ_raw$Date_stamp <- as.Date(AQ_raw$Date_stamp)

AQ_raw <- AQ_raw |>
  filter(year(Date_stamp) >= 2019, year(Date_stamp) <= 2025)

#  Aggregate hourly → daily, averaging across sensors where >1 active ──
AQ_daily <- AQ_raw |>
  dplyr::group_by(Date_stamp) |>
  dplyr::summarise(
    pm25 = mean(pm2.5_alt.3.4, na.rm = TRUE),
    .groups = 'drop'
  )

#  Derive time variables 
AQ_daily <- AQ_daily |>
  mutate(
    Year       = as.factor(year(Date_stamp)),
    Month      = month(Date_stamp, label = TRUE, abbr = FALSE),
    Julian_Day = yday(Date_stamp),
    Season     = case_when(
      month(Date_stamp) %in% c(1, 2, 3)   ~ 'Winter',
      month(Date_stamp) %in% c(4, 5, 6)   ~ 'Spring',
      month(Date_stamp) %in% c(7, 8, 9)   ~ 'Summer',
      month(Date_stamp) %in% c(10, 11, 12) ~ 'Fall'
    ),
    Season = ordered(Season, levels = c('Winter', 'Spring', 'Summer', 'Fall'))
  )


#  Descriptive stats: mean, median, SD, min, max per Season × Year 
desc_stats <- AQ_daily |>
  dplyr::group_by(Year, Season) |>
  dplyr::summarise(
    n      = n(),
    mean   = round(mean(pm25, na.rm = TRUE), 2),
    median = round(median(pm25, na.rm = TRUE), 2),
    sd     = round(sd(pm25, na.rm = TRUE), 2),
    min    = round(min(pm25, na.rm = TRUE), 2),
    max    = round(max(pm25, na.rm = TRUE), 2),
    .groups = 'drop'
  )

write_csv(desc_stats, 'air_quality/pm25_seasonal_descriptive_stats.csv')
print(desc_stats, n = Inf)


# Colour palette: one colour per year 
year_levels <- levels(AQ_daily$Year)
year_cols   <- setNames(
  colorRampPalette(c('#2166ac', '#74add1', '#abd9e9', '#fee090', '#f46d43', '#d73027', '#a50026'))(length(year_levels)),
  year_levels
)


# Plot 1: Day-of-year time series, coloured by year 
sdd_lines <- data.frame(
  xintercept = c(160.43, 165.14, 167.72),
  label      = factor(c('2023', '2024', '2025'), levels = year_levels)
)

p1 <- ggplot(
    dplyr::filter(AQ_daily, Year != 2019),
    aes(x = Julian_Day, y = pm25, colour = Year)
  ) +
  geom_vline(
    data = sdd_lines,
    aes(xintercept = xintercept, colour = label),
    linewidth = 0.6, linetype = 'dashed', show.legend = FALSE
  ) +
  geom_smooth(method = 'loess', span = 0.15, se = FALSE, linewidth = 0.8) +
  scale_colour_manual(values = year_cols, breaks = year_levels) +
  scale_x_continuous(
    breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
    labels = month.abb,
    expand = c(0, 0)
  ) +
  labs(
    x      = 'Day of Year',
    y      = expression(PM[2.5] ~ '(' * mu * g~m^{-3} * ')'),
    colour = 'Year'
  ) +
  guides(
    colour = guide_legend(position = 'right', ncol = 1, title.position = 'top')
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank()
  )
p1


ggsave('air_quality/pm25_seasonal_timeseries.png', p1, width = 10, height = 5, dpi = 300)





#  Plot 2: Grouped bar chart: mean PM2.5 by year, bars coloured by season 


season_cols <- c(
  'Winter' = '#fee090',
  'Spring' = '#f46d43',
  'Summer' = '#d73027',
  'Fall'   = '#a50026'
)

p2 <- ggplot(dplyr::filter(desc_stats, Year != 2019), aes(x = Year, y = mean, fill = Season)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           colour = 'black', linewidth = 0.3) +
  scale_fill_manual(values = season_cols) +
  labs(
    x    = 'Year',
    y    = expression('Mean' ~ PM[2.5] ~ '(' * mu * g~m^{-3} * ')'),
    fill = 'Season',
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position  = 'right'
  )

ggsave('air_quality/pm25_seasonal_barplot.png', p2, width = 10, height = 5, dpi = 300)


# patchwork p1 and p2
library(patchwork)
p_combined <- p1 / p2 + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 'bold', size = 14))
p_combined

ggsave('air_quality/pm25_grouped_timeseries_barplot.png', p_combined, width = 10, height = 10, dpi = 300)



