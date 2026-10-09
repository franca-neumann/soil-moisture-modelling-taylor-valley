# the aim is to plot predictor developement for aspa 131 trend pixels
# plot all predictors for trend pixels, each individual pixel? 
library(terra)
library(sf)
library(purrr)
AOA_season_stack <- rast("data/processed_data/AOA_masks/season_wise_masks.tif")
AOA_season_stack[AOA_season_stack==1] <- NA
AOA_overall_mask <- rast("data/processed_data/AOA_masks/overall_mask.tif")
AOA_overall_mask[AOA_overall_mask==1] <- NA


aspa_131 <- st_read("data/raw_data/aspa_131/APAD_ID_036/131.shp")
aspa_131 <- st_transform(aspa_131, crs = "EPSG:32758")
aspa_131_v <- vect(aspa_131)

mask_pred_to_aspa <- function(x, season, aspa = aspa_131_v){
  message("processing season", season )
  x$stack <- crop(x$stack, aspa)
  x$stack <- mask(x$stack, aspa)
  return(x)
}

prediction_list_aspa <- imap(
  prediction_list,
  ~ mask_pred_to_aspa(.x, season = .y)
)

mask_event_to_aspa <- function(x, season, aspa = aspa_131_v){
  x <- crop(x, aspa)
  x <- mask(x, aspa)
  return(x)
}

AOA_season_stack_aspa <- mask_event_to_aspa(AOA_season_stack)
AOA_overall_mask_aspa <- mask_event_to_aspa(AOA_overall_mask)


load_predictor_stacks <- function(path) {
  
  files <- list.files(
    path = path,
    pattern = "^predictors_\\d{8}\\.tif$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  dates <- sub(
    "^predictors_(\\d{8})\\.tif$",
    "\\1",
    basename(files)
  )
  
  files |>
    purrr::set_names(dates) |>
    purrr::map(terra::rast)
  
}

all_stacks_AOI <- load_predictor_stacks(
  "D:/antarctica/data/processed_data/snow_stacks"
)

# all in the same EPSG? 
all(
  purrr::map_lgl(
    all_stacks_AOI,
    ~terra::same.crs(.x, "EPSG:32758")
  )
) #yes

all_stacks_aspa <- map(all_stacks_AOI, mask_event_to_aspa)

#load trend and mask to aspa
season_mean_trend <- rast("data/processed_data/season_mean_trend_analysis_zyp_zhang/season_mean_trend_zyp_zhang.tif")
season_mean_trend_aspa <- mask_event_to_aspa(season_mean_trend)
plot(season_mean_trend_aspa)
trend_aspa <- season_mean_trend_aspa$practical_relevance
plot(trend_aspa)

slope_aspa <- season_mean_trend_aspa$slope_total

library(terra)
library(dplyr)
library(purrr)


##########################################################################

library(terra)
library(dplyr)
library(purrr)
library(tidyr)
library(lubridate)
library(ggplot2)

# relevant pixels

cells <- which(values(trend_aspa) == 1)

xy <- terra::xyFromCell(trend_aspa, cells)

pixel_vect <- terra::vect(
  xy,
  type = "points",
  crs = terra::crs(trend_aspa)
)

pixel_vect$pixel_id <- seq_along(cells)
pixel_vect$cell <- cells

# slope for pixels
pixel_vect$sen_slope <- values(slope_aspa)[cells]


## predictor values

predictor_df <- purrr::imap_dfr(
  all_stacks_AOI,
  function(r, date) {
    
    # resample if not matching
    if (!isTRUE(all.equal(
      terra::res(r),
      terra::res(trend_aspa)
    ))) {
      
      r <- terra::resample(
        r,
        trend_aspa,
        method = "average"
      )
    }
    
    vals <- terra::extract(
      r,
      pixel_vect
    )
    
    vals |>
      mutate(
        pixel_id = pixel_vect$pixel_id,
        cell = pixel_vect$cell,
        sen_slope = pixel_vect$sen_slope,
        date = as.Date(date, "%Y%m%d")
      ) |>
      select(
        pixel_id,
        cell,
        sen_slope,
        date,
        everything(),
        -ID
      )
    
  }
)


## seson mean snow dist class

snow_season <- predictor_df |>
  filter(month(date) %in% c(11, 12, 1, 2)) |>
  mutate(
    season = if_else(
      month(date) >= 11,
      paste0(
        substr(year(date), 3, 4),
        "|",
        substr(year(date) + 1, 3, 4)
      ),
      paste0(
        substr(year(date) - 1, 3, 4),
        "|",
        substr(year(date), 3, 4)
      )
    )
  ) |>
  filter(!season %in% c("17|18", "19|20")) |>
  group_by(
    pixel_id,
    cell,
    sen_slope,
    season
  ) |>
  summarise(
    snow_dist_class = median(
      snow_dist_class,
      na.rm = TRUE
    ),
    .groups = "drop"
  )



snow_season$season <- factor(
  snow_season$season,
  levels = sort(unique(snow_season$season))
)

snow_season$season_num <- seq_along(
  levels(snow_season$season)
)[snow_season$season]


# plot
windowsFonts(
  Calibri = windowsFont("Calibri")
)


library(RColorBrewer)

slope_cols <- colorRampPalette(
  brewer.pal(11, "BrBG")
)(100)


slope_cols_neg <- slope_cols[0:51]


ggplot(
  snow_season,
  aes(
    x = season_num,
    y = snow_dist_class,
    group = pixel_id,
    colour = sen_slope
  )
) +
  geom_line(
    alpha = 0.6,
    linewidth = 0.6,
    linetype = "dashed"
  ) +
  geom_point(
    size = 1.5
  ) +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    se = TRUE,
    colour = "black",
    linewidth = 1.2
  ) +
  scale_x_continuous(
    breaks = seq_along(levels(snow_season$season)),
    labels = levels(snow_season$season)
  ) +
  scale_y_continuous(
    breaks = seq(
      floor(min(snow_season$snow_dist_class, na.rm = TRUE)),
      ceiling(max(snow_season$snow_dist_class, na.rm = TRUE)),
      by = 1
    )
  ) +
  scale_colour_gradientn(
    colours = slope_cols_neg,
    limits = c(
      min(snow_season$sen_slope, na.rm = TRUE),
      0
    ),
    name = "Sen's slope"
  ) +
  theme_classic(
    base_family = "Calibri"
  ) +
  labs(
    x = "Austral summer season",
    y = "Median snow distance class"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )
ggsave(
  filename = "plots/snow_distance_class_trend_pixels.png",
  width = 160,
  height = 110,
  units = "mm",
  dpi = 600,
  bg = "white"
)

#########################################################################
# load and check snow layers
load_snow_stacks <- function(path) {
  
  files <- list.files(
    path = path,
    pattern = "^snow_\\d{8}\\.tif$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  dates <- sub(
    "^snow_(\\d{8})\\.tif$",
    "\\1",
    basename(files)
  )
  
  files |>
    purrr::set_names(dates) |>
    purrr::map(terra::rast)
  
}

all_stacks_AOI_snow <- load_snow_stacks(
  "D:/antarctica/data/processed_data/snow_stacks"
)

# all in the same EPSG? 
all(
  purrr::map_lgl(
    all_stacks_AOI_snow,
    ~terra::same.crs(.x, "EPSG:32758")
  )
) #no

all_stacks_AOI_snow <- purrr::map(
  all_stacks_AOI_snow,
  ~ if (!terra::same.crs(.x, "EPSG:32758")) {
    terra::project(.x, "EPSG:32758")
  } else {
    .x
  }
)


all_stacks_aspa_snow <- map(all_stacks_AOI_snow, mask_event_to_aspa)

#mask to AOI
glaciers <- st_read("C:/Users/mini4/Desktop/Uni/Master/antarctica/data/raw_data/Antarctic_GIS/MDV_GLACIERS_STREAMS_POLYGONS/MDV_GLACIERS.shp")
lakes_ponds <- st_read("C:/Users/mini4/Desktop/Uni/Master/antarctica/data/raw_data/GIS layers/MCMLTER_GIS_Export_Layers/USGS_Map_1970/Lakes_and_Poonds_epsg32758.gpkg")
glaciers <- st_transform(glaciers, crs = st_crs(lakes_ponds))
AOI_mask <- st_union(glaciers, lakes_ponds)
plot(AOI_mask)
AOI_mask <- st_transform(AOI_mask, crs = "EPSG:32758")

all_stacks_aspa_snow_mask_trend <- purrr::map(
  all_stacks_aspa_snow,
  ~ {
    trend_res <- terra::resample(
      trend_aspa,
      .x,
      method = "near"
    )
    
    terra::ifel(
      trend_res == 1,
      2,
      .x
    )
  }
)

walk(all_stacks_aspa_snow_mask_trend, ~plot(.x))
