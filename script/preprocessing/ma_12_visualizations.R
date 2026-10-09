#visualizations
library(terra)
library(sf)
library(purrr)
library(tmap)
library(lubridate)
library(stringr)
library(av)

season_folders <- list.files("data/processed_data/predictions/", full.names = T)
season_folders <- season_folders[-3] # exclude seasons which are not fulfilling the requirements (17|18) (19|20)
season_folders <- season_folders[-4]
season_folders
season_folders_name <- list.files("data/processed_data/predictions/", full.names = F)
season_folders_name <- season_folders_name[-3] # exclude seasons which are not fulfilling the requirements (17|18) (19|20)
season_folders_name <- season_folders_name[-4]
extract_date <- function(x) {
  as.Date(str_extract(x, "\\d{8}"), format = "%Y%m%d")
}

load_predictions <- function(season_path) {
  
  prediction_dir <- file.path(season_path, "prediction_AOA")
  
  files <- list.files(prediction_dir, full.names = TRUE)
  
  dates <- extract_date(files)
  
  keep <- !is.na(dates)
  files <- files[keep]
  dates <- dates[keep]
  
  keep_first <- !duplicated(dates)
  files <- files[keep_first]
  dates <- dates[keep_first]
  
  rasters <- map(files, rast)
  
  # sensor reference
  s2_idx <- which(str_detect(files, "S2"))[1]
  ps_idx <- which(str_detect(files, "PS"))[1]
  
  ref_idx <- if (!is.na(s2_idx)) {
    s2_idx
  } else if (!is.na(ps_idx)) {
    ps_idx
  } else {
    stop("No S2 or PS reference found in: ", season_path)
  }
  
  ref <- rasters[[ref_idx]]
  
  
  rasters <- map(rasters, ~{
    if (!compareGeom(.x, ref, stopOnError = FALSE)) {
      resample(.x, ref, method = "bilinear")
    } else {
      .x
    }
  })
  
  s <- rast(rasters)
  
  names(s) <- basename(files)
  
  list(
    stack = s,
    dates = dates,
    ref_type = if (!is.na(s2_idx)) "S2" else "PS"
  )
}

prediction_list <- map(season_folders, load_predictions)
names(prediction_list) <- season_folders_name

glaciers <- st_read("data/raw_data/Antarctic_GIS/MDV_GLACIERS_STREAMS_POLYGONS/MDV_GLACIERS.shp")
lakes_ponds <- st_read("data/raw_data/GIS layers/MCMLTER_GIS_Export_Layers/USGS_Map_1970/Lakes_and_Poonds_epsg32758.gpkg")
lakes_ponds <- st_transform(lakes_ponds, crs = st_crs(glaciers))
AOI_mask <- st_union(glaciers, lakes_ponds)
AOI_mask <- st_transform(AOI_mask, crs("EPSG:32758"))
str(AOI_mask)

# one color scale for all prediction situations
mins <- c()
maxs <- c()

for(season in names(prediction_list)) {
  
  r <- prediction_list[[season]]$stack
  
  mm <- global(
    r,
    c("min", "max"),
    na.rm = TRUE
  )
  
  mins <- c(mins, mm[,1])
  maxs <- c(maxs, mm[,2])
}

global_min <- min(mins, na.rm = TRUE)
global_max <- max(maxs, na.rm = TRUE)

global_min #0.1027383
global_max #0.3574558

# limit AOI to study area

 

AOI_extend <- st_read("C:/Users/mini4/Desktop/Uni/Master/antarctica/data/processed_data/AOI_Taylor.geojson")
AOI_extend <- st_transform(AOI_extend, crs = st_crs(AOI_mask))

AOI <- st_intersection(AOI_mask, AOI_extend)

#load AOA 
AOA_season_stack <- rast("data/processed_data/AOA_masks/season_wise_masks.tif")
AOA_season_stack[AOA_season_stack==1] <- NA
AOA_overall_mask <- rast("data/processed_data/AOA_masks/overall_mask.tif")
AOA_overall_mask[AOA_overall_mask==1] <- NA
# prepare and build frames
tmap_mode("plot")

#outdir <- "data/processed_data/vis_frames/predictions/"
outdir <- "data/processed_data/vis_frames/Taylor_Valley_AOA_masked/predictions/"
dir.create(outdir)
#AOI_vect <- vect(AOI_mask)
AOI_vect <- vect(AOI)

frame_id <- 1

for(season in names(prediction_list)) {
  
  stack <- prediction_list[[season]]$stack
  dates <- prediction_list[[season]]$dates
  
  for(i in seq_len(nlyr(stack))) {
    
    cat("Processing", season, "-", i, "/", nlyr(stack), "\n")
    
    # AOA-Layer
    if (season %in% names(AOA_season_stack)) {
      AOA_mask <- AOA_season_stack[[season]]
    } else {
      AOA_mask <- AOA_overall_mask
    }
    
    sm <- stack[[i]]
    
    AOI_raster <- rasterize(
      AOI_vect,
      sm,
      field = 1,
      background = NA
    )
    
    AOI_raster <- mask(AOI_raster, AOI_extend)
    
    title_text <- paste0(
      "Season ",
      season,
      "\n",
      format(dates[i], "%d.%m.%Y")
    )
    
    tm <-
      
      # out of AOI grey 
      tm_shape(AOI_raster) +
      tm_raster(
        col.scale = tm_scale_categorical(
          values = "grey80"
        ),
        col.legend = tm_legend(
          show = FALSE
        )
      ) +
      tm_add_legend(
        type = "polygons",
        labels = "Out of AOI",
        fill = "grey80"
      ) +
      # prediction
      
      tm_shape(sm) +
      tm_raster(
        col.scale = tm_scale_continuous(
          values = "blues",
          limits = c(global_min, global_max)
        ),
        col_alpha = 1,
        
        # no missing in legend
        col.legend = tm_legend(
          title = "VWC",
          show_na = FALSE
        ),
        
      ) +
      # Out of AOA
      tm_shape(AOA_mask) +
      tm_raster(
        col.scale = tm_scale_categorical(
          values = c("0" = "grey95")
        ),
        col.legend = tm_legend(show = FALSE)
      ) +tm_add_legend(
        type = "polygons",
        labels = "Out of AOA",
        fill = "grey95"
      ) +
      # No data legend
      tm_add_legend(
        type = "polygons",
        labels = "No data",
        fill = "white"
      ) +
      
      # title 
      tm_title(title_text) +
      
      tm_layout(
        text.fontfamily = "Calibri",
        legend.outside = TRUE,
        frame = FALSE
      )
    
    outfile <- file.path(
      outdir,
      sprintf("frame_%04d.png", frame_id)
    )
    
    tmap_save(
      tm,
      filename = outfile,
      width = 1800,
      height = 1400,
      dpi = 300
    )
    
    frame_id <- frame_id + 1
  }
}

#mp4 video

pngs <- list.files(
  outdir,
  pattern = "\\.png$",
  full.names = TRUE
)


pngs <- sort(pngs)

av_encode_video(
  input = pngs,
  output = file.path(outdir, "VWC_all_timesteps.mp4"),
  framerate = 2,
  vfilter = "scale=trunc(iw/2)*2:trunc(ih/2)*2"
)
rm(pngs)
##################################################################################
#same for season_means

season_mean_path <- list.files("data/processed_data/season_mean_stacks/", full.names = T)
season_mean_names <- list.files("data/processed_data/season_mean_stacks/", full.names = F)
season_mean_names <- substring(season_mean_names, first = 1, last = 5)
season_mean_list <- map(season_mean_path, rast)
names(season_mean_list) <- season_mean_names
#str(season_mean_list)
#plot(season_mean_list[[1]]$mean)
#plot(season_mean_list[[1]]$qa_flag)

outdir_mean <- "data/processed_data/vis_frames/Taylor_Valley_AOA_masked/season_means/"
frame_id <- 1

for(i in seq_along(season_mean_list)) {
  
  mean <- season_mean_list[[i]]$mean
  season <- names(season_mean_list[i])
  
  #matching AOA layer
  if (season %in% names(AOA_season_stack)) {
    AOA_mask <- AOA_season_stack[[season]]
  } else {
    AOA_mask <- AOA_overall_mask
  }
  
  cat("Processing", season)
  
  AOI_raster <- rasterize(
    AOI_vect,
    mean,
    field = 1,
    background = NA
  )
  
  AOI_raster <- mask(AOI_raster, AOI_extend)
  
  title_text <- paste0(
    "Season ",
    season)
  
  
  tm <-
    
    # out of AOI grey 
    tm_shape(AOI_raster) +
    tm_raster(
      col.scale = tm_scale_categorical(
        values = "grey80"
      ),
      col.legend = tm_legend(
        show = FALSE,
      )
    ) +
    tm_add_legend(
      type = "polygons",
      labels = "Out of AOI",
      fill = "grey80"
    ) +
    # prediction
    
    tm_shape(mean) +
    tm_raster(
      col.scale = tm_scale_continuous(
        values = "blues",
        limits = c(global_min, global_max)
      ),
      col_alpha = 1,
      
      # no missing in legend
      col.legend = tm_legend(
        title = "VWC",
        show_na = FALSE
      ),
      
    ) +
    tm_shape(AOA_mask) +
    tm_raster(
      col.scale = tm_scale_categorical(
        values = c("0" = "grey95")
      ),
      col.legend = tm_legend(show = FALSE)
    ) +
    tm_add_legend(
      type = "polygons",
      labels = "Out of AOA",
      fill = "grey95"
    ) +
    # out of QA legend
    tm_add_legend(
      type = "polygons",
      labels = "Out of QA",
      fill = "white"
    ) +
    
    # title 
    tm_title(title_text) +
    
    tm_layout(
      text.fontfamily = "Calibri",
      legend.outside = TRUE,
      frame = FALSE
    )
  
  outfile <- file.path(
    outdir_mean,
    sprintf("frame_%04d.png", frame_id)
  )
  
  tmap_save(
    tm,
    filename = outfile,
    width = 1800,
    height = 1400,
    dpi = 300
  )
  
  frame_id <- frame_id + 1
  
}

#mp4 video

pngs_mean <- list.files(
  outdir_mean,
  pattern = "\\.png$",
  full.names = TRUE
)


pngs_mean <- sort(pngs_mean)

av_encode_video(
  input = pngs_mean,
  output = file.path(outdir_mean, "VWC_season_means.mp4"),
  framerate = 2,
  vfilter = "scale=trunc(iw/2)*2:trunc(ih/2)*2"
)

#####################################################################################
# season trend visualization
#label information
#0 = no trend
#1 = neg relevant
#2 = pos relevant
#3 = sig irrelevant

season_trend_path <- list.files("data/processed_data/season_trend_analysis_zyp_zhang/", full.names = T)
season_ov_path <- list.files("data/processed_data/season_mean_trend_analysis_zyp_zhang/", full.names = T)
season_trend_path <- c(season_trend_path, season_ov_path)
season_trend_names <- list.files("data/processed_data/season_trend_analysis_zyp_zhang/", full.names = F)
season_trend_list_zyp <- map(season_trend_path, rast)

season_trend_names <- substring(season_trend_names, first = 1, last = 5)
season_ov_name <- "Season means"
season_trend_names <-c(season_trend_names, season_ov_name)
names(season_trend_list_zyp) <- season_trend_names
#plot(season_trend_list_zyp[[1]])
str(season_trend_list_zyp)

outdir_trend <- "data/processed_data/vis_frames/Taylor_Valley_AOA_masked_zhang/season_trends/"
#dir.create(outdir_trend)

# --------------------------------------------------
# globale Skalen
# --------------------------------------------------

all_slope <- c()
all_tau <- c()

for(season in names(season_trend_list_zyp)) {
  
  r <- season_trend_list_zyp[[season]]
  
  all_slope <- c(
    all_slope,
    values(r$slope_total)
  )
  
  all_tau <- c(
    all_tau,
    values(r$tau)
  )
}

global_slope <- range(all_slope, na.rm = TRUE)
global_tau   <- range(all_tau, na.rm = TRUE)


# maps
trend_cols <- c(
  "grey45",   # no trend
  "#D08B3E",  # negative trend
  "#4A7FB8",  # positive trend
  "grey25"    # irrelevant trend
)

windowsFonts(
  Calibri = windowsFont("Calibri")
)
tmap_options(text.fontfamily = "Calibri")

# as function
create_trend_maps <- function(
    season,
    r,
    AOA_mask,
    AOI_vect,
    AOI_extend,
    global_slope,
    global_tau,
    trend_cols,
    outdir_trend,
    height_value = 1800,
    width_value = 2200
) {
  
  cat("Processing", season, "\n")
  
  trend <- r$practical_relevance
  slope <- r$slope_total
  tau   <- r$tau
  valid_mask <- r$valid_mask
  
  # stable
  
  valid_mask[valid_mask == 1] <- NA
  valid_mask[valid_mask == 0] <- 1
  
  
  # AOI
  
  AOI_raster <- rasterize(
    AOI_vect,
    trend,
    field = 1,
    background = NA
  )
  
  AOI_raster <- mask(
    AOI_raster,
    AOI_extend
  )
  
  
  # slope, tau, varying
  
  trend_class <- r$practical_relevance
  
  keep <- trend_class %in% c(1, 2)
  
  slope[!keep] <- NA
  tau[!keep] <- NA
  
  varying_mask <- trend_class %in% c(0, 3)
  
  varying_mask[varying_mask == 0] <- NA
  varying_mask[varying_mask == 1] <- 1
  
  
  # trend 
  
  trend <- as.factor(trend)
  
  levels(trend) <- data.frame(
    value = c(0,1,2,3),
    label = c(
      "Varying conditions",
      "Negative trend",
      "Positive trend",
      "Irrelevant trend"
    )
  )
  
  
  tm_trend <-
    
    tm_shape(trend) +
    tm_raster(
      col.scale = tm_scale_categorical(
        values = trend_cols
      ),
      col.legend = tm_legend(
        title = paste0(
          season,
          "\nTrend"
        ),
        show_na = FALSE
      )
    ) +
    
    tm_shape(valid_mask) +
    tm_raster(
      col.scale = tm_scale_categorical(
        values = "grey60"
      ),
      col.legend = tm_legend(show = FALSE)
    ) +
    
    tm_shape(AOI_raster) +
    tm_raster(
      col.scale = tm_scale_categorical(
        values = "grey80"
      ),
      col.legend = tm_legend(show = FALSE)
    ) +
    
    tm_shape(AOA_mask) +
    tm_raster(
      col.scale = tm_scale_categorical(
        values = c("0" = "grey95")
      ),
      col.legend = tm_legend(show = FALSE)
    ) +
    
    tm_add_legend(
      type = "polygons",
      labels = "Stable conditions",
      fill = "grey60"
    ) +
    
    tm_add_legend(
      type = "polygons",
      labels = "Out of AOI",
      fill = "grey80"
    ) +
    
    tm_add_legend(
      type = "polygons",
      labels = "Out of AOA",
      fill = "grey95"
    ) +
    
    tm_add_legend(
      type = "polygons",
      labels = "Out of QA",
      fill = "white"
    ) +
    
    tm_layout(
      text.fontfamily = "Calibri",
      legend.outside = TRUE,
      legend.outside.position = "right"
    )
  
  
  # functions for overlays
  
  add_overlays <- function(map) {
    
    map +
      tm_shape(varying_mask) +
      tm_raster(
        col.scale = tm_scale_categorical(
          values = "grey45"
        ),
        col.legend = tm_legend(show = FALSE)
      ) +
      
      tm_shape(valid_mask) +
      tm_raster(
        col.scale = tm_scale_categorical(
          values = "grey60"
        ),
        col.legend = tm_legend(show = FALSE)
      ) +
      
      tm_shape(AOI_raster) +
      tm_raster(
        col.scale = tm_scale_categorical(
          values = "grey80"
        ),
        col.legend = tm_legend(show = FALSE)
      ) +
      
      tm_shape(AOA_mask) +
      tm_raster(
        col.scale = tm_scale_categorical(
          values = c("0" = "grey95")
        ),
        col.legend = tm_legend(show = FALSE)
      ) +
      
      tm_add_legend(
        type = "polygons",
        labels = "Varying conditions",
        fill = "grey45"
      ) +
      
      tm_add_legend(
        type = "polygons",
        labels = "Stable conditions",
        fill = "grey60"
      ) +
      
      tm_add_legend(
        type = "polygons",
        labels = "Out of AOI",
        fill = "grey80"
      ) +
      
      tm_add_legend(
        type = "polygons",
        labels = "Out of AOA",
        fill = "grey95"
      ) +
      
      tm_add_legend(
        type = "polygons",
        labels = "Out of QA",
        fill = "white"
      ) +
      
      tm_layout(
        text.fontfamily = "Calibri",
        legend.outside = TRUE,
        legend.outside.position = "right"
      )
  }
  
  
  # slope
  
  tm_slope <-
    
    tm_shape(slope) +
    tm_raster(
      col.scale = tm_scale_continuous(
        values = "brewer.br_bg",
        limits = global_slope
      ),
      col.legend = tm_legend(
        title = paste0(
          season,
          "\nSen's slope"
        ),
        reverse = TRUE,
        show_na = FALSE
      )
    )
  
  tm_slope <- add_overlays(tm_slope)
  
  
  # tau
  
  tm_tau <-
    
    tm_shape(tau) +
    tm_raster(
      col.scale = tm_scale_continuous(
        values = "viridis",
        limits = global_tau
      ),
      col.legend = tm_legend(
        title = paste0(
          season,
          "\nKendall's tau"
        ),
        reverse = TRUE,
        show_na = FALSE
      )
    )
  
  tm_tau <- add_overlays(tm_tau)
  
  
  # save
  
  files <- list(
    trend = file.path(
      outdir_trend,
      paste0("trend_", season, ".png")
    ),
    slope = file.path(
      outdir_trend,
      paste0("slope_", season, ".png")
    ),
    tau = file.path(
      outdir_trend,
      paste0("tau_", season, ".png")
    )
  )
  
  
  tmap_save(tm_trend, files$trend,
            width = width_value, height = height_value, dpi = 300)
  
  tmap_save(tm_slope, files$slope,
            width = width_value, height = height_value, dpi = 300)
  
  tmap_save(tm_tau, files$tau,
            width = width_value, height = height_value, dpi = 300)
  
  
  return(
    list(
      season = season,
      maps = list(
        trend = tm_trend,
        slope = tm_slope,
        tau = tm_tau
      ),
      files = files
    )
  )
}


trend_results <- purrr::imap(
  season_trend_list_zyp,
  ~ create_trend_maps(
    season = .y,
    r = .x,
    AOA_mask = if(.y %in% names(AOA_season_stack)) {
      AOA_season_stack[[.y]]
    } else {
      AOA_overall_mask
    },
    AOI_vect = AOI_vect,
    AOI_extend = AOI_extend,
    global_slope = global_slope,
    global_tau = global_tau,
    trend_cols = trend_cols,
    outdir_trend = outdir_trend
  )
)

maps <- purrr::map(
  trend_results,
  "maps"
) |> 
  unlist(recursive = FALSE)

trend_files <- purrr::map_chr(
  trend_results,
  ~ .x$files$trend
)

overview <- tmap_arrange(
  maps,
  ncol = 3
)

tmap_save(
  overview,
  filename = file.path(
    outdir_trend,
    "trend_overview.png"
  ),
  width = 12,
  height = 24,
  units = "in",
  dpi = 300
)



av_encode_video(
  input = trend_files,
  output = file.path(
    outdir_trend,
    "trend_animation.mp4"
  ),
  framerate = 1
)

#################################################################

#event effects for all anomalies, foehn days marked
foehn_all <- read.csv("data/processed_data/foehn_all_steinhoff_15min.csv")
event_path <- list.files("data/processed_data/all_dates_event_effects_zhang/", full.names = T)
event_names <- basename(event_path)
event_names <- substring(event_names, first = 7, last = 16)
library(stringr)

event_names <- str_extract(event_names, "\\d{4}-\\d{2}-\\d{2}")
event_list <- map(event_path, rast)
names(event_list) <- event_names
event_list <- purrr::map(event_list, ~ .x[["lyr1"]])
outdir_events <- "data/processed_data/vis_frames/Taylor_Valley_AOA_masked_zhang/event_effects_all_foehn_marked_final/"
dir.create(outdir_events)

#global symmetric scale with 0 as center

max_abs <- 0

for(i in seq_along(event_list)) {
  
  r <- event_list[[i]]
  
  rng <- terra::global(r, fun = "range", na.rm = TRUE)
  
  rmin <- rng[1,1]
  rmax <- rng[1,2]
  
  max_abs <- max(
    max_abs,
    abs(rmin),
    abs(rmax)
  )
}

global_range <- c(-max_abs, max_abs)

brks <- pretty(c(-max_abs, max_abs), n = 8)



########################################################
# as function
create_event_map <- function(
    event_name,
    r,
    AOA_mask,
    AOI_vect,
    AOI_extend,
    foehn_all,
    global_range,
    brks,
    outdir_events,
    frame_id
) {
  
  cat("Processing", event_name, "\n")
  
  
  
  AOI_raster <- rasterize(
    AOI_vect,
    r,
    field = 1,
    background = NA
  )
  
  AOI_raster <- mask(
    AOI_raster,
    AOI_extend
  )
  
  
  # mark foehn
  
  foehn_stations <- unique(
    foehn_all$station[
      foehn_all$date == event_name &
        foehn_all$foehn_day
    ]
  )
  
  if(length(foehn_stations) > 0) {
    
    foehn_text <- paste(
      "Foehn day:",
      paste(
        foehn_stations,
        collapse = ", "
      )
    )
    
  } else {
    
    foehn_text <- "No foehn"
    
  }
  
  
  title_text <- paste(
    format(
      as.Date(event_name),
      "%d.%m.%Y"
    ),
    "\n",
    foehn_text
  )
  
  
  # map
  
  tm <-
    
    tm_shape(r) +
    tm_raster(
      col.scale = tm_scale_continuous(
        values = "-vik",
        limits = global_range
      ),
      col.legend = tm_legend(
        title = expression(Delta * "VWC"),
        show_na = FALSE,
        at = brks,
        reverse = TRUE
      )
    ) +
    
    # AOI
    tm_shape(AOI_raster) +
    tm_raster(
      col.scale = tm_scale_categorical(
        values = "grey80"
      ),
      col.legend = tm_legend(
        show = FALSE
      )
    ) +
    
    # AOA
    tm_shape(AOA_mask) +
    tm_raster(
      col.scale = tm_scale_categorical(
        values = c(
          "0" = "grey95"
        )
      ),
      col.legend = tm_legend(
        show = FALSE
      )
    ) +
    
    
    tm_add_legend(
      type = "polygons",
      labels = "Out of AOI",
      fill = "grey80"
    ) +
    
    tm_add_legend(
      type = "polygons",
      labels = "Out of AOA",
      fill = "grey95"
    ) +
    
    tm_add_legend(
      type = "polygons",
      labels = "Out of QA",
      fill = "white"
    ) +
    
    tm_title(title_text) +
    
    tm_layout(
      text.fontfamily = "Calibri",
      legend.outside = TRUE,
      legend.outside.position = "right",
      frame = FALSE
    )
  
  
  
  outfile <- file.path(
    outdir_events,
    sprintf(
      "frame_%03d.png",
      frame_id
    )
  )
  
  tmap_save(
    tm,
    filename = outfile,
    width = 1800,
    height = 1400,
    dpi = 300
  )
  
  
  return(outfile)
}

frame_files <- purrr::imap_chr(
  event_list,
  function(r, event_name) {
    
    event_date <- as.Date(event_name)
    
    event_season <- paste0(
      substr(format(event_date, "%Y"), 3, 4),
      "_",
      substr(format(event_date + months(6), "%Y"), 3, 4)
    )
    
    
    AOA_mask <- if(event_season %in% names(AOA_season_stack)) {
      
      AOA_season_stack[[event_season]]
      
    } else {
      
      AOA_overall_mask
      
    }
    
    
    create_event_map(
      event_name = event_name,
      r = r,
      AOA_mask = AOA_mask,
      AOI_vect = AOI_vect,
      AOI_extend = AOI_extend,
      foehn_all = foehn_all,
      global_range = global_range,
      brks = brks,
      outdir_events = outdir_events,
      frame_id = match(event_name, names(event_list))
    )
    
  }
)

# video

av_encode_video(
  input = frame_files,
  output = file.path(
    outdir_events,
    "event_effects_2.mp4"
  ),
  framerate = 2
)

################################################################################

#for aspa 131
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

# one color scale for all prediction situations; same as for all predictions
mins <- c()
maxs <- c()

for(season in names(prediction_list)) {
  
  r <- prediction_list[[season]]$stack
  
  mm <- global(
    r,
    c("min", "max"),
    na.rm = TRUE
  )
  if (all(is.na(mm))) next
  mins <- c(mins, mm[,1])
  maxs <- c(maxs, mm[,2])
}

global_min <- min(mins, na.rm = TRUE)
global_max <- max(maxs, na.rm = TRUE)

global_min
global_max

# prepare and build frames
tmap_mode("plot")

outdir_pred_aspa <- "data/processed_data/vis_frames/aspa_131/aspa_131_AOA/predictions_aspa131/"
dir.create(outdir_pred_aspa)

frame_id <- 1

for(season in names(prediction_list_aspa)) {
  
  stack <- prediction_list_aspa[[season]]$stack
  dates <- prediction_list_aspa[[season]]$dates
  
  for(i in seq_len(nlyr(stack))) {
    
    cat("Processing", season, "-", i, "/", nlyr(stack), "\n")
    
    # AOA-Layer
    if (season %in% names(AOA_season_stack_aspa)) {
      AOA_mask <- AOA_season_stack_aspa[[season]]
    } else {
      AOA_mask <- AOA_overall_mask_aspa
    }
    
    sm <- stack[[i]]
    
    AOI_raster <- rasterize(
      AOI_vect,
      sm,
      field = 1,
      background = NA
    )
    
    AOI_raster <- mask(AOI_raster, AOI_extend)
    
    title_text <- paste0(
      "Season ",
      season,
      "\n",
      format(dates[i], "%d.%m.%Y")
    )
    
    tm <-
      
      # out of AOI grey 
      tm_shape(AOI_raster) +
      tm_raster(
        col.scale = tm_scale_categorical(
          values = "grey80"
        ),
        col.legend = tm_legend(
          show = FALSE
        )
      ) +
      tm_add_legend(
        type = "polygons",
        labels = "Out of AOI",
        fill = "grey80"
      ) +
      # prediction
      
      tm_shape(sm) +
      tm_raster(
        col.scale = tm_scale_continuous(
          values = "blues",
          limits = c(global_min, global_max)
        ),
        col_alpha = 1,
        
        # no missing in legend
        col.legend = tm_legend(
          title = "VWC",
          show_na = FALSE
        ),
        
      ) +
      # Out of AOA
      tm_shape(AOA_mask) +
      tm_raster(
        col.scale = tm_scale_categorical(
          values = c("0" = "grey95")
        ),
        col.legend = tm_legend(show = FALSE)
      ) +tm_add_legend(
        type = "polygons",
        labels = "Out of AOA",
        fill = "grey95"
      ) +
      # No data legend
      tm_add_legend(
        type = "polygons",
        labels = "No data",
        fill = "white"
      ) +
      
      # title 
      tm_title(title_text) +
      
      tm_layout(
        text.fontfamily = "Calibri",
        legend.outside = TRUE,
        frame = FALSE
      )
    
    outfile <- file.path(
      outdir_pred_aspa,
      sprintf("frame_%04d.png", frame_id)
    )
    
    tmap_save(
      tm,
      filename = outfile,
      width = 1800,
      height = 1400,
      dpi = 300
    )
    
    frame_id <- frame_id + 1
  }
}

#mp4 video

pngs <- list.files(
  outdir_pred_aspa,
  pattern = "\\.png$",
  full.names = TRUE
)


pngs <- sort(pngs)

av_encode_video(
  input = pngs,
  output = file.path(outdir_pred_aspa, "VWC_all_timesteps_aspa.mp4"),
  framerate = 2,
  vfilter = "scale=trunc(iw/2)*2:trunc(ih/2)*2"
)

################################################################################
# event list
mask_event_to_aspa <- function(x, season, aspa = aspa_131_v){
  x <- crop(x, aspa)
  x <- mask(x, aspa)
  return(x)
}

AOA_season_stack_aspa <- mask_event_to_aspa(AOA_season_stack)
AOA_overall_mask_aspa <- mask_event_to_aspa(AOA_overall_mask)
event_list_aspa <- map(event_list, mask_event_to_aspa)


names(event_list_aspa) <- event_names

outdir_events_aspa <- "data/processed_data/vis_frames/aspa_131/aspa_131_AOA_zhang/event_effects_all_foehn_marked_aspa131_final/"
dir.create(outdir_events_aspa)

#frame files aspa, same global min max as for whole valley
frame_files_aspa <- purrr::imap_chr(
  event_list_aspa,
  function(r, event_name) {
    
    event_date <- as.Date(event_name)
    
    event_season <- paste0(
      substr(format(event_date, "%Y"), 3, 4),
      "_",
      substr(format(event_date + months(6), "%Y"), 3, 4)
    )
    
    
    AOA_mask <- if(event_season %in% names(AOA_season_stack)) {
      
      AOA_season_stack_aspa[[event_season]]
      
    } else {
      
      AOA_overall_mask_aspa
      
    }
    
    
    create_event_map(
      event_name = event_name,
      r = r,
      AOA_mask = AOA_mask,
      AOI_vect = AOI_vect,
      AOI_extend = AOI_extend,
      foehn_all = foehn_all,
      global_range = global_range,
      brks = brks,
      outdir_events = outdir_events_aspa,
      frame_id = match(event_name, names(event_list))
    )
    
  }
)


# video

av_encode_video(
  input = frame_files_aspa,
  output = file.path(
    outdir_events_aspa,
    "event_effects_2.mp4"
  ),
  framerate = 2
)
##############################################################
#season means

season_mean_list_aspa <- map(season_mean_list, mask_event_to_aspa)
names(season_mean_list_aspa) <- season_mean_names

outdir_mean_aspa <- "data/processed_data/vis_frames/aspa_131/aspa_131_AOA/season_means_aspa131/"
dir.create(outdir_mean_aspa)
frame_id <- 1

for(i in seq_along(season_mean_list_aspa)) {
  
  mean <- season_mean_list_aspa[[i]]$mean
  season <- names(season_mean_list_aspa[i])
  
  #matching AOA layer
  if (season %in% names(AOA_season_stack_aspa)) {
    AOA_mask <- AOA_season_stack_aspa[[season]]
  } else {
    AOA_mask <- AOA_overall_mask_aspa
  }
  
  cat("Processing", season)
  
  AOI_raster <- rasterize(
    AOI_vect,
    mean,
    field = 1,
    background = NA
  )
  
  AOI_raster <- mask(AOI_raster, AOI_extend)
  
  title_text <- paste0(
    "Season ",
    season)
  
  
  tm <-
    
    # out of AOI grey 
    tm_shape(AOI_raster) +
    tm_raster(
      col.scale = tm_scale_categorical(
        values = "grey80"
      ),
      col.legend = tm_legend(
        show = FALSE,
      )
    ) +
    tm_add_legend(
      type = "polygons",
      labels = "Out of AOI",
      fill = "grey80"
    ) +
    # prediction
    
    tm_shape(mean) +
    tm_raster(
      col.scale = tm_scale_continuous(
        values = "blues",
        limits = c(global_min, global_max)
      ),
      col_alpha = 1,
      
      # no missing in legend
      col.legend = tm_legend(
        title = "VWC",
        show_na = FALSE
      ),
      
    ) +
    tm_shape(AOA_mask) +
    tm_raster(
      col.scale = tm_scale_categorical(
        values = c("0" = "grey95")
      ),
      col.legend = tm_legend(show = FALSE)
    ) +
    tm_add_legend(
      type = "polygons",
      labels = "Out of AOA",
      fill = "grey95"
    ) +
    # out of QA legend
    tm_add_legend(
      type = "polygons",
      labels = "Out of QA",
      fill = "white"
    ) +
    
    # title 
    tm_title(title_text) +
    
    tm_layout(
      text.fontfamily = "Calibri",
      legend.outside = TRUE,
      frame = FALSE
    )
      
  outfile <- file.path(
    outdir_mean_aspa,
    sprintf("frame_%04d.png", frame_id)
  )
  
  tmap_save(
    tm,
    filename = outfile,
    width = 1800,
    height = 1400,
    dpi = 300
  )
  
  frame_id <- frame_id + 1
  
}


#mp4 video

pngs_mean <- list.files(
  outdir_mean_aspa,
  pattern = "\\.png$",
  full.names = TRUE
)


pngs_mean <- sort(pngs_mean)

av_encode_video(
  input = pngs_mean,
  output = file.path(outdir_mean_aspa, "VWC_season_means.mp4"),
  framerate = 2,
  vfilter = "scale=trunc(iw/2)*2:trunc(ih/2)*2"
)

#################################################################
#season_mean_trends
season_trend_list_zyp_aspa <- map(season_trend_list_zyp, mask_event_to_aspa)
names(season_trend_list_zyp_aspa) <- season_trend_names
#plot(season_trend_list_zyp[[1]])
outdir_trend_aspa <- "data/processed_data/vis_frames/aspa_131/aspa_131_AOA_zhang/season_trends_aspa131/"
dir.create(outdir_trend_aspa)

# same scale as for whole valley for comparability 
trend_results_aspa <- purrr::imap(
  season_trend_list_zyp_aspa,
  ~ create_trend_maps(
    season = .y,
    r = .x,
    AOA_mask = if(.y %in% names(AOA_season_stack_aspa)) {
      AOA_season_stack_aspa[[.y]]
    } else {
      AOA_overall_mask_aspa
    },
    AOI_vect = AOI_vect,
    AOI_extend = AOI_extend,
    global_slope = global_slope,
    global_tau = global_tau,
    trend_cols = trend_cols,
    outdir_trend = outdir_trend_aspa,
    width_value = 2000,
    height_value = 1400
  )
)


trend_files_aspa <- purrr::map_chr(
  trend_results,
  ~ .x$files$trend
)

av_encode_video(
  input = trend_files_aspa,
  output = file.path(
    outdir_trend_aspa,
    "trend_animation.mp4"
  ),
  framerate = 1
)


####################################################################################
#visualisation of stable conditions- mean of season means
str(season_mean_list)
season_mean_list_r <- map(season_mean_list, ~resample(.x, season_mean_list[[1]]))
mean_layers <- map(season_mean_list_r, function(x) x$mean)
season_mean_stack <- rast(mean_layers)
season_mean_mean <- mean(season_mean_stack, na.rm =T)
plot(season_mean_mean)

season_means_trend_layer <- rast("data/processed_data/season_mean_trend_analysis_zyp_zhang/season_mean_trend_zyp_zhang.tif")
summary(season_means_trend_layer$practical_relevance)

mean_plot <- ifel(
  is.na(season_means_trend_layer$practical_relevance),
  season_mean_mean,
  NA
)

unstable <- season_means_trend_layer$valid_mask
unstable[unstable == 0] <- NA


season_mean_mean_p <- 
  
  tm_shape(mean_plot) +
  tm_raster(
    col.scale = tm_scale_continuous(
      values = "blues",
      limits = c(global_min, global_max)
    ),
    col.legend = tm_legend(
      title = "VWC",
      show_na = FALSE
    )
  ) +
  
  # Unstable conditions
  tm_shape(unstable) +
  tm_raster(
    col.scale = tm_scale_categorical(
      values = "grey60"
    ),
    col.legend = tm_legend(
      show = FALSE
    )
  ) +
  tm_add_legend(
    type = "polygons",
    labels = "Unstable conditions",
    fill = "grey60"
  ) +
  
  # Out of AOI
  tm_shape(AOI_raster) +
  tm_raster(
    col.scale = tm_scale_categorical(
      values = "grey80"
    ),
    col.legend = tm_legend(
      show = FALSE
    )
  ) +
  tm_add_legend(
    type = "polygons",
    labels = "Out of AOI",
    fill = "grey80"
  ) +
  
  # Out of AOA 
  tm_shape(AOA_overall_mask) +
  tm_raster(
    col.scale = tm_scale_categorical(
      values = c("0" = "grey95")
    ),
    col.legend = tm_legend(
      show = FALSE
    )
  ) +
  tm_add_legend(
    type = "polygons",
    labels = "Out of AOA",
    fill = "grey95"
  ) +
  
  # Out of QA
  tm_add_legend(
    type = "polygons",
    labels = "Out of QA",
    fill = "white"
  )+
  
  tm_layout(
    text.fontfamily = "Calibri")


tmap_save(
  season_mean_mean_p,
  filename = file.path(
    "data/processed_data/vis_frames/Taylor_Valley_AOA_masked/season_mean_mean.png"
  ),
  width = 1800,
  height = 1400,
  dpi = 300
)


# stable conditions for ASPA 131
season_mean_list_aspa_r <- map(season_mean_list_aspa, ~resample(.x, season_mean_list_aspa[[1]]))
mean_layers_aspa <- map(season_mean_list_aspa_r, function(x) x$mean)
season_mean_stack_aspa <- rast(mean_layers_aspa)
season_mean_mean_aspa <- mean(season_mean_stack_aspa, na.rm =T)
plot(season_mean_mean_aspa)

ov_season_mean_aspa <- mask_event_to_aspa(season_means_trend_layer)

mean_plot_aspa <- ifel(
  is.na(ov_season_mean_aspa$practical_relevance),
  season_mean_mean_aspa,
  NA
)

unstable_aspa <- ov_season_mean_aspa$valid_mask
unstable_aspa[unstable_aspa == 0] <- NA

AOI_raster <- rasterize(
  AOI_vect,
  prediction_list[1]$`15_16`$stack$S2_20151218_pred_AOA.tif,
  field = 1,
  background = NA
)

season_mean_mean_aspa_p <- 
  
  tm_shape(mean_plot_aspa) +
  tm_raster(
    col.scale = tm_scale_continuous(
      values = "blues",
      limits = c(global_min, global_max)
    ),
    col.legend = tm_legend(
      title = "VWC",
      show_na = FALSE
    )
  ) +
  
  # Unstable conditions
  tm_shape(unstable_aspa) +
  tm_raster(
    col.scale = tm_scale_categorical(
      values = "grey60"
    ),
    col.legend = tm_legend(
      show = FALSE
    )
  ) +
  tm_add_legend(
    type = "polygons",
    labels = "Unstable conditions",
    fill = "grey60"
  ) +
  
  # Out of AOI
  tm_shape(AOI_raster) +
  tm_raster(
    col.scale = tm_scale_categorical(
      values = "grey80"
    ),
    col.legend = tm_legend(
      show = FALSE
    )
  ) +
  tm_add_legend(
    type = "polygons",
    labels = "Out of AOI",
    fill = "grey80"
  ) +
  
  # Out of AOA 
  tm_shape(AOA_overall_mask_aspa) +
  tm_raster(
    col.scale = tm_scale_categorical(
      values = c("0" = "grey95")
    ),
    col.legend = tm_legend(
      show = FALSE
    )
  ) +
  tm_add_legend(
    type = "polygons",
    labels = "Out of AOA",
    fill = "grey95"
  ) +
  
  # Out of QA
  tm_add_legend(
    type = "polygons",
    labels = "Out of QA",
    fill = "white"
  )+
  tm_layout(
    text.fontfamily = "Calibri")


tmap_save(
  season_mean_mean_aspa_p,
  filename = file.path("data/processed_data/vis_frames/aspa_131/aspa_131_AOA/season_trends_aspa131/season_mean_mean_aspa131.png"),
  width = 1800,
  height = 1400,
  dpi = 300
)
#############################################################################
#overall overview map 
LTER_met_station_locations <- st_read("data/raw_data/GIS layers/MCMLTER_GIS_Export_Layers/LTER/met_station_locations.shp")
LTER_met_station_locations <- st_transform(LTER_met_station_locations, "EPSG:32758")

stations <- c("EXEM", "FRLM", "HOEM")
LTER_stations_taylor <- LTER_met_station_locations %>% filter(METLOCID %in% stations)

#own rgb
rgb_dir <- "data/raw_data/S2/15_16/S2A_MSIL2A_20160114T203522_N0500_R085_T57CWP_20231014T144817.SAFE/GRANULE/L2A_T57CWP_A002942_20160114T203520/IMG_DATA/R10m/"
b02 <- rast(paste0(rgb_dir,"T57CWP_20160114T203522_B02_10m.jp2"))
b03 <- rast(paste0(rgb_dir,"T57CWP_20160114T203522_B03_10m.jp2"))
b04 <- rast(paste0(rgb_dir,"T57CWP_20160114T203522_B04_10m.jp2"))
rgb <- c(b04,b03,b02)
names(rgb) <- c("red", "green", "blue")
AOI_rgb <- st_transform(AOI_extend, crs(rgb))
AOI_rgb <- st_buffer(AOI_rgb, 100)
rgb <- crop(rgb, AOI_rgb)

rgb_d <- rgb/1000
#plot(LTER_stations_taylor)
TMS_location <- st_read("data/processed_data/TMS_locations.gpkg")
TMS_location <- st_transform(TMS_location, crs(rgb))
qs <- global(
  rgb_d,
  quantile,
  probs = c(0.05, 0.95),
  na.rm = TRUE
)

qs
rgb_stretch <- rgb_d

for (i in 1:3) {
  
  low  <- qs$X5.[i]
  high <- qs$X95.[i]
  
  rgb_stretch[[i]] <- clamp(
    (rgb_d[[i]] - low) / (high - low),
    lower = 0,
    upper = 1
  )
}

rgb_stretch <- clamp(rgb_stretch, 0, 1)
gamma <- 0.6   # smaller = brighter

rgb_plot <- rgb_stretch ^ gamma

tm_shape(rgb_plot) +
  tm_rgb(max.value = 1)

LTER_stations_taylor$station <- factor(
  LTER_stations_taylor$METLOCID,
  levels=c("HOEM","EXEM","FRLM")
)
LTER_stations_taylor$station <- factor(
  LTER_stations_taylor$METLOCID,
  levels = c("BOYM", "EXEM", "FRLM"),
  labels = c(
    "Lake Bonney",
    "Explorer's Cove",
    "Lake Fryxell"
  )
)

LTER_stations_taylor <- st_transform(LTER_stations_taylor, crs(rgb))
aspa_131_rgb <- st_transform(aspa_131, crs(rgb))

AOI_bbox <- st_as_sfc(st_bbox(AOI_extend))



overview_map <-
  
  tm_shape(rgb_plot) +
  tm_rgb(max.value = 1) +
  
  tm_shape(AOI_rgb) +
  tm_borders(
    col="#4A7FB8",
    lwd=4
  )


#Taylor Valley overview map
#streams
streams <- st_read("data/raw_data/GIS layers/MCMLTER_GIS_Export_Layers/streams_monitored.gpkg")
mapview::mapview(streams)
streams <- st_transform(streams, st_crs(AOI_extend))
streams <- st_intersection(streams, AOI_extend)
streams <- st_transform(streams, crs(rgb))


mapview::mapview(streams2)
streams <- st_transform(streams, st_crs(AOI_extend))
streams <- st_intersection(streams, AOI_extend)
streams <- st_transform(streams, crs(rgb))


#ALMS network
alms_locations <- st_read("data/raw_data/VWC_measurements_goosef/alms_staions.gpkg")
alms_locations <- st_transform(alms_locations, crs(rgb))
mapview::mapview(alms_locations)
alms_locations


col_aoi      <- "grey90"
col_aspa     <- "lightgreen"

col_hoare    <- "#FEE08B"
col_explorer <- "#FDAE61"
col_fryxell  <- "#F46D43"

col_tms      <- "lightblue"

#streams

stream_names <- as.character(unique(streams$STREAM_NAM))

stream_cols <- colorRampPalette(
  c(
    "#08306B",
    "#08519C",
    "#2171B5",
    "#4292C6",
    "#6BAED6",
    "#9ECAE1"
  )
)(length(stream_names))

names(stream_cols) <- stream_names

#ALMS stations
alms_cols <- c(
  "ALMS1" = "#6A3D9A",
  "ALMS3" = "#88419D",
  "ALMS4" = "#9E9AC8",
  "ALMS6" = "#CBC9E2"
)

#LTER  met stations

st_hoare <- LTER_stations_taylor[
  LTER_stations_taylor$METLOCID == "HOEM", ]

st_explorer <- LTER_stations_taylor[
  LTER_stations_taylor$METLOCID == "EXEM", ]

st_fryxell <- LTER_stations_taylor[
  LTER_stations_taylor$METLOCID == "FRLM", ]


overview_map <-
  
  tm_shape(rgb_plot) +
  tm_rgb(max.value = 1)+
  
  
#AOI

tm_shape(AOI_extend) +
  tm_borders(
    col = col_aoi,
    lwd = 3
  ) +
  
#ASPA 131

tm_shape(aspa_131_rgb) +
  tm_borders(
    col = col_aspa,
    lwd = 3
  ) +
  
#TMS network

tm_shape(TMS_location) +
  tm_borders(
    col = col_tms,
    lwd = 2
  ) +
  
#lter lake hoare

tm_shape(st_hoare) +
  tm_symbols(
    col = col_hoare,
    shape = 21,
    size = 0.5,
    border.col = "black"
  ) +
  
# lter explorers cove

tm_shape(st_explorer) +
  tm_symbols(
    col = col_explorer,
    shape = 21,
    size = 0.5,
    border.col = "black"
  ) +
  
#lter lake fryxell

tm_shape(st_fryxell) +
  tm_symbols(
    col = col_fryxell,
    shape = 21,
    size = 0.5,
    border.col = "black"
  ) +
  
#legend

tm_add_legend(
  type = "lines",
  labels = "Study area",
  col = col_aoi,
  lwd = 2
) +
  
  tm_add_legend(
    type = "lines",
    labels = "ASPA 131",
    col = col_aspa,
    lwd = 2
  ) +
  
  tm_add_legend(
    type = "lines",
    labels = "TMS sensor network",
    col = col_tms,
    lwd = 2
  ) +
  
  tm_add_legend(
    type = "symbols",
    labels = "Lake Hoare LTER station",
    fill = col_hoare,
    col = "black",
    shape = 21,
    size = 0.5
  ) +
  
  tm_add_legend(
    type = "symbols",
    labels = "Explorer's Cove LTER station",
    fill = col_explorer,
    col = "black",
    shape = 21,
    size = 0.5
  ) +
  
  tm_add_legend(
    type = "symbols",
    labels = "Lake Fryxell LTER station",
    fill = col_fryxell,
    col = "black",
    shape = 21,
    size = 0.5
  ) +
  tm_credits(
    text = paste(
      "Basemap:Copernicus Sentinel data [14 Jan 2016],",
      "downloaded from Copernicus Browser"
    ),
    position = c("LEFT", "BOTTOM"),
    size = 0.3,
    bg.color = "white",
    bg.alpha = 0.7
  ) +
  
tm_compass(
  type = "arrow",
  position = c("right", "top"),
  size = 1,
  show.labels = F,
  bg.color = "lightgrey",
  bg.alpha = 0.3,
) +
  

tm_scalebar(
  position = c("left", "bottom"),
  bg.color = "lightgrey",
  bg.alpha = 0.3,
  breaks = c(0,1,2)
) +
  
 

tm_layout(
  frame = FALSE,
  legend.outside = TRUE,
  legend.position = tm_pos_out("right", "center"),
  legend.bg.color = "white",
  legend.bg.alpha = 0.95
)
overview_map
################################################

#map with streams and ALMS network in addition


col_aoi      <- "grey90"
col_aspa     <- "lightgreen"

col_hoare    <- "#FEE08B"
col_explorer <- "#FDAE61"
col_fryxell  <- "#F46D43"

col_tms      <- "#007C83"


#LTER stations

st_hoare <- LTER_stations_taylor[
  LTER_stations_taylor$METLOCID == "HOEM", ]

st_explorer <- LTER_stations_taylor[
  LTER_stations_taylor$METLOCID == "EXEM", ]

st_fryxell <- LTER_stations_taylor[
  LTER_stations_taylor$METLOCID == "FRLM", ]


#streams

stream_names <- as.character(unique(streams$STREAM_NAM))

stream_cols <- colorRampPalette(
  c(
    "#08306B",
    "#08519C",
    "#2171B5",
    "#4292C6",
    "#6BAED6",
    "#9ECAE1"
  )
)(length(stream_names))

names(stream_cols) <- stream_names


#ALMS

alms_names <- as.character(alms_locations$station)

alms_cols <- c(
  "ALMS1" = "#54278F",
  "ALMS3" = "#756BB1",
  "ALMS4" = "#9E9AC8",
  "ALMS6" = "#CBC9E2"
)


#map
overview_map <-
  
  tm_shape(rgb_plot) +
  tm_rgb(
    col.scale = tm_scale_rgb(max_color_value = 1)
  ) +
  
# AOI

tm_shape(AOI_extend) +
  tm_borders(
    col = col_aoi,
    lwd = 3
  ) +
  
#ASPA
  
tm_shape(aspa_131_rgb) +
  tm_borders(
    col = col_aspa,
    lwd = 3
  )


#streams

for (i in seq_along(stream_names)) {
  
  stream_i <- streams[
    as.character(streams$STREAM_NAM) == stream_names[i],
  ]
  
  overview_map <-
    overview_map +
    tm_shape(stream_i) +
    tm_lines(
      col = stream_cols[i],
      lwd = 2
    )
}


#TMS

overview_map <-
  overview_map +
  
  tm_shape(TMS_location) +
  tm_borders(
    col = col_tms,
    lwd = 2
  ) +
  
  


tm_shape(st_hoare) +
  tm_symbols(
    col = col_hoare,
    shape = 21,
    size = 0.5,
    border.col = "black"
  ) +
  
  

tm_shape(st_explorer) +
  tm_symbols(
    col = col_explorer,
    shape = 21,
    size = 0.5,
    border.col = "black"
  ) +
  
  


tm_shape(st_fryxell) +
  tm_symbols(
    col = col_fryxell,
    shape = 21,
    size = 0.5,
    border.col = "black"
  )




for (i in seq_along(alms_names)) {
  
  alms_i <- alms_locations[
    as.character(alms_locations$station) == alms_names[i],
  ]
  
  overview_map <-
    overview_map +
    tm_shape(alms_i) +
    tm_symbols(
      col = alms_cols[alms_names[i]],
      shape = 21,
      size = 0.35,
      border.col = "black"
    )
}




overview_map <-
  
  overview_map +
  


tm_add_legend(
  type = "lines",
  labels = "Study area",
  col = col_aoi,
  lwd = 2
) +
  
  

tm_add_legend(
  type = "lines",
  labels = "ASPA 131",
  col = col_aspa,
  lwd = 2
) +
  


tm_add_legend(
  type = "lines",
  labels = "TMS sensor network",
  col = col_tms,
  lwd = 2
) +
  


tm_add_legend(
  type = "lines",
  labels = stream_names,
  col = unname(stream_cols[stream_names]),
  lwd = 2
) +
  


tm_add_legend(
  type = "symbols",
  labels = "Lake Hoare LTER station",
  fill = col_hoare,
  col = "black",
  shape = 21,
  size = 0.5
) +
  
 

tm_add_legend(
  type = "symbols",
  labels = "Explorer's Cove LTER station",
  fill = col_explorer,
  col = "black",
  shape = 21,
  size = 0.5
) +
  
  

tm_add_legend(
  type = "symbols",
  labels = "Lake Fryxell LTER station",
  fill = col_fryxell,
  col = "black",
  shape = 21,
  size = 0.5
) +
  


tm_add_legend(
  type = "symbols",
  labels = c("ALMS 1", "ALMS 3", "ALMS 4", "ALMS 6"),
  fill = unname(alms_cols[alms_names]),
  col = "black",
  shape = 21,
  size = 0.35
) +
  
  
 # credits

tm_credits(
  text = paste(
    "Basemap:Copernicus Sentinel data [14 Jan 2016],",
    "downloaded from Copernicus Browser"
  ),
  position = c("LEFT", "BOTTOM"),
  size = 0.3,
  bg.color = "white",
  bg.alpha = 0.7
) +
  
  


tm_compass(
  type = "arrow",
  position = c("right", "top"),
  size = 1,
  show.labels = F,
  bg.color = "lightgrey",
  bg.alpha = 0.3
) +
  
  
 

tm_scalebar(
  position = c("left", "bottom"),
  bg.color = "lightgrey",
  bg.alpha = 0.3,
  breaks = c(0, 1, 2)
) +
  
  

tm_layout(
  frame = FALSE,
  fontfamily = "Calibri",
  legend.outside = TRUE,
  legend.position = tm_pos_out("right", "center"),
  legend.bg.color = "white",
  legend.bg.alpha = 0.95
)




overview_map
#################################################

tmap_save(
  overview_map,
  filename = file.path("plots/overview_map2.png"),
  width = 1800,
  height = 900,
  dpi = 300
)

antarctica <- rnaturalearth::ne_countries(
  scale="medium",
  returnclass="sf"
)

antarctica <- antarctica |>
  filter(admin=="Antarctica")

antarctica <- st_transform(antarctica,32757)
AOI_center <- st_centroid(st_union(AOI_rgb))
inset <-
  
  tm_shape(antarctica)+
  tm_polygons(
    col="grey80",
    border.col="grey30"
  )+
  
  tm_shape(AOI_center)+
  tm_symbols(
    fill="orange",
    shape=21,
    size=1.5
  )+
  
  tm_layout(frame=TRUE)

inset

tmap_save(
  inset,
  filename = file.path("plots/antarctica_position.png"),
  width = 1800,
  height = 1800,
  dpi = 300
)

# plot LPD

load_lpd <- function(season_path) {
  
  prediction_dir <- file.path(season_path, "LPD")
  
  files <- list.files(prediction_dir, full.names = TRUE)
  
  dates <- extract_date(files)
  
  keep <- !is.na(dates)
  files <- files[keep]
  dates <- dates[keep]
  
  keep_first <- !duplicated(dates)
  files <- files[keep_first]
  dates <- dates[keep_first]
  
  rasters <- map(files, rast)
  
  
  s2_idx <- which(str_detect(files, "S2"))[1]
  ps_idx <- which(str_detect(files, "PS"))[1]
  
  ref_idx <- if (!is.na(s2_idx)) {
    s2_idx
  } else if (!is.na(ps_idx)) {
    ps_idx
  } else {
    stop("No S2 or PS reference found in: ", season_path)
  }
  
  ref <- rasters[[ref_idx]]
  
  
  rasters <- map(rasters, ~{
    if (!compareGeom(.x, ref, stopOnError = FALSE)) {
      resample(.x, ref, method = "bilinear")
    } else {
      .x
    }
  })
  
  s <- rast(rasters)
  
  names(s) <- basename(files)
  
  list(
    stack = s,
    dates = dates,
    ref_type = if (!is.na(s2_idx)) "S2" else "PS"
  )
}

lpd_list <- map(season_folders, load_lpd)
names(lpd_list) <- season_folders_name

mins <- c()
maxs <- c()

for(season in names(lpd_list)) {
  
  r <- lpd_list[[season]]$stack
  
  mm <- global(
    r,
    c("min", "max"),
    na.rm = TRUE
  )
  
  mins <- c(mins, mm[,1])
  maxs <- c(maxs, mm[,2])
}

global_min <- min(mins, na.rm = TRUE)
global_max <- max(maxs, na.rm = TRUE)

global_min #0
global_max #216

#outdir <- "data/processed_data/vis_frames/predictions/"
outdir <- "data/processed_data/vis_frames/LPD/"
dir.create(outdir)


frame_id <- 1

for(season in names(lpd_list)) {
  
  stack <- lpd_list[[season]]$stack
  dates <- lpd_list[[season]]$dates
  
  for(i in seq_len(nlyr(stack))) {
    
    cat("Processing", season, "-", i, "/", nlyr(stack), "\n")
    
    
    sm <- stack[[i]]
    
    
    
    tm <-AOI_raster <- rasterize(
      AOI_vect,
      sm,
      field = 1,
      background = NA
    )
    
    AOI_raster <- mask(AOI_raster, AOI_extend)
    
    title_text <- paste0(
      "Season ",
      season,
      "\n",
      format(dates[i], "%d.%m.%Y")
    )
    
    tm <-
      
      # out of AOI grey 
      tm_shape(AOI_raster) +
      tm_raster(
        col.scale = tm_scale_categorical(
          values = "grey80"
        ),
        col.legend = tm_legend(
          show = FALSE
        )
      ) +
      tm_add_legend(
        type = "polygons",
        labels = "Out of AOI",
        fill = "grey80"
      ) +
      
      
      # lpd
      
      tm_shape(sm) +
      tm_raster(
        col.scale = tm_scale_continuous(
          values = "viridis",
          limits = c(global_min, global_max)
        ),
        col_alpha = 1,
        
        # no missing in legend
        col.legend = tm_legend(
          title = "LPD",
          show_na = FALSE
        ),
        
      ) +
      
      # No data legend
      tm_add_legend(
        type = "polygons",
        labels = "No data",
        fill = "white"
      ) +
      
      # title 
      tm_title(title_text) +
      
      tm_layout(
        text.fontfamily = "Calibri",
        legend.outside = TRUE,
        frame = FALSE
      )
    
    outfile <- file.path(
      outdir,
      sprintf("frame_%04d.png", frame_id)
    )
    
    tmap_save(
      tm,
      filename = outfile,
      width = 1800,
      height = 1400,
      dpi = 300
    )
    
    frame_id <- frame_id + 1
  }
}

#mp4 video

pngs <- list.files(
  outdir,
  pattern = "\\.png$",
  full.names = TRUE
)


pngs <- sort(pngs)

av_encode_video(
  input = pngs,
  output = file.path(outdir, "LPD_all_timesteps.mp4"),
  framerate = 2,
  vfilter = "scale=trunc(iw/2)*2:trunc(ih/2)*2"
)
rm(pngs)