library(terra)
library(CAST)
library(sf)
library(dplyr)
# The aim is to generate an AOI mask, that includes all out of AOI pixels
# aka all pixels that have been out of AOI at least once

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

load_AOA <- function(season_path) {
  
  prediction_dir <- file.path(season_path, "AOA")
  
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

AOA_list <- map(season_folders, load_AOA)
names(AOA_list) <- season_folders_name

ref <- AOA_list[[1]]$stack[[1]]
season_masks <- vector("list", length(AOA_list))

for(i in seq_along(AOA_list)) {
  
  x <- AOA_list[[i]]$stack
  
  # resample to S2 resolution, assigning the majority value if in question
  if (!compareGeom(x, ref, stopOnError = FALSE)) {
    x <- resample(x, ref, method = "mode")
  }
  
  # AOI handling
  all_na <- app(is.na(x), fun = all)
  
  # out of AOA if the pixel was out of AOÍ at least once during the season
  zero_any <- app(x == 0, fun = any, na.rm = TRUE)
  
  # season mask
  season_masks[[i]] <- ifel(
    all_na,
    NA,
    ifel(zero_any, 0, 1)
  )
}


season_stack <- rast(season_masks)
plot(season_stack)
names(season_stack)<- season_folders_name
writeRaster(season_stack, "data/processed_data/AOA_masks/season_wise_masks.tif")

all_na <- app(is.na(season_stack), fun = all)
zero_any <- app(season_stack == 0, fun = any, na.rm = TRUE)

overall_mask <- ifel(
  all_na,
  NA,
  ifel(zero_any, 0, 1)
)

plot(overall_mask)
writeRaster(overall_mask, "data/processed_data/AOA_masks/overall_mask.tif")
