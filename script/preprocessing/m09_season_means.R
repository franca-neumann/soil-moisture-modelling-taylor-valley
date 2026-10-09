# calculate season wise prediction means

library(terra)
library(purrr)
library(stringr)
library(lubridate)

season_folders <- list.files("data/processed_data/predictions/", full.names = T)
season_folders <- season_folders[-3] # exclude seasons which are not fulfilling the requirements (17|18) (19|20)
season_folders <- season_folders[-4]
season_folders

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

calc_season_mean_qa <- function(pred_obj) {
  
  pred_stack <- pred_obj$stack
  dates <- pred_obj$dates
  months <- format(dates, "%m")
  
  out <- app(pred_stack, fun = function(x) {
    
    valid <- !is.na(x)
    
    n_valid <- sum(valid)
    
    n_months <- if (n_valid > 0) {
      length(unique(months[valid]))
    } else {
      0
    }
    
    time_span <- if (n_valid > 1) {
      as.numeric(max(dates[valid]) - min(dates[valid]))
    } else {
      0
    }
    
    qa <- as.numeric(
      n_valid >= 8 &&
        n_months >= 2 &&
        time_span >= 40
    )
    
    mean_val <- if (qa == 1) {
      mean(x, na.rm = TRUE)
    } else {
      NA_real_
    }
    
    c(
      mean = mean_val,
      n_valid = n_valid,
      n_months = n_months,
      span_days = time_span,
      qa_flag = qa
    )
  })
  
  names(out) <- c(
    "mean",
    "n_valid",
    "n_months",
    "span_days",
    "qa_flag"
  )
  
  out
}


prediction_list <- map(season_folders, load_predictions)

season_means <- map(prediction_list, calc_season_mean_qa)

names(season_means) <- basename(season_folders)

for (i in seq_along(season_means)) {
  
  out_path <- file.path(
    "data/processed_data/season_mean_stacks",
    paste0(
      names(season_means)[i],
      "_inkl_foehn.tif"
    )
  )
  
  writeRaster(
    season_means[[i]],
    out_path,
    overwrite = TRUE
  )
}

