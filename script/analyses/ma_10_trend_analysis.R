# Mann- Kendall trend analysis
# pixel wise trend for seasonal means
library(terra)
library(purrr)
library(trend)
library(zyp)
library(lubridate)
library(dplyr)
library(readxl)
library(stringr)
library(tibble)
library(sf)
# load list of season means, combine mean layers and pr  project to the S2 spatial resolution

season_mean_path <- list.files("data/processed_data/season_mean_stacks/", full.names = T)
season_mean_names <- list.files("data/processed_data/season_mean_stacks/", full.names = F)
season_mean_names <- substring(season_mean_names, first = 1, last = 5)
season_mean_stack_list <- map(season_mean_path, rast)
season_mean_list_mean_only <- map(season_mean_stack_list, ~ .x[["mean"]])
names(season_mean_list_mean_only) <- season_mean_names

S2_reference <- season_mean_list_mean_only[[1]]

resample_to_S2 <- function(raster, reference = S2_reference) {
  
  if (!compareGeom(raster, reference, stopOnError = FALSE)) {
    
    raster <- resample(raster, reference)
  }
  
  return(raster)
}

season_mean_list_s2geom <- map(season_mean_list_mean_only, resample_to_S2)

season_mean_stack <- rast(season_mean_list_s2geom)
plot(season_mean_stack)


##############################################################################
#trend analysis zyp; adjusted for uneven time steps
trend_analysis_zyp <- function(
    prediction_obj,
    season_name = "unknown_season",
    alpha = 0.01,
    precision_threshold = 0.01,
    min_n = 8,
    zyp_method = "zhang" #"yuepilon" is alternative, but not recomended. 
) {
  
  library(terra)
  library(zyp)
  
  message("Calculating trend analysis for season: ", season_name)
  
  # input
  
  x_stack <- prediction_obj$stack
  dates <- as.Date(prediction_obj$dates)
  
  if (length(dates) != nlyr(x_stack)) {
    stop("Number of dates does not match number of raster layers.")
  }
  
  first_date <- min(dates)
  
  date_num <- as.numeric(
    dates - first_date
  )
  
  # SD layer
  
  sd_of_mean <- app(
    x_stack,
    function(x) sd(x, na.rm = TRUE)
  )
  
  names(sd_of_mean) <- "sd_of_mean"
  
  # valid filter
  
  valid_mask <- sd_of_mean > precision_threshold
  names(valid_mask) <- "valid_mask"
  
  x_stack_masked <- mask(
    x_stack,
    valid_mask,
    maskvalue = FALSE
  )
  

  n_obs <- app(
    x_stack_masked,
    function(x) sum(!is.na(x))
  )
  
  names(n_obs) <- "n_obs"
  
  # zyp trend analysis
  
  zyp_results <- app(
    x_stack_masked,
    fun = function(x) {
      
      valid <- !is.na(x)
      x_valid <- x[valid]
      
      if (length(x_valid) < min_n) {
        return(rep(NA_real_, 11))
      }
      
      dates_valid <- date_num[valid]
      
      duration_days <-
        max(dates_valid) - min(dates_valid)
      
      if (duration_days <= 0) {
        return(rep(NA_real_, 11))
      }
      
      z <- tryCatch(
        zyp::zyp.trend.vector(
          y = x_valid,
          x = dates_valid,
          method = zyp_method,
          conf.intervals = TRUE
        ),
        error = function(e) NULL
      )
      
      if (is.null(z)) {
        return(rep(NA_real_, 11))
      }
      
      # extract/prepare values
      
      trend <- as.numeric(z[["trend"]])
      tau <- as.numeric(z[["tau"]])
      p_val <- as.numeric(z[["sig"]])
      autocor <- as.numeric(z[["autocor"]])
      
      ci_lower <- as.numeric(z[["lbound"]])
      ci_upper <- as.numeric(z[["ubound"]])
      intercept <- as.numeric(z[["intercept"]])
      
      
      slope_total <- trend * duration_days
      
      ci_lower_total <-
        ci_lower * duration_days
      
      ci_upper_total <-
        ci_upper * duration_days
      
      # practical relevance
      practical_class <- 0
      
      if (!is.na(p_val) &&
          p_val < alpha) {
        
        if (ci_upper_total <
            -precision_threshold) {
          
          practical_class <- 1
          
        } else if (
          ci_lower_total >
          precision_threshold
        ) {
          
          practical_class <- 2
          
        } else {
          
          practical_class <- 3
        }
      }
      
      c(
        tau = tau,
        p_value = p_val,
        sen_slope = trend,
        slope_total = slope_total,
        ci_lower_total = ci_lower_total,
        ci_upper_total = ci_upper_total,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        autocor = autocor,
        practical_relevance =
          practical_class,
        intercept = intercept
      )
    }
  )
  
  names(zyp_results) <- c(
    "tau",
    "p_value",
    "sen_slope",
    "slope_total",
    "ci_lower_total",
    "ci_upper_total",
    "ci_lower",
    "ci_upper",
    "autocor",
    "practical_relevance",
    "intercept"
  )
  
  # significance filter
  
  sig_mask <-
    zyp_results[["p_value"]] < alpha
  
  names(sig_mask) <- "significant"
  
  # output
  
  out <- c(
    sd_of_mean,
    valid_mask,
    n_obs,
    zyp_results,
    sig_mask
  )
  
  return(out)
}

######
#label information
#0 = no trend
#1 = neg relevant
#2 = pos relevant
#3 = sig irrelevant


#season_mean_zyp
#assigning pseaudo dates, because its necessary for the function 
season_dates <- c(
  "2016-01-15",
  "2017-01-15",
  "2019-01-15",
  "2021-01-15",
  "2022-01-15",
  "2023-01-15",
  "2024-01-15",
  "2025-01-15",
  "2026-01-15"
)

season_mean_obj <- list(
  stack = season_mean_stack,
  dates = season_dates
)
# min_n = 8
season_mean_zyp <- trend_analysis_zyp(
  season_mean_obj,
  season_name = "mean"
)

writeRaster(season_mean_zyp, 
            "data/processed_data/season_mean_trend_analysis_zyp_zhang/season_mean_trend_zyp_zhang.tif",
            datatype = "FLT4S",
            overwrite = TRUE)
#zyp based trend analysis, out of AOA not excluded yet, but will be during visualisation and for foehn analysis
season_folders <- list.files("data/processed_data/predictions/", full.names = T)
season_folders <- season_folders[-3] # exclude seasons wich are not fulfilling the requirements (not enough observations, to short temporal coverage) (17|18) (19|20)
season_folders <- season_folders[-4]
season_folders

extract_date <- function(x) {
  as.Date(str_extract(x, "\\d{8}"), format = "%Y%m%d")
}


load_predictions_raw <- function(season_path) {
  
  prediction_dir <- file.path(season_path, "prediction")
  
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

prediction_list_raw <- map(season_folders, load_predictions_raw)
names(prediction_list_raw) <- basename(season_folders)

season_trend_list_raw <- map2(prediction_list_raw, names(prediction_list_raw), trend_analysis_zyp)
names(season_trend_list_raw) <- basename(season_folders)

##
 

for (i in seq_along(season_trend_list_raw)) {
  
  out_path <- file.path(
    "data/processed_data/season_trend_analysis_zyp_zhang_safe",
    paste0(
      names(season_trend_list_raw)[i],
      "_season_trend_analysis_zhang.tif"
    )
  )
  
  writeRaster(
    season_trend_list_raw[[i]],
    out_path,
    datatype = "FLT4S",
    overwrite = TRUE
  )
  
}

str(season_trend_list_raw)
season_trend_list_raw_path <- list.files("data/processed_data/season_trend_analysis_zyp_zhang/", full.names = T)
season_trend_list_raw <- map(season_trend_list_raw_path, rast)
trend_names <- substring(basename(season_trend_list_raw_path), first = 1, last = 5)
season_trend_list_raw <- map(season_trend_list_raw_path, rast)
names(season_trend_list_raw) <- trend_names 
##############################################################################
# additional trend metrics

calculate_trend_error <- function(
    season_obj,
    season_name,
    trend_obj,
    save_pixel_rasters = TRUE
){
  
  library(terra)
  
  obs_stack <- season_obj$stack
  dates     <- as.Date(season_obj$dates)
  
  message("processing season: ", season_name)
  # trend pixel
  
  trend_mask <-
    trend_obj$practical_relevance %in% c(1, 2)
  
 
  keep_layers <- rep(TRUE, length(dates)) #necessary, because first date is base of trend estimation, therefore trend error == 0 by definition
  keep_layers[1] <- FALSE
  
  if(sum(keep_layers) < 2){
    
    return(NULL)
    
  }
  
  dates_use <- dates[keep_layers]
  
  
  
  intercept <- obs_stack[[1]]
  slope     <- trend_obj$sen_slope
  
  date_num <- as.numeric(
    dates_use - min(dates)
  )
  
  # pixel wise metrics
  
  pixel_sum <- NULL
  pixel_sumsq <- NULL
  n_used <- 0
  
  time_mean_error <- c()
  time_median_error <- c()
  err_list <- list()
  n_used <- 0
  for(i in which(keep_layers)){
    
    pred <- intercept + slope * as.numeric(dates[i] - min(dates))
    
    err <- abs(
      obs_stack[[i]] - pred
    )
    
    err <- mask(
      err,
      trend_mask,
      maskvalue = FALSE
    )
    
    err_list[[length(err_list) + 1]] <- err
    # mean error 
    
    time_mean_error <- c(
      time_mean_error,
      terra::global(err, fun = "mean", na.rm = TRUE)[1,1]
    )
    
    time_median_error <- c(
      time_median_error,
      terra::global(
        err,
        fun = function(x, ...) median(x, na.rm = TRUE)
      )[1, 1]
    )
    
    

    
    if(is.null(pixel_sum)){
      
      pixel_sum   <- err
      pixel_sumsq <- err^2
      
    } else {
      
      pixel_sum   <- pixel_sum + err
      pixel_sumsq <- pixel_sumsq + err^2
      
    }
    
    n_used <- n_used + 1
    
    terra::tmpFiles(remove = TRUE)
    
  }
  
  
  # pixelwise mean/sd
  
  pixel_mean_error <- pixel_sum / n_used
  
  pixel_sd_error <- sqrt(
    (pixel_sumsq / n_used) -
      pixel_mean_error^2
  )
  
  pixel_median_error <- app(
    rast(err_list),
    median,
    na.rm = TRUE,
    filename = paste0(
      "data/processed_data/trend_error_stats_zyp/",
      season_name,
      "_pixel_median_error.tif"
    ),
    overwrite = TRUE
  )
  # overall metrics
  
  overall_mean_error <- mean(
    time_mean_error,
    na.rm = TRUE
  )
  
  overall_sd_error <- sd(
    time_mean_error,
    na.rm = TRUE
  )
  
  overall_median_error <- median(
    time_median_error,
    na.rm = TRUE
  )
  
  overall_median_sd_error <- sd(
    time_median_error,
    na.rm = TRUE
  )
  # output
  
  result <- list(
    
    pixel_mean_error =
      if(save_pixel_rasters) pixel_mean_error else NULL,
    
    pixel_sd_error =
      if(save_pixel_rasters) pixel_sd_error else NULL,
    
    pixel_median_error =
      if(save_pixel_rasters) pixel_median_error else NULL,
    
    time_mean_error = time_mean_error,
    
    overall_mean_error = overall_mean_error,
    
    overall_sd_error = overall_sd_error,
    
    time_median_error = time_median_error,
    
    overall_median_error = overall_median_error,
    
    overall_median_sd_error = overall_median_sd_error,
    
    n_timesteps = n_used
    
  )
  
  
  return(result)
}

#run
prediction_error_stats <- purrr::imap(
  prediction_list_raw,
  ~ calculate_trend_error(
    season_obj  = .x,
    trend_obj   = season_trend_list_raw[[.y]],
    season_name = .y
  )
)

prediction_error_stats

for(season in names(prediction_error_stats)) {
  
  x <- prediction_error_stats[[season]]
  
  writeRaster(
    x$pixel_mean_error,
    paste0(
      "data/processed_data/trend_error_stats_zyp_zhang/",
      season,
      "_pixel_mean_error.tif"
    ),
    overwrite = TRUE
  )
  
  writeRaster(
    x$pixel_median_error,
    paste0(
      "data/processed_data/trend_error_stats_zyp_zhang/",
      season,
      "_pixel_median_error.tif"
    ),
    overwrite = TRUE
  )
  
  writeRaster(
    x$pixel_sd_error,
    paste0(
      "data/processed_data/trend_error_stats_zyp_zhang/",
      season,
      "_pixel_sd_error.tif"
    ),
    overwrite = TRUE
  )
}



prediction_error_stats_meta <- purrr::map(
  prediction_error_stats,
  ~ list(
    time_mean_error = .x$time_mean_error,
    overall_mean_error = .x$overall_mean_error,
    overall_sd_error = .x$overall_sd_error,
    time_median_error = .x$time_median_error,
    overall_median_error = .x$overall_median_error,
    overall_median_sd_error = .x$overall_median_sd_error,
    n_timesteps = .x$n_timesteps
  )
)

saveRDS(
  prediction_error_stats_meta,
  "data/processed_data/trend_error_stats_zyp_zhang/prediction_error_stats_meta.rds"
)

#reload
prediction_error_stats <- readRDS(
  "data/processed_data/trend_error_stats_zyp_zhang/prediction_error_stats_meta.rds"
)


for(season in names(prediction_error_stats)) {
  
  prediction_error_stats[[season]]$pixel_mean_error <- rast(
    file.path(
      "data/processed_data/trend_error_stats_zyp_zhang",
      paste0(season, "_pixel_mean_error.tif")
    )
  )
  
  prediction_error_stats[[season]]$pixel_sd_error <- rast(
    file.path(
      "data/processed_data/trend_error_stats_zyp_zhang",
      paste0(season, "_pixel_sd_error.tif")
    )
  )
  
  prediction_error_stats[[season]]$pixel_median_error <- rast(
    file.path(
      "data/processed_data/trend_error_stats_zyp_zhang",
      paste0(season, "_pixel_median_error.tif")
    )
  )
}

prediction_error_stats

# is there a correlation between kendalls'tau and the median pixel error
tau <- values(season_trend_list_raw[["25_26"]]$tau)
err <- values(prediction_error_stats[["25_26"]]$pixel_median_error)

df <- data.frame(
  tau = tau,
  err = err
)

df <- na.omit(df)

cor.test(
  abs(df$tau),
  df$median,
  method = "spearman"
)
#not really a correlation found, for no time step. Therefore I might follow a different
#approach; one considering a filtering and using the season_mean instead, another 
# taking the pixelwise median error into consideration  

df$tau_bin <- cut(
  df$tau,
  breaks = seq(-1, 1, by = 0.1)
)

agg <- aggregate(
  median ~ tau_bin,
  data = df,
  median
)

plot(
  agg$median,
  type = "b",
  xaxt = "n"
)

axis(
  1,
  at = seq_len(nrow(agg)),
  labels = agg$tau_bin,
  las = 2
)



###################################################
#trend comparison based on calculated season means within QA
compare_trend_vs_mean <- function(
    season_obj,
    trend_obj,
    trend_median_error,
    season_mean
){
  
  library(terra)
  
  obs_stack <- season_obj$stack
  dates     <- as.Date(season_obj$dates)
  
  trend_mask <-
    trend_obj$practical_relevance %in% c(1, 2)
  
  keep_layers <- rep(TRUE, length(dates)) #as first was excluded for trend error, this needs to be done here as well
  keep_layers[1] <- FALSE
  
  if(sum(keep_layers) < 2){
    return(NULL)
  }
  
  # median mean error
  
  err_list <- vector("list", sum(keep_layers))
  
  k <- 1
  
  for(i in which(keep_layers)){
    
    err <- abs(
      obs_stack[[i]] - season_mean
    )
    
    err <- mask(
      err,
      trend_mask,
      maskvalue = FALSE
    )
    
    err_list[[k]] <- err
    k <- k + 1
  }
  
  mean_model_median_error <- app(
    rast(err_list),
    median,
    na.rm = TRUE
  )
  
  rm(err_list)
  gc()
  
  # overall trend error
  
  trend_err_all <- mask(
    trend_median_error,
    trend_mask,
    maskvalue = FALSE
  )
  
  mean_err_all <- mask(
    mean_model_median_error,
    trend_mask,
    maskvalue = FALSE
  )
  
  improvement_all <- mean_err_all - trend_err_all
  
  # output
  
  list(
    
    mean_model_median_error = mean_model_median_error,
    
    improvement_all = improvement_all,
    
    median_improvement_all =
      global(
        improvement_all,
        median,
        na.rm = TRUE
      )[1,1],
    
    mean_improvement_all =
      global(
        improvement_all,
        mean,
        na.rm = TRUE
      )[1,1],
    
    proportion_mean_better_all =
      global(
        mean_err_all < trend_err_all,
        mean,
        na.rm = TRUE
      )[1,1]
    
  )
}


comparison <- purrr::imap(
  prediction_list_raw,
  ~ compare_trend_vs_mean(
    season_obj = .x,
    trend_obj = season_trend_list_raw[[.y]],
    trend_median_error =
      prediction_error_stats[[.y]]$pixel_median_error,
    season_mean =
      season_mean_list_mean_only[[.y]]
  )
)


###################################################
# save comparison
for(season in names(comparison)) {
  
  x <- comparison[[season]]
  
  writeRaster(
    x$mean_model_median_error,
    paste0(
      "data/processed_data/comparison_zhang/",
      season,
      "_mean_model_median_error.tif"
    ),
    overwrite = TRUE
  )
  
  writeRaster(
    x$improvement,
    paste0(
      "data/processed_data/comparison_zhang/",
      season,
      "_improvement_threshold.tif"
    ),
    overwrite = TRUE
  )
  
  writeRaster(
    x$improvement_all,
    paste0(
      "data/processed_data/comparison_zhang/",
      season,
      "_improvement_all.tif"
    ),
    overwrite = TRUE
  )
}

comparison_meta <- purrr::map(
  comparison,
  ~ list(
    median_improvement =
      .x$median_improvement,
    
    mean_improvement =
      .x$mean_improvement,
    
    proportion_mean_better =
      .x$proportion_mean_better,
    
    median_improvement_all =
      .x$median_improvement_all,
    
    mean_improvement_all =
      .x$mean_improvement_all,
    
    proportion_mean_better_all =
      .x$proportion_mean_better_all
  )
)

saveRDS(
  comparison_meta,
  "data/processed_data/comparison_zhang/comparison_meta.rds"
)

#relaod comparison
comparison <- list()

for(season in names(comparison_meta)) {
  
  comparison[[season]] <- list(
    
    mean_model_median_error =
      rast(
        paste0(
          "data/processed_data/comparison/",
          season,
          "_mean_model_median_error.tif"
        )
      ),
    
    improvement =
      rast(
        paste0(
          "data/processed_data/comparison/",
          season,
          "_improvement_threshold.tif"
        )
      ),
    
    improvement_all =
      rast(
        paste0(
          "data/processed_data/comparison/",
          season,
          "_improvement_all.tif"
        )
      )
  )
  
  comparison[[season]] <- c(
    comparison[[season]],
    comparison_meta[[season]]
  )
}

#median tau and slope

season_trend_path <- list.files("data/processed_data/season_trend_analysis_zyp_zhang/", full.names = T)
season_ov_path <- list.files("data/processed_data/season_mean_trend_analysis_zyp_zhang/", full.names = T)
season_trend_path <- c(season_trend_path, season_ov_path)
season_trend_names <- list.files("data/processed_data/season_trend_analysis_zyp_zhang/", full.names = F)
season_trend_list_zyp <- map(season_trend_path, rast)

season_trend_names <- substring(season_trend_names, first = 1, last = 5)
season_ov_name <- "Season means"
season_trend_names <-c(season_trend_names, season_ov_name)
names(season_trend_list_zyp) <- season_trend_names
season_trend_list_zyp


season_means_trend_layer <- season_trend_list_zyp[["Season means"]]
plot(season_means_trend_layer)

season_means_practical_relevance <- season_means_trend_layer$practical_relevance
season_means_slope_total <- season_means_trend_layer$slope_total

seson_means_slope_total <- mask(season_means_trend_layer$slope_total, season_means_practical_relevance, maskvalue = c(1,2), inverse = T)
slope_total_values <- values(season_means_slope_total, na.rm = T)
slope_pos <- slope_total_values[slope_total_values[, "slope_total"] > 0.01]
median(slope_pos) #0.02314247
anyNA(slope_pos)#F
n_slope_pos <- length(slope_pos)
slope_neg <- slope_total_values[slope_total_values[, "slope_total"] < -0.01]
median(slope_neg) #-0.01923276
anyNA(slope_neg)#F
n_slope_neg <- length(slope_neg)

n_slope_total <- n_slope_pos + n_slope_neg

perc_s_p <- n_slope_pos/n_slope_total
perc_s_p # 0.1657959

perc_s_n <- n_slope_neg/n_slope_total
perc_s_n # 0.8342041


season_means_tau <- mask(season_means_trend_layer$tau, season_means_practical_relevance, maskvalue = c(1,2), inverse = T)
tau_values <- values(season_means_tau, na.rm = T)
tau_pos <- tau_values[tau_values[,"tau"]>0]
min(tau_pos)
median(tau_pos)
tau_pos

tau_neg <- tau_values[tau_values[,"tau"]<0]
max(tau_neg)

# VWC condition areas
#load AOA 
AOA_season_stack <- rast("data/processed_data/AOA_masks/season_wise_masks.tif")
AOA_season_stack[AOA_season_stack==1] <- NA
AOA_overall_mask <- rast("data/processed_data/AOA_masks/overall_mask.tif")
AOA_overall_mask[AOA_overall_mask==1] <- NA

glaciers <- st_read("data/raw_data/Antarctic_GIS/MDV_GLACIERS_STREAMS_POLYGONS/MDV_GLACIERS.shp")
lakes_ponds <- st_read("data/raw_data/GIS layers/MCMLTER_GIS_Export_Layers/USGS_Map_1970/Lakes_and_Poonds_epsg32758.gpkg")
lakes_ponds <- st_transform(lakes_ponds, crs = st_crs(glaciers))
AOI_mask <- st_union(glaciers, lakes_ponds)
AOI_mask <- st_transform(AOI_mask, crs("EPSG:32758"))
str(AOI_mask)

# limit AOI to study area
AOI_extend <- st_read("C:/Users/mini4/Desktop/Uni/Master/antarctica/data/processed_data/AOI_Taylor.geojson")
AOI_extend <- st_transform(AOI_extend, crs = st_crs(AOI_mask))
AOI <- st_intersection(AOI_mask, AOI_extend)

r <- season_means_trend_layer

AOA_overall_mask


r <- season_means_trend_layer

# AOI

AOI_outer <- terra::rasterize(
  terra::vect(AOI_extend),
  r[[1]],
  field = 1,
  background = NA
)

AOI_exclude <- terra::rasterize(
  terra::vect(AOI),
  r[[1]],
  field = 1,
  background = NA
)

# AOA

AOA <- AOA_overall_mask

practical <- r[["practical_relevance"]]
valid     <- r[["valid_mask"]]


geo_mask <-
  !is.na(AOI_outer) &
  is.na(AOI_exclude) &
  is.na(AOA)

# trend based classification

final_class <- ifel(
  geo_mask & valid == 0,
  1,   # Stable
  ifel(
    geo_mask & valid == 1 & practical %in% c(0, 3),
    2,   # Varying conditions
    ifel(
      geo_mask & valid == 1 & practical == 1,
      3,   # Drying trend
      ifel(
        geo_mask & valid == 1 & practical == 2,
        4,   # Wetting trend
        NA
      )
    )
  )
)

names(final_class) <- "class"

# area

cell_area <- terra::cellSize(
  r[[1]],
  unit = "m"
)

#area per class

df <- data.frame(
  class   = terra::values(final_class)[, 1],
  area_m2 = terra::values(cell_area)[, 1]
) |>
  filter(
    !is.na(class),
    !is.na(area_m2)
  )


valid_AOI_mask <-
  !is.na(AOI_outer) &
  is.na(AOI_exclude)

valid_AOI_area <- cell_area
valid_AOI_area[!valid_AOI_mask] <- NA

total_valid_AOI_area <- terra::global(
  valid_AOI_area,
  "sum",
  na.rm = TRUE
)[[1]]


df <- data.frame(
  class   = terra::values(final_class)[, 1],
  area_m2 = terra::values(cell_area)[, 1]
) |>
  filter(
    !is.na(class),
    !is.na(area_m2)
  )

area_summary <- df |>
  mutate(
    class = factor(
      class,
      levels = 1:4,
      labels = c(
        "Stable",
        "Varying conditions",
        "Drying trend",
        "Wetting trend"
      )
    )
  ) |>
  group_by(class) |>
  summarise(
    area_m2 = sum(area_m2),
    .groups = "drop"
  ) |>
  mutate(
    area_km2 = area_m2 / 1e6,
    percentage_AOI = area_m2 / total_valid_AOI_area * 100
  )

area_summary

# control 

sum(area_summary$percentage_AOI)

table(
  values(final_class),
  useNA = "ifany"
)

# slope and tau
## slope
slope_total <- season_means_trend_layer$slope_total
tau         <- season_means_trend_layer$tau


# Drying trend = final_class 3
slope_drying <- mask(
  slope_total,
  final_class,
  maskvalues = c(3),
  inverse = T
)

# Wetting trend = final_class 4
slope_wetting <- mask(
  slope_total,
  final_class,
  maskvalues = c(4),
  inverse = T
)

slope_drying_values <- values(
  slope_drying,
  na.rm = TRUE
)[, "slope_total"]

slope_wetting_values <- values(
  slope_wetting,
  na.rm = TRUE
)[, "slope_total"]


slope_stats <- tibble(
  class = c(
    "Drying trend",
    "Wetting trend"
  ),
  n = c(
    length(slope_drying_values),
    length(slope_wetting_values)
  ),
  median_slope = c(
    median(slope_drying_values),
    median(slope_wetting_values)
  ),
  min_slope = c(
    min(slope_drying_values),
    min(slope_wetting_values)
  ),
  max_slope = c(
    max(slope_drying_values),
    max(slope_wetting_values)
  )
)

slope_stats

##tau
# Drying trend = final_class 3
tau_drying <- mask(
  tau,
  final_class,
  maskvalues = c(3),
  inverse = T
)

# Wetting trend = final_class 4
tau_wetting <- mask(
  tau,
  final_class,
  maskvalues = c(4),
  inverse = T
)

tau_drying_values <- values(
  tau_drying,
  na.rm = TRUE
)[, "tau"]

tau_wetting_values <- values(
  tau_wetting,
  na.rm = TRUE
)[, "tau"]


tau_stats <- tibble(
  class = c(
    "Drying trend",
    "Wetting trend"
  ),
  n = c(
    length(tau_drying_values),
    length(tau_wetting_values)
  ),
  median_tau = c(
    median(tau_drying_values),
    median(tau_wetting_values)
  ),
  min_tau = c(
    min(tau_drying_values),
    min(tau_wetting_values)
  ),
  max_tau = c(
    max(tau_drying_values),
    max(tau_wetting_values)
  )
)

tau_stats
