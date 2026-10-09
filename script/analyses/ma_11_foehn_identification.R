# write foehn identification algorithm
library(lubridate)
library(purrr)
library(dplyr)
library(stringr)
library(tidyr)
library(zyp)
library(sf)
library(terra)
#met_data_frlm <- read.csv("D:/antarctica/data/raw_data/met_data_LTER/mcmlter-clim-frlm_hourly-2026-02-28.csv", header = T)
#hourly
#met_data_path <- list.files("D:/antarctica/data/raw_data/met_data_LTER/", pattern = "hourly", full.names = T)
#met_stations <- gsub("mcmlter-clim-|_hourly.*", "", basename(met_data_path))
# 15 min
met_data_path <- list.files("D:/antarctica/data/raw_data/met_data_LTER/", pattern = "15min", full.names = T)
met_stations <- gsub("mcmlter-clim-|_15min.*", "", basename(met_data_path))
#met_data_frlm$date <- as.Date(met_data_frlm$date_time)
#names(met_data_frlm)

##########################################
#prepare met_data

filter_season <- function(met_data) {
  met_data %>%
    mutate(date = as.Date(date_time),
           year = year(date),
           month = month(date),
           season_year = if_else(month >= 10, year + 1L, year)) %>%
    filter(
      season_year >= 2015,
      month %in% c(11, 12, 1, 2)
    ) %>%
    dplyr::select(-year)
}

get_tage_pos_5 <- function(met_data) {
  unique(
    met_data$date[
      met_data$soiltemp1_5cm_degc > 0
    ]
  )
}

#met_ref <- filter_season(met_data)
#tage_pos_5 <- get_tage_pos_5(met_ref)

filter_positive_days <- function(met_data, tage_pos) {
  met_data %>%
    filter(date %in% tage_pos)
}

#met_filtered <- filter_positive_days(
#  met_ref,
#  tage_pos_5
#)
#############################################

#foehn identification for 15 min data as input (as suggested by speirs/ steihoff)
foehn_identification_steinhoff <- function(met_data){
  
  met_data <- filter_season(met_data)
  
  met_data %>%
    
    mutate(
      date_time = ymd_hms(date_time),
      date = as.Date(date_time)
    ) %>%
    
    arrange(date_time) %>%
    
    mutate(
      
      # change per hour (4x15min)
      
      d_temp =
        airtemp_3m_degc -
        lag(airtemp_3m_degc, 4),
      
      d_rh =
        rhh2o_3m_pct -
        lag(rhh2o_3m_pct, 4),
      
      # air temperature
      
      temp_trigger =
        d_temp >= 1 |
        airtemp_3m_degc > 0,
      
      # relative humidity
      
      rh_trigger =
        d_rh <= -5 |
        rhh2o_3m_pct < 30,
      
      # wind speed and direction
      
      wind_ok =
        wspdmax_ms > 5,
      
      dir_ok =
        between(wdir_deg, 180, 315)
    ) %>%
    
    mutate(
      
      # 6h window 
      
      temp_6h =
        zoo::rollapply(
          temp_trigger,
          width = 24,
          FUN = any,
          fill = FALSE,
          align = "right"
        ),
      
      rh_6h =
        zoo::rollapply(
          rh_trigger,
          width = 24,
          FUN = any,
          fill = FALSE,
          align = "right"
        ),
      
      wind_6h =
        zoo::rollapply(
          wind_ok,
          width = 24,
          FUN = function(x)
            sum(x, na.rm = TRUE) >= 20,
          fill = FALSE,
          align = "right"
        ),
      
      dir_6h =
        zoo::rollapply(
          dir_ok,
          width = 24,
          FUN = function(x)
            sum(x, na.rm = TRUE) >= 20,
          fill = FALSE,
          align = "right"
        )
    ) %>%
    
    mutate(
      foehn =
        temp_6h &
        rh_6h &
        wind_6h &
        dir_6h
    ) %>%
    
    group_by(date) %>%
    
    mutate(
      foehn_day = any(foehn, na.rm = TRUE)
    ) %>%
    
    ungroup()
}

###################################################################
#foehn_frlm <- foehn_identification_spiers(met_data_frlm)
#unique(foehn_frlm[foehn_frlm$foehn_day == T,]$date)
# for all met stations
met_data_list <- map(met_data_path, ~read.csv(.x, header = T))
names(met_data_list) <- met_stations

foehn_data_list <- map(met_data_list, foehn_identification_steinhoff)
foehn_data_list <- map2(
  foehn_data_list,
  met_stations,
  \(df, st) {
    df$station <- st
    df
  }
)
foehn_all <- bind_rows(foehn_data_list)
write.csv(foehn_all, "data/processed_data/foehn_all_steinhoff_15min.csv")
foehn_all <- read.csv("data/processed_data/foehn_all_steinhoff_15min.csv")
#write_xlsx(foehn_all, "data/processed_data/foehn_all_steinhoff_15min.xslx")
#foehn_all <- read_excel("data/processed_data/foehn_all_steinhoff_15min.xslx")
#test my foehn identification approach with Martes data
#foehn_frlm_marte <- read.csv("data/raw_data/foehn_dataset_Marte/fryxell_foehn_results.csv", header = T)
#intersect(as.character(foehn_all[foehn_all$station == "frlm" & foehn_all$foehn == T,]$date_time),
# as.character(foehn_frlm_marte[foehn_frlm_marte$foehn_Speirs == 1,]$date_time))

#ft_marte <- foehn_frlm_marte %>%
# filter(foehn_Speirs == 1,
#       year(as.Date(date_time)) >= 2015)

#ft_marte <- filter_season(ft_marte)

#setdiff(as.character(unique(date(as.Date(ft_marte$date_time)))),as.character(unique(date(as.Date(foehn_all[foehn_all$station == "frlm" & foehn_all$foehn, ]$date_time)))))
#setdiff(as.character(unique(date(as.Date(foehn_all[foehn_all$station == "frlm" & foehn_all$foehn, ]$date_time)))),as.character(unique(date(as.Date(ft_marte$date_time)))))
#intersect(as.character(unique(date(as.Date(foehn_all[foehn_all$station == "frlm" & foehn_all$foehn, ]$date_time)))),as.character(unique(date(as.Date(ft_marte$date_time)))))
#dates where I have available data
season_folders <- list.files("data/processed_data/predictions/", full.names = T)
season_folders <- season_folders[-3] # exclude seasons wich are not fulfilling the requirements (17|18) (19|20)
season_folders <- season_folders[-4]

extract_date <- function(x) {
  as.Date(str_extract(x, "\\d{8}"), format = "%Y%m%d")
}

prediction_dir <- file.path(season_folders, "prediction")
files <- list.files(prediction_dir, full.names = TRUE)
dates <- extract_date(files)
dates <- as.character(dates)

# for frlm station
#foehn_dates <- unique(
# as.character(
#  foehn_frlm$date[foehn_frlm$foehn_day]
#))
#intersect(foehn_dates, dates)

#for all stations
foehn_dates_all <- unique(
  as.character(
    foehn_all$date[foehn_all$foehn_day]
  )
)

intersect(foehn_dates_all, dates)

foehn_dates_stations <- foehn_all |>
  dplyr::filter(foehn_day) |>
  dplyr::distinct(
    date,
    station
  )
write_xlsx(foehn_dates_stations, "data/processed_data/foehn_dates_stations_steinhoff_15min.xlsx")
library(readxl)
foehn_dates_stations <- read_xlsx("data/processed_data/foehn_dates_stations_steinhoff_15min.xlsx")  
foehn_dates_rel <- foehn_dates_stations %>%
  filter(station %in% c("frlm", "hoem", "exem"))
foehn_dates_all <- unique(foehn_dates_rel$date)
available_events <- as.Date(intersect(as.Date(all_dates), as.Date(foehn_dates_all)),origin = "1970-01-01")
all_dates <- prediction_list |>
  purrr::map("dates") |>
  unlist(use.names = FALSE) |>
  as.Date(origin = "1970-01-01")
#############################################################
#foehn analysis

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

###################################################
#event effects based on difference to expected mean in cases where we have a significant 
# within season trend
trend_path <- list.files("data/processed_data/season_trend_analysis_zyp_zhang/", full.names = T)
#trend_path <- trend_path[-10] #exclude season_mean_calc
trend_names <- substring(basename(trend_path), first = 1, last = 5)
trend_list <- map(trend_path, rast)
names(trend_list) <- trend_names

# attach improvemet_all to trend_list
for(season in names(trend_list)) {
  
  trend_list[[season]]$improvement_all <- rast(
    paste0(
      "data/processed_data/comparison_zhang/",
      season,
      "_improvement_all.tif"
    )
  )
  
}
# attach pixel_median error to trend_list
for(season in names(trend_list)) {
  
  trend_list[[season]]$pixel_median_error <- rast(
    paste0(
      "data/processed_data/trend_error_stats_zyp_zhang/",
      season,
      "_pixel_median_error.tif"
    )
  )
  
}

#load season mean stacks

season_mean_path <- list.files("data/processed_data/season_mean_stacks/", full.names = T)
season_mean_names <- list.files("data/processed_data/season_mean_stacks/", full.names = F)
season_mean_names <- substring(season_mean_names, first = 1, last = 5)
season_mean_list <- map(season_mean_path, rast)
names(season_mean_list) <- season_mean_names

#########################################
#write directly version
calc_event_effect_trend_dependent_write <- function(
    event_date,
    prediction_list,
    trend_list,
    foehn_dates
) {
  
  library(terra)
  library(purrr)
  
  event_date <- as.Date(event_date)
  
  # identify season
  
  season_name <- names(prediction_list)[
    map_lgl(
      prediction_list,
      ~ event_date %in% .x$dates
    )
  ]
  
  if (length(season_name) == 0) {
    return(NULL)
  }
  
  pred_obj <- prediction_list[[season_name]]
  
  stack <- pred_obj$stack
  dates <- as.Date(pred_obj$dates)
  
  trend_rast <- trend_list[[season_name]]
  
  # occurring event?
  
  event_idx <- which(dates == event_date)
  
  if (length(event_idx) != 1) {
    return(NULL)
  }
  
  
  keep_idx <- which(
    !(dates %in% as.Date(foehn_dates)) &
      dates != event_date
  )
  
  if (length(keep_idx) < 2) {
    return(NULL)
  }
  
  season_mean <- app(
    stack[[keep_idx]],
    mean,
    na.rm = TRUE
  )
  
  season_sd <- app(
    stack[[keep_idx]],
    sd,
    na.rm = TRUE
  )
  
  pixel_median_error <- trend_rast$pixel_median_error #referred to as trend error in text
  
  event_layer <- stack[[event_idx]]
  
  
  
  first_date <- min(dates)
  
  days_since_start <-
    as.numeric(event_date - first_date)
  
  trend_reference <-
    stack[[1]] +
    trend_rast[["sen_slope"]] *
    days_since_start
  
  # only use trend if practical relevant (class 1,2)
  
  trend_mask <-((
    trend_rast[["practical_relevance"]] == 1 |
      trend_rast[["practical_relevance"]] == 2)
    &
      (trend_rast$improvement_all > 0.01))
  
  
  reference <- ifel(
    trend_mask,
    trend_reference,
    season_mean
  )
  
  
  event_effect <- event_layer - reference
  
 #relevance thresholds
  
  event_effect <- ifel(
    
    trend_mask,
    
    ifel(
      abs(event_effect) < 0.01 |
        abs(event_effect) < season_sd |
        abs(event_effect) < pixel_median_error,
      0,
      event_effect
    ),
    
    ifel(
      abs(event_effect) < 0.01 |
        abs(event_effect) < season_sd,
      0,
      event_effect
    )
    
  )
  
  
  event_name <- paste0(season_name,"_", as.character(event_date), "_event_effect.tif")
  event_path <- paste0("data/processed_data/all_dates_event_effects_zhang/", event_name)
  writeRaster(event_effect, event_path )
  gc()
}

prediction_dates <- sort(unique(as.Date(dates)))

walk(
  prediction_dates,
  ~ calc_event_effect_trend_dependent_write(
    event_date = .x,
    prediction_list = prediction_list,
    trend_list = trend_list,
    foehn_dates = foehn_dates_all
  )
)

######################################################################
#improved version
#write directly version
calc_event_effect_trend_dependent_write_improved <- function(
    event_date,
    prediction_list,
    trend_list,
    season_mean_list,
    foehn_dates
) {
  
  library(terra)
  library(purrr)
  
  event_date <- as.Date(event_date)
  
  # identify season
  
  season_name <- names(prediction_list)[
    map_lgl(
      prediction_list,
      ~ event_date %in% .x$dates
    )
  ]
  
  if (length(season_name) == 0) {
    return(NULL)
  }
  
  pred_obj <- prediction_list[[season_name]]
  
  stack <- pred_obj$stack
  dates <- as.Date(pred_obj$dates)
  
  trend_rast <- trend_list[[season_name]]
  season_mean <- season_mean_list[[season_name]]
  # occurring event?
  
  event_idx <- which(dates == event_date)
  
  if (length(event_idx) != 1) {
    return(NULL)
  }
  
  
  #check if there are enough non-foehn days to perform a meaningful analysis
  keep_idx <- which(
    !(dates %in% as.Date(foehn_dates)) &
      dates != event_date
  )
  
  if (length(keep_idx) < 2) {
    return(NULL)
  }
  
  season_sd <- app(
    stack,
    sd,
    na.rm = TRUE
)
  
  pixel_median_error <- trend_rast$pixel_median_error #refered to as trend error in text
 
  
  event_layer <- stack[[event_idx]]
  
  # trend based VWC
  
  first_date <- min(dates)
  
  days_since_start <-
    as.numeric(event_date - first_date)
  
  trend_reference <-
    stack[[1]] +
    trend_rast[["sen_slope"]] *
    days_since_start
  
  # use trend only when pratical relevant (class 1,2)
  
  trend_mask <-((
    trend_rast[["practical_relevance"]] == 1 |
      trend_rast[["practical_relevance"]] == 2)
    &
      (trend_rast$improvement_all > 0.01))
  
 
  
  reference <- ifel(
    trend_mask,
    trend_reference,
    season_mean
  )
  
  # VWC anomaly (= event effect)
  
  event_effect <- event_layer - reference
  
  # relevance threshold
  
  event_effect <- ifel(
    
    trend_mask,
    
    ifel(
      abs(event_effect) < 0.01 |
        abs(event_effect) < season_sd |
        abs(event_effect) < pixel_median_error,
      0,
      event_effect
    ),
    
    ifel(
      abs(event_effect) < 0.01 |
        abs(event_effect) < season_sd,
      0,
      event_effect
    )
    
  )
  
  
  event_name <- paste0(season_name,"_", as.character(event_date), "_event_effect.tif")
  event_path <- paste0("data/processed_data/all_dates_event_effects_zhang/", event_name)
  writeRaster(event_effect, event_path )
  gc()
}

prediction_dates <- sort(unique(as.Date(dates)))

# check ratio of foehn dates to non-foehn dates
library(dplyr)
library(lubridate)


foehn_dates <- as.Date(foehn_dates_matching)


all_dates <- as.Date(prediction_dates)

season_df <- tibble(
  date = all_dates
) %>%
  mutate(
    year = year(date),
    season = if_else(month(date) >= 11, year(date), year(date) - 1)
  )


season_df <- season_df %>%
  mutate(
    foehn = date %in% foehn_dates
  )

# foehn ratio per season
season_ratio <- season_df %>%
  group_by(season) %>%
  summarise(
    n_prediction_days = n(),
    n_foehn_days = sum(foehn),
    n_non_foehn_days = sum(!foehn),
    foehn_ratio = n_foehn_days / n_prediction_days,
    non_foehn_ratio = n_non_foehn_days / n_prediction_days,
    .groups = "drop"
  )

season_ratio

dir.create("D:/antarctica/terra_tmp", showWarnings = T)

terra::terraOptions(
  tempdir = "D:/antarctica/terra_tmp"
)

terra::terraOptions()

walk(
  prediction_dates,
  ~ calc_event_effect_trend_dependent_write_improved(
    event_date = .x,
    prediction_list = prediction_list,
    trend_list = trend_list,
    season_mean_list = season_mean_list,
    foehn_dates = foehn_dates_all
  )
)


walk(
  prediction_dates[81:length(prediction_dates)],
  ~ calc_event_effect_trend_dependent_write_improved(
    event_date = .x,
    prediction_list = prediction_list,
    trend_list = trend_list,
    season_mean_list = season_mean_list,
    foehn_dates = foehn_dates_all
  )
)


