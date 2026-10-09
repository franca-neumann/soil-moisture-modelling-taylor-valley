# event effect statistics

library(terra)
library(lubridate)
library(emmeans)
library(glmmTMB)
library(fitdistrplus)
library(caret)
library(purrr)
library(dplyr)
library(tidyr)
library(sf)
library(stringr)
library(DHARMa)



all_event_path <- list.files("data/processed_data/all_dates_event_effects_zhang/", full.names = T)
dates <- substring(basename(all_event_path), first =  7, last =  16)
dates <- as.character(dates)
all_event_list <- map(all_event_path, rast)
names(all_event_list) <- dates
all_event_list <- purrr::map(all_event_list, ~ .x[["lyr1"]])


foehn_all <- read.csv("data/processed_data/foehn_all_steinhoff_15min.csv")    

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

prediction_dates <- sort(unique(as.Date(dates)))

##################################################
#frames for analysis on McMLTER weatherstation data, this analysis was not used for the final thesis
# daily mean VWC

daily_df <- foehn_all %>%
  
  mutate(
    date_time = parse_date_time(
      date_time,
      orders = c("ymd HMS", "ymd HM", "ymd")
    ),
    date = as.Date(date_time)
  ) %>%
  
  filter(
    !is.na(date),
    !is.na(soilvwc_5cm_m3m3)
  ) %>%
  
 # filter to dates with available predictions
  
  group_by(
    station,
    season_year,
    date
  ) %>%
  
  summarise(
    vwc = mean(soilvwc_5cm_m3m3, na.rm = TRUE),
    foehn = any(foehn_day, na.rm = TRUE),
    .groups = "drop"
  )

# season wise trend

event_df <- daily_df %>%
  
  group_by(
    station,
    season_year
  ) %>%
  
  group_modify(~{
    
    dat <- arrange(.x, date)
    
    if(nrow(dat) < 8)
      return(tibble())
    
    
    season_mean <- mean(
      dat$vwc[!dat$foehn],
      na.rm = TRUE
    )
    
    season_sd <- sd(
      dat$vwc[!dat$foehn],
      na.rm = TRUE
    )
    
    
    first_date <- min(dat$date)
    
    dat$date_num <- as.numeric(
      dat$date - first_date
    )
    
    
    z <- tryCatch(
      zyp.trend.vector(
        y = dat$vwc,
        x = dat$date_num,
        method = "zhang",
        conf.intervals = TRUE
      ),
      error = function(e) NULL
    )
    
    if(is.null(z)) {
      
      dat$reference <- season_mean
      
    } else {
      
      trend <- unname(z[["trend"]])
      intercept <- unname(z[["intercept"]])
      p_value <- unname(z[["sig"]])
      
      # prectical relevance
      
      trend_total <-
        abs(trend) *
        max(dat$date_num)
      
      relevant <- trend_total > season_sd
      
      first_value <- dat$vwc[1]
      
      trend_pred <-
        first_value +
        trend * dat$date_num
      
      trend_error <- median(
        abs(
          dat$vwc[!dat$foehn] -
            trend_pred[!dat$foehn]
        ),
        na.rm = TRUE
      )
      
      mean_error <- median(
        abs(
          dat$vwc[!dat$foehn] -
            season_mean
        ),
        na.rm = TRUE
      )
      
      improvement_all <-
        trend_error - mean_error
      
      use_trend <-
        !is.na(p_value) &&
        p_value < 0.05 &&
        relevant&&
        improvement_all < 0.01
      
      if(use_trend){
        
        dat$reference <- trend_pred
        
        
      } else {
        
        dat$reference <- season_mean
        
      }
      
      dat$trend <- trend
      dat$intercept <- intercept
      dat$p_value <- p_value
      dat$trend_relevant <- relevant
      dat$use_trend <- use_trend
    }
    
    # VWC anomaly
    
    dat$event_effect <-
      dat$vwc  - dat$reference
    
    dat$event_effect[
      abs(dat$event_effect) < 0.01
    ] <- 0
    
    dat$event_effect[
      abs(dat$event_effect) < season_sd
    ] <- 0
    
    dat$season_mean <- season_mean
    dat$season_sd <- season_sd
    
    dat
    
  }) %>%
  
  ungroup()

###############################################################################
#mod on real (McMLTER network) data (LTER met stations)

hist(event_df$vwc)
descdist(event_df$vwc)

# beta distribution

event_df$season_factor <- as.factor(event_df$season_year)
event_df$foehn_factor <- as.factor(event_df$foehn)
event_df$vwc_rd <- round(event_df$vwc, 2)

event_df$doy <- lubridate::yday(event_df$date)

event_df <- event_df %>%
  arrange(station, date) %>%
  group_by(station) %>%
  mutate(
    vwc_p1 = lead(vwc, n = 1)
  ) %>%
  ungroup()
event_df$month <- month(event_df$date)
event_df$month <- as.factor(event_df$month)

mod_vwc <- glmmTMB(vwc ~ foehn*season_factor + foehn*month +station ,
                   data = event_df, 
                   dispformula = ~ foehn*season_factor,
                   family = beta_family(link = "logit"))
summary(mod_vwc)

mod_vwc_p1 <- glmmTMB(vwc_p1 ~ foehn*season_factor + station ,
                      data = event_df, 
                      dispformula = ~ foehn*season_factor,
                      family = beta_family(link = "logit"))
summary(mod_vwc_p1)

car::Anova(mod_vwc, type = "III")

car::Anova(mod_vwc_p1, type = "III")

emm <- emmeans(
  mod_vwc,
  ~ foehn | season_factor,
  type = "response"
)

pairs(emm)
summary(emm)

contr <- pairs(
  emm,
  reverse = TRUE
)

#
emm_p1 <- emmeans(
  mod_vwc_p1,
  ~ foehn | season_factor,
  type = "response"
)

pairs(emm_p1)
summary(emm_p1)

contr <- pairs(
  emm_p1,
  reverse = TRUE
)
#
contr_df <- as.data.frame(
  summary(
    confint(
      pairs(emm, reverse = TRUE)
    )
  )
)

str(contr_df)
cld_df <- multcomp::cld(emm, Letters = letters, adjust = "sidak")


ggplot(cld_df, aes(x = season_factor, y = response, color = foehn, group = foehn)) +
  
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    position = position_dodge(width = 0.3),
    width = 0.1
  ) +
  
  geom_text(
    aes(label = .group),
    position = position_dodge(width = 0.3),
    vjust = -1
  ) +
  
  theme_bw() +
  labs(
    x = "Season",
    y = "Estimated mean VWC",
    color = "Foehn"
  )


simres <- simulateResiduals(mod_vwc)
testDispersion(simres)
testUniformity(simres)

plotResiduals(simres,
              fitted(mod_vwc))
testDispersion(simres, alternative = "greater")

#event effect model
event_df$event_occurs <- event_df$event_effect > 0
event_df$event_occurs <- as.factor(event_df$event_occurs)

mod_occ <- glmmTMB(
  event_occurs ~ foehn * season_factor + foehn*month + station,
  data = event_df,
  family = binomial()
)

summary(mod_occ)
car::Anova(mod_occ, type = "III")
simresocc <- simulateResiduals(mod_occ)
testUniformity(simresocc)
testDispersion(simresocc)

emm_station <- emmeans(
  mod_occ,
  ~ station,
  type = "response"
)

summary(emm_station)
pairs(emm_station, reverse = T)
contrast(pairs(emm_station))

#
event_size_df <- event_df[event_df$event_effect != 0,]
mod_size <- glmmTMB(
  event_effect ~ foehn * season_factor + station,
  data = event_size_df,
  family = gaussian()
)
summary(mod_size)
car::Anova(mod_size, type = "III")
simressize <- simulateResiduals(mod_size)
testUniformity(simressize)
testDispersion(simressize)

###############################################################################

# models on prediction data

season_mean_path <- list.files("data/processed_data/season_mean_stacks/", full.names = T)
season_names <- substring(basename(season_mean_path), first = 1, last = 5)

season_mean_list <- map(season_mean_path, rast)
names(season_mean_list) <- season_names

date_lookup <- imap_dfr(
  prediction_list,
  ~ tibble(
    date = .x$dates,
    season = .y
  )
)

#random sample
# sample with pixel id

#mask to AOA

mask_event_to_aoa <- function(x,  aoa = AOA_overall_mask){
  if (!compareGeom(x, aoa, stopOnError = FALSE)) {
    x <- resample(x, aoa, method = "average")
  }
  x <- mask(x, aoa, maskvalue = 0)
  return(x)
}

AOA_season_stack <- rast("data/processed_data/AOA_masks/season_wise_masks.tif")
AOA_season_stack[AOA_season_stack==1] <- NA
AOA_overall_mask <- rast("data/processed_data/AOA_masks/overall_mask.tif")
AOA_overall_mask[AOA_overall_mask==1] <- NA

all_event_list <- map(all_event_list, mask_event_to_aoa)
set.seed(42)

sample_pts <- spatSample(
  all_event_list[[1]],
  size = 5000,
  method = "random",
  na.rm = TRUE,
  as.points = TRUE
)

# pixel-id
sample_pts$pixel_id <- seq_len(nrow(sample_pts))
#sample_pts_sf <- st_as_sf(sample_pts)
#st_write(sample_pts_sf, "data/processed_data/sample_pts.gpkg")
# extract values

pixel_df <- purrr::imap_dfr(
  all_event_list,
  \(r, date_name){
    
    ## anomalies
    event_vals <- terra::extract(
      r,
      sample_pts,
      ID = FALSE
    )[[1]]
    
    ## add season
    season <- names(prediction_list)[
      purrr::map_lgl(
        prediction_list,
        \(x) date_name %in% x$dates
      )
    ]
    
    n
    layer <- match(
      date_name,
      prediction_list[[season]]$dates
    )
    
    
    vwc_vals <- terra::extract(
      prediction_list[[season]]$stack[[layer]],
      sample_pts,
      ID = FALSE
    )[[1]]
    
    tibble(
      date = as.Date(date_name),
      pixel_id = sample_pts$pixel_id,
      event_effect = event_vals,
      vwc = vwc_vals
    )
  }
)



pixel_df <- pixel_df |>
  mutate(
    year = lubridate::year(date),
    month = lubridate::month(date),
    
    season_year = case_when(
      month %in% c(11, 12) ~ year + 1,
      month %in% c(1, 2) ~ year
    ),
    
    season_factor = factor(season_year),
    
    season = case_when(
      season_year == 2016 ~ "15_16",
      season_year == 2017 ~ "16_17",
      season_year == 2019 ~ "18_19",
      season_year == 2021 ~ "20_21",
      season_year == 2022 ~ "21_22",
      season_year == 2023 ~ "22_23",
      season_year == 2024 ~ "23_24",
      season_year == 2025 ~ "24_25",
      season_year == 2026 ~ "25_26"
    )
  )


# season quartiles

vwc_class_list <- purrr::map(
  season_mean_list,
  \(r){
    
    q <- quantile(
      values(r$mean),
      probs = c(0.25, 0.5, 0.75),
      na.rm = TRUE
    )
    
    classify(
      r$mean,
      matrix(
        c(
          -Inf, q[1], 1,
          q[1], q[2], 2,
          q[2], q[3], 3,
          q[3], Inf, 4
        ),
        ncol = 3,
        byrow = TRUE
      )
    )
  }
)

names(vwc_class_list) <- season_names

# global quartiles

global_vals <- unlist(
  lapply(
    season_mean_list,
    \(x) values(x$mean, mat = FALSE)
  ),
  use.names = FALSE
)

global_q <- quantile(
  global_vals,
  probs = c(0.25, 0.5, 0.75),
  na.rm = TRUE
)

#global_q
#25%       50%       75% 
#0.1358552 0.1518939 0.1757852 

#min(global_vals, na.rm = T)
#[1] 0.1089491
#max(global_vals, na.rm = T)
#[1] 0.3304235

rm(global_vals)
gc()
##
global_vals_all <- unlist(
  lapply(
    prediction_list,
    \(x) values(x$stack, mat = FALSE)
  ),
  use.names = FALSE
)

global_q_all <- quantile(
  global_vals_all,
  probs = c(0.25, 0.5, 0.75),
  na.rm = TRUE
)

#global_q_all
#25%       50%       75% 
#0.1327533 0.1506227 0.1742195 

# min(global_vals_all, na.rm = T)
#[1] 0.1041551
# max(global_vals_all, na.rm = T)
#[1] 0.3574558
rm(global_vals_all)

# global vwc class

global_vwc_class_list <- purrr::map(
  season_mean_list,
  \(r)
  classify(
    r$mean,
    matrix(
      c(
        -Inf,        global_q[1], 1,
        global_q[1], global_q[2], 2,
        global_q[2], global_q[3], 3,
        global_q[3], Inf,         4
      ),
      ncol = 3,
      byrow = TRUE
    )
  )
)

names(global_vwc_class_list) <- season_names

#attach vwc and vwc class

pixel_df <- pixel_df |>
  group_by(season) |>
  group_modify(
    \(d, key){
      
      cls <- terra::extract(
        vwc_class_list[[key$season]],
        sample_pts,
        ID = FALSE
      )[[1]]
      
      global_cls <- terra::extract(
        global_vwc_class_list[[key$season]],
        sample_pts,
        ID = FALSE
      )[[1]]
      
      
      d$vwc_class <- factor(
        cls[d$pixel_id],
        levels = 1:4,
        labels = c(
          "dry",
          "intermediate_dry",
          "intermediate_wet",
          "wet"
        )
      )
      
      d$global_vwc_class <- factor(
        global_cls[d$pixel_id],
        levels = 1:4,
        labels = c(
          "dry",
          "intermediate_dry",
          "intermediate_wet",
          "wet"
        )
      )
      
      d
    }
  ) |>
  ungroup()
#############################################################
# add foehn

foehn_daily <- foehn_all %>%
  dplyr::group_by(station, season_year, date) %>%
  dplyr::summarise(
    foehn_day = dplyr::first(foehn_day),
    .groups = "drop"
  )

foehn_daily <- foehn_daily %>%
  dplyr::arrange(station, date) %>%
  dplyr::group_by(station) %>%
  dplyr::mutate(
    event = cumsum(
      foehn_day != dplyr::lag(foehn_day, default = FALSE)
    )
  ) %>%
  dplyr::group_by(station, event) %>%
  dplyr::mutate(
    foehn_duration = ifelse(first(foehn_day), dplyr::n(), 0L)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(-event)

foehn_all <- foehn_all %>%
  dplyr::left_join(
    foehn_daily %>%
      dplyr::select(station, date, foehn_duration),
    by = c("station", "date")
  )

foehn_info_df <- foehn_daily %>%
  dplyr::filter(station %in% c("boym", "exem", "frlm")) %>%
  dplyr::group_by(season_year, date) %>%
  dplyr::summarise(
    foehn = sum(foehn_day) >= 2,
    foehn_duration = ifelse(
      sum(foehn_day) >= 2,
      max(foehn_duration[foehn_day]),
      0L
    ),
    .groups = "drop"
  )


pixel_df$date <- as.character(pixel_df$date)

pixel_df <- pixel_df |>
  left_join(
    foehn_info_df,
    by = "date"
  )


pixel_df <- pixel_df |>
  mutate(
    event_occurs = event_effect != 0
  )

pixel_df$month <- month(pixel_df$date)
pixel_df$month <- as.factor(pixel_df$month)

#write pixel_df
str(pixel_df)
write.csv(pixel_df, "data/processed_data/statistic_models_final/pixel_df_mod_on_prediction_data_final.csv")
pixel_df <- read.csv("data/processed_data/statistic_models_final/pixel_df_mod_on_prediction_data_final.csv")
############################################################################################
#models with months (m) as additional fixed effect
# mod_occ md

build_anova_df <- function(model, t = "III"){
  anova <- car::Anova(model, type = t)
  
  anova_df <- data.frame(
    Term = rownames(anova),
    anova,
    row.names = NULL
  )
  
  anova_df <- anova_df %>%
    filter(Term != "(Intercept)") %>%
    mutate(
      Chisq = round(Chisq, 2),
      p = case_when(
        Pr..Chisq. < 0.001 ~ "<0.001",
        TRUE ~ sprintf("%.3f", Pr..Chisq.)
      )
    ) %>%
    dplyr::select(
      Term,
      Chisq,
      Df,
      p
    )
  return(anova_df)
}

#we don't have any foehn in November
pixel_12_df <- pixel_df[pixel_df$month != 11,]
pixel_12_df$month <- as.factor(pixel_12_df$month)
pixel_12_df$month <- droplevels(pixel_12_df$month)
pixel_12_df$season_factor <- as.factor(pixel_12_df$season_factor)

sp_mod_occ_md <- glmmTMB(
  event_occurs ~ foehn * season_factor + foehn * global_vwc_class + foehn*month +(1|pixel_id),
  dispformula = ~ season_factor,
  data = pixel_12_df,
  family = binomial()
)

dir.create("data/processed_data/statistic_models_final_2")
saveRDS(sp_mod_occ_md, "data/processed_data/statistic_models_final_2/sp_mod_occ_md.RDS")
sp_mod_occ_md <- readRDS("data/processed_data/statistic_models_final/sp_mod_occ_md.RDS")

sp_mod_occ_simres <- simulateResiduals(sp_mod_occ_md)
saveRDS(sp_mod_occ_simres, "data/processed_data/statistic_models_final_2/sp_mod_occ_simres.RDS")

testUniformity(sp_mod_occ_simres)
testDispersion(sp_mod_occ_simres)

summary(sp_mod_occ_md) 

car::Anova(sp_mod_occ_md, type = "III")

sp_mod_occ_anova_df <- build_anova_df(sp_mod_occ_md)
write.csv(sp_mod_occ_anova_df, "data/processed_data/statistic_models_final_2/sp_mod_occ_anova_df.csv", row.names = F)

sub_12_df <- pixel_12_df |>
  slice_sample(n = 50000)

#mod_checkoccmd <- glmmTMB(
#  event_occurs ~ foehn * season_factor+ foehn * global_vwc_class + foehn*month +(1|pixel_id),
#  dispformula = ~ season_factor,
#  data = sub_12_df,
 # family = binomial()
#)
#sp_mod_occ_simresmd <- simulateResiduals(mod_checkoccmd)
#testUniformity(sp_mod_occ_simresmd) #not sig
#testDispersion(sp_mod_occ_simresmd) #not sig

# mod_size wg
#filter to event effect != 0
size_12_df <- pixel_12_df %>%
  filter(event_effect != 0)

sp_mod_size_md <- glmmTMB(
  event_effect ~ foehn * season_factor + foehn*global_vwc_class +foehn*month + (1|pixel_id),
  dispformula = ~ season_factor,
  data = size_12_df,
  family = gaussian()
)
diagnose(sp_mod_size_md)


saveRDS(sp_mod_size_md, "data/processed_data/statistic_models_final_2/sp_mod_size_md.RDS")
sp_mod_size_md <- readRDS("data/processed_data/statistic_models_fina_2/sp_mod_size_md.RDS")
sp_mod_size_simres <- simulateResiduals(sp_mod_size_md)
saveRDS(sp_mod_size_simres, "data/processed_data/statistic_models_final_2/sp_mod_size_simres.RDS")

testUniformity(sp_mod_size_simres)
testDispersion(sp_mod_size_simres)

summary(sp_mod_size_md)

car::Anova(sp_mod_size_md, type = "III")

sp_mod_size_anova_df <- build_anova_df(sp_mod_size_md)
write.csv(sp_mod_size_anova_df, "data/processed_data/statistic_models_final_2/sp_mod_size_anova_df.csv", row.names = F)

#sub_pixel_12_df <- size_12_df |>
 # slice_sample(n = 50000)

#mod_checksizemd <- glmmTMB(
#  event_effect ~ foehn * season_factor+ foehn * global_vwc_class + foehn*month + (1|pixel_id),
#  dispformula = ~ season_factor,
#  data = sub_pixel_12_df,
#  family = gaussian()
#)

#sp_mod_size_simresmd <- simulateResiduals(mod_checksizemd)
#testUniformity(sp_mod_size_simresmd) #ok
#testDispersion(sp_mod_size_simresmd) #ok

# vwc month(m) duration(d)

vwc_mod_pixel_md <- glmmTMB(vwc ~ foehn*season_factor + foehn*month + (1|pixel_id),
                            dispformula = ~ season_factor,
                            data = pixel_12_df,
                            family = beta_family()
)

summary(vwc_mod_pixel_md)
car::Anova(vwc_mod_pixel_md, type = "III")
saveRDS(vwc_mod_pixel_md, "data/processed_data/statistic_models_final_2/sp_mod_pixel_m.RDS")
#vwc_mod_pixel_md <- readRDS("data/processed_data/statistic_models_final_2/sp_mod_pixel_m.RDS")

sp_mod_vwc_simres <- simulateResiduals(vwc_mod_pixel_md)
saveRDS(sp_mod_vwc_simres, "data/processed_data/statistic_models_final_2/sp_mod_vwc_simres.RDS")
sp_mod_vwc_simres <- readRDS("data/processed_data/statistic_models_final_2/sp_mod_vwc_simres.RDS")

testUniformity(sp_mod_vwc_simres)
testDispersion(sp_mod_vwc_simres)
testOutliers(sp_mod_vwc_simres)

sp_mod_vwc_anova_df <- build_anova_df(vwc_mod_pixel_md)
write.csv(sp_mod_vwc_anova_df, "data/processed_data/statistic_models_final_2/sp_mod_vwc_anova_df.csv", row.names = F)


diagnose(vwc_mod_pixel_md)
#mod_checkvwcmd <- glmmTMB(
#  vwc ~ foehn * season_factor+ foehn*month +(1|pixel_id),
#  dispformula = ~ season_factor,
#  data = sub_12_df,
#  family = beta_family()
#)
#sp_mod_vwc_simresmd <- simulateResiduals(mod_checkvwcmd)
#testUniformity(sp_mod_vwc_simresmd) #not sig
#testDispersion(sp_mod_vwc_simresmd) #not sig

##########################################################################################
#prepare plots
#anomaly occurrence
#month
emm_occ_month_group <- emmeans(
  sp_mod_occ_md,
  ~ month | foehn,
  type = "response"
)

emm_occ_month_group_df <- summary(emm_occ_month_group) |> as_tibble()

delta_occ_month_group <- emm_occ_month_group_df |>
  dplyr::group_by(foehn) |>
  dplyr::summarise(
    `12 - 1` = prob[month == "12"] - prob[month == "1"],
    `12 - 2` = prob[month == "12"] - prob[month == "2"],
    `1 - 2`  = prob[month == "1"]  - prob[month == "2"],
    .groups = "drop"
  )
delta_occ_month_group
pairs(emm_occ_month_group)
# season
emm_occ_season_group <- emmeans(
  sp_mod_occ_md,
  ~ season_factor | foehn,
  type = "response"
)

emm_occ_season_group_df <- summary(emm_occ_season_group) |> as_tibble()

delta_occ_season_group <- emm_occ_season_group_df |>
  dplyr::select(season_factor, foehn, prob) |>
  tidyr::pivot_wider(
    names_from = season_factor,
    values_from = prob
  ) |>
  dplyr::mutate(
    `2016 - 2017` = `2016` - `2017`,
    `2016 - 2019` = `2016` - `2019`,
    `2016 - 2021` = `2016` - `2021`,
    `2016 - 2022` = `2016` - `2022`,
    `2016 - 2023` = `2016` - `2023`,
    `2016 - 2024` = `2016` - `2024`,
    `2016 - 2025` = `2016` - `2025`,
    `2016 - 2026` = `2016` - `2026`,
    
    `2017 - 2019` = `2017` - `2019`,
    `2017 - 2021` = `2017` - `2021`,
    `2017 - 2022` = `2017` - `2022`,
    `2017 - 2023` = `2017` - `2023`,
    `2017 - 2024` = `2017` - `2024`,
    `2017 - 2025` = `2017` - `2025`,
    `2017 - 2026` = `2017` - `2026`,
    
    `2019 - 2021` = `2019` - `2021`,
    `2019 - 2022` = `2019` - `2022`,
    `2019 - 2023` = `2019` - `2023`,
    `2019 - 2024` = `2019` - `2024`,
    `2019 - 2025` = `2019` - `2025`,
    `2019 - 2026` = `2019` - `2026`,
    
    `2021 - 2022` = `2021` - `2022`,
    `2021 - 2023` = `2021` - `2023`,
    `2021 - 2024` = `2021` - `2024`,
    `2021 - 2025` = `2021` - `2025`,
    `2021 - 2026` = `2021` - `2026`,
    
    `2022 - 2023` = `2022` - `2023`,
    `2022 - 2024` = `2022` - `2024`,
    `2022 - 2025` = `2022` - `2025`,
    `2022 - 2026` = `2022` - `2026`,
    
    `2023 - 2024` = `2023` - `2024`,
    `2023 - 2025` = `2023` - `2025`,
    `2023 - 2026` = `2023` - `2026`,
    
    `2024 - 2025` = `2024` - `2025`,
    `2024 - 2026` = `2024` - `2026`,
    
    `2025 - 2026` = `2025` - `2026`
  )


delta_occ_season_group
pairs(emm_occ_season_group)

# letters occ foehn T
library(multcompView)

pair_df_true <- summary(pairs(emm_occ_season_group)) |>
  dplyr::filter(foehn == TRUE) |>
  dplyr::mutate(
    contrast = contrast |>
      gsub("season_factor", "", x = _) |>
      gsub(" / ", " - ", x = _)
  ) |>
  dplyr::left_join(
    delta_occ_season_group |>
      dplyr::filter(foehn == TRUE) |>
      tidyr::pivot_longer(
        cols = -foehn,
        names_to = "contrast",
        values_to = "delta"
      ),
    by = "contrast"
  ) |>
  dplyr::mutate(
    different = p.value < 0.05 & abs(delta) >= 0.01
  )
sum(is.na(pair_df_true$delta))

x_true <- pair_df_true$different

names(x_true) <- pair_df_true$contrast |>
  gsub(" - ", "-", x = _)

letters_occ_true <- multcompView::multcompLetters(x_true)$Letters

#letters occ foehn F
# occ season FALSE

pair_df_occ_false <- summary(pairs(emm_occ_season_group)) |>
  dplyr::filter(foehn == FALSE) |>
  dplyr::mutate(
    contrast = contrast |>
      gsub("season_factor", "", x = _) |>
      gsub(" / ", " - ", x = _)
  ) |>
  dplyr::left_join(
    delta_occ_season_group |>
      dplyr::filter(foehn == FALSE) |>
      tidyr::pivot_longer(
        cols = -foehn,
        names_to = "contrast",
        values_to = "delta"
      ),
    by = "contrast"
  ) |>
  dplyr::mutate(
    different = p.value < 0.05 & abs(delta) >= 0.01
  )


x_occ_false <- pair_df_occ_false$different

names(x_occ_false) <- pair_df_occ_false$contrast |>
  gsub(" - ", "-", x = _)

str(multcompView::multcompLetters(x_occ_false)$Letters)
letters_occ_false <-multcompView::multcompLetters(x_occ_false)$Letters

# VWC class
emm_occ_vwc_group <- emmeans(
  sp_mod_occ_md,
  ~ global_vwc_class | foehn,
  type = "response"
)

emm_occ_vwc_group_df <- summary(emm_occ_vwc_group) |> as_tibble()

delta_occ_vwc_group <- emm_occ_vwc_group_df |>
  dplyr::group_by(foehn) |>
  dplyr::summarise(
    `dry - intermediate_dry` = prob[global_vwc_class == "dry"] - prob[global_vwc_class == "intermediate_dry"],
    `dry - intermediate_wet` = prob[global_vwc_class == "dry"] - prob[global_vwc_class == "intermediate_wet"],
    `dry - wet` = prob[global_vwc_class == "dry"] - prob[global_vwc_class == "wet"],
    `intermediate_dry - wet` = prob[global_vwc_class == "intermediate_dry"] - prob[global_vwc_class == "wet"],
    `intermediate_dry - intermediate_wet` = prob[global_vwc_class == "intermediate_dry"] - prob[global_vwc_class == "intermediate_wet"],
    `intermediate_wet - wet` = prob[global_vwc_class == "intermediate_wet"] - prob[global_vwc_class == "wet"],
    .groups = "drop"
  )


delta_occ_vwc_group

delta_occ_season_group |>
  dplyr::select(foehn, dplyr::contains(" - ")) |>
  tidyr::pivot_longer(
    cols = -foehn,
    names_to = "comparison",
    values_to = "delta"
  ) |>
  dplyr::filter(abs(delta) < 0.01)

pairs(emm_occ_vwc_group)
# anomaly size
## month
emm_size_month_group <- emmeans(
  sp_mod_size_md,
  ~ month | foehn,
  type = "response"
)

emm_size_month_group_df <- summary(emm_size_month_group) |> as_tibble()

delta_size_month_group <- emm_size_month_group_df |>
  dplyr::group_by(foehn) |>
  dplyr::summarise(
    `12 - 1` = emmean[month == "12"] - emmean[month == "1"],
    `12 - 2` = emmean[month == "12"] - emmean[month == "2"],
    `1 - 2`  = emmean[month == "1"]  - emmean[month == "2"],
    .groups = "drop"
  )

delta_size_month_group
#alternative
#month_contrasts <- contrast(
#  emm_size_month_group,
#  method = list(
#    "12 - 1" = c(-1, 0, 1),
#    "12 - 2" = c(0, -1, 1),
#    "1 - 2"  = c(1, -1, 0)
#  )
#)

summary(month_contrasts)
pairs(emm_size_month_group)
## season 
emm_size_season_group <- emmeans(
  sp_mod_size_md,
  ~ season_factor | foehn,
  type = "response"
)

emm_size_season_group_df <- summary(emm_size_season_group) |> as_tibble()


delta_size_season_group <- emm_size_season_group_df |>
  dplyr::select(season_factor, foehn, emmean) |>
  tidyr::pivot_wider(
    names_from = season_factor,
    values_from = emmean
  ) |>
  dplyr::mutate(
    `2016 - 2017` = `2016` - `2017`,
    `2016 - 2019` = `2016` - `2019`,
    `2016 - 2021` = `2016` - `2021`,
    `2016 - 2022` = `2016` - `2022`,
    `2016 - 2023` = `2016` - `2023`,
    `2016 - 2024` = `2016` - `2024`,
    `2016 - 2025` = `2016` - `2025`,
    `2016 - 2026` = `2016` - `2026`,
    
    `2017 - 2019` = `2017` - `2019`,
    `2017 - 2021` = `2017` - `2021`,
    `2017 - 2022` = `2017` - `2022`,
    `2017 - 2023` = `2017` - `2023`,
    `2017 - 2024` = `2017` - `2024`,
    `2017 - 2025` = `2017` - `2025`,
    `2017 - 2026` = `2017` - `2026`,
    
    `2019 - 2021` = `2019` - `2021`,
    `2019 - 2022` = `2019` - `2022`,
    `2019 - 2023` = `2019` - `2023`,
    `2019 - 2024` = `2019` - `2024`,
    `2019 - 2025` = `2019` - `2025`,
    `2019 - 2026` = `2019` - `2026`,
    
    `2021 - 2022` = `2021` - `2022`,
    `2021 - 2023` = `2021` - `2023`,
    `2021 - 2024` = `2021` - `2024`,
    `2021 - 2025` = `2021` - `2025`,
    `2021 - 2026` = `2021` - `2026`,
    
    `2022 - 2023` = `2022` - `2023`,
    `2022 - 2024` = `2022` - `2024`,
    `2022 - 2025` = `2022` - `2025`,
    `2022 - 2026` = `2022` - `2026`,
    
    `2023 - 2024` = `2023` - `2024`,
    `2023 - 2025` = `2023` - `2025`,
    `2023 - 2026` = `2023` - `2026`,
    
    `2024 - 2025` = `2024` - `2025`,
    `2024 - 2026` = `2024` - `2026`,
    
    `2025 - 2026` = `2025` - `2026`
  )
delta_size_season_group

# size season letters TRUE

pair_df_size_true <- summary(pairs(emm_size_season_group)) |>
  dplyr::filter(foehn == TRUE) |>
  dplyr::mutate(
    contrast = contrast |>
      gsub("season_factor", "", x = _) |>
      gsub(" / ", " - ", x = _)
  ) |>
  dplyr::left_join(
    delta_size_season_group |>
      dplyr::filter(foehn == TRUE) |>
      tidyr::pivot_longer(
        cols = -foehn,
        names_to = "contrast",
        values_to = "delta"
      ),
    by = "contrast"
  ) |>
  dplyr::mutate(
    different = p.value < 0.05 & abs(delta) >= 0.01
  )


x_size_true <- pair_df_size_true$different

names(x_size_true) <- pair_df_size_true$contrast |>
  gsub(" - ", "-", x = _)

letters_size_true <- multcompView::multcompLetters(x_size_true)$Letters

# size season letters FALSE

pair_df_size_false <- summary(pairs(emm_size_season_group)) |>
  dplyr::filter(foehn == FALSE) |>
  dplyr::mutate(
    contrast = contrast |>
      gsub("season_factor", "", x = _) |>
      gsub(" / ", " - ", x = _)
  ) |>
  dplyr::left_join(
    delta_size_season_group |>
      dplyr::filter(foehn == FALSE) |>
      tidyr::pivot_longer(
        cols = -foehn,
        names_to = "contrast",
        values_to = "delta"
      ),
    by = "contrast"
  ) |>
  dplyr::mutate(
    different = p.value < 0.05 & abs(delta) >= 0.01
  )


x_size_false <- pair_df_size_false$different

names(x_size_false) <- pair_df_size_false$contrast |>
  gsub(" - ", "-", x = _)

letters_size_false <- multcompView::multcompLetters(x_size_false)$Letters

## VWC class
emm_size_vwc_group <- emmeans(
  sp_mod_size_md,
  ~ global_vwc_class | foehn,
  type = "response"
)

emm_size_vwc_group_df <- summary(emm_size_vwc_group) |> as_tibble()

delta_size_vwc_group <- emm_size_vwc_group_df |>
  dplyr::group_by(foehn) |>
  dplyr::summarise(
    `dry - intermediate_dry` = emmean[global_vwc_class == "dry"] - emmean[global_vwc_class == "intermediate_dry"],
    `dry - intermediate_wet` = emmean[global_vwc_class == "dry"] - emmean[global_vwc_class == "intermediate_wet"],
    `dry - wet` = emmean[global_vwc_class == "dry"] - emmean[global_vwc_class == "wet"],
    `intermediate_dry - wet` = emmean[global_vwc_class == "intermediate_dry"] - emmean[global_vwc_class == "wet"],
    `intermediate_dry - intermediate_wet` = emmean[global_vwc_class == "intermediate_dry"] - emmean[global_vwc_class == "intermediate_wet"],
    `intermediate_wet - wet` = emmean[global_vwc_class == "intermediate_wet"] - emmean[global_vwc_class == "wet"],
    .groups = "drop"
  )

delta_size_vwc_group
pairs(emm_size_vwc_group)
############################################################
# interaction plot function
plot_emm_interaction <- function(
    emm,
    x,
    y,
    ylab,
    xlab = x,
    ylim = NULL,
    LCL = "asymp.LCL",
    UCL = "asymp.UCL"
) {
  
  df <- summary(emm) |>
    as_tibble() |>
    mutate(
      foehn = as.character(foehn),
      foehn = recode(
        foehn,
        "TRUE" = "Foehn",
        "FALSE" = "No foehn"
      )
    )
  
  
  sig <- pairs(emm, reverse = TRUE) |>
    summary() |>
    as_tibble() |>
    mutate(
      stars = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE ~ ""
      )
    )
  
  
  wide <- df |>
    dplyr::select(all_of(x), foehn, all_of(y)) |>
    pivot_wider(
      names_from = foehn,
      values_from = all_of(y)
    )
  
  
  if(all(c("Foehn","No foehn") %in% names(wide))){
    
    wide <- wide |>
      mutate(
        diff = Foehn - `No foehn`,
        relevant = abs(diff) > 0.01
      )
    
  } else {
    
    wide <- wide |>
      mutate(
        relevant = FALSE
      )
  }
  
  
  df[[x]] <- as.factor(df[[x]])
  sig[[x]] <- as.factor(sig[[x]])
  wide[[x]] <- as.factor(wide[[x]])
  
  
  # sort
  if(x == "month"){
    
    df[[x]] <- factor(
      df[[x]],
      levels = c("12","1","2")
    )
    
    sig[[x]] <- factor(
      sig[[x]],
      levels = c("12","1","2")
    )
    
  }
  
  
  if(x == "season_factor"){
    
    season_levels <- c(
      "2016","2017","2019",
      "2021","2022","2023",
      "2024","2025","2026"
    )
    
    df[[x]] <- factor(
      df[[x]],
      levels = season_levels
    )
    
    sig[[x]] <- factor(
      sig[[x]],
      levels = season_levels
    )
  }
  
  
  if(is.null(ylim)){
    
    ymin <- min(df[[LCL]], na.rm = TRUE)
    ymax <- max(df[[UCL]], na.rm = TRUE)
    
  } else {
    
    ymin <- ylim[1]
    ymax <- ylim[2]
  }
  
  
  anno <- left_join(
    sig,
    wide,
    by = x
  ) |>
    mutate(
      star_colour = ifelse(
        relevant,
        "black",
        NA   # use "grey60" if significant but irrelevant stars should be plottted
      ),
      y = ymax + (ymax-ymin)*0.08
    )
  
  
  foehn_cols_new <- c(
    "No foehn" = "#0072B2",
    "Foehn" = "#E69F00"
  )
  
  
  p <- ggplot(
    df,
    aes(
      x = .data[[x]],
      y = .data[[y]]
    )
  ) +
    
    geom_point(
      aes(colour = foehn),
      size = 2.5,
      position = position_dodge(width = 0.4)
    ) +
    
    geom_errorbar(
      aes(
        ymin = .data[[LCL]],
        ymax = .data[[UCL]],
        colour = foehn
      ),
      width = 0.15,
      position = position_dodge(width = 0.4)
    )+
    
    geom_text(
      data = anno,
      aes(
        x = .data[[x]],
        y = y,
        label = stars
      ),
      inherit.aes = FALSE,
      colour = anno$star_colour,
      size = 5,
      show.legend = FALSE
    ) +
    
    scale_colour_manual(
      values = foehn_cols_new
    ) +
    
    coord_cartesian(
      ylim = c(
        ymin,
        ymax + (ymax-ymin)*0.15
      )
    ) +
    
    scale_y_continuous(
      expand = c(0,0)
    )+
    
    theme_bw() +
    theme(
      text = element_text(
        family = "Calibri",
        size = 14
      ),
      axis.title = element_text(size = 16),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(),
      panel.grid.minor.y = element_line(),
      strip.text = element_blank()
    ) +
    
    labs(
      x = xlab,
      y = ylab,
      colour = NULL
    )
  
  
  
  
  if(x == "season_factor"){
    
    p <- p +
      scale_x_discrete(
        labels = c(
          "2016"="15|16",
          "2017"="16|17",
          "2019"="18|19",
          "2021"="20|21",
          "2022"="21|22",
          "2023"="22|23",
          "2024"="23|24",
          "2025"="24|25",
          "2026"="25|26"
        )
      )
  }
  # abbreviation, for intermediate
  if (x == "global_vwc_class") {
    
    p <- p +
      scale_x_discrete(
        labels = c(
          "dry"="dry",
          "intermediate_dry"="int. dry",
          "intermediate_wet"="int. wet",
          "wet"="wet"
        )
      ) 
    
  }
  
  
  if(x == "month"){
    
    p <- p +
      scale_x_discrete(
        limits = c("12","1","2"),
        labels = c(
          "12" = "Dec",
          "1"  = "Jan",
          "2"  = "Feb"
        )
      )
  }
  
  

  if(x %in% c("season_factor", "global_vwc_class")){
    
    p <- p +
      theme(
        axis.text.x = element_text(
          angle = 30,
          hjust = 1,
          vjust = 1
        )
      )
  }
  
  return(p)
}

#initialize Calibri
windowsFonts(
  Calibri = windowsFont("Calibri")
)
#interaction plots
#occurence
ylim_occ <- c(
  0,
  max(
    emm_occ_month_group_df$asymp.UCL,
    emm_occ_vwc_group_df$asymp.UCL,
    emm_occ_season_group_df$asymp.UCL,
    na.rm = TRUE
  )
)


emm_month_occ_md <- emmeans(
  sp_mod_occ_md,
  ~ foehn | month,
  type = "response"
)


p_occ_month_md <- plot_emm_interaction(
  emm_month_occ_md,
  "month",
  "prob",
  "Estimated occurrence probability",
  ylim = ylim_occ,
  xlab = "Month"
)
p_occ_month_md 


# occ glob_vwc

emm_gvwc_class_occ_md <- emmeans(
  sp_mod_occ_md,
  ~ foehn | global_vwc_class,
  type = "response"
)


p_occ_vwc_md <- plot_emm_interaction(
  emm_gvwc_class_occ_md,
  "global_vwc_class",
  "prob",
  "Estimated occurrence probability",
  ylim = ylim_occ,
  xlab = "VWC class"
)

p_occ_vwc_md

# occ season

emm_season_occ_md <- emmeans(
  sp_mod_occ_md,
  ~ foehn | season_factor,
  type = "response"
)

p_occ_season_md <- plot_emm_interaction(
  emm_season_occ_md,
  "season_factor",
  "prob",
  "Estimated occurrence probability",
  ylim = ylim_occ,
  xlab = "Austral summer season"
)

p_occ_season_md

#(p_occ_month|p_occ_season_md)

ylim_size <- c(
  min(
    emm_size_month_group_df$lower.CL,
    emm_size_vwc_group_df$lower.CL,
    emm_size_season_group_df$lower.CL
  )-0.01,
  max(
    emm_size_month_group_df$upper.CL,
    emm_size_vwc_group_df$upper.CL,
    emm_size_season_group_df$upper.CL
  )
)

#ylim_size <- c(0,
#    max(
#      emm_size_month_group_df$asymp.UCL,
#      emm_size_vwc_group_df$asymp.UCL,
#     emm_size_season_group_df$asymp.UCL
#   )
#)


emm_month_size_md <- emmeans(
  sp_mod_size_md,
  ~ foehn | month,
  type = "response"
)


p_size_month_md <- plot_emm_interaction(
  emm_month_size_md,
  "month",
  "emmean",
  "Estimated mean event size",
  ylim = ylim_size,
  xlab = "Month",
  UCL = "upper.CL",
  LCL = "lower.CL"
)
p_size_month_md


# size gwc class
emm_gvwc_class_size_md <- emmeans(
  sp_mod_size_md,
  ~ foehn | global_vwc_class,
  type = "response"
)

p_size_vwc_md <- plot_emm_interaction(
  emm_gvwc_class_size_md,
  "global_vwc_class",
  "emmean",
  "Estimated mean event size",
  ylim = ylim_size,
  xlab = "VWC class",
  UCL = "upper.CL",
  LCL = "lower.CL"
)

p_size_vwc_md

# size season 
emm_season_size_md <- emmeans(
  sp_mod_size_md,
  ~ foehn | season_factor,
  type = "response"
)

p_size_season_md <- plot_emm_interaction(
  emm_season_size_md,
  "season_factor",
  "emmean",
  "Estimated mean event size",
  ylim = ylim_size,
  xlab = "Austral summer season",
  UCL = "upper.CL",
  LCL = "lower.CL"
)

p_size_season_md

#(p_size_month_md|p_size_season_md)

############################################################
#foehn_not foehn plots
library(ggplot2)
library(dplyr)
library(tibble)

plot_emm_letters <- function(
    emm_df,
    letters_df,
    x,
    y,
    ylab,
    xlab = x,
    ylim = NULL,
    UCL = "asymp.UCL",
    LCL = "asymp.LCL"
    
) {
  
  df <- emm_df |>
    as_tibble() |>
    mutate(
      foehn = as.character(foehn),
      foehn = recode(
        foehn,
        "TRUE" = "Foehn",
        "FALSE" = "No foehn"
      )
    )
  
  
  letters_df <- letters_df |>
    mutate(
      foehn = as.character(foehn),
      foehn = recode(
        foehn,
        "TRUE" = "Foehn",
        "FALSE" = "No foehn"
      )
    )
  
  
  
  df[[x]] <- as.factor(df[[x]])
  
  if (x == "month") {
    
    df[[x]] <- factor(
      df[[x]],
      levels = c("12","1","2")
    )
    
    letters_df[[x]] <- factor(
      letters_df[[x]],
      levels = c("12","1","2")
    )
    
  }
  
  
  if (x == "season_factor") {
    
    season_levels <- c(
      "2016","2017","2019",
      "2021","2022","2023",
      "2024","2025","2026"
    )
    
    df[[x]] <- factor(
      df[[x]],
      levels = season_levels
    )
    
    letters_df[[x]] <- factor(
      letters_df[[x]],
      levels = season_levels
    )
  }
  
  
  # add letters
  df <- df |>
    left_join(
      letters_df,
      by = c("foehn", x)
    )
  
  
  
  if (is.null(ylim)) {
    
    ymin <- min(df[[LCL]], na.rm = TRUE)
    ymax <- max(df[[UCL]], na.rm = TRUE)
    
  } else {
    
    ymin <- ylim[1]
    ymax <- ylim[2]
    
  }
  
  
  # letters at same height
  anno <- df |>
    group_by(foehn) |>
    summarise(
      y = max(.data[[LCL]], na.rm = TRUE) + 
        (ymax - ymin) * 0.12,
      .groups = "drop"
    ) |>
    right_join(
      df |> 
        dplyr::select(foehn, all_of(x), letter) |>
        distinct(),
      by = "foehn"
    )
  
  foehn_cols_new <- c(
    "No foehn" = "#0072B2",
    "Foehn" = "#E69F00"
  )
  
  
  p <- ggplot(
    df,
    aes(
      x = .data[[x]],
      y = .data[[y]],
      colour = foehn
    )
  ) +
    
    geom_point(
      size = 2.5,
      position = position_dodge(width = 0.35)
    ) +
    
    geom_errorbar(
      aes(
        ymin = .data[[LCL]],
        ymax = .data[[UCL]]
      ),
      width = 0.15,
      position = position_dodge(width = 0.35)
    ) +
    
    geom_text(
      data = anno,
      aes(
        x = .data[[x]],
        y = y,
        label = letter
      ),
      inherit.aes = FALSE,
      colour = "black",
      size = 3.5,
      family = "Calibri",
      vjust = 0,
      show.legend = FALSE
    )+
    
    facet_wrap(
      ~foehn,
      ncol = 1
    ) +
    
    scale_colour_manual(
      values = foehn_cols_new
    ) +
    coord_cartesian(
      ylim = c(
        ymin,
        ymax + (ymax-ymin)*0.25
      )
    ) +
    
    scale_y_continuous(
      expand = c(0,0)
    )+
    theme_bw() +
    
    theme(
      text = element_text(
        family = "Calibri",
        size = 14
      ),
      
      axis.title = element_text(
        size = 16
      ),
      
      axis.text = element_text(
        size = 14
      ),
      
      legend.text = element_text(
        size = 14
      ),
      
      legend.title = element_blank(),
      
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      
      panel.grid.major.y = element_line(),
      panel.grid.minor.y = element_line(),
      
      strip.text = element_blank(),
      
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5
      ),
      panel.spacing.y = unit(0.7, "lines") #set to 0.3 for occurrence plots
    ) +
    
    labs(
      x = xlab,
      y = ylab,
      colour = NULL
    )
  
 
  if (x == "season_factor") {
    
    p <- p +
      scale_x_discrete(
        labels = c(
          "2016"="15|16",
          "2017"="16|17",
          "2019"="18|19",
          "2021"="20|21",
          "2022"="21|22",
          "2023"="22|23",
          "2024"="23|24",
          "2025"="24|25",
          "2026"="25|26"
        )
      ) 
    
  }
  
  # abbreviation, for intermediate
  if (x == "global_vwc_class") {
    
    p <- p +
      scale_x_discrete(
        labels = c(
          "dry"="dry",
          "intermediate_dry"="int. dry",
          "intermediate_wet"="int. wet",
          "wet"="wet"
        )
      ) 
    
  }
 
  if (x == "month") {
    
    p <- p +
      scale_x_discrete(
        limits = c("12","1","2"),
        labels = c(
          "12" = "Dec",
          "1"  = "Jan",
          "2"  = "Feb"
        )
      )
  }
  
  if(x %in% c("season_factor", "global_vwc_class")){
    
    p <- p +
      theme(
        axis.text.x = element_text(
          angle = 30,
          hjust = 1,
          vjust = 1
        )
      )
  }
  
  return(p)
}

##############################################################################
# letter objects
letters_occ_month <- bind_rows(
  tibble(
    foehn = "FALSE",
    month = c(12,1,2),
    letter = c("a","b","c")
  ),
  tibble(
    foehn = "TRUE",
    month = c(12,1,2),
    letter = c("a","b","c")
  )
)

letters_occ_vwc <- tibble(
  foehn = factor(
    rep(c(FALSE, TRUE), each = 4),
    levels = c(FALSE, TRUE)
  ),
  global_vwc_class = rep(
    c("dry",
      "intermediate_dry",
      "intermediate_wet",
      "wet"),
    2
  ),
  letter = c(
    "a", "b", "c", "b",
    "a", "b", "b", "b"
  )
)

letters_occ_season <- bind_rows(
  tibble(
    foehn = "FALSE",
    season_factor = names(letters_occ_false),
    letter = unname(letters_occ_false)
  ),
  tibble(
    foehn = "TRUE",
    season_factor = names(letters_occ_true),
    letter = unname(letters_occ_true)
  )
)


letters_size_season <- bind_rows(
  tibble(
    foehn = "FALSE",
    season_factor = names(letters_size_false),
    letter = unname(letters_size_false)
  ),
  tibble(
    foehn = "TRUE",
    season_factor = names(letters_size_true),
    letter = unname(letters_size_true)
  )
)


letters_size_month <- tibble(
  foehn = factor(
    rep(c(FALSE, TRUE), each = 3),
    levels = c(FALSE, TRUE)
  ),
  month = rep(c("12", "1", "2"), 2),
  letter = c(
    "a", "b", "c",
    "a", "a", "b"
  )
)

letters_size_vwc <- tibble(
  foehn = factor(
    rep(c(FALSE, TRUE), each = 4),
    levels = c(FALSE, TRUE)
  ),
  global_vwc_class = rep(
    c("dry",
      "intermediate_dry",
      "intermediate_wet",
      "wet"),
    2
  ),
  letter = c(
    "a", "a", "a", "a",
    "a", "a", "a", "a"
  )
)

#######################
p_occ_month_letters <- plot_emm_letters(
  emm_occ_month_group_df,
  letters_occ_month,
  "month",
  "prob",
  "Estimated occurrence probability",
  ylim = ylim_occ,
  xlab = "Month"
)

p_occ_season_letters <- plot_emm_letters(
  emm_occ_season_group_df,
  letters_occ_season,
  "season_factor",
  "prob",
  "Estimated occurrence probability",
  ylim = ylim_occ,
  xlab = "Austral summer season"
)

p_occ_vwc_letters <- plot_emm_letters(
  emm_occ_vwc_group_df,
  letters_occ_vwc,
  "global_vwc_class",
  "prob",
  "Estimated occurrence probability",
  ylim = ylim_occ,
  xlab = "VWC class"
)

# size

p_size_month_letters <- plot_emm_letters(
  emm_size_month_group_df,
  letters_size_month,
  "month",
  "emmean",
  "Estimated mean event size",
  ylim = ylim_size,
  xlab = "Month",
  LCL = "lower.CL",
  UCL = "upper.CL"
)

p_size_month_letters

p_size_season_letters <- plot_emm_letters(
  emm_size_season_group_df,
  letters_size_season,
  "season_factor",
  "emmean",
  "Estimated mean event size",
  ylim = ylim_size,
  xlab = "Austral summer season",
  LCL = "lower.CL",
  UCL = "upper.CL"
)

p_size_vwc_letters <- plot_emm_letters(
  emm_size_vwc_group_df,
  letters_size_vwc,
  "global_vwc_class",
  "emmean",
  "Estimated mean event size",
  ylim = ylim_size,
  xlab = "VWC class",
  LCL = "lower.CL",
  UCL = "upper.CL"
)
###############################################################################
#patchwork plots
library(patchwork)
clean_for_grid <- function(p, keep_legend = FALSE){
  
  p <- p +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      
      axis.text.x = element_text(
        family = "Calibri",
        size = 12,
        angle = 45,
        hjust = 1
      ),
      
      axis.text.y = element_text(
        family = "Calibri",
        size = 12
      ),
      
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  if(keep_legend){
    
    p <- p +
      theme(
        legend.position = "bottom",
        legend.text = element_text(
          family = "Calibri",
          size = 12
        )
      )
    
  } else {
    
    p <- p +
      guides(
        colour = "none"
      ) +
      theme(
        legend.position = "none"
      )
  }
  
  return(p)
}

occ_month_md_clean <- clean_for_grid(
  p_occ_month_md,
  keep_legend = TRUE
)

occ_month_letters_clean <- clean_for_grid(
  p_occ_month_letters,
  keep_legend = FALSE
)

occ_season_md_clean <- clean_for_grid(
  p_occ_season_md,
  keep_legend = FALSE
)

occ_season_letters_clean <- clean_for_grid(
  p_occ_season_letters,
  keep_legend = FALSE
)

occ_vwc_md_clean <- clean_for_grid(
  p_occ_vwc_md,
  keep_legend = FALSE
)

occ_vwc_letters_clean <- clean_for_grid(
  p_occ_vwc_letters,
  keep_legend = FALSE
)

occ_combined <-
  
  (
    occ_month_md_clean |
      occ_month_letters_clean
  ) /
  (
    occ_season_md_clean |
      occ_season_letters_clean
  ) /
  (
    occ_vwc_md_clean |
      occ_vwc_letters_clean
  ) +
  
  plot_layout(
    guides = "collect"
  ) &
  
  theme(
    legend.position = "bottom"
  )
occ_combined

library(grid)

occ_final <-
  
  wrap_elements(
    textGrob(
      "Estimated anomaly occurrence probability",
      rot = 90,
      gp = gpar(
        fontsize = 16,
        fontfamily = "Calibri"
      )
    )
  ) +
  occ_combined +
  
  plot_layout(
    widths = c(0.05, 1)
  )

occ_final
# occ_comb_no vwc


occ_combined_red <-
  
  (
    occ_month_md_clean |
      occ_month_letters_clean
  ) /
  (
    occ_season_md_clean |
      occ_season_letters_clean
  )  +
  
  plot_layout(
    guides = "collect"
  ) &
  
  theme(
    legend.position = "bottom"
  )
occ_combined_red

occ_final_red <-
  
  wrap_elements(
    textGrob(
      "Estimated anomaly occurrence probability",
      rot = 90,
      gp = gpar(
        fontsize = 16,
        fontfamily = "Calibri"
      )
    )
  ) +
  occ_combined_red +
  
  plot_layout(
    widths = c(0.05, 1)
  )

occ_final_red

#size
size_month_md_clean <- clean_for_grid(
  p_size_month_md,
  keep_legend = TRUE
)

size_month_letters_clean <- clean_for_grid(
  p_size_month_letters,
  keep_legend = FALSE
)

size_season_md_clean <- clean_for_grid(
  p_size_season_md,
  keep_legend = FALSE
)

size_season_letters_clean <- clean_for_grid(
  p_size_season_letters,
  keep_legend = FALSE
)

size_vwc_md_clean <- clean_for_grid(
  p_size_vwc_md,
  keep_legend = FALSE
)

size_vwc_letters_clean <- clean_for_grid(
  p_size_vwc_letters,
  keep_legend = FALSE
)


size_combined <-
  
  (
    size_month_md_clean |
      size_month_letters_clean
  ) /
  (
    size_season_md_clean |
      size_season_letters_clean
  ) /
  (
    size_vwc_md_clean |
      size_vwc_letters_clean
  ) +
  
  plot_layout(
    guides = "collect"
  ) &
  
  theme(
    legend.position = "bottom"
  )

size_combined

size_final <-
  
  wrap_elements(
    textGrob(
      "Estimated mean anomaly size",
      rot = 90,
      gp = gpar(
        fontsize = 16,
        fontfamily = "Calibri"
      )
    )
  ) +
  size_combined +
  
  plot_layout(
    widths = c(0.05, 1)
  )


size_final
# size final reduced
size_combined_red <-
  
  (
    size_month_md_clean |
      size_month_letters_clean
  ) /
  (
    size_season_md_clean |
      size_season_letters_clean
  ) +
  
  plot_layout(
    guides = "collect"
  ) &
  
  theme(
    legend.position = "bottom"
  )


size_combined_red


size_final_red <-
  
  wrap_elements(
    textGrob(
      "Estimated mean anomaly size",
      rot = 90,
      gp = gpar(
        fontsize = 16,
        fontfamily = "Calibri"
      )
    )
  ) +
  size_combined_red +
  
  plot_layout(
    widths = c(0.05, 1)
  )


size_final_red

#save
library(Cairo)
#final full plots
CairoPNG(
  filename = "plots/occurrence_probability_pixel_data_black_only_AOA_final_3.png",
  width = 4200,
  height = 5250,
  res = 600
)

print(occ_final)

dev.off()


CairoPNG(
  filename = "plots/event_size_pixel_data_black_only_AOA_final_3.png",
  width = 4200,
  height = 5250,
  res = 600
)

print(size_final)

dev.off()
#final reduced plots
CairoPNG(
  filename = "plots/occurrence_probability_pixel_data_black_only_red_AOA_final_2.png",
  width = 4200,
  height = 4200,
  res = 600
)

print(occ_final_red)

dev.off()


CairoPNG(
  filename = "plots/event_size_pixel_data_black_only_red_AOA_final_2.png",
  width = 4200,
  height = 4200,
  res = 600
)

print(size_final_red)

dev.off()

####################
#plots for vwc model 
# month
emm_vwc_month_group <- emmeans(
  vwc_mod_pixel_md,
  ~ month | foehn,
  type = "response"
)

emm_vwc_month_group_df <- summary(emm_vwc_month_group) |> as_tibble()

delta_vwc_month_group <- emm_vwc_month_group_df |>
  dplyr::group_by(foehn) |>
  dplyr::summarise(
    `12 - 1` = response[month == "12"] - response[month == "1"],
    `12 - 2` = response[month == "12"] - response[month == "2"],
    `1 - 2`  = response[month == "1"]  - response[month == "2"],
    .groups = "drop"
  )
delta_vwc_month_group
pairs(emm_vwc_month_group)
#vwc season
emm_vwc_season_group <- emmeans(
  vwc_mod_pixel_md,
  ~ season_factor | foehn,
  type = "response"
)

emm_vwc_season_group_df <- summary(emm_vwc_season_group) |> as_tibble()


delta_vwc_season_group <- emm_vwc_season_group_df |>
  dplyr::select(season_factor, foehn, response) |>
  tidyr::pivot_wider(
    names_from = season_factor,
    values_from = response
  ) |>
  dplyr::mutate(
    `2016 - 2017` = `2016` - `2017`,
    `2016 - 2019` = `2016` - `2019`,
    `2016 - 2021` = `2016` - `2021`,
    `2016 - 2022` = `2016` - `2022`,
    `2016 - 2023` = `2016` - `2023`,
    `2016 - 2024` = `2016` - `2024`,
    `2016 - 2025` = `2016` - `2025`,
    `2016 - 2026` = `2016` - `2026`,
    
    `2017 - 2019` = `2017` - `2019`,
    `2017 - 2021` = `2017` - `2021`,
    `2017 - 2022` = `2017` - `2022`,
    `2017 - 2023` = `2017` - `2023`,
    `2017 - 2024` = `2017` - `2024`,
    `2017 - 2025` = `2017` - `2025`,
    `2017 - 2026` = `2017` - `2026`,
    
    `2019 - 2021` = `2019` - `2021`,
    `2019 - 2022` = `2019` - `2022`,
    `2019 - 2023` = `2019` - `2023`,
    `2019 - 2024` = `2019` - `2024`,
    `2019 - 2025` = `2019` - `2025`,
    `2019 - 2026` = `2019` - `2026`,
    
    `2021 - 2022` = `2021` - `2022`,
    `2021 - 2023` = `2021` - `2023`,
    `2021 - 2024` = `2021` - `2024`,
    `2021 - 2025` = `2021` - `2025`,
    `2021 - 2026` = `2021` - `2026`,
    
    `2022 - 2023` = `2022` - `2023`,
    `2022 - 2024` = `2022` - `2024`,
    `2022 - 2025` = `2022` - `2025`,
    `2022 - 2026` = `2022` - `2026`,
    
    `2023 - 2024` = `2023` - `2024`,
    `2023 - 2025` = `2023` - `2025`,
    `2023 - 2026` = `2023` - `2026`,
    
    `2024 - 2025` = `2024` - `2025`,
    `2024 - 2026` = `2024` - `2026`,
    
    `2025 - 2026` = `2025` - `2026`
  )
delta_vwc_season_group

# vwc season TRUE

pair_df_vwc_true <- summary(pairs(emm_vwc_season_group)) |>
  dplyr::filter(foehn == TRUE) |>
  dplyr::mutate(
    contrast = contrast |>
      gsub("season_factor", "", x = _) |>
      gsub(" / ", " - ", x = _)
  ) |>
  dplyr::left_join(
    delta_vwc_season_group |>
      dplyr::filter(foehn == TRUE) |>
      tidyr::pivot_longer(
        cols = -foehn,
        names_to = "contrast",
        values_to = "delta"
      ),
    by = "contrast"
  ) |>
  dplyr::mutate(
    different = p.value < 0.05 & abs(delta) >= 0.01
  )


x_vwc_true <- pair_df_vwc_true$different

names(x_vwc_true) <- pair_df_vwc_true$contrast |>
  gsub(" - ", "-", x = _)

letters_vwc_true <- multcompView::multcompLetters(x_vwc_true)$Letters



# vwc season FALSE

pair_df_vwc_false <- summary(pairs(emm_vwc_season_group)) |>
  dplyr::filter(foehn == FALSE) |>
  dplyr::mutate(
    contrast = contrast |>
      gsub("season_factor", "", x = _) |>
      gsub(" / ", " - ", x = _)
  ) |>
  dplyr::left_join(
    delta_vwc_season_group |>
      dplyr::filter(foehn == FALSE) |>
      tidyr::pivot_longer(
        cols = -foehn,
        names_to = "contrast",
        values_to = "delta"
      ),
    by = "contrast"
  ) |>
  dplyr::mutate(
    different = p.value < 0.05 & abs(delta) >= 0.01
  )


x_vwc_false <- pair_df_vwc_false$different

names(x_vwc_false) <- pair_df_vwc_false$contrast |>
  gsub(" - ", "-", x = _)

letters_vwc_false <- multcompView::multcompLetters(x_vwc_false)$Letters

#letters

letters_vwc_season <- bind_rows(
  tibble(
    foehn = "FALSE",
    season_factor = names(letters_vwc_false),
    letter = unname(letters_vwc_false)
  ),
  tibble(
    foehn = "TRUE",
    season_factor = names(letters_vwc_true),
    letter = unname(letters_vwc_true)
  )
)


letters_vwc_month <- tibble(
  foehn = factor(
    rep(c(FALSE, TRUE), each = 3),
    levels = c(FALSE, TRUE)
  ),
  month = rep(c("12", "1", "2"), 2),
  letter = c(
    "a", "ab", "c",
    "a", "ab", "c"
  )
)


ylim_vwc <- c(0,
              max(
                emm_vwc_month_group_df$asymp.UCL,
                emm_vwc_season_group_df$asymp.UCL
              )
)

#letter plots
p_vwc_month_letters <- plot_emm_letters(
  emm_vwc_month_group_df,
  letters_vwc_month,
  "month",
  "response",
  "Estimated mean VWC",
  ylim = ylim_vwc,
  xlab = "Month"
)

p_vwc_month_letters

p_vwc_season_letters <- plot_emm_letters(
  emm_vwc_season_group_df,
  letters_vwc_season,
  "season_factor",
  "response",
  "Estimated mean VWC",
  ylim = ylim_vwc,
  xlab = "Austral summer season"
)

p_vwc_season_letters

#################################################
#plots interaction vwc

emm_month_vwc_md <- emmeans(
  vwc_mod_pixel_md,
  ~ foehn | month,
  type = "response"
)


p_vwc_month_md <- plot_emm_interaction(
  emm_month_vwc_md,
  "month",
  "response",
  "Estimated mean VWC",
  ylim = ylim_vwc,
  xlab = "Month"
)
p_vwc_month_md


emm_season_vwc_md <- emmeans(
  vwc_mod_pixel_md,
  ~ foehn | season_factor,
  type = "response"
)

p_vwc_season_md <- plot_emm_interaction(
  emm_season_vwc_md,
  "season_factor",
  "response",
  "Estimated mean VWC",
  ylim = ylim_vwc,
  xlab = "Austral summer season"
)

p_vwc_season_md

# combined plot VWC
vwc_month_md_clean <- clean_for_grid(
  p_vwc_month_md,
  keep_legend = TRUE
)

vwc_month_letters_clean <- clean_for_grid(
  p_vwc_month_letters,
  keep_legend = FALSE
)

vwc_season_md_clean <- clean_for_grid(
  p_vwc_season_md,
  keep_legend = FALSE
)

vwc_season_letters_clean <- clean_for_grid(
  p_vwc_season_letters,
  keep_legend = FALSE
)


vwc_combined <-
  (
    vwc_month_md_clean | vwc_month_letters_clean
  ) /
  (
    vwc_season_md_clean | vwc_season_letters_clean
  ) +
  plot_layout(
    heights = c(1, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )
vwc_combined


vwc_final <-
  
  wrap_elements(
    textGrob(
      "Estimated mean VWC",
      rot = 90,
      gp = gpar(
        fontvwc = 16,
        fontfamily = "Calibri"
      )
    )
  ) +
  vwc_combined +
  
  plot_layout(
    widths = c(0.05, 1)
  )


vwc_final

#save
CairoPNG(
  filename = "plots/vwc_pixel_data_only_black_AOA_final.png",
  width = 4200,
  height = 4200,
  res = 600
)

print(vwc_final)

dev.off()

#############################################################################################


# AlMS data/ Goosef data
ALMS_1 <- read.csv("data/raw_data/VWC_measurements_goosef/knb-lter-mcm.4020.1/mcmlter-soil-alms01_f6-vwc-20220622.csv")
#head(ALMS_1)
ALMS_3 <-read.csv("data/raw_data/VWC_measurements_goosef/knb-lter-mcm.4022.1/mcmlter-soil-alms03_vg-vwc-20220622.csv")
#head(ALMS_3)
ALMS_4 <-read.csv("data/raw_data/VWC_measurements_goosef/knb-lter-mcm.4023.1/mcmlter-soil-alms04_gc-vwc-20220622.csv")
#head(ALMS_4)
ALMS_6 <-read.csv("data/raw_data/VWC_measurements_goosef/knb-lter-mcm.4024.1/mcmlter-soil-alms06_wtb-vwc-20220622.csv")
#head(ALMS_6)
min(ALMS_1$date_time)
max(ALMS_1$date_time)
#ALMS 2 is not located within our AOI

ALMS_full <- left_join(ALMS_1, ALMS_3, by = c("date_time"))

library(sf)
# coordinates from web information
# Station coordinates (WGS84)
ALMS_stations_df <- data.frame(
  station = c("ALMS1", "ALMS3", "ALMS4", "ALMS6"),
  longitude = c(
    163.2530667,
    163.2547333,
    163.0576667,
    162.9785833
  ),
  latitude = c(
    -77.6083333,
    -77.6091667,
    -77.6241833,
    -77.6423333
  )
)

# Convert to sf object (WGS84 geographic coordinates)
alms_stations_sf <- st_as_sf(
  ALMS_stations_df,
  coords = c("longitude", "latitude"),
  crs = 4326
)

st_write(alms_stations_sf, "data/raw_data/VWC_measurements_goosef/alms_staions.gpkg")
alms_stations_sf <- st_read("data/raw_data/VWC_measurements_goosef/alms_staions.gpkg")
#mapview::mapview(alms_stations_sf)
alms_all <- bind_rows(ALMS_1, ALMS_3, ALMS_4, ALMS_6)

alms_sites <- alms_all %>%
  pivot_longer(
    cols = matches("^[a-e]_(5|15)cm$"),
    names_to = c("plot", "depth"),
    names_pattern = "([a-e])_(5|15)cm",
    values_to = "vwc"
  ) %>%
  mutate(
    site_id = paste0(station, "_", plot),
    depth = paste0("vwc_", depth)
  ) %>%
  select(date_time, station, site_id, depth, vwc) %>%
  pivot_wider(
    names_from = depth,
    values_from = vwc
  )

head(alms_sites)

# calculate daily mean
alms_daily_df <- alms_sites %>%
  mutate(
    date_time = parse_date_time(
      date_time,
      orders = c("mdy HM", "mdy HMS")
    ),
    date = as.Date(date_time)
  ) %>%
  filter(
    !is.na(date)
  ) %>%
  group_by(
    station,
    site_id,
    date
  ) %>%
  summarise(
    vwc_5 = mean(vwc_5, na.rm = TRUE),
    vwc_15 = mean(vwc_15, na.rm = TRUE),
    .groups = "drop"
  )
head(alms_daily_df)

quantile(alms_daily_df$vwc_5, na.rm =  T)
#    25%         50%         75%        
#0.030000000 0.042020833 0.064000000 
#median(alms_daily_df$vwc_5, na.rm = T)
# only for prediction dates

alms_prediction_df <- alms_daily_df %>%
  filter(date %in% as.Date(prediction_dates))
quantile(alms_prediction_df$vwc_5, na.rm = T)
#0%        25%        50%        75%       100% 
#0.02770833 0.04124632 0.07593655 0.19449169 0.40681250 
#compare to prediction data in this area
alms_buffer_sf <- st_buffer(alms_stations_sf, 50)
pred_VWC_alms_50 <- map(prediction_list, ~terra::extract(.x$stack, alms_buffer_sf))
pred_VWC_alms_50_vals <- pred_VWC_alms_50 |>
  map(\(x) {
    x |>
      dplyr::select(-ID) |>
      unlist(use.names = FALSE)
  }) |>
  unlist(use.names = FALSE)

quantile(pred_VWC_alms_50_vals, na.rm = T)
#0%       25%       50%       75%      100% 
#0.1070685 0.1499152 0.1670341 0.1920400 0.3527313 
#median(alms_prediction_df$vwc_5, na.rm = T)
#[1] 0.07593655 
run_alms_event_analysis <- function(data, response){
  
  response <- rlang::ensym(response)
  data <- data %>%
    mutate(
      month = lubridate::month(date),
      year = lubridate::year(date),
      season_year = if_else(month >= 11, year + 1L, year)
    )
  data %>%
    group_by(site_id, station, season_year) %>%
    group_modify(~{
      
      dat <- arrange(.x, date)
      
      dat$vwc <- dplyr::pull(dat, !!response)
      
      dat <- dat[!is.na(dat$vwc), ]
      
      if(nrow(dat) < 8)
        return(tibble())
      
      season_mean <- mean(dat$vwc)
      
      season_sd <- sd(dat$vwc)
      
      first_date <- min(dat$date)
      
      dat$date_num <- as.numeric(dat$date - first_date)
      
      z <- tryCatch(
        zyp.trend.vector(
          y = dat$vwc,
          x = dat$date_num,
          method = "zhang",
          conf.intervals = TRUE
        ),
        error = function(e) NULL
      )
      
      if(is.null(z)){
        
        dat$reference <- season_mean
        
      } else{
        
        trend <- unname(z[["trend"]])
        intercept <- unname(z[["intercept"]])
        p_value <- unname(z[["sig"]])
        
        trend_total <- abs(trend) * max(dat$date_num)
        
        relevant <- trend_total > season_sd
        
        first_value <- dat$vwc[1]
        
        trend_pred <- first_value +
          trend * dat$date_num
        
        trend_error <- median(
          abs(dat$vwc - trend_pred),
          na.rm = TRUE
        )
        
        mean_error <- median(
          abs(dat$vwc - season_mean),
          na.rm = TRUE
        )
        
       
        improvement_all <- mean_error - trend_error
        
        use_trend <-
          !is.na(p_value) &&
          p_value < 0.05 &&
          relevant &&
          improvement_all > 0.03
        
        dat$reference <-
          if(use_trend) trend_pred else season_mean
        
        dat$trend <- trend
        dat$intercept <- intercept
        dat$p_value <- p_value
        dat$trend_relevant <- relevant
        dat$trend_error <- trend_error
        dat$mean_error <- mean_error
        dat$improvement_all <- improvement_all
        dat$use_trend <- use_trend
      }
      
      dat$event_effect <- dat$vwc - dat$reference
      
      dat$event_effect[
        abs(dat$event_effect) < 0.03
      ] <- 0
      
      dat$event_effect[
        abs(dat$event_effect) < season_sd
      ] <- 0
      
      dat$season_mean <- season_mean
      dat$season_sd <- season_sd
      
      dat
      
    }) %>%
    ungroup()
}

alms_event_df_5cm <- run_alms_event_analysis(
  alms_prediction_df,
  vwc_5
)

#alms_event_df_15cm <- run_alms_event_analysis(
#  alms_prediction_df,
#  vwc_15
#)

alms_event_df_5cm$date <- as.character(alms_event_df_5cm$date)

alms_event_df_5cm <- alms_event_df_5cm |>
  left_join(
    foehn_info_df,
    by = "date"
  )


alms_event_df_5cm <- alms_event_df_5cm |>
  mutate(
    event_occurs = event_effect != 0
  )

alms_event_df_5cm$month <- as.factor(alms_event_df_5cm$month)

alms_event_df_5cm$season_year <- alms_event_df_5cm$season_year.x

# calculate new quartiles because there could be an error due to calibration offsets. 

alms_q_5 <- quantile(
  alms_event_df_5cm$vwc,
  probs = c(0.25, 0.50, 0.75),
  na.rm = TRUE
)

alms_event_df_5cm <- alms_event_df_5cm %>%
  mutate(
    global_vwc_class = cut(
      vwc,
      breaks = c(-Inf, alms_q_5, Inf),
      labels = c(
        "dry",
        "intermediate_dry",
        "intermediate_wet",
        "wet"
      ),
      include.lowest = TRUE,
      ordered_result = TRUE
    )
  )

class_check_season_5 <- alms_event_df_5cm %>%
  count(season_year, global_vwc_class) %>%
  complete(
    season_year,
    vwc_class = levels(alms_event_df_5cm$vwc_class),
    fill = list(n = 0)
  ) %>%
  arrange(season_year, vwc_class)

class_check_season_5
alms_q_5
global_q

ggplot(
  alms_event_df_5cm,
  aes(
    x = date,
    y = vwc,
    colour = site_id,
    group = site_id
  )
) +
  geom_line() +
  geom_point(size = 1) +
  facet_grid(
    season_year ~ station,
    scales = "free_x"
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "VWC (5 cm)",
    colour = "Site"
  )
ggplot(
  alms_event_df_5cm,
  aes(
    x = site_id,
    y = vwc
  )
) +
  geom_boxplot(outlier.size = 0.8) +
  facet_grid(
    season_year ~ station
  ) +
  theme_classic() +
  labs(
    x = "Site",
    y = "VWC (5 cm)"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )
# on first sight there are no obvious outliers. let's check the model results and 
# go back to this if neccessary

write.csv(alms_event_df_5cm, "data/processed_data/statistic_models_final/ALMS_test/alms_event_df_5cm.csv")
alms_event_df_5cm <- read.csv("data/processed_data/statistic_models_final/ALMS_test/alms_event_df_5cm.csv")
#models for ALMS data

#we don't have any foehn in November
alms_pixel_12_df <- alms_event_df_5cm[alms_event_df_5cm$month != 11,]
alms_pixel_12_df$month <- as.factor(alms_pixel_12_df$month)
alms_pixel_12_df$month <- droplevels(alms_pixel_12_df$month)
alms_pixel_12_df$season_factor <- as.factor(alms_pixel_12_df$season_year)

alms_mod_occ <- glmmTMB(
  event_occurs ~ foehn +season_factor + global_vwc_class + foehn*month +(1|site_id),
  dispformula = ~ season_factor,
  data = alms_pixel_12_df,
  family = binomial()
)

#saveRDS(sp_mod_occ_md, "data/processed_data/statistic_models_final/sp_mod_occ_md.RDS")

summary(alms_mod_occ) 

car::Anova(alms_mod_occ, type = "III")

alms_anova_occ_df <- build_anova_df(alms_mod_occ)
write.csv(alms_anova_occ_df, "data/processed_data/statistic_models_final/ALMS_test/alms_anova_occ_df.csv", row.names = F)

alms_mod_occ_simresmd <- simulateResiduals(alms_mod_occ)
testUniformity(alms_mod_occ_simresmd) #not sig
testDispersion(alms_mod_occ_simresmd) #not sig

#emm_foehn_occ <- emmeans(
#  alms_mod_occ,
#  ~ foehn,
#  type = "response"
#)

#emm_foehn_occ
#pairs(emm_foehn)
# mod_size wg
#filter to event effect != 0
alms_size_12_df <- alms_pixel_12_df %>%
  filter(event_effect != 0)


#alms_size_12_df <- alms_size_12_df[alms_size_12_df$season_factor != 2017,]
#alms_size_12_df$season_factor <- droplevels(alms_size_12_df$season_factor)

alms_mod_size <- glmmTMB(
  event_effect ~ foehn + season_factor + global_vwc_class+ month + (1|site_id),
  dispformula = ~ season_factor,
  data = alms_size_12_df,
  family = gaussian()
)

#saveRDS(sp_mod_size_md, "data/processed_data/statistic_models_final/sp_mod_size_md.RDS")

summary(alms_mod_size)

car::Anova(alms_mod_size, type = "II")

alms_anova_size_df <- build_anova_df(alms_mod_size, t= "II")
write.csv(alms_anova_size_df, "data/processed_data/statistic_models_final/ALMS_test/alms_anova_size_df.csv", row.names = F)

alms_mod_size_simresmd <- simulateResiduals(alms_mod_size)
testUniformity(alms_mod_size_simresmd) #ok
testDispersion(alms_mod_size_simresmd) #ok

# vwc month(m) 

alms_vwc_mod <- glmmTMB(vwc_5 ~ foehn*season_factor+month  + (1|site_id),
                            dispformula = ~ season_factor,
                            data = alms_pixel_12_df,
                            family = beta_family()
)

summary(alms_vwc_mod)
car::Anova(alms_vwc_mod, type = "III")

alms_anova_vwc_df <- build_anova_df(alms_vwc_mod)
write.csv(alms_anova_vwc_df, "data/processed_data/statistic_models_final/ALMS_test/alms_anova_vwc_df.csv", row.names = F)


#saveRDS(vwc_mod_pixel_md, "data/processed_data/statistic_models_final/sp_mod_pixel_m.RDS")

alms_mod_vwc_simresmd <- simulateResiduals(alms_vwc_mod)
testUniformity(alms_mod_vwc_simresmd) #not sig
testDispersion(alms_mod_vwc_simresmd) #not sig

#emm_foehn <- emmeans(
#  alms_vwc_mod,
#  ~ foehn,
#  type = "response"
#)

#emm_foehn
#pairs(emm_foehn) #0.01; sensitivity of this sensor is 0.03-no difference!  

#as foehn as main effect is not significant in any case while other main effects are sig, we just report
#it but don't plot this. 

############################################################################################
#plots

#plot main effect
library(multcomp)
plot_emm_main_effect <- function(
    emm,
    letters_df,
    x,
    y,
    ylab,
    xlab = x,
    ylim = NULL,
    LCL = "lower.CL",
    UCL = "upper.CL"
) {
  
  df <- emm |>
    as.data.frame() |>
    tibble::as_tibble()
  
  letters_df <- letters_df |>
    as.data.frame() |>
    tibble::as_tibble()
  
  # x as character for joining
  df <- df |>
    dplyr::mutate(
      x_value = as.character(.data[[x]])
    )
  
  letters_df <- letters_df |>
    dplyr::mutate(
      x_value = as.character(.data[[x]])
    )
  
  # Join grouping letters
  df <- df |>
    dplyr::left_join(
      letters_df |>
        dplyr::select(
          x_value,
          .group
        ),
      by = "x_value"
    )
  
  # sort
  if (x == "month") {
    
    df <- df |>
      dplyr::mutate(
        x_value = factor(
          x_value,
          levels = c("12", "1", "2")
        )
      )
    
  } else if (x == "season_factor") {
    
    df <- df |>
      dplyr::mutate(
        x_value = factor(
          x_value,
          levels = c(
            "2016",
            "2017",
            "2019",
            "2021",
            "2022",
            "2023",
            "2024",
            "2025",
            "2026"
          )
        )
      )
    
  } else if (x == "global_vwc_class") {
    
    df <- df |>
      dplyr::mutate(
        x_value = factor(
          x_value,
          levels = c(
            "dry",
            "intermediate_dry",
            "intermediate_wet",
            "wet"
          )
        )
      )
  }
  
  # Position of letters
  y_range <- diff(
    range(
      c(
        df[[LCL]],
        df[[UCL]]
      ),
      na.rm = TRUE
    )
  )
  
  if (y_range == 0) {
    y_range <- 0.01
  }
  
  letters_y <- max(
    df[[UCL]],
    na.rm = TRUE
  ) + 0.10 * y_range
  
  df <- df |>
    dplyr::mutate(
      letters_y = letters_y
    )
  
  # Y limits
  if (is.null(ylim)) {
    
    ymin <- min(
      df[[LCL]],
      na.rm = TRUE
    )
    
    ymax <- max(
      df[[UCL]],
      na.rm = TRUE
    )
    
    ylim_use <- c(
      ymin,
      ymax + 0.25 * (ymax - ymin)
    )
    
  } else {
    
    ylim_use <- ylim
  }
  
  # X-axis labels
  x_labels <- waiver()
  
  if (x == "month") {
    
    x_labels <- c(
      "12" = "Dec",
      "1" = "Jan",
      "2" = "Feb"
    )
    
  } else if (x == "season_factor") {
    
    x_labels <- c(
      "2016" = "15|16",
      "2017" = "16|17",
      "2019" = "18|19",
      "2021" = "20|21",
      "2022" = "21|22",
      "2023" = "22|23",
      "2024" = "23|24",
      "2025" = "24|25",
      "2026" = "25|26"
    )
    
  } else if (x == "global_vwc_class") {
    
    x_labels <- c(
      "dry" = "dry",
      "intermediate_dry" = "i_dry",
      "intermediate_wet" = "i_wet",
      "wet" = "wet"
    )
  }
  
  ggplot(
    df,
    aes(
      x = x_value,
      y = .data[[y]]
    )
  ) +
    geom_errorbar(
      aes(
        ymin = .data[[LCL]],
        ymax = .data[[UCL]]
      ),
      width = 0.12,
      linewidth = 0.6
    ) +
    geom_point(
      size = 3
    ) +
    geom_text(
      aes(
        y = letters_y,
        label = .group
      ),
      size = 5,
      family = "Calibri",
      vjust = 0
    ) +
    scale_x_discrete(
      labels = x_labels
    ) +
    coord_cartesian(
      ylim = ylim_use,
      clip = "off"
    ) +
    labs(
      x = xlab,
      y = ylab
    ) +
    theme_bw() +
    theme(
      text = element_text(
        family = "Calibri",
        size = 14
      ),
      axis.text.x = element_text(
        size = 14
      ),
      axis.text.y = element_text(
        size = 14
      ),
      axis.title = element_text(
        size = 14
      ),
      plot.margin = margin(
        5.5, 5.5, 15, 5.5
      )
    )
}

# alms vwc

emm_alms_vwc_season <- emmeans(
  alms_vwc_mod,
  ~ foehn | season_factor,
  type = "response"
)

emm_alms_vwc_season_df <- summary(
  emm_alms_vwc_season
) |>
  as_tibble()

# Foehn effect within each season
delta_alms_vwc_season <- emm_alms_vwc_season_df |>
  dplyr::group_by(season_factor) |>
  dplyr::reframe(
    `Foehn - No foehn` =
      response[foehn == TRUE] -
      response[foehn == FALSE],
    .groups = "drop"
  )

delta_alms_vwc_season

# Pairwise comparisons
pairs_alms_vwc_season <- pairs(
  emm_alms_vwc_season,
  reverse = TRUE
)

pairs_alms_vwc_season


#month sig main effect

emm_alms_vwc_month <- emmeans(
  alms_vwc_mod,
  ~ month,
  type = "response"
)

emm_alms_vwc_month_df <- summary(
  emm_alms_vwc_month
) |>
  as_tibble()

# Differences between months
delta_alms_vwc_month <- emm_alms_vwc_month_df |>
  dplyr::reframe(
    `12 - 1` =
      response[month == "12"] -
      response[month == "1"],
    `12 - 2` =
      response[month == "12"] -
      response[month == "2"],
    `1 - 2` =
      response[month == "1"] -
      response[month == "2"]
  )

delta_alms_vwc_month

# Pairwise comparisons
pairs_alms_vwc_month <- pairs(
  emm_alms_vwc_month
)

pairs_alms_vwc_month


# alms occurrence

emm_alms_occ_season <- emmeans(
  alms_mod_occ,
  ~ season_factor,
  type = "response"
)

emm_alms_occ_season_df <- summary(
  emm_alms_occ_season
) |>
  as_tibble()

# Differences between seasons
delta_alms_occ_season <- emm_alms_occ_season_df |>
  dplyr::reframe(
    `2016 - 2017` =
      prob[season_factor == "2016"] -
      prob[season_factor == "2017"],
    `2017 - 2019` =
      prob[season_factor == "2017"] -
      prob[season_factor == "2019"],
    `2019 - 2021` =
      prob[season_factor == "2019"] -
      prob[season_factor == "2021"]
    
  )

delta_alms_occ_season

pairs_alms_occ_season <- pairs(
  emm_alms_occ_season
)

pairs_alms_occ_season
#all other not sig

#alms occurrence vwc

emm_alms_occ_vwc <- emmeans(
  alms_mod_occ,
  ~ global_vwc_class,
  type = "response"
)

emm_alms_occ_vwc_df <- summary(
  emm_alms_occ_vwc
) |>
  as_tibble()

# Pairwise differences
delta_alms_occ_vwc <- emm_alms_occ_vwc_df |>
  dplyr::reframe(
    `dry - intermediate-dry` =
      prob[global_vwc_class == "dry"] -
      prob[global_vwc_class == "intermediate_dry"],
    `dry - intermediate-wet` =
      prob[global_vwc_class == "dry"] -
      prob[global_vwc_class == "intermediate_wet"],
    `dry - wet` =
      prob[global_vwc_class == "dry"] -
      prob[global_vwc_class == "wet"],
    `intermediate-dry - intermediate-wet` =
      prob[global_vwc_class == "intermediate_dry"] -
      prob[global_vwc_class == "intermediate_wet"],
    `intermediate-dry - wet` =
      prob[global_vwc_class == "intermediate_dry"] -
      prob[global_vwc_class == "wet"],
    `intermediate-wet - wet` =
      prob[global_vwc_class == "intermediate_wet"] -
      prob[global_vwc_class == "wet"]
  )

delta_alms_occ_vwc

pairs_alms_occ_vwc <- pairs(
  emm_alms_occ_vwc
)

pairs_alms_occ_vwc


# alms occ month

emm_alms_occ_month <- emmeans(
  alms_mod_occ,
  ~ month,
  type = "response"
)

emm_alms_occ_month_df <- summary(
  emm_alms_occ_month
) |>
  as_tibble()

delta_alms_occ_month <- emm_alms_occ_month_df |>
  dplyr::reframe(
    `12 - 1` =
      prob[month == "12"] -
      prob[month == "1"],
    `12 - 2` =
      prob[month == "12"] -
      prob[month == "2"],
    `1 - 2` =
      prob[month == "1"] -
      prob[month == "2"]
  )

delta_alms_occ_month

pairs_alms_occ_month <- pairs(
  emm_alms_occ_month
)

pairs_alms_occ_month


# alms anomaly magnitude

# ALMS size season

emm_alms_size_season <- emmeans(
  alms_mod_size,
  ~ season_factor
)

emm_alms_size_season_df <- summary(
  emm_alms_size_season
) |>
  as_tibble()

delta_alms_size_season <- emm_alms_size_season_df |>
  dplyr::reframe(
    `2016 - 2017` =
      emmean[season_factor == "2016"] -
      emmean[season_factor == "2017"],
    `2017 - 2021` =
      emmean[season_factor == "2017"] -
      emmean[season_factor == "2021"],
    `2019 - 2021` =
      emmean[season_factor == "2019"] -
      emmean[season_factor == "2021"]
  )

delta_alms_size_season

pairs_alms_size_season <- pairs(
  emm_alms_size_season
)

pairs_alms_size_season

# alms size vwc class

emm_alms_size_vwc <- emmeans(
  alms_mod_size,
  ~ global_vwc_class
)

emm_alms_size_vwc_df <- summary(
  emm_alms_size_vwc
) |>
  as_tibble()

delta_alms_size_vwc <- emm_alms_size_vwc_df |>
  dplyr::reframe(
    `dry - intermediate-dry` =
      emmean[global_vwc_class == "dry"] -
      emmean[global_vwc_class == "intermediate_dry"],
    `dry - intermediate-wet` =
      emmean[global_vwc_class == "dry"] -
      emmean[global_vwc_class == "intermediate_wet"],
    `dry - wet` =
      emmean[global_vwc_class == "dry"] -
      emmean[global_vwc_class == "wet"],
    `intermediate-dry - intermediate-wet` =
      emmean[global_vwc_class == "intermediate_dry"] -
      emmean[global_vwc_class == "intermediate_wet"],
    `intermediate-dry - wet` =
      emmean[global_vwc_class == "intermediate_dry"] -
      emmean[global_vwc_class == "wet"],
    `intermediate-wet - wet` =
      emmean[global_vwc_class == "intermediate_wet"] -
      emmean[global_vwc_class == "wet"]
  )

delta_alms_size_vwc

pairs_alms_size_vwc <- pairs(
  emm_alms_size_vwc
)

pairs_alms_size_vwc


# alms size month

emm_alms_size_month <- emmeans(
  alms_mod_size,
  ~ month
)

emm_alms_size_month_df <- summary(
  emm_alms_size_month
) |>
  as_tibble()

delta_alms_size_month <- emm_alms_size_month_df |>
  dplyr::reframe(
    `12 - 1` =
      emmean[month == "12"] -
      emmean[month == "1"],
    `12 - 2` =
      emmean[month == "12"] -
      emmean[month == "2"],
    `1 - 2` =
      emmean[month == "1"] -
      emmean[month == "2"]
  )

delta_alms_size_month

pairs_alms_size_month <- pairs(
  emm_alms_size_month
)

pairs_alms_size_month


# ALMS plots
ylim_alms_vwc_season <- range(
  c(
    emm_alms_vwc_season_df$asymp.LCL,
    emm_alms_vwc_season_df$asymp.UCL
  ),
  na.rm = TRUE
)

p_alms_vwc_season <- plot_emm_interaction(
  emm_alms_vwc_season,
  x = "season_factor",
  y = "response",
  ylab = "Estimated VWC",
  xlab = "Season",
  ylim = ylim_alms_vwc_season,
  LCL = "asymp.LCL",
  UCL = "asymp.UCL"
)


letters_alms_vwc_month <- tibble::tibble(
  month = factor(c("12", "1", "2"), levels = c("12", "1", "2")),
  .group = c("a", "a", "a")
)

p_alms_vwc_month <- plot_emm_main_effect(
  emm = emm_alms_vwc_month,
  letters_df = letters_alms_vwc_month,
  x = "month",
  y = "response",
  ylab = "Estimated VWC",
  xlab = "Month",
  LCL = "asymp.LCL",
  UCL = "asymp.UCL"
)

#occ

letters_alms_occ_season <- tibble::tibble(
  season_factor = factor(
    c("2016", "2017", "2019", "2021"),
    levels = c("2016", "2017", "2019", "2021")
  ),
  .group = c("a", "a", "a", "a")
)

p_alms_occ_season <- plot_emm_main_effect(
  emm = emm_alms_occ_season,
  letters_df = letters_alms_occ_season,
  x = "season_factor",
  y = "prob",
  ylab = "Estimated anomaly occurrence probability",
  xlab = "Season",
  LCL = "asymp.LCL",
  UCL = "asymp.UCL"
)


letters_alms_occ_vwc_class <- tibble::tibble(
  global_vwc_class = factor(
    c("dry", "intermediate_dry", "intermediate_wet", "wet"),
    levels = c("dry", "intermediate_dry", "intermediate_wet", "wet")
  ),
  .group = c("a", "a", "a", "a")
)

p_alms_occ_vwc <- plot_emm_main_effect(
  emm = emm_alms_occ_vwc,
  letters_df = letters_alms_occ_vwc_class,
  x = "global_vwc_class",
  y = "prob",
  ylab = "Estimated anomaly occurrence probability",
  xlab = "Global VWC class",
  LCL = "asymp.LCL",
  UCL = "asymp.UCL"
)


letters_alms_occ_month <- tibble::tibble(
  month = factor(c("12", "1", "2"), levels = c("12", "1", "2")),
  .group = c("a", "b", "ab")
)

p_alms_occ_month <- plot_emm_main_effect(
  emm = emm_alms_occ_month,
  letters_df = letters_alms_occ_month,
  x = "month",
  y = "prob",
  ylab = "Estimated anomaly occurrence probability",
  xlab = "Month",
  LCL = "asymp.LCL",
  UCL = "asymp.UCL"
)

#size
letters_alms_size_season <- tibble::tibble(
  season_factor = factor(
    c("2016", "2017", "2019", "2021"),
    levels = c("2016", "2017", "2019", "2021")
  ),
  .group = c("a", "a", "a", "a")
)

p_alms_size_season <- plot_emm_main_effect(
  emm = emm_alms_size_season,
  letters_df = letters_alms_size_season,
  x = "season_factor",
  y = "emmean",
  ylab = "Estimated VWC anomaly size",
  xlab = "Season",
  LCL = "lower.CL",
  UCL = "upper.CL"
)


letters_alms_size_vwc_class <- tibble::tibble(
  global_vwc_class = factor(
    c("dry", "intermediate_dry", "intermediate_wet", "wet"),
    levels = c("dry", "intermediate_dry", "intermediate_wet", "wet")
  ),
  .group = c("a", "ab", "c", "d")
)

p_alms_size_vwc <- plot_emm_main_effect(
  emm = emm_alms_size_vwc,
  letters_df = letters_alms_size_vwc_class,
  x = "global_vwc_class",
  y = "emmean",
  ylab = "Estimated VWC anomaly size",
  xlab = "Global VWC class",
  LCL = "lower.CL",
  UCL = "upper.CL"
)



letters_alms_size_month <- tibble::tibble(
  month = factor(c("12", "1", "2"), levels = c("12", "1", "2")),
  .group = c("a", "b", "b")
)

p_alms_size_month <- plot_emm_main_effect(
  emm = emm_alms_size_month,
  letters_df = letters_alms_size_month,
  x = "month",
  y = "emmean",
  ylab = "Estimated VWC anomaly size",
  xlab = "Month",
  LCL = "lower.CL",
  UCL = "upper.CL"
)


p_alms_vwc_season
p_alms_vwc_month

p_alms_occ_season
p_alms_occ_vwc
p_alms_occ_month

p_alms_size_season
p_alms_size_vwc
p_alms_size_month

#save

ggsave(
  "plots/alms/alms_vwc_season.png",
  p_alms_vwc_season,
  width = 8.3,
  height = 6,
  dpi = 300
)

ggsave(
  "plots/alms/alms_vwc_month.png",
  p_alms_vwc_month,
  width = 8.3,
  height = 6,
  dpi = 300
)

p_alms_occ_season_grid <- clean_for_grid(
  p_alms_occ_season
)

p_alms_occ_vwc_grid <- clean_for_grid(
  p_alms_occ_vwc
)

p_alms_occ_month_grid <- clean_for_grid(
  p_alms_occ_month
)


# combined plots

p_alms_occ_combined <- (
  p_alms_occ_season_grid /
    p_alms_occ_vwc_grid /
    p_alms_occ_month_grid
)



p_alms_occ_combined <- patchwork::wrap_plots(
  p_alms_occ_month_grid,
  p_alms_occ_season_grid,
  p_alms_occ_vwc_grid,
  ncol = 1
) +
  patchwork::plot_annotation(
    theme = theme(
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    )
  )

occ_alms_combined <-
  
  (
    p_alms_occ_month_grid
  ) /
  (
    p_alms_occ_season_grid
  ) /
  (
    p_alms_occ_vwc_grid
  ) +
  
  plot_layout(
    guides = "collect"
  ) &
  
  theme(
    legend.position = "bottom"
  )


occ_alms_final <-
  
  wrap_elements(
    textGrob(
      "Estimated anomaly occurrence probability",
      rot = 90,
      gp = gpar(
        fontsize = 16,
        fontfamily = "Calibri"
      )
    )
  ) +
  occ_alms_combined +
  
  plot_layout(
    widths = c(0.05, 1)
  )
occ_alms_final



p_alms_size_season_grid <- clean_for_grid(
  p_alms_size_season
)

p_alms_size_vwc_grid <- clean_for_grid(
  p_alms_size_vwc
)

p_alms_size_month_grid <- clean_for_grid(
  p_alms_size_month
)


size_alms_combined <-
  
  (
    p_alms_size_month_grid
  ) /
  (
    p_alms_size_season_grid
  ) /
  (
    p_alms_size_vwc_grid
  ) +
  
  plot_layout(
    guides = "collect"
  ) &
  
  theme(
    legend.position = "bottom"
  )


size_alms_final <-
  
  wrap_elements(
    textGrob(
      "Estimated mean anomaly size",
      rot = 90,
      gp = gpar(
        fontsize = 16,
        fontfamily = "Calibri"
      )
    )
  ) +
  size_alms_combined +
  
  plot_layout(
    widths = c(0.05, 1)
  )
size_alms_final



ggsave(
  "plots/alms/alms_occ_combined.png",
  occ_alms_final,
  width = 8.3,
  height = 10,
  dpi = 300
)

ggsave(
  "plots/alms/alms_size_combined.png",
  size_alms_final,
  width = 8.3,
  height = 10,
  dpi = 300
)

##########################################################################################
############################################################################################
# prediction data model with event effect > 0.03
#we don't have any foehn in November
pixel_12_df <- pixel_df[pixel_df$month != 11,]

pixel_12_df$month <- as.factor(pixel_12_df$month)
pixel_12_df$month <- droplevels(pixel_12_df$month)

pixel_12_df$season_factor <- as.factor(pixel_12_df$season_factor)

pixel_12_df_03 <- pixel_12_df
pixel_12_df_03$event_effect[abs(pixel_12_df_03$event_effect) <= 0.03] <- 0

sp_03_mod_occ_md <- glmmTMB(
  event_occurs ~ foehn * season_factor + foehn * global_vwc_class + foehn*month +(1|pixel_id),
  dispformula = ~ season_factor,
  data = pixel_12_df_03,
  family = binomial()
)

saveRDS(sp_03_mod_occ_md, "data/processed_data/statistic_models_final/models_event_effect03/sp_03_mod_occ_md.RDS")
#sp_03_mod_occ_md <- readRDS("data/processed_data/statistic_models_final/models_event_effect03/sp_03_mod_occ_md.RDS")

sp_mod_occ_simres_03 <- simulateResiduals(sp_03_mod_occ_md)
saveRDS(sp_mod_occ_simres_03, "data/processed_data/statistic_models_final/models_event_effect03/sp_mod_occ_simres_03.RDS")

testUniformity(sp_mod_occ_simres_03)
testDispersion(sp_mod_occ_simres_03)

summary(sp_03_mod_occ_md) 

car::Anova(sp_03_mod_occ_md, type = "III")
anova_df_occ_03 <- build_anova_df(sp_03_mod_occ_md)
write.csv(anova_df_occ_03, "data/processed_data/statistic_models_final/models_event_effect03/anova_df_occ_03.csv", row.names = F)

# mod_size wg
#filter to event effect != 0
size_12_df_03 <- pixel_12_df_03 %>%
  filter(event_effect != 0)

sp_03_mod_size_md <- glmmTMB(
  event_effect ~ foehn * season_factor + foehn*global_vwc_class +foehn*month + (1|pixel_id),
  dispformula = ~ season_factor,
  data = size_12_df_03,
  family = gaussian()
)
diagnose(sp_03_mod_size_md)


saveRDS(sp_03_mod_size_md, "data/processed_data/statistic_models_final/models_event_effect03/sp_03_mod_size_md.RDS")
#sp_03_mod_size_md <- readRDS("data/processed_data/statistic_models_final/models_event_effect03/sp_03_mod_size_md.RDS")

sp_mod_size_simres_03 <- simulateResiduals(sp_03_mod_size_md)
saveRDS(sp_mod_size_simres_03, "data/processed_data/statistic_models_final/models_event_effect03/sp_mod_size_simres_03.RDS")

testUniformity(sp_mod_size_simres_03)
testDispersion(sp_mod_size_simres_03)

summary(sp_03_mod_size_md)

car::Anova(sp_03_mod_size_md, type = "III")
anova_df_size_03 <- build_anova_df(sp_03_mod_size_md)
write.csv(anova_df_size_03, "data/processed_data/statistic_models_final/models_event_effect03/anova_df_size_03.csv", row.names = F)


# vwc month(m) 
# not relevant, same as for non filtered event effect data, was vwc was not affected
vwc_mod_pixel_md_03 <- glmmTMB(vwc ~ foehn*season_factor + foehn*month + (1|pixel_id),
                            dispformula = ~ season_factor,
                            data = pixel_12_df_03,
                            family = beta_family()
)

summary(vwc_mod_pixel_md_03)
car::Anova(vwc_mod_pixel_md_03, type = "III")
saveRDS(vwc_mod_pixel_md_03, "data/processed_data/statistic_models_final/models_event_effect03/sp_03_mod_pixel_m.RDS")

diagnose(vwc_mod_pixel_md_03)
mod_checkvwcmd_03 <- glmmTMB(
  vwc ~ foehn * season_factor+ foehn*month +(1|pixel_id),
  dispformula = ~ season_factor,
  data = sub_12_df_03,
  family = beta_family()
)
sp_mod_vwc_simresmd_03 <- simulateResiduals(mod_checkvwcmd_03)
testUniformity(sp_mod_vwc_simresmd_03) #not sig
testDispersion(sp_mod_vwc_simresmd_03) #not sig

##########################################################################################
#prepare plots
#month
emm_occ_month_group_03 <- emmeans(
  sp_03_mod_occ_md,
  ~ month | foehn,
  type = "response"
)

emm_occ_month_group_df_03 <- summary(emm_occ_month_group_03) |> as_tibble()

delta_occ_month_group_03 <- emm_occ_month_group_df_03 |>
  dplyr::group_by(foehn) |>
  dplyr::reframe(
    `12 - 1` = prob[month == "12"] - prob[month == "1"],
    `12 - 2` = prob[month == "12"] - prob[month == "2"],
    `1 - 2`  = prob[month == "1"]  - prob[month == "2"],
    .groups = "drop"
  )
delta_occ_month_group_03
pairs(emm_occ_month_group_03)
# season
emm_occ_season_group_03 <- emmeans(
  sp_03_mod_occ_md,
  ~ season_factor | foehn,
  type = "response"
)

emm_occ_season_group_df_03 <- summary(emm_occ_season_group_03) |> as_tibble()

delta_occ_season_group_03 <- emm_occ_season_group_df_03 |>
  dplyr::select(season_factor, foehn, prob) |>
  tidyr::pivot_wider(
    names_from = season_factor,
    values_from = prob
  ) |>
  dplyr::mutate(
    `2016 - 2017` = `2016` - `2017`,
    `2016 - 2019` = `2016` - `2019`,
    `2016 - 2021` = `2016` - `2021`,
    `2016 - 2022` = `2016` - `2022`,
    `2016 - 2023` = `2016` - `2023`,
    `2016 - 2024` = `2016` - `2024`,
    `2016 - 2025` = `2016` - `2025`,
    `2016 - 2026` = `2016` - `2026`,
    
    `2017 - 2019` = `2017` - `2019`,
    `2017 - 2021` = `2017` - `2021`,
    `2017 - 2022` = `2017` - `2022`,
    `2017 - 2023` = `2017` - `2023`,
    `2017 - 2024` = `2017` - `2024`,
    `2017 - 2025` = `2017` - `2025`,
    `2017 - 2026` = `2017` - `2026`,
    
    `2019 - 2021` = `2019` - `2021`,
    `2019 - 2022` = `2019` - `2022`,
    `2019 - 2023` = `2019` - `2023`,
    `2019 - 2024` = `2019` - `2024`,
    `2019 - 2025` = `2019` - `2025`,
    `2019 - 2026` = `2019` - `2026`,
    
    `2021 - 2022` = `2021` - `2022`,
    `2021 - 2023` = `2021` - `2023`,
    `2021 - 2024` = `2021` - `2024`,
    `2021 - 2025` = `2021` - `2025`,
    `2021 - 2026` = `2021` - `2026`,
    
    `2022 - 2023` = `2022` - `2023`,
    `2022 - 2024` = `2022` - `2024`,
    `2022 - 2025` = `2022` - `2025`,
    `2022 - 2026` = `2022` - `2026`,
    
    `2023 - 2024` = `2023` - `2024`,
    `2023 - 2025` = `2023` - `2025`,
    `2023 - 2026` = `2023` - `2026`,
    
    `2024 - 2025` = `2024` - `2025`,
    `2024 - 2026` = `2024` - `2026`,
    
    `2025 - 2026` = `2025` - `2026`
  )

delta_occ_season_group_03
pairs(emm_occ_season_group_03)

# letters occ foehn T
library(multcompView)

pair_df_true_03 <- summary(pairs(emm_occ_season_group_03)) |>
  dplyr::filter(foehn == TRUE) |>
  dplyr::mutate(
    contrast = contrast |>
      gsub("season_factor", "", x = _) |>
      gsub(" / ", " - ", x = _)
  ) |>
  dplyr::left_join(
    delta_occ_season_group_03 |>
      dplyr::filter(foehn == TRUE) |>
      tidyr::pivot_longer(
        cols = -foehn,
        names_to = "contrast",
        values_to = "delta"
      ),
    by = "contrast"
  ) |>
  dplyr::mutate(
    different = p.value < 0.05 & abs(delta) >= 0.01
  )
sum(is.na(pair_df_true_03$delta))

x_true_03 <- pair_df_true_03$different

names(x_true_03) <- pair_df_true_03$contrast |>
  gsub(" - ", "-", x = _)

letters_occ_true_03 <- multcompView::multcompLetters(x_true_03)$Letters

#letters occ foehn F
# occ season FALSE

pair_df_occ_false_03 <- summary(pairs(emm_occ_season_group_03)) |>
  dplyr::filter(foehn == FALSE) |>
  dplyr::mutate(
    contrast = contrast |>
      gsub("season_factor", "", x = _) |>
      gsub(" / ", " - ", x = _)
  ) |>
  dplyr::left_join(
    delta_occ_season_group_03 |>
      dplyr::filter(foehn == FALSE) |>
      tidyr::pivot_longer(
        cols = -foehn,
        names_to = "contrast",
        values_to = "delta"
      ),
    by = "contrast"
  ) |>
  dplyr::mutate(
    different = p.value < 0.05 & abs(delta) >= 0.01
  )


x_occ_false_03 <- pair_df_occ_false_03$different

names(x_occ_false_03) <- pair_df_occ_false_03$contrast |>
  gsub(" - ", "-", x = _)

str(multcompView::multcompLetters(x_occ_false_03)$Letters)
letters_occ_false_03 <-multcompView::multcompLetters(x_occ_false_03)$Letters

# VWC class
emm_occ_vwc_group_03 <- emmeans(
  sp_03_mod_occ_md,
  ~ global_vwc_class | foehn,
  type = "response"
)

emm_occ_vwc_group_df_03 <- summary(emm_occ_vwc_group_03) |> as_tibble()

delta_occ_vwc_group_03 <- emm_occ_vwc_group_df_03 |>
  dplyr::group_by(foehn) |>
  dplyr::summarise(
    `dry - intermediate_dry` = prob[global_vwc_class == "dry"] - prob[global_vwc_class == "intermediate_dry"],
    `dry - intermediate_wet` = prob[global_vwc_class == "dry"] - prob[global_vwc_class == "intermediate_wet"],
    `dry - wet` = prob[global_vwc_class == "dry"] - prob[global_vwc_class == "wet"],
    `intermediate_dry - wet` = prob[global_vwc_class == "intermediate_dry"] - prob[global_vwc_class == "wet"],
    `intermediate_dry - intermediate_wet` = prob[global_vwc_class == "intermediate_dry"] - prob[global_vwc_class == "intermediate_wet"],
    `intermediate_wet - wet` = prob[global_vwc_class == "intermediate_wet"] - prob[global_vwc_class == "wet"],
    .groups = "drop"
  )


delta_occ_vwc_group_03

pairs(emm_occ_vwc_group_03)
#size
#month
emm_size_month_group_03 <- emmeans(
  sp_03_mod_size_md,
  ~ month | foehn,
  type = "response"
)

emm_size_month_group_df_03 <- summary(emm_size_month_group_03) |> as_tibble()

delta_size_month_group_03 <- emm_size_month_group_df_03 |>
  dplyr::group_by(foehn) |>
  dplyr::summarise(
    `12 - 1` = emmean[month == "12"] - emmean[month == "1"],
    `12 - 2` = emmean[month == "12"] - emmean[month == "2"],
    `1 - 2`  = emmean[month == "1"]  - emmean[month == "2"],
    .groups = "drop"
  )

delta_size_month_group_03
pairs(emm_size_month_group_03)
#season
emm_size_season_group_03 <- emmeans(
  sp_03_mod_size_md,
  ~ season_factor | foehn,
  type = "response"
)

emm_size_season_group_df_03 <- summary(emm_size_season_group_03) |> as_tibble()

delta_size_season_group_03 <- emm_size_season_group_df_03 |>
  dplyr::select(season_factor, foehn, emmean) |>
  tidyr::pivot_wider(
    names_from = season_factor,
    values_from = emmean
  ) |>
  dplyr::mutate(
    `2016 - 2017` = `2016` - `2017`,
    `2016 - 2019` = `2016` - `2019`,
    `2016 - 2021` = `2016` - `2021`,
    `2016 - 2022` = `2016` - `2022`,
    `2016 - 2023` = `2016` - `2023`,
    `2016 - 2024` = `2016` - `2024`,
    `2016 - 2025` = `2016` - `2025`,
    `2016 - 2026` = `2016` - `2026`,
    
    `2017 - 2019` = `2017` - `2019`,
    `2017 - 2021` = `2017` - `2021`,
    `2017 - 2022` = `2017` - `2022`,
    `2017 - 2023` = `2017` - `2023`,
    `2017 - 2024` = `2017` - `2024`,
    `2017 - 2025` = `2017` - `2025`,
    `2017 - 2026` = `2017` - `2026`,
    
    `2019 - 2021` = `2019` - `2021`,
    `2019 - 2022` = `2019` - `2022`,
    `2019 - 2023` = `2019` - `2023`,
    `2019 - 2024` = `2019` - `2024`,
    `2019 - 2025` = `2019` - `2025`,
    `2019 - 2026` = `2019` - `2026`,
    
    `2021 - 2022` = `2021` - `2022`,
    `2021 - 2023` = `2021` - `2023`,
    `2021 - 2024` = `2021` - `2024`,
    `2021 - 2025` = `2021` - `2025`,
    `2021 - 2026` = `2021` - `2026`,
    
    `2022 - 2023` = `2022` - `2023`,
    `2022 - 2024` = `2022` - `2024`,
    `2022 - 2025` = `2022` - `2025`,
    `2022 - 2026` = `2022` - `2026`,
    
    `2023 - 2024` = `2023` - `2024`,
    `2023 - 2025` = `2023` - `2025`,
    `2023 - 2026` = `2023` - `2026`,
    
    `2024 - 2025` = `2024` - `2025`,
    `2024 - 2026` = `2024` - `2026`,
    
    `2025 - 2026` = `2025` - `2026`
  )
delta_size_season_group_03

# size season TRUE

pair_df_size_true_03 <- summary(pairs(emm_size_season_group_03)) |>
  dplyr::filter(foehn == TRUE) |>
  dplyr::mutate(
    contrast = contrast |>
      gsub("season_factor", "", x = _) |>
      gsub(" / ", " - ", x = _)
  ) |>
  dplyr::left_join(
    delta_size_season_group_03 |>
      dplyr::filter(foehn == TRUE) |>
      tidyr::pivot_longer(
        cols = -foehn,
        names_to = "contrast",
        values_to = "delta"
      ),
    by = "contrast"
  ) |>
  dplyr::mutate(
    different = p.value < 0.05 & abs(delta) >= 0.01
  )


x_size_true_03 <- pair_df_size_true_03$different

names(x_size_true_03) <- pair_df_size_true_03$contrast |>
  gsub(" - ", "-", x = _)

letters_size_true_03 <- multcompView::multcompLetters(x_size_true_03)$Letters

# size season FALSE

pair_df_size_false_03 <- summary(pairs(emm_size_season_group_03)) |>
  dplyr::filter(foehn == FALSE) |>
  dplyr::mutate(
    contrast = contrast |>
      gsub("season_factor", "", x = _) |>
      gsub(" / ", " - ", x = _)
  ) |>
  dplyr::left_join(
    delta_size_season_group_03 |>
      dplyr::filter(foehn == FALSE) |>
      tidyr::pivot_longer(
        cols = -foehn,
        names_to = "contrast",
        values_to = "delta"
      ),
    by = "contrast"
  ) |>
  dplyr::mutate(
    different = p.value < 0.05 & abs(delta) >= 0.01
  )


x_size_false_03 <- pair_df_size_false_03$different

names(x_size_false_03) <- pair_df_size_false_03$contrast |>
  gsub(" - ", "-", x = _)

letters_size_false_03 <- multcompView::multcompLetters(x_size_false_03)$Letters

# VWC class
emm_size_vwc_group_03 <- emmeans(
  sp_03_mod_size_md,
  ~ global_vwc_class | foehn,
  type = "response"
)

emm_size_vwc_group_df_03 <- summary(emm_size_vwc_group_03) |> as_tibble()

delta_size_vwc_group_03 <- emm_size_vwc_group_df_03 |>
  dplyr::group_by(foehn) |>
  dplyr::summarise(
    `dry - intermediate_dry` = emmean[global_vwc_class == "dry"] - emmean[global_vwc_class == "intermediate_dry"],
    `dry - intermediate_wet` = emmean[global_vwc_class == "dry"] - emmean[global_vwc_class == "intermediate_wet"],
    `dry - wet` = emmean[global_vwc_class == "dry"] - emmean[global_vwc_class == "wet"],
    `intermediate_dry - wet` = emmean[global_vwc_class == "intermediate_dry"] - emmean[global_vwc_class == "wet"],
    `intermediate_dry - intermediate_wet` = emmean[global_vwc_class == "intermediate_dry"] - emmean[global_vwc_class == "intermediate_wet"],
    `intermediate_wet - wet` = emmean[global_vwc_class == "intermediate_wet"] - emmean[global_vwc_class == "wet"],
    .groups = "drop"
  )

delta_size_vwc_group_03
pairs(emm_size_vwc_group_03)
############################################################

#interaction plots
#occurence
ylim_occ_03 <- c(
  0,
  max(
    emm_occ_month_group_df_03$asymp.UCL,
    emm_occ_vwc_group_df_03$asymp.UCL,
    emm_occ_season_group_df_03$asymp.UCL,
    na.rm = TRUE
  )
)


emm_month_occ_md_03 <- emmeans(
  sp_03_mod_occ_md,
  ~ foehn | month,
  type = "response"
)


p_occ_month_md_03 <- plot_emm_interaction(
  emm_month_occ_md_03,
  "month",
  "prob",
  "Estimated occurrence probability",
  ylim = ylim_occ_03,
  xlab = "Month"
)
p_occ_month_md_03 


# occ glob_vwc

emm_gvwc_class_occ_md_03 <- emmeans(
  sp_03_mod_occ_md,
  ~ foehn | global_vwc_class,
  type = "response"
)


p_occ_vwc_md_03 <- plot_emm_interaction(
  emm_gvwc_class_occ_md_03,
  "global_vwc_class",
  "prob",
  "Estimated occurrence probability",
  ylim = ylim_occ_03,
  xlab = "VWC class"
)

p_occ_vwc_md_03

# occ season

emm_season_occ_md_03 <- emmeans(
  sp_03_mod_occ_md,
  ~ foehn | season_factor,
  type = "response"
)

p_occ_season_md_03 <- plot_emm_interaction(
  emm_season_occ_md_03,
  "season_factor",
  "prob",
  "Estimated occurrence probability",
  ylim = ylim_occ_03,
  xlab = "Austral summer season"
)

p_occ_season_md_03

#(p_occ_month|p_occ_season_md)

ylim_size_03 <- c(
  min(
    emm_size_month_group_df_03$lower.CL,
    emm_size_vwc_group_df_03$lower.CL,
    emm_size_season_group_df_03$lower.CL
  )-0.01,
  max(
    emm_size_month_group_df_03$upper.CL,
    emm_size_vwc_group_df_03$upper.CL,
    emm_size_season_group_df_03$upper.CL
  )
)


emm_month_size_md_03 <- emmeans(
  sp_03_mod_size_md,
  ~ foehn | month,
  type = "response"
)


p_size_month_md_03 <- plot_emm_interaction(
  emm_month_size_md_03,
  "month",
  "emmean",
  "Estimated mean event size",
  ylim = ylim_size_03,
  xlab = "Month",
  UCL = "upper.CL",
  LCL = "lower.CL"
)
p_size_month_md_03

# size season 
emm_season_size_md_03 <- emmeans(
  sp_03_mod_size_md,
  ~ foehn | season_factor,
  type = "response"
)

p_size_season_md_03 <- plot_emm_interaction(
  emm_season_size_md_03,
  "season_factor",
  "emmean",
  "Estimated mean event size",
  ylim = ylim_size_03,
  xlab = "Austral summer season",
  UCL = "upper.CL",
  LCL = "lower.CL"
)

p_size_season_md_03

emm_gvwc_class_size_md_03 <- emmeans(
  sp_03_mod_size_md,
  ~ foehn | global_vwc_class,
  type = "response"
)

p_size_vwc_md_03 <- plot_emm_interaction(
  emm_gvwc_class_size_md_03,
  "global_vwc_class",
  "emmean",
  "Estimated mean event size",
  ylim = ylim_size_03,
  xlab = "VWC class",
  UCL = "upper.CL",
  LCL = "lower.CL"
)
#(p_size_month_md|p_size_season_md)

############################################################
#foehn_not foehn plots

# letter objects
letters_occ_month_03 <- bind_rows(
  tibble(
    foehn = "FALSE",
    month = c(12,1,2),
    letter = c("a","b","c")
  ),
  tibble(
    foehn = "TRUE",
    month = c(12,1,2),
    letter = c("a","b","c")
  )
)

letters_occ_vwc_03 <- tibble(
  foehn = factor(
    rep(c(FALSE, TRUE), each = 4),
    levels = c(FALSE, TRUE)
  ),
  global_vwc_class = rep(
    c("dry",
      "intermediate_dry",
      "intermediate_wet",
      "wet"),
    2
  ),
  letter = c(
    "a", "b", "b", "b",
    "a", "b", "b", "b"
  )
)

letters_occ_season_03 <- bind_rows(
  tibble(
    foehn = "FALSE",
    season_factor = names(letters_occ_false_03),
    letter = unname(letters_occ_false_03)
  ),
  tibble(
    foehn = "TRUE",
    season_factor = names(letters_occ_true_03),
    letter = unname(letters_occ_true_03)
  )
)


letters_size_season_03 <- bind_rows(
  tibble(
    foehn = "FALSE",
    season_factor = names(letters_size_false_03),
    letter = unname(letters_size_false_03)
  ),
  tibble(
    foehn = "TRUE",
    season_factor = names(letters_size_true_03),
    letter = unname(letters_size_true_03)
  )
)


letters_size_month_03 <- tibble(
  foehn = factor(
    rep(c(FALSE, TRUE), each = 3),
    levels = c(FALSE, TRUE)
  ),
  month = rep(c("12", "1", "2"), 2),
  letter = c(
    "a", "b", "b",
    "a", "b", "ab"
  )
)

letters_size_vwc_03 <- tibble(
  foehn = factor(
    rep(c(FALSE, TRUE), each = 4),
    levels = c(FALSE, TRUE)
  ),
  global_vwc_class = rep(
    c("dry",
      "intermediate_dry",
      "intermediate_wet",
      "wet"),
    2
  ),
  letter = c(
    "a", "a", "a", "a",
    "a", "b", "b", "b"
  )
)


#######################
p_occ_month_letters_03 <- plot_emm_letters(
  emm_occ_month_group_df_03,
  letters_occ_month_03,
  "month",
  "prob",
  "Estimated occurrence probability",
  ylim = ylim_occ_03,
  xlab = "Month"
)

p_occ_season_letters_03 <- plot_emm_letters(
  emm_occ_season_group_df_03,
  letters_occ_season_03,
  "season_factor",
  "prob",
  "Estimated occurrence probability",
  ylim = ylim_occ_03,
  xlab = "Austral summer season"
)

p_occ_vwc_letters_03 <- plot_emm_letters(
  emm_occ_vwc_group_df_03,
  letters_occ_vwc_03,
  "global_vwc_class",
  "prob",
  "Estimated occurrence probability",
  ylim = ylim_occ_03,
  xlab = "VWC class"
)



# size

p_size_month_letters_03 <- plot_emm_letters(
  emm_size_month_group_df_03,
  letters_size_month_03,
  "month",
  "emmean",
  "Estimated mean event size",
  ylim = ylim_size_03,
  xlab = "Month",
  LCL = "lower.CL",
  UCL = "upper.CL"
)

p_size_season_letters_03 <- plot_emm_letters(
  emm_size_season_group_df_03,
  letters_size_season_03,
  "season_factor",
  "emmean",
  "Estimated mean event size",
  ylim = ylim_size_03,
  xlab = "Austral summer season",
  LCL = "lower.CL",
  UCL = "upper.CL"
)


p_size_vwc_letters_03 <- plot_emm_letters(
  emm_size_vwc_group_df_03,
  letters_size_vwc_03,
  "global_vwc_class",
  "emmean",
  "Estimated mean event size",
  ylim = ylim_size_03,
  xlab = "VWC class",
  LCL = "lower.CL",
  UCL = "upper.CL"
)
###############################################################################
#patchwork plots

occ_month_md_clean_03 <- clean_for_grid(
  p_occ_month_md_03,
  keep_legend = TRUE
)

occ_month_letters_clean_03 <- clean_for_grid(
  p_occ_month_letters_03,
  keep_legend = FALSE
)

occ_season_md_clean_03 <- clean_for_grid(
  p_occ_season_md_03,
  keep_legend = FALSE
)

occ_season_letters_clean_03 <- clean_for_grid(
  p_occ_season_letters_03,
  keep_legend = FALSE
)

occ_vwc_md_clean_03 <- clean_for_grid(
  p_occ_vwc_md_03,
  keep_legend = FALSE
)

occ_vwc_letters_clean_03 <- clean_for_grid(
  p_occ_vwc_letters_03,
  keep_legend = FALSE
)

occ_combined_03 <-
  
  (
    occ_month_md_clean_03 |
      occ_month_letters_clean_03
  ) /
  (
    occ_season_md_clean_03 |
      occ_season_letters_clean_03
  ) /
  (
    occ_vwc_md_clean_03 |
      occ_vwc_letters_clean_03
  ) +
  
  plot_layout(
    guides = "collect"
  ) &
  
  theme(
    legend.position = "bottom"
  )
occ_combined_03

library(grid)

occ_final_03 <-
  
  wrap_elements(
    textGrob(
      "Estimated anomaly occurrence probability",
      rot = 90,
      gp = gpar(
        fontsize = 16,
        fontfamily = "Calibri"
      )
    )
  ) +
  occ_combined_03 +
  
  plot_layout(
    widths = c(0.05, 1)
  )

occ_final_03
# occ_comb_no vwc


occ_combined_red_03 <-
  
  (
    occ_month_md_clean_03 |
      occ_month_letters_clean_03
  ) /
  (
    occ_season_md_clean_03 |
      occ_season_letters_clean_03
  )  +
  
  plot_layout(
    guides = "collect"
  ) &
  
  theme(
    legend.position = "bottom"
  )
occ_combined_red_03

occ_final_red_03 <-
  
  wrap_elements(
    textGrob(
      "Estimated anomaly occurrence probability",
      rot = 90,
      gp = gpar(
        fontsize = 16,
        fontfamily = "Calibri"
      )
    )
  ) +
  occ_combined_red_03 +
  
  plot_layout(
    widths = c(0.05, 1)
  )

occ_final_red_03

#size
size_month_md_clean_03 <- clean_for_grid(
  p_size_month_md_03,
  keep_legend = TRUE
)

size_month_letters_clean_03 <- clean_for_grid(
  p_size_month_letters_03,
  keep_legend = FALSE
)

size_season_md_clean_03 <- clean_for_grid(
  p_size_season_md_03,
  keep_legend = FALSE
)

size_season_letters_clean_03 <- clean_for_grid(
  p_size_season_letters_03,
  keep_legend = FALSE
)

size_vwc_md_clean_03 <- clean_for_grid(
  p_size_vwc_md_03,
  keep_legend = FALSE
)

size_vwc_letters_clean_03 <- clean_for_grid(
  p_size_vwc_letters_03,
  keep_legend = FALSE
)


size_combined_03 <-
  
  (
    size_month_md_clean_03 |
      size_month_letters_clean_03
  ) /
  (
    size_season_md_clean_03 |
      size_season_letters_clean_03
  ) /
  (
    size_vwc_md_clean_03 |
      size_vwc_letters_clean_03
  ) +
  
  plot_layout(
    guides = "collect"
  ) &
  
  theme(
    legend.position = "bottom"
  )

size_combined_03

size_final_03 <-
  
  wrap_elements(
    textGrob(
      "Estimated mean anomaly size",
      rot = 90,
      gp = gpar(
        fontsize = 16,
        fontfamily = "Calibri"
      )
    )
  ) +
  size_combined_03 +
  
  plot_layout(
    widths = c(0.05, 1)
  )


size_final_03
# size final reduced
size_combined_red_03 <-
  
  (
    size_month_md_clean_03 |
      size_month_letters_clean_03
  ) /
  (
    size_season_md_clean_03 |
      size_season_letters_clean_03
  ) +
  
  plot_layout(
    guides = "collect"
  ) &
  
  theme(
    legend.position = "bottom"
  )


size_combined_red_03


size_final_red_03 <-
  
  wrap_elements(
    textGrob(
      "Estimated mean anomaly size",
      rot = 90,
      gp = gpar(
        fontsize = 16,
        fontfamily = "Calibri"
      )
    )
  ) +
  size_combined_red_03 +
  
  plot_layout(
    widths = c(0.05, 1)
  )


size_final_red_03

#save
library(Cairo)
#final full plots
CairoPNG(
  filename = "plots/statistic_models_final_03/occurrence_probability_pixel_data_black_only_AOA_final_2_03.png",
  width = 4200,
  height = 5250,
  res = 600
)

print(occ_final_03)

dev.off()

CairoPNG(
  filename = "plots/statistic_models_final_03/event_size_pixel_data_black_only_AOA_final_2_03.png",
  width = 4200,
  height = 5250,
  res = 600
)

print(size_final_03)

dev.off()

#final reduced plots
CairoPNG(
  filename = "plots/statistic_models_final_03/occurrence_probability_pixel_data_black_only_red_AOA_final_2_03.png",
  width = 4200,
  height = 4200,
  res = 600
)

print(occ_final_red_03)

dev.off()


CairoPNG(
  filename = "plots/statistic_models_final_03/event_size_pixel_data_black_only_red_AOA_final_2_03.png",
  width = 4200,
  height = 4200,
  res = 600
)

print(size_final_red_03)

dev.off()

