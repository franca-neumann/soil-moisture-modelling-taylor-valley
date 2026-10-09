#streams significance test. 
#The idea is that the streams should show significantly higher wetness than streamless areas.
#Ofcourse, there are other reasons for wet areas as well, mainly snowpatches. But I assume that, beause in 
#general Taylor Valley has more dry areas than wet areas that we can see a significan ce because all treams should
#be systematically wetter.

#load all streams and combine them wihle keeping the stream name for stream_id random effect
#extract stream VWC predictions while keeping the date/layer information

#extract all not stream conditions while keeping the date/layer information
#combine both

#GLMM distribution; beta familiy, VWC~stream_info + (1|stream_name) + (1|date), data=, familiy = binom familiy

library(terra)
library(sf)
library(purrr)
library(dplyr)
library(tidyr)
library(glmmTMB)
library(car)
library(DHARMa)

pred_paths <- list.files("D:/antarctica/data/processed_data/model_k4_p5/pred_in_AOI_and_AOA/", full.names = T)
names_pred <- list.files("D:/antarctica/data/processed_data/model_k4_p5/pred_in_AOI_and_AOA/", full.names = F)
names <- substr(names_pred, start = 16, stop = 23)
pred_list <- map(pred_paths, rast)
pred_stack <- rast(pred_list)
names(pred_stack) <- names

stream_path <- list.files("C:/Users/mini4/Desktop/Uni/Master/antarctica/data/raw_data/mcmlter-gis-watershed_shapefiles-20231020/mcmlter-gis-watershed-shapefiles/MDV_streams/", pattern = ".shp", full.names =  T)
stream_names <- list.files("C:/Users/mini4/Desktop/Uni/Master/antarctica/data/raw_data/mcmlter-gis-watershed_shapefiles-20231020/mcmlter-gis-watershed-shapefiles/MDV_streams/", pattern = ".shp")
stream_list <- map(stream_path, st_read)
stream_list <- map(stream_list, ~st_transform(.x, crs = crs(pred_stack)))
streams <- bind_rows(stream_list)

streams$stream_names <- stream_names

streams_vect <- vect(streams) 

streams_vwc <- terra::extract(pred_stack, streams_vect, bind=F)

stream_long <- pivot_longer(
  streams_vwc,
  cols = -ID,
  names_to = "date",
  values_to = "VWC"
)



masked_pred_stack <- mask(pred_stack, mask = streams_vect, inverse = T)
not_stream_vwc <- values(masked_pred_stack, na.rm = T)

not_stream_df <- as.data.frame(not_stream_vwc)
not_stream_df$stream_info <- 0  #not_stream_pixel

not_stream_long <- pivot_longer(
  not_stream_df,
  cols = -stream_info,
  names_to = "date",
  values_to = "VWC"
)

not_stream_long$ID <- 0 


set.seed(42)

# subsample stream pixel
stream_sample <- stream_long %>%
  group_by(ID, date) %>%
  sample_n(size = min(1000, n())) %>%
  ungroup()

# n stream pixel per date
n_fluss_per_date <- stream_sample %>%
  group_by(date) %>%
  summarise(n_fluss = n(), .groups = "drop")



compare_data <- bind_rows(stream_long, dry_non_stream) %>%
  mutate(stream_info = ifelse(row_number() <= nrow(stream_long), 1, 0))  # 1 = Stream, 0 = Dry

compare_data$ID <- as.factor(compare_data$ID)      
compare_data$date <- as.factor(compare_data$date)



#############################
#dry_non_stream_sample <- not_stream_long %>%
 # group_by(date) %>%
  #filter(VWC <= quantile(VWC, 0.5, na.rm = TRUE)) %>%  #0.5 quantile as upper bouandary for comparison with dry conditions
  #inner_join(n_fluss_per_date, by = "date") %>%         
  #group_modify(~ {
  #  n_to_sample <- unique(.x$n_fluss)                  
  #  sample_n(.x, size = n_to_sample*1, replace = FALSE) 
 # }) %>%
#  ungroup()

######same for mean instead of median
dry_non_stream_sample <- not_stream_long %>%
  group_by(date) %>%
  filter(VWC <= mean(VWC, na.rm = TRUE)) %>%  #mean VWC as upper boundary
  inner_join(n_fluss_per_date, by = "date") %>%         
  group_modify(~ {
    n_to_sample <- unique(.x$n_fluss)                  
    sample_n(.x, size = n_to_sample*1, replace = FALSE) 
  }) %>%
  ungroup()
######
stream_dry_data <- bind_rows(stream_sample, dry_non_stream_sample) %>%
  mutate(
    stream_info = ifelse(row_number() <= nrow(stream_sample), 1, 0), 
    ID = as.factor(ID),
    date = as.factor(date)
  )

stream_dry_data$stream_info_factor <- as.factor(stream_dry_data$stream_info)


model_stream_dry <- glmmTMB(
  VWC ~ stream_info_factor + (1|date) + (1|ID),
  data = stream_dry_data,
  dispformula = ~ stream_info, #allows for different dispersion per stream_info
  family = beta_family()
) 


car::Anova(model_stream_dry, type = "II")
summary(model_stream_dry)
###check model assumptions
sim_res <- DHARMa::simulateResiduals(
  fittedModel = model_stream_dry)

plot(sim_res)
#testCategorical(sim_res, stream_dry_data$stream_info_factor)
testUniformity(sim_res) #D = 0.05
testDispersion(sim_res) #dispersion = 1.16
#testZeroInflation(sim_res)
#testOutliers(sim_res)
#################

mean(stream_long$VWC, na.rm = T)

for (i in seq_along(stream_names)) {
  m <- mean(stream_long[stream_long$ID == i, ]$VWC, na.rm = TRUE)
  cat(stream_names[[i]], ":", m, "\n")
}

#####
#################
#test if the difference between the mean stream VWC and non- stream VWC is bigger than
#the sensor sensitivity
> # Intercept = No stream
  #from summary (model)
  intercept <- -1.85420
stream_beta <- 0.21427
vwc_non_stream <- exp(intercept) / (1 + exp(intercept))
vwc_stream     <- exp(intercept + stream_beta) / (1 + exp(intercept + stream_beta))
vwc_non_stream  # mean no stream VWC

wc_stream      # mean stream VWC

vwc_stream - vwc_non_stream  # has to be > 0.01 (sensitivity)
