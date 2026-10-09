setwd("D:/antarctica")
setwd("C:/Users/mini4/Desktop/Uni/Master/antarctica")
library(terra)
library(lubridate)
library(sf)
full_stack <- rast("data/processed_data/S1/subset_0_of_S1A_IW_GRDH_1SSH_20200616T124804_20200616T124829_033042_03D3DB_D5E2_Orb_Cal_Stack_TC.tif")
summer_stack <- rast("data/processed_data/S1/summer_ref_stack_Stack_TC.tif")
winter_sub <- full_stack[[-2]]
winter_mean <- mean(winter_sub)
winter_mean

ratio <- summer/winter_mean
#hist(values(ratio), na.rm = T)

ratio_25_01_2022 <- summer_stack$summer_ref_stack_Stack_TC_1 / winter_mean
ratio_db_25_01_2022 <- 10 * log10(ratio_25_01_2022)
snow_patch_25_01 <- ratio_db_25_01_2022 < -2.
plot(snow_patch_25_01)
snow_patch_i_25_01 <- ifel(snow_patch_25_01, 1, 0)

ratio_13_01_2022 <- summer_stack$summer_ref_stack_Stack_TC_2 / winter_mean
ratio_db_13_01_2022 <- 10 * log10(ratio_13_01_2022)
snow_patch_13_01 <- ratio_db_13_01_2022 < -2.0
plot(snow_patch_13_01)
snow_patch_i_13_01 <- ifel(snow_patch_13_01, 1, 0)
snow_reference <- rast("C:/Users/mini4/Desktop/Uni/Master/antarctica/data/processed_data/snow_stack.tif")

r_13_01 <- project(snow_reference[[3]], snow_patch_13_01, method = "near")
plot(r_13_01)

r_25_01 <- project(snow_reference[[7]], snow_patch_i_25_01, method = "near")
plot(r_25_01)


glaciers <- st_read("C:/Users/mini4/Desktop/Uni/Master/antarctica/data/raw_data/Antarctic_GIS/MDV_GLACIERS_STREAMS_POLYGONS/MDV_GLACIERS.shp")
lakes_ponds <- st_read("C:/Users/mini4/Desktop/Uni/Master/antarctica/data/raw_data/GIS layers/MCMLTER_GIS_Export_Layers/USGS_Map_1970/Lakes_and_Poonds_epsg32758.gpkg")

lakes_ponds <- st_transform(lakes_ponds, crs = st_crs(glaciers))


AOI_mask <- st_union(glaciers, lakes_ponds)
#plot(AOI_mask)

AOI_extend <- st_read("C:/Users/mini4/Desktop/Uni/Master/antarctica/data/processed_data/AOI_Taylor.geojson")
AOI_extend <- st_transform(AOI_extend, crs = st_crs(AOI_mask))

r_25_01_m <- mask(r_25_01, mask = AOI_mask, inverse = T)
r_13_01_m <- mask(r_13_01, mask = AOI_mask, inverse = T)



sp_25_01_m <- mask(snow_patch_i_25_01, mask = AOI_mask, inverse = T)
sp_25_01_m_c <- mask(sp_25_01_m, mask = AOI_extend)
sp_13_01_m <- mask(snow_patch_i_13_01, mask = AOI_mask, inverse = T)
sp_13_01_m_c <- mask(sp_13_01_m, mask = AOI_extend)

plot(sp_25_01_m_c)
plot(r_25_01_m_c)

plot(sp_13_01_m_c)
plot(r_13_01_m)

compareGeom(sp_13_01_m_c, r_13_01_m_c)#check needed 
library(MLmetrics)
v_p_13_01 <- values(sp_13_01_m_c)
v_r_13_01 <- values(r_13_01_m)
table(v_p_13_01, v_r_13_01)

v_p_25_01 <- values(sp_25_01_m_c)
v_r_25_01 <- values(r_25_01_m)
table(v_p_25_01, v_r_25_01)

#Foehn
foehn_marte <- read.csv("data/raw_data/data_and_info_Eva_Bendix/foehn_dataset_Marte/fryxell_foehn_results.csv")
head(foehn_marte)
s_w_foehn <- foehn_marte[foehn_marte$foehn_Speirs == 1 & foehn_marte$foehn_Wiesenekker == 1,]
s_w_foehn$date <- date(s_w_foehn$date_time)
unique(s_w_foehn$date)

s_w_foehn_15 <- s_w_foehn[
  s_w_foehn$date >= as.Date("2015-12-18"),]
unique(s_w_foehn_15$date)
sort(unique(s_w_foehn_15$date))

library(dplyr)
result <- s_w_foehn_15 %>%
  group_by(date) %>%
  summarise(
    Speirs_hours = sum(foehn_Speirs) / 4,
    Wiesenekker_hours = sum(foehn_Wiesenekker) / 4
  )
result
subset(result, date == as.Date("2021-01-22"))
