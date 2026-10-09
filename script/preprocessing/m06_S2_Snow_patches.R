
library(terra)
library(diptest)
library(purrr)
library(sf)
library(caret)
library(dplyr)
library(tibble)
############
Taylor_outline <- st_read("C:/Users/mini4/Desktop/Uni/Master/antarctica/data/processed_data/AOI_Taylor.geojson")
Taylor_outline <- st_transform(Taylor_outline, crs(snow_prob_test))
#snow_prob_test <- rast("D:/antarctica/data/raw_data/S2/SP_test/S2B_MSIL2A_20211228T203619_N0500_R085_T57CWP_20221227T125854.SAFE/GRANULE/L2A_T57CWP_A025136_20211228T203622/QI_DATA/MSK_SNWPRB_20m.jp2")
#snow_prob_test <- crop(snow_prob_test, Taylor_outline)
#plot(snow_prob_test)

#cloud_prob_test <- rast("D:/antarctica/data/raw_data/S2/SP_test/S2B_MSIL2A_20211228T203619_N0500_R085_T57CWP_20221227T125854.SAFE/GRANULE/L2A_T57CWP_A025136_20211228T203622/QI_DATA/MSK_CLDPRB_20m.jp2")
#cloud_prob_test <- crop(cloud_prob_test, Taylor_outline)
#plot(cloud_prob_test)

#cls <- rast("D:/antarctica/data/raw_data/S2/SP_test/S2B_MSIL2A_20211228T203619_N0500_R085_T57CWP_20221227T125854.SAFE/GRANULE/L2A_T57CWP_A025136_20211228T203622/IMG_DATA/R20m/T57CWP_20211228T203619_SCL_20m.jp2")
#cls <- crop(cls, Taylor_outline)
#plot(cls)
#cloud_cls <- ifel(cls %in% c(8, 9), cls, NA)
#cloud_cls <- mask(cloud_cls, Taylor_outline)
#plot(cloud_cls)

AOI_sample <- st_read("D:/antarctica/data/processed_data/AOI_BST_validation.gpkg") 
AOI_sample <- st_transform(AOI_sample, crs(blue_20211228))

#blue_20211228 <- rast("D:/antarctica/data/raw_data/S2/SP_test/S2B_MSIL2A_20211228T203619_N0500_R085_T57CWP_20221227T125854.SAFE/GRANULE/L2A_T57CWP_A025136_20211228T203622/IMG_DATA/R10m/T57CWP_20211228T203619_B02_10m.jp2")
#blue_20211228$B03 <- rast("D:/antarctica/data/raw_data/S2/SP_test/S2B_MSIL2A_20211228T203619_N0500_R085_T57CWP_20221227T125854.SAFE/GRANULE/L2A_T57CWP_A025136_20211228T203622/IMG_DATA/R10m/T57CWP_20211228T203619_B03_10m.jp2")
#B11 <- rast("D:/antarctica/data/raw_data/S2/SP_test/S2B_MSIL2A_20211228T203619_N0500_R085_T57CWP_20221227T125854.SAFE/GRANULE/L2A_T57CWP_A025136_20211228T203622/IMG_DATA/R20m/T57CWP_20211228T203619_B11_20m.jp2")
#B11 <- crop(B11, Taylor_outline)

#blue_20211228 <- crop(blue_20211228, Taylor_outline)
#B11 <- resample(B11, blue_20211228)
#blue_20211228$B11 <- B11

#blue_20211228$NDSI <- (blue_20211228$B03 - blue_20211228$B11)/(blue_20211228$B03 + blue_20211228$B11)

#plot(blue_20211228) 
#plot(blue_20211228$NDSI)
#blue_20211228$snow_NDSI <- ifel(blue_20211228$NDSI > 0.39, 1, 0)
#plot(blue_20211228$snow_NDSI)

#msk <- ifel(cloud_prob_test > 50, NA, 1)
#msk <- resample(msk, blue_20211228)
#plot(msk)
#blue_cm <- mask(blue_20211228, msk)
#plot(blue_cm$snow_NDSI)

#AOI
glaciers <- st_read("C:/Users/mini4/Desktop/Uni/Master/antarctica/data/raw_data/Antarctic_GIS/MDV_GLACIERS_STREAMS_POLYGONS/MDV_GLACIERS.shp")
lakes_ponds <- st_read("C:/Users/mini4/Desktop/Uni/Master/antarctica/data/raw_data/GIS layers/MCMLTER_GIS_Export_Layers/USGS_Map_1970/Lakes_and_Poonds_epsg32758.gpkg")
glaciers <- st_transform(glaciers, crs = st_crs(lakes_ponds))
AOI_mask <- st_union(glaciers, lakes_ponds)
plot(AOI_mask)
AOI_mask <- st_transform(AOI_mask, crs = crs(blue_20211228))

##########################

#Blue Snow Threshold (BST)

blue_snow_threshold <- function(raster, blue, alpha = 0.05) {
  
  values_blue <- values(raster[[blue]])
  values_blue <- values_blue[!is.na(values_blue)]
  
  mean_blue <- mean(values_blue)
  
  # high snow image? mean_blue > 0.7
  if (mean_blue > 0.7) {
    t <- 0.7
    message("mean_blue > 0.7; Threshold t = 0.7")
    
  } else {
    # Dip-test bimodality
    dip_result <- dip.test(values_blue)
    
    if (dip_result$p.value < alpha) {
      # Bimodal; density and smoothing
      dens <- density(values_blue, bw = "nrd0") 
      y <- dens$y
      x <- dens$x
      local_min <- which(diff(sign(diff(y))) == 2) + 1
      t_candidates <- x[local_min]
      
      # first min > mean_blue
      t_candidates <- t_candidates[t_candidates > mean_blue]
      
      if (length(t_candidates) > 0) {
        t <- t_candidates[1]
        message("Bimodal distribution; t based on first min of smoothed density histogram")
      } else {
        # fallback: mean
        t <- mean_blue
        message("Bimodal but no local min found; t based on mean")
      }
      
    } else {
      # unimodal -> mean
      t <- mean_blue
      message("Unimodal distribution; t based on mean")
    }
  }
  
  # snow layer
  snow_layer <- ifel(raster[[blue]] > t, 1, 0)
  raster[["snow"]] <- snow_layer
  
  cat("Final Threshold t =", t, "\n")
  return(raster)
}


  

#vt1 <- rast("data/processed_data/Planet_scope_merged/merged_20211226.tif")
#t1$blue_norm <- t1$blue/10000
#s1 <- blue_snow_threshold(t1, blue = "blue_norm")
#plot(s1$snow)
#blue_cm$blue_norm <- blue_cm$T57CWP_20211228T203619_B02_10m/10000
#s1 <- blue_snow_threshold(blue_cm, blue = "blue_norm")
#plot(s1$snow)
#s2 <- blue_snow_threshold(blue_20211228, blue = "NDSI")
# plot(s2$snow)

###
#bm <- mask(blue_20211228, Taylor_outline)
#bm_cm <- mask(bm, msk)
#s3 <- blue_snow_threshold(bm_cm, blue = "NDSI")
#plot(s3$snow)

#bm$blue_norm <- bm$T57CWP_20211228T203619_B02_10m/10000
#s4 <- blue_snow_threshold(bm_cm, blue = "blue_norm")
#plot(s4$snow)


###################################################################
#load all Sentine-2 scenes from a base direction 
base_dir <- "D:/antarctica/data/raw_data/S2/SP_test"

safe_dirs <- list.dirs(
  base_dir,
  recursive = FALSE,
  full.names = TRUE
)

load_s2_scene <- function(safe_path) {
  
  b02 <- rast(list.files(safe_path, "B02_10m.jp2$", recursive=TRUE, full.names=TRUE))
  b03 <- rast(list.files(safe_path, "B03_10m.jp2$", recursive=TRUE, full.names=TRUE))
  b11 <- rast(list.files(safe_path, "B11_20m.jp2$", recursive=TRUE, full.names=TRUE))
  scl <- rast(list.files(safe_path, "SCL_20m.jp2$", recursive=TRUE, full.names=TRUE))
  scl <- resample(scl, b02)
  
  cls <- ifel(scl %in% c(8, 9), scl, NA) #medium and high probability for clouds
  csh <- ifel(scl == 3, scl, NA ) #cloud shadow
  si <- ifel(scl == 11, scl, NA) #snow or ice
  
  b11 <- resample(b11, b02)
  
  s2 <- c(b02, b03, b11, cls, csh, si, scl)
  names(s2) <- c("blue", "green", "swir", "cls", "csh", "si", "scl")
  return(s2)
}


##
#writing directory for writing the stratified point samples   
setwd("D:/antarctica")
out_dir <- file.path(
  getwd(),
  "data",
  "processed_data",
  "BST_validation",
  "Individ_pointsample"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


#function to generate stratified random samples; 3 classes both snow/no snow and 
#NDSI and NDSI_bst indicate different categories
one_step_indiv_strat_rs <- function(safe_path, AOI, outline){
  
  #load S2
  s2 <- load_s2_scene(safe_path)
 
   outline <- st_transform(outline, crs(s2))
  AOI <- st_transform(AOI, crs(s2))
  
  s2 <- crop(s2, outline)
  cls <- s2$cls
  s2 <- mask(s2, cls, inverse = T)
  s2 <- mask(s2, outline)
  #NDSI
  s2$NDSI <- (s2$green - s2$swir) / (s2$green + s2$swir)
  
  #snow classification based on common threshold
  snow_ndsi <- ifel(s2$NDSI > 0.4, 1, 0) #from Sentinel documentation
  names(snow_ndsi) <- "snow"
  
  #BST with NDSI as "blue" band (NDSI shows bimodal distribution as well, 
  #0.4 is known to work well for global applications but not neccesary at the local scale)
  #as for all timesteps a bimodal doistribution was found, we can use the blue_snow_threshold function 
  #for this try out and if it is successful, an own NDSI_snow_threshold function 
  #will be implemented afterwards, that then will be capable of dealing with non bimodal 
  #scearios for the NDSI as well.
  ndsi_r <- s2[["NDSI"]]
  names(ndsi_r) <- "NDSI"
  
  bst_ndsi <- blue_snow_threshold(ndsi_r, blue = "NDSI")
  
  #BST normal
  bst_r <- s2[["blue"]]
  names(bst_r) <- "blue"
  bst_r$blue_norm <- bst_r$blue/10000
  
  bst_blue <- blue_snow_threshold(bst_r, blue = "blue_norm")
  
  #4 class raster (TP/FP/FN/TN) there will be only three classes, depending on bst_ndsi >/< 0.4
  class_raster <- snow_ndsi * 2 + bst_ndsi[["snow"]]
  names(class_raster) <- "class"
  
  #AOI mask to retrieve the sample only from the middle of the valley (excluding 
  #the glaciers and lakes, aim is to find the snow patches for the validation approach 
  #as this is what we are mainly interested in)
  class_raster <- mask(class_raster, AOI)
  
  #stratified random sampling
  set.seed(123)
  
  samples <- spatSample(
    class_raster,
    size = 20,
    method = "stratified",
    na.rm = TRUE,
    as.points = TRUE
  )
  
  samples <- st_as_sf(samples)
  
  #extract values
  ref  <- terra::extract(snow_ndsi, samples)
  pred <- terra::extract(bst_ndsi[["snow"]], samples)
  pred2 <- terra::extract(bst_blue[["snow"]], samples)
  
  samples$snow_ndsi <- ref$snow
  samples$snow_bst_ndsi  <- pred$snow
  samples$snow_bst_blue  <- pred2$snow
  
  #datetime
  name <- stringr::str_extract(safe_path, "\\d{8}")
  samples$datetime <- lubridate::ymd(name)
  
  #output
 
  out_file <- file.path(
  out_dir,
  paste0("complete_strat_rd_", name, ".gpkg")
)

write_sf(samples, out_file, append = FALSE)

  message("wrote ", name)
  
  gc()
}
######


 #batch
safe_dirs <- list.dirs(
  "D:/antarctica/data/raw_data/S2/SP_test",
  recursive = FALSE,
  full.names = TRUE
)

safe_dirs <- safe_dirs[grepl("\\.SAFE$", safe_dirs)]

one_step_indiv_strat_rs(safe_dirs[1], AOI_sample, Taylor_outline)
purrr::map(safe_dirs, ~ one_step_indiv_strat_rs(.x, AOI_sample, Taylor_outline))

##################################################
#generate merged and cropped to taylor S2 scenes
out_dir2 <- file.path(
  getwd(),
  "data",
  "processed_data",
  "BST_validation",
  "S2_merged"
)

dir.create(out_dir2, recursive = TRUE, showWarnings = FALSE)

#write merged and cropped scenes for generating ground truth data
load_write_s2_scene <- function(safe_path, outline) {
  
  b02 <- rast(list.files(safe_path, "B02_10m.jp2$", recursive=TRUE, full.names=TRUE))
  b03 <- rast(list.files(safe_path, "B03_10m.jp2$", recursive=TRUE, full.names=TRUE))
  b11 <- rast(list.files(safe_path, "B11_20m.jp2$", recursive=TRUE, full.names=TRUE))
  b04 <- rast(list.files(safe_path, "B04_10m.jp2$", recursive=TRUE, full.names=TRUE))
  
  b11 <- resample(b11, b02)
 
  s2 <- c(b02, b03, b11, b04)
  names(s2) <- c("blue", "green", "swir", "red")
  
  outline <- st_transform(outline, crs(s2))
  
  s2 <- crop(s2, outline)
  
  name <- stringr::str_extract(safe_path, "\\d{8}")
  
  out_file <- file.path(
    out_dir2,
    paste0("sen2", name, ".tif")
  )
  
  writeRaster(s2, out_file, overwrite = T)
  
  message("wrote ", name)
  
  gc()
}

purrr::map(safe_dirs, ~ load_write_s2_scene(.x, Taylor_outline))
##################################################
##################################################################

#after manually assigning the ground truth column, 
# calculate the F score now

Bst_path <- list.files("data/processed_data/BST_validation/Individ_pointsample/", pattern = ".gpkg$", full.names = T)
Bst_names <- list.files("data/processed_data/BST_validation/Individ_pointsample/", pattern = ".gpkg$",full.names = F)
bst_all <- map(Bst_path, st_read)
target_crs <- st_crs(bst_all[[1]])
bst_all <- lapply(bst_all, st_transform, crs = target_crs)
bst_list <- bst_all
Bst_nsdi_b4 <- bst_all[[7]]
bst_all <- bst_all[-7]

##################
#function to calculate the F1 scores of the different snow identification approaches
indiv_f_scores <- function(test_data){
  
  dt <- unique(test_data$datetime)[1]
  
  out <- data.frame(
    datetime = dt,
    
    F1_NDSI = F_meas(
      data = factor(test_data$snow_ndsi),
      reference = factor(test_data$snow_gt)
    ),
    
    F1_NDSI_BST = F_meas(
      data = factor(test_data$snow_bst_ndsi),
      reference = factor(test_data$snow_gt)
    ),
    
    F1_BLUE_BST = F_meas(
      data = factor(test_data$snow_bst_blue),
      reference = factor(test_data$snow_gt)
    )
  )
  
  return(out)
}


fscore_results <- purrr::map_dfr(bst_list, indiv_f_scores)
fscore_results #F1 scores for individual dates

bst_all <- do.call(rbind, bst_all)


F1_NDSI <- F_meas(data = factor(bst_all$snow_ndsi), 
                    reference = factor(bst_all$snow_gt))
F1_NDSI #0.74

F1_NDSI_bst <- F_meas(data = factor(bst_all$snow_bst_ndsi), 
                  reference = factor(bst_all$snow_gt))
F1_NDSI_bst #0.9

F1_orig_BST <- F_meas(data = factor(bst_all$snow_bst_blue), 
                      reference = factor(bst_all$snow_gt))
F1_orig_BST #0.86 

F1_outl_ndsi <- F_meas(data = factor(Bst_nsdi_b4$snow_ndsi), 
                  reference = factor(Bst_nsdi_b4$snow_gt))
F1_outl_ndsi #0.98

F1_outl_bst_ndsi <- F_meas(data = factor(Bst_nsdi_b4$snow_bst_ndsi), 
                  reference = factor(Bst_nsdi_b4$snow_gt))
F1_outl_bst_ndsi #0.69


# the F1 score for the blue snow threshold for PlanetScope was approx. 0.9
# the individual dates' F1 score do have a broader range than this

####
#NDSI snow Threshold; if no bimodal distribution was identified (with effect size D at least 0.001)
#the common TH value of 0.4 will be used as the fall back. 
NDSI_snow_threshold <- function(raster, blue, alpha = 0.05, D = 0.001, fall_back = 0.4) {
  
  values_blue <- values(raster[[blue]])
  values_blue <- values_blue[!is.na(values_blue)]
  
  mean_blue <- mean(values_blue)
  
  {
    dip_result <- dip.test(values_blue)
    
    if (dip_result$p.value < alpha && dip_result$statistic > D) {
      
      dens <- density(values_blue, bw = "nrd0")
      y <- dens$y
      x <- dens$x
      
      local_min <- which(diff(sign(diff(y))) == 2) + 1
      t_candidates <- x[local_min]
      t_candidates <- t_candidates[t_candidates > mean_blue]
      
      if (length(t_candidates) > 0) {
        t <- t_candidates[1]
        message("Bimodal distribution; t based on first min of smoothed density histogram, t = ",t)
      } else {
        t <- fall_back
        message("Bimodal but no local min found; t set to fall back value, t = ", t)
      }
      
    } else {
      
      t <- fall_back
      message("Unimodal distribution; t set to fall back value, t = ",t)
    }
  }
  
  snow_layer <- ifel(raster[[blue]] > t, 1, 0)
  raster[["snow"]] <- snow_layer
  
  cat("Final Threshold t =", t, "\n")
  
  return(raster)
}

##
#generate additional test data for the S2 Scene classification layer and a refined BST NDSI approach.
add_test_data <- function(safe_path, samples, AOI, outline, out_dir){
  
  # load S2
  s2 <- load_s2_scene(safe_path)
  
  outline <- st_transform(outline, crs(s2))
  AOI <- st_transform(AOI, crs(s2))
  
  s2 <- crop(s2, outline)
  cls <- s2$cls
  csh <- ifel(s2$scl == 3, 1, NA)
  s2 <- mask(s2, cls, inverse = T)
  s2 <- mask(s2, outline)
  
 
  
  # NDSI
  s2$NDSI <- (s2$green - s2$swir) / (s2$green + s2$swir)
  th_data <- mask(s2, csh, inverse = T) #exclude cloud shadow influenced data from the data used to find the threshlod
  # NDSI common Th
  snow_ndsi <- ifel(s2$NDSI > 0.4, 1, 0) #from Sentinel documentation
  names(snow_ndsi) <- "snow"
  
  # NDSI snow threshold
  ndsi_r <- th_data[["NDSI"]]
  names(ndsi_r) <- "NDSI"
  
  ndsi_full <- s2[["NDSI"]]
  names(ndsi_full) <- "NDSI"
  
  bst_ndsi <- NDSI_snow_threshold(ndsi_r, ndsi_full, blue = "NDSI")
  
  # BST normal
  bst_r <- th_data[["blue"]]
  names(bst_r) <- "blue"
  bst_r$blue_norm <- bst_r$blue/10000
 
  bst_full <- s2[["blue"]]
  names(bst_full) <- "blue"
  bst_full$blue_norm <- bst_full$blue/10000
  
  bst_blue <- blue_snow_threshold(bst_r, bst_full, blue = "blue_norm")
  
  # class raster (TP/FP/FN/TN)
  class_raster <- snow_ndsi[["snow"]] * 2 + bst_ndsi[["snow"]]
  names(class_raster) <- "class"
  
  # AOI mask
  class_raster <- mask(class_raster, AOI)
  
  # S2 snow_ice
  s2_snow_ice <- ifel(s2$scl == 11, 1, 0)
  names(s2_snow_ice) <- "snow"
  
  # Extract values
  pred <- terra::extract(bst_ndsi[["snow"]], samples)
  pred2 <- terra::extract(bst_blue[["snow"]], samples)
  pred3 <- terra::extract(s2_snow_ice, samples)
  pred4 <- terra::extract(class_raster, samples)
  
  samples$snow_bst_ndsi_ref  <- pred$snow
  samples$snow_bst_blue_ref  <- pred2$snow
  samples$snow_s2 <- pred3$snow
  samples$new_classes <- pred4$class
  # datetime
  name <- stringr::str_extract(safe_path, "\\d{8}")
 
  
  # output
  
  out_file <- file.path(
    out_dir,
    paste0("complete_strat_rd_ref", name, ".gpkg")
  )
  
  write_sf(samples, out_file, append = FALSE)
  
  message("wrote ", name)
  
  gc()
}
##############################################################

#NDSI snow Threshold; one raster only 
NDSI_snow_threshold_1r <- function(raster, blue, alpha = 0.05, D = 0.001, fall_back = 0.4) {
  
  values_blue <- values(raster[[blue]])
  values_blue <- values_blue[!is.na(values_blue)]
  
  mean_blue <- mean(values_blue)
  
  {
    dip_result <- dip.test(values_blue)
    
    if (dip_result$p.value < alpha && dip_result$statistic > D) {
      
      dens <- density(values_blue, bw = "nrd0")
      y <- dens$y
      x <- dens$x
      
      local_min <- which(diff(sign(diff(y))) == 2) + 1
      t_candidates <- x[local_min]
      t_candidates <- t_candidates[t_candidates > mean_blue]
      
      if (length(t_candidates) > 0) {
        t <- t_candidates[1]
        message("Bimodal distribution; t based on first min of smoothed density histogram, t = ",t)
      } else {
        t <- fall_back
        message("Bimodal but no local min found; t set to fall back value, t = ", t)
      }
      
    } else {
      
      t <- fall_back
      message("Unimodal distribution; t set to fall back value, t = ",t)
    }
  }
  
  snow_layer <- ifel(raster[[blue]] > t, 1, 0)
  raster[["snow"]] <- snow_layer
  
  cat("Final Threshold t =", t, "\n")
  
  return(raster)
}

##
#final function for additional test data
add_test_data_1r <- function(safe_path, samples, AOI, outline, out_dir){
  
  # S2
  s2 <- load_s2_scene(safe_path)
  
  outline <- st_transform(outline, crs(s2))
  AOI <- st_transform(AOI, crs(s2))
  
  s2 <- crop(s2, outline)
  cls <- s2$cls
  s2 <- mask(s2, cls, inverse = T)
  s2 <- mask(s2, outline)
  
  
  # NDSI
  s2$NDSI <- (s2$green - s2$swir) / (s2$green + s2$swir)
 
  snow_ndsi <- ifel(s2$NDSI > 0.4, 1, 0) #from Sentinel documentation
  names(snow_ndsi) <- "snow"
  
  # NDSI snow threshold 
  ndsi_r <- s2[["NDSI"]]
  names(ndsi_r) <- "NDSI"
  
  bst_ndsi <- NDSI_snow_threshold_1r(ndsi_r, blue = "NDSI")
  
  # BST 
  bst_r <- s2[["blue"]]
  names(bst_r) <- "blue"
  bst_r$blue_norm <- bst_r$blue/10000
  
  bst_full <- s2[["blue"]]
  names(bst_full) <- "blue"
  bst_full$blue_norm <- bst_full$blue/10000
  
  bst_blue <- blue_snow_threshold(bst_r, bst_full, blue = "blue_norm")
  
  # class raster (TP/FP/FN/TN)
  class_raster <- snow_ndsi[["snow"]] * 2 + bst_ndsi[["snow"]]
  names(class_raster) <- "class"
  
  # AOI mask
  class_raster <- mask(class_raster, AOI)
  
  # S2 snow_ice
  s2_snow_ice <- ifel(s2$scl == 11, 1, 0)
  names(s2_snow_ice) <- "snow"
 
  # extract values
  pred <- terra::extract(bst_ndsi[["snow"]], samples)
  pred2 <- terra::extract(bst_blue[["snow"]], samples)
  pred3 <- terra::extract(s2_snow_ice, samples)
  pred4 <- terra::extract(class_raster, samples)
  
  samples$snow_bst_ndsi_ref  <- pred$snow
  samples$snow_bst_blue_ref  <- pred2$snow
  samples$snow_s2 <- pred3$snow
  samples$new_classes <- pred4$class
  
  # datetime
  name <- stringr::str_extract(safe_path, "\\d{8}")
  
  # output
  out_file <- file.path(
    out_dir,
    paste0("complete_strat_rd_ref_1r", name, ".gpkg")
  )
  
  write_sf(samples, out_file, append = FALSE)
  
  message("wrote ", name)
  
  gc()
}

#
out_dir5 <- file.path(
  getwd(),
  "data",
  "processed_data",
  "BST_validation",
  "add_samples_excls_fallback_value_1r"
)

dir.create(out_dir5, recursive = TRUE, showWarnings = FALSE)

purrr::map2(
  safe_match,
  bst_list,
  ~ add_test_data_1r(.x, .y, AOI_sample, Taylor_outline, out_dir5)
)

##############################################################

safe_dirs2 <- safe_dirs[-7]
out_dir3 <- file.path(
  getwd(),
  "data",
  "processed_data",
  "BST_validation",
  "add_samples_t2"
)

dir.create(out_dir3, recursive = TRUE, showWarnings = FALSE)

safe_dates <- stringr::str_extract(safe_dirs, "\\d{8}")
sample_dates <- stringr::str_extract(Bst_path, "\\d{8}")
idx <- match(sample_dates, safe_dates)
safe_match <- safe_dirs[idx]

purrr::map2(
  safe_match,
  bst_list,
  ~ add_test_data(.x, .y, AOI_sample, Taylor_outline, out_dir3)
)

##
out_dir4 <- file.path(
  getwd(),
  "data",
  "processed_data",
  "BST_validation",
  "add_samples_excls_fallback_value"
)

dir.create(out_dir4, recursive = TRUE, showWarnings = FALSE)

purrr::map2(
  safe_match,
  bst_list,
  ~ add_test_data(.x, .y, AOI_sample, Taylor_outline, out_dir4)
)

#####
Bst_path2 <- list.files("data/processed_data/BST_validation/add_samples/", pattern = ".gpkg$", full.names = T)
Bst_names2 <- list.files("data/processed_data/BST_validation/add_samples/", pattern = ".gpkg$",full.names = F)
bst_all2 <- map(Bst_path2, st_read)
target_crs <- st_crs(bst_all2[[1]])
bst_all2 <- lapply(bst_all, st_transform, crs = target_crs)
bst_list <- bst_all
Bst_nsdi_b4 <- bst_all[[7]]
bst_all <- bst_all[-7]

###
Bst_path3 <- list.files("data/processed_data/BST_validation/add_samples_t2/", pattern = ".gpkg$", full.names = T)
Bst_names3 <- list.files("data/processed_data/BST_validation/add_samples_t2/", pattern = ".gpkg$",full.names = F)
bst_all3 <- map(Bst_path3, st_read)
target_crs <- st_crs(bst_all3[[1]])
bst_all3 <- lapply(bst_all3, st_transform, crs = target_crs)
bst_list3 <- bst_all3
bst_all3 <- do.call(rbind, bst_all3)

#
Bst_path4 <- list.files("data/processed_data/BST_validation/add_samples_excls_fallback_value/", pattern = ".gpkg$", full.names = T)
Bst_names4 <- list.files("data/processed_data/BST_validation/add_samples_excls_fallback_value/", pattern = ".gpkg$",full.names = F)
bst_all4 <- map(Bst_path4, st_read)
target_crs <- st_crs(bst_all4[[1]])
bst_all4 <- lapply(bst_all4, st_transform, crs = target_crs)
bst_list4 <- bst_all4
bst_all4 <- do.call(rbind, bst_all4)

#

Bst_path5 <- list.files("data/processed_data/BST_validation/add_samples_excls_fallback_value_1r/", pattern = ".gpkg$", full.names = T)
Bst_names5 <- list.files("data/processed_data/BST_validation/add_samples_excls_fallback_value_1r/", pattern = ".gpkg$",full.names = F)
bst_all5 <- map(Bst_path5, st_read)
target_crs <- st_crs(bst_all5[[1]])
bst_all5 <- lapply(bst_all5, st_transform, crs = target_crs)
bst_list5 <- bst_all5
bst_all5 <- do.call(rbind, bst_all5)
##
#function for datewise F scores, including the adjusted BST and the S2 snow classification from SCL layer
indiv_f_scores2 <- function(test_data){
  
  dt <- unique(test_data$datetime)[1]
  
  out <- data.frame(
    datetime = dt,
    
    F1_NDSI = F_meas(
      data = factor(test_data$snow_ndsi),
      reference = factor(test_data$snow_gt)
    ),
    
    F1_NDSI_BST = F_meas(
      data = factor(test_data$snow_bst_ndsi),
      reference = factor(test_data$snow_gt)
    ),
    
    F1_NDSI_BST_ref = F_meas(
      data = factor(test_data$snow_bst_ndsi_ref, levels = c(0,1)),
      reference = factor(test_data$snow_gt)
    ),
    
    F1_BLUE_BST = F_meas(
      data = factor(test_data$snow_bst_blue),
      reference = factor(test_data$snow_gt)
    ),
    
    F1_BLUE_BST_ref = F_meas(
      data = factor(test_data$snow_bst_blue_ref, levels = c(0,1)),
      reference = factor(test_data$snow_gt)
    ),
    F1_snow = F_meas(
      data = factor(test_data$snow_s2, levels = c(0,1)),
      reference = factor(test_data$snow_gt)
    )
    
  )
  
  return(out)
}

bst_all2_4 <- bst_all2[4]
bst_all2 <- bst_all2[-4]
fscore_results_2 <- purrr::map_dfr(bst_all2, indiv_f_scores2)
fscore_results_3 <- purrr::map_dfr(bst_list3, indiv_f_scores2)
fscore_results_4 <- purrr::map_dfr(bst_list4, indiv_f_scores2)
fscore_results_5 <- purrr::map_dfr(bst_list5, indiv_f_scores2)

fscore_results_4
fscore_results_5

#overall F scores
indiv_f_scores2(bst_all4)
indiv_f_scores2(bst_all5)

###

# identification of the probability of an F score to be the best for an individual date
methods <- c(
  "F1_NDSI",
  "F1_NDSI_BST_ref",
  "F1_BLUE_BST_ref",
  "F1_snow"
)

best_results <- map_dfr(1:nrow(fscore_results_5), function(i) {
  
  vals <- unlist(fscore_results_5[i, methods])
  best <- names(vals)[vals == max(vals)]
  
  tibble(
    datetime = fscore_results_5$datetime[i],
    best_method = best
  )
})

best_summary <- tibble(best_method = methods) %>%
  left_join(
    best_results %>% count(best_method, name = "wins"),
    by = "best_method"
  ) %>%
  mutate(
    wins = ifelse(is.na(wins), 0L, wins),
    percent = round(100 * wins / nrow(fscore_results_5), 1)
  ) %>%
  arrange(desc(wins))

best_summary
##


#load S2 scene for AOI
load_s2_taylor <- function(safe_path, outline, AOI){
  
s2 <- load_s2_scene(safe_path)
outline <- st_transform(outline, crs(s2))
AOI <- st_transform(AOI, crs(s2))
s2 <- crop(s2, outline)
cls <- s2$cls
s2 <- mask(s2, cls, inverse = T)
s2 <- mask(s2, outline)
s2 <- mask(s2, AOI, inverse = T)
s2$NDSI <- (s2$green - s2$swir) / (s2$green + s2$swir)
csh <- s2$csh
s2$NDSI_no_shadow <- mask(s2$NDSI, csh, inverse = T)
return(s2)
}

#s2_scenes_talyor <- purrr::map(safe_dirs, ~ load_s2_taylor(.x, Taylor_outline, AOI_sample))
#s2_220128 <- load_s2_taylor(safe_dirs[6],Taylor_outline, AOI_mask)
#hist(values(s2_220128$NDSI_no_shadow))
#hist(values(s2_220128$NDSI))
#scl_da <- ifel(s2_220128$scl == 2, 1, NA)
#s2_csh_da <- mask(s2_220128$NDSI_no_shadow, scl_da, inverse = T)
#hist(values(s2_csh_da))
#dip_result <- dip.test(values(s2_csh_da))

##
#s2_221122 <- load_s2_taylor(safe_dirs[8],Taylor_outline, AOI_mask)
#hist(values(s2_221122$NDSI_no_shadow))
#hist(values(s2_221122$NDSI))
#scl_da_11 <- ifel(s2_221122$scl == 2, 1, NA)
#s2_11_csh_da <- mask(s2_221122$NDSI_no_shadow, scl_da_11, inverse = T)
#hist(values(s2_11_csh_da))
#dip_result <- dip.test(values(s2_11_csh_da))

#s2_211228 <- load_s2_taylor(safe_dirs[1],Taylor_outline, AOI_mask)
#hist(values(s2_211228$NDSI_no_shadow))
#hist(values(s2_211228$NDSI))
#scl_da_28 <- ifel(s2_211228$scl == 2, 1, NA)
#s2_28_csh_da <- mask(s2_211228$NDSI_no_shadow, scl_da_28, inverse = T)
#hist(values(s2_28_csh_da))
#dip_result <- dip.test(values(s2_28_csh_da))
####


#so; calculate snowpatch layer for all scenes, but together in one stack, add time
#argument and calculate distance to closest snowpatch
load_s2_for_NDSI_snow_th <- function(safe_path, outline){
  
    b03 <- rast(list.files(safe_path, "B03_10m.jp2$", recursive=TRUE, full.names=TRUE))
    b11 <- rast(list.files(safe_path, "B11_20m.jp2$", recursive=TRUE, full.names=TRUE))
    b11 <- resample(b11, b03)
    scl <- rast(list.files(safe_path, "SCL_20m.jp2$", recursive=TRUE, full.names=TRUE))
    scl <- resample(scl, b03)
    
    cls <- ifel(scl %in% c(8, 9), scl, NA) #medium and high probability for clouds
    
    s2 <- c(b03, b11, cls)
    names(s2) <- c( "green", "swir", "cls")
  
    outline <- st_transform(outline, crs(s2))
    s2 <- crop(s2, outline)
    cls <- s2$cls
    s2 <- mask(s2, cls, inverse = T)
    s2 <- mask(s2, outline)
    s2$NDSI <- (s2$green - s2$swir) / (s2$green + s2$swir)
    return(s2)
}


get_snow_layer_NDSI_ST <- function(safe_path, outline){
  s2 <- load_s2_for_NDSI_snow_th(safe_path, outline)
  
  bst_ndsi <- NDSI_snow_threshold_1r(ndsi_r, blue = "NDSI")
  snow <- bst_ndsi[["snow"]]
  name <- stringr::str_extract(safe_path, "\\d{8}")
  time(snow) <- ymd(name)
  message(name, " done")
  rm(s2)
  return(snow)
}
##########################################################################
NDSI_snow_threshold_t <- function(raster, blue,
                                  alpha = 0.05,
                                  D = 0.001,
                                  fall_back = 0.4) {
  
  values_blue <- values(raster[[blue]])
  values_blue <- values_blue[!is.na(values_blue)]
  
  mean_blue <- mean(values_blue)
  
  dip_result <- dip.test(values_blue)
  
  if (dip_result$p.value < alpha &&
      dip_result$statistic > D) {
    
    dens <- density(values_blue, bw = "nrd0")
    y <- dens$y
    x <- dens$x
    
    local_min <- which(diff(sign(diff(y))) == 2) + 1
    t_candidates <- x[local_min]
    t_candidates <- t_candidates[t_candidates > mean_blue]
    
    if (length(t_candidates) > 0) {
      
      t <- t_candidates[1]
      
      message(
        "Bimodal distribution; first local minimum used. t = ",
        round(t, 3)
      )
      
    } else {
      
      t <- fall_back
      
      message(
        "Bimodal distribution, but no local minimum found; fallback used. t = ", t)
    }
    
  } else {
    
    t <- fall_back
    
    message(
      "Unimodal distribution; fallback used. t = ", t)
  }
  
  raster[["snow"]] <- ifel(raster[[blue]] > t, 1, 0)
  
  message("Final threshold = ", round(t, 3))
  
  list(
    raster = raster,
    threshold = t
  )
}


get_snow_layer_NDSI_ST <- function(safe_path, outline){
  
  s2 <- load_s2_for_NDSI_snow_th(safe_path, outline)
  
  res <- NDSI_snow_threshold_t(s2, blue = "NDSI")
  
  snow <- res$raster[["snow"]]
  th   <- res$threshold
  
  name <- stringr::str_extract(safe_path, "\\d{8}")
  date <- lubridate::ymd(name)
  
  time(snow) <- date
  message("calculated ", date)
  
  tibble::tibble(
    date = date,
    threshold = th,
    snow = list(snow)
  )
}



base_dir <- "D:/antarctica/data/raw_data/S2/16_17/"

zip_files <- list.files(base_dir, pattern = "\\.zip$", full.names = TRUE)
purrr::walk(zip_files, unzip, exdir = base_dir)

#S2 season 15_16
safe_dirs_15_16 <- list.files(
  base_dir,
  pattern = "\\.SAFE$",
  full.names = TRUE
)

res_15_16 <- map_dfr(safe_dirs_15_16, get_snow_layer_NDSI_ST,
        outline = Taylor_outline)
#write Th table
df_out <- res_15_16 %>%
  dplyr::select(date, threshold)

write.csv(
  df_out,
  file = file.path("D:/antarctica/data/processed_data/snow_stacks/Thresholds_NDSI/", "thresholds_2015_2016.csv"),
  row.names = FALSE
)

snow_list_15_16 <- res_15_16$snow
#save snow layers
out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/S2_15_16/snow_stack/"

purrr::walk2(
  res_15_16$date,
  res_15_16$snow,
  function(d, r) {
    
    name <- format(d, "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)

snow_dist_list_15_16 <- map(snow_list_15_16, ~gridDist(.x, target = 1))

out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/S2_15_16/snow_dist_stack/"

purrr::walk2(
  res_15_16$date,
  snow_dist_list_15_16,
  function(d, r) {
    
    name <- format(d, "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_dist_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)

#S2 season 16_17
base_dir <- "D:/antarctica/data/raw_data/S2/16_17/"
safe_dirs_16_17 <- list.files(
  base_dir,
  pattern = "\\.SAFE$",
  full.names = TRUE
)

res_16_17 <- map_dfr(safe_dirs_16_17, get_snow_layer_NDSI_ST,
                     outline = Taylor_outline)
#write Th table
df_out <- res_16_17 %>%
  dplyr::select(date, threshold)

write.csv(
  df_out,
  file = file.path("D:/antarctica/data/processed_data/snow_stacks/Thresholds_NDSI/", "thresholds_2016_2017.csv"),
  row.names = FALSE
)

snow_list_16_17 <- res_16_17$snow
#save snow layers
out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/S2_16_17/snow_stack/"

purrr::walk2(
  res_16_17$date,
  res_16_17$snow,
  function(d, r) {
    
    name <- format(d, "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)

snow_dist_list_16_17 <- map(snow_list_16_17, ~gridDist(.x, target = 1))

out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/S2_16_17/snow_dist_stack/"

purrr::walk2(
  res_16_17$date,
  snow_dist_list_16_17,
  function(d, r) {
    
    name <- format(d, "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_dist_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)

#S2 season 17_18
base_dir <- "D:/antarctica/data/raw_data/S2/17_18/"

zip_files <- list.files(base_dir, pattern = "\\.zip$", full.names = TRUE)
purrr::walk(zip_files, unzip, exdir = base_dir)

safe_dirs_17_18 <- list.files(
  base_dir,
  pattern = "\\.SAFE$",
  full.names = TRUE
)

res_17_18 <- map_dfr(safe_dirs_17_18, get_snow_layer_NDSI_ST,
                     outline = Taylor_outline)
#write Th table
df_out <- res_17_18 %>%
  dplyr::select(date, threshold)

write.csv(
  df_out,
  file = file.path("D:/antarctica/data/processed_data/snow_stacks/Thresholds_NDSI/", "thresholds_2017_2018.csv"),
  row.names = FALSE
)

snow_list_17_18 <- res_17_18$snow
#save snow layers
out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/S2_17_18/snow_stack/"

purrr::walk2(
  res_17_18$date,
  res_17_18$snow,
  function(d, r) {
    
    name <- format(d, "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)

snow_dist_list_17_18 <- map(snow_list_17_18, ~gridDist(.x, target = 1))

out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/S2_17_18/snow_dist_stack/"

purrr::walk2(
  res_17_18$date,
  snow_dist_list_17_18,
  function(d, r) {
    
    name <- format(d, "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_dist_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)

#S2 season 18_19
base_dir <- "D:/antarctica/data/raw_data/S2/18_19/"

zip_files <- list.files(base_dir, pattern = "\\.zip$", full.names = TRUE)
purrr::walk(zip_files, unzip, exdir = base_dir)

safe_dirs_18_19 <- list.files(
  base_dir,
  pattern = "\\.SAFE$",
  full.names = TRUE
)

res_18_19 <- map_dfr(safe_dirs_18_19, get_snow_layer_NDSI_ST,
                     outline = Taylor_outline)
#write Th table
df_out <- res_18_19 %>%
  dplyr::select(date, threshold)

write.csv(
  df_out,
  file = file.path("D:/antarctica/data/processed_data/snow_stacks/Thresholds_NDSI/", "thresholds_2018_2019.csv"),
  row.names = FALSE
)

snow_list_18_19 <- res_18_19$snow
#save snow layers
out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/S2_18_19/snow_stack/"

purrr::walk2(
  res_18_19$date,
  res_18_19$snow,
  function(d, r) {
    
    name <- format(d, "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)

snow_dist_list_18_19 <- map(snow_list_18_19, ~gridDist(.x, target = 1))

out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/S2_18_19/snow_dist_stack/"

purrr::walk2(
  res_18_19$date,
  snow_dist_list_18_19,
  function(d, r) {
    
    name <- format(d, "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_dist_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)

#S2 season 19_20
base_dir <- "D:/antarctica/data/raw_data/S2/19_20/"

zip_files <- list.files(base_dir, pattern = "\\.zip$", full.names = TRUE)
purrr::walk(zip_files, unzip, exdir = base_dir)

safe_dirs_19_20 <- list.files(
  base_dir,
  pattern = "\\.SAFE$",
  full.names = TRUE
)

res_19_20 <- map_dfr(safe_dirs_19_20, get_snow_layer_NDSI_ST,
                     outline = Taylor_outline)
#write Th table
df_out <- res_19_20 %>%
  dplyr::select(date, threshold)

write.csv(
  df_out,
  file = file.path("D:/antarctica/data/processed_data/snow_stacks/Thresholds_NDSI/", "thresholds_2019_2020.csv"),
  row.names = FALSE
)

snow_list_19_20 <- res_19_20$snow
#save snow layers
out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/S2_19_20/snow_stack/"

purrr::walk2(
  res_19_20$date,
  res_19_20$snow,
  function(d, r) {
    
    name <- format(d, "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)

snow_dist_list_19_20 <- map(snow_list_19_20, ~gridDist(.x, target = 1))

out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/S2_19_20/snow_dist_stack/"

purrr::walk2(
  res_19_20$date,
  snow_dist_list_19_20,
  function(d, r) {
    
    name <- format(d, "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_dist_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)

#S2 season 20_21
base_dir <- "D:/antarctica/data/raw_data/S2/20_21/"

zip_files <- list.files(base_dir, pattern = "\\.zip$", full.names = TRUE)
purrr::walk(zip_files, unzip, exdir = base_dir)

safe_dirs_20_21 <- list.files(
  base_dir,
  pattern = "\\.SAFE$",
  full.names = TRUE
)

res_20_21 <- map_dfr(safe_dirs_20_21, get_snow_layer_NDSI_ST,
                     outline = Taylor_outline)
#write Th table
df_out <- res_20_21 %>%
  dplyr::select(date, threshold)

write.csv(
  df_out,
  file = file.path("D:/antarctica/data/processed_data/snow_stacks/Thresholds_NDSI/", "thresholds_2020_2021.csv"),
  row.names = FALSE
)

snow_list_20_21 <- res_20_21$snow
#save snow layers
out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/S2_20_21/snow_stack/"

purrr::walk2(
  res_20_21$date,
  res_20_21$snow,
  function(d, r) {
    
    name <- format(d, "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)

snow_dist_list_20_21 <- map(snow_list_20_21, ~gridDist(.x, target = 1))

out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/S2_20_21/snow_dist_stack/"

purrr::walk2(
  res_20_21$date,
  snow_dist_list_20_21,
  function(d, r) {
    
    name <- format(d, "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_dist_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)

#S2 season 24_25
base_dir <- "D:/antarctica/data/raw_data/S2/24_25/"

zip_files <- list.files(base_dir, pattern = "\\.zip$", full.names = TRUE)
purrr::walk(zip_files, unzip, exdir = base_dir)

safe_dirs_24_25 <- list.files(
  base_dir,
  pattern = "\\.SAFE$",
  full.names = TRUE
)

res_24_25 <- map_dfr(safe_dirs_24_25, get_snow_layer_NDSI_ST,
                     outline = Taylor_outline)
#write Th table
df_out <- res_24_25 %>%
  dplyr::select(date, threshold)

write.csv(
  df_out,
  file = file.path("D:/antarctica/data/processed_data/snow_stacks/Thresholds_NDSI/", "thresholds_2024_2025.csv"),
  row.names = FALSE
)

snow_list_24_25 <- res_24_25$snow
#save snow layers
out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/S2_24_25/snow_stack/"

purrr::walk2(
  res_24_25$date,
  res_24_25$snow,
  function(d, r) {
    
    name <- format(d, "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)

snow_dist_list_24_25 <- map(snow_list_24_25, ~gridDist(.x, target = 1))

out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/S2_24_25/snow_dist_stack/"

purrr::walk2(
  res_24_25$date,
  snow_dist_list_24_25,
  function(d, r) {
    
    name <- format(d, "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_dist_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)

#S2 season 25_26
base_dir <- "D:/antarctica/data/raw_data/S2/25_26/"

zip_files <- list.files(base_dir, pattern = "\\.zip$", full.names = TRUE)
purrr::walk(zip_files, unzip, exdir = base_dir)

safe_dirs_25_26 <- list.files(
  base_dir,
  pattern = "\\.SAFE$",
  full.names = TRUE
)

res_25_26 <- map_dfr(safe_dirs_25_26, get_snow_layer_NDSI_ST,
                     outline = Taylor_outline)
#write Th table
df_out <- res_25_26 %>%
  dplyr::select(date, threshold)

write.csv(
  df_out,
  file = file.path("D:/antarctica/data/processed_data/snow_stacks/Thresholds_NDSI/", "thresholds_2025_2026.csv"),
  row.names = FALSE
)

snow_list_25_26 <- res_25_26$snow
#save snow layers
out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/S2_25_26/snow_stack/"

purrr::walk2(
  res_25_26$date,
  res_25_26$snow,
  function(d, r) {
    
    name <- format(d, "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)

snow_dist_list_25_26 <- map(snow_list_25_26, ~gridDist(.x, target = 1))

out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/S2_25_26/snow_dist_stack/"

purrr::walk2(
  res_25_26$date,
  snow_dist_list_25_26,
  function(d, r) {
    
    name <- format(d, "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_dist_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)
#########################################################################
#snow dist PlanetScope


raster_path_taylor23_24 <- list.files("D:/antarctica/data/raw_data/taylor23_24_psscene_analytic_sr_udm2/", pattern = "composite.tif$" ,full.names = T)
raster_names_taylor23_24 <- list.files("D:/antarctica/data/raw_data/taylor23_24_psscene_analytic_sr_udm2/", pattern = "composite.tif$" ,full.names = F)
dates_taylor_23_24 <- unique(substring(raster_names_taylor23_24, first = 1, last = 10))

get_snow_layer_BST <- function(raster_path, raster_name){
  
  r <- rast(raster_path)
  
  # blue band = first layer
  raster1 <- r[[1]]
  
  raster1 <- raster1 / 10000
  names(raster1) <- "blue_norm"
  
  raster1 <- blue_snow_threshold(raster1, blue = "blue_norm")
  
  snow <- raster1[["snow"]]
  
  name <- substring(raster_name, first = 1, last = 10)
  time(snow) <- ymd(name)
  
  message(name, " done")
  
  return(snow)
}

#merge scenes from the same date if we had more than one scene per day
df <- data.frame(
  raster_path = raster_path_taylor23_24,
  raster_name = raster_names_taylor23_24,
  date = substr(raster_names_taylor23_24, 1, 10)
)

df <- df |>
  split(~ date) |>
  lapply(function(x){
    
    if(nrow(x) == 1) return(x)
    
    tmp <- tempfile(fileext = ".tif")
    
    writeRaster(
      do.call(merge, lapply(x$raster_path, rast)),
      tmp,
      overwrite = TRUE
    )
    
    x[1, "raster_path"] <- tmp
    x[1, ]
  }) |>
  do.call(rbind, args = _)


#apply to list
snow_list_PS_23_24 <- map2(df$raster_path,df$raster_name, get_snow_layer_BST)

snow_dist_list_PS_23_24 <- map(snow_list_PS_23_24, ~gridDist(.x, target = 1))

#save snow layers
out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/PS_23_24/snow_stack/"

purrr::walk2(
 dates_taylor_23_24,
  snow_list_PS_23_24,
  function(d, r) {
    
    name <- format(as.Date(d), "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)


out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/PS_23_24/snow_dist_stack/"

purrr::walk2(
  dates_taylor_23_24,
  snow_dist_list_PS_23_24,
  function(d, r) {
    
    name <- format(as.Date(d), "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_dist_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)

##21_23 rest
raster_path_taylor21_23 <- list.files("D:/antarctica/data/raw_data/taylor21_23_rest_psscene_analytic_sr_udm2/", pattern = "composite.tif$" ,full.names = T)
raster_names_taylor21_23 <- list.files("D:/antarctica/data/raw_data/taylor21_23_rest_psscene_analytic_sr_udm2/", pattern = "composite.tif$" ,full.names = F)
dates_taylor_21_23 <- unique(substring(raster_names_taylor21_23, first = 1, last = 10))

#merge scenes from the same date if we had more than one scene per day
df <- data.frame(
  raster_path = raster_path_taylor21_23,
  raster_name = raster_names_taylor21_23,
  date = substr(raster_names_taylor21_23, 1, 10)
)

df <- df |>
  split(~ date) |>
  lapply(function(x){
    
    if(nrow(x) == 1) return(x)
    
    tmp <- tempfile(fileext = ".tif")
    
    writeRaster(
      do.call(merge, lapply(x$raster_path, rast)),
      tmp,
      overwrite = TRUE
    )
    
    x[1, "raster_path"] <- tmp
    x[1, ]
  }) |>
  do.call(rbind, args = _)


#apply to list
snow_list_PS_21_23 <- map2(df$raster_path,df$raster_name, get_snow_layer_BST)

snow_dist_list_PS_21_23 <- map(snow_list_PS_21_23, ~gridDist(.x, target = 1))

#save snow layers
out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/PS_21_23_rest/snow_stack/"

purrr::walk2(
  dates_taylor_21_23,
  snow_list_PS_21_23,
  function(d, r) {
    
    name <- format(as.Date(d), "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)


out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/PS_21_23_rest/snow_dist_stack/"

purrr::walk2(
  dates_taylor_21_23,
  snow_dist_list_PS_21_23,
  function(d, r) {
    
    name <- format(as.Date(d), "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_dist_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)

## PS classic addition
raster_path_taylorclassic <- list.files("D:/antarctica/data/raw_data/PS_classic_17_18_18_19_first21_psscene_analytic_sr_udm2/", pattern = "composite.tif$" ,full.names = T)
raster_names_taylorclassic <- list.files("D:/antarctica/data/raw_data/PS_classic_17_18_18_19_first21_psscene_analytic_sr_udm2/", pattern = "composite.tif$" ,full.names = F)
dates_taylor_classic <- unique(substring(raster_names_taylorclassic, first = 1, last = 10))

#merge scenes from the same date if we had more than one scene per day
df <- data.frame(
  raster_path = raster_path_taylorclassic,
  raster_name = raster_names_taylorclassic,
  date = substr(raster_names_taylorclassic, 1, 10)
)

df <- df |>
  split(~ date) |>
  lapply(function(x){
    
    if(nrow(x) == 1) return(x)
    
    tmp <- tempfile(fileext = ".tif")
    
    writeRaster(
      do.call(merge, lapply(x$raster_path, rast)),
      tmp,
      overwrite = TRUE
    )
    
    x[1, "raster_path"] <- tmp
    x[1, ]
  }) |>
  do.call(rbind, args = _)


#apply to list
snow_list_PS_classic <- map2(df$raster_path,df$raster_name, get_snow_layer_BST)

snow_dist_list_PS_classic <- map(snow_list_PS_classic, ~gridDist(.x, target = 1))

#save snow layers
out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/PS_classic/snow_stack/"

purrr::walk2(
  dates_taylor_classic,
  snow_list_PS_classic,
  function(d, r) {
    
    name <- format(as.Date(d), "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)


out_dir_raster <- "D:/antarctica/data/processed_data/snow_stacks/PS_classic/snow_dist_stack/"

purrr::walk2(
  dates_taylor_classic,
  snow_dist_list_PS_classic,
  function(d, r) {
    
    name <- format(as.Date(d), "%Y%m%d") 
    out_file <- file.path(
      out_dir_raster,
      paste0("snow_dist_", name, ".tif")
    )
    
    writeRaster(r, out_file, overwrite = TRUE)
  }
)

