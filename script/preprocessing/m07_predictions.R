library(terra)
library(CAST)
library(caret)
library(sf)
library(lubridate)

model_k4_p5 <- readRDS("data/processed_data/model_k4_p5/model_k4_p5_mtry2.RDS")
tdi_k4_p5 = trainDI(model_k4_p5, verbose = FALSE)

glaciers <- st_read("data/raw_data/Antarctic_GIS/MDV_GLACIERS_STREAMS_POLYGONS/MDV_GLACIERS.shp")
lakes_ponds <- st_read("data/raw_data/GIS layers/MCMLTER_GIS_Export_Layers/USGS_Map_1970/Lakes_and_Poonds_epsg32758.gpkg")
glaciers <- st_transform(glaciers, crs = st_crs(lakes_ponds))
AOI_mask <- st_union(glaciers, lakes_ponds)
plot(AOI_mask)
AOI_mask <- st_transform(AOI_mask, crs = "EPSG:32758")

Planet_outline <- st_read("C:/Users/mini4/Desktop/Uni/Master/antarctica/data/processed_data/Planet_small_outline.geojson")
Planet_outline <- st_transform(Planet_outline, crs = "EPSG:32758")

###############
predictors <- c( "daily_insol", "snow_dist_class", 
                 "daily_glob_rad", "TWI", "air_temp")

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

#####################
#LPD, AOA stepwise
outdir_LPD <- "data/processed_data/model_k4_p5/LPD_tif"
dir.create(outdir_LPD, showWarnings = T)

outdir_AOA <- "data/processed_data/model_k4_p5/AOA_tif"
dir.create(outdir_AOA, showWarnings = T)

calc_write_LPD <- function(pred, model_daily, tdi_daily, nm){
  message(nm)
  r <- aoa(pred, model_daily, trainDI = tdi_daily, LPD = T)
  LPD <- r$LPD$LPD
  AOA <- r$AOA$AOA
  outpath_LPD <- file.path(outdir_LPD, paste0("lpd_",nm,".tif"))
  writeRaster(LPD, outpath_LPD, overwrite = T)
  outpath_AOA <- file.path(outdir_AOA, paste0("aoa_",nm,".tif"))
  writeRaster(AOA, outpath_AOA, overwrite = T)
  rm(LPD)
  gc()
  message("wrote ", nm)
}

purrr::walk2(
  all_stacks_AOI,names(all_stacks_AOI), ~
    calc_write_LPD(.x,
                   model_daily = model_daily,
                   tdi_daily   = tdi_daily,
                   nm = .y)
)



#mask all_stack_AOI to indiv aoa
aoa_path <- list.files("data/processed_data/model_k4_p5/AOA_tif", full.names = T)
aoa_list <- map(aoa_path, rast)
names(aoa_list) <- dates

all_stacks_AOI_pred <- map(all_stacks_AOI,~ .x[[predictors]])
all_stacks_AOI_AOA <- map2(all_stacks_AOI_pred, aoa_list, ~ mask(.x, .y$AOA, maskvalues = 0)) 

traindat_sf_sel_pred <- traindat_sf %>%
  dplyr::select(
    daily_insol,
    snow_dist_class,
    daily_glob_rad,
    TWI,
    air_temp,
    geom  
  )



#for whole list 
dist_list <- map(all_stacks_AOI_AOA, ~ geodist(traindat_sf_sel_pred,.x,cvfolds=trainids$indexOut, type = "feature"))
plot(dist_list[[1]])+ scale_x_log10(labels=round)
names(dist_list)

saveRDS(dist_list,"data/processed_data/model_k4_p5/dist_list.RDS")
#save all feature space distance plots- not working automatically 

out_dir_fd <- "data/processed_data/model_k4_p5/feature_dist_plots/"
dir.create(out_dir_fd, showWarnings = T)


for (nm in names(dist_list)) {
  file_name <- file.path(out_dir_fd, paste0(nm, "_fd.png"))
  png(filename = file_name, width = 800, height = 600)
  x <- dist_list[[nm]]
  print(plot(x))
  dev.off()
}

#predictions
prediction_list <- map(all_stacks_AOI, ~predict(.x, model_daily, na.rm = T))
#save predictions
outdir <- "data/processed_data/model_k4_p5/pred_in_AOI/"
dir.create(outdir, showWarnings = T)

for (nm in names(prediction_list)) {
  
  r <- prediction_list[[nm]]
  
  outpath <- file.path(outdir, paste0("pred_in_AOI_", nm, ".tif"))
  
  writeRaster(r, outpath, overwrite = TRUE)
}


#pred in AOI and AOA
#predictions in AOI and AOrast#predictions in AOI and AOA
pred_aoa_list <- map2(prediction_list, aoa_list, ~ {
  x <- .x
  x[.y$AOA == 0] <- NA
  x
})

names(pred_aoa_list)
#save predictions in AOI and AOA
outdir <- "data/processed_data/model_k4_p5/pred_in_AOI_and_AOA/"
dir.create(outdir, showWarnings = T)

for (nm in names(pred_aoa_list)) {
  
  r <- pred_aoa_list[[nm]]
  
  outpath <- file.path(outdir, paste0("pred_in_AOI_AOA", nm, ".tif"))
  
  writeRaster(r, outpath, overwrite = TRUE)
}


###############################################

process_prediction_stack <- function(
    stack_path,
    model,
    tdi,
    AOI,
    outline,
    out_base = "data/processed_data/predictions"
) {
  
  file_name <- basename(stack_path)
  
  message("Processing: ", file_name)
  
  # date and sensor
  
  date_chr <- stringr::str_extract(
    file_name,
    "\\d{8}"
  )
  
  date_obj <- lubridate::ymd(
    date_chr
  )
  
  sensor <- dplyr::case_when(
    grepl("S2", stack_path) ~ "S2",
    grepl("PS", stack_path) ~ "PS",
    TRUE ~ "UNK"
  )
  
  # season
  
  yr <- lubridate::year(date_obj)
  mo <- lubridate::month(date_obj)
  
  season <- if (mo >= 11) {
    
    paste0(
      substr(yr, 3, 4),
      "_",
      substr(yr + 1, 3, 4)
    )
    
  } else {
    
    paste0(
      substr(yr - 1, 3, 4),
      "_",
      substr(yr, 3, 4)
    )
  }
  
  # outdir
  
  outdir_pred <- file.path(
    out_base,
    season,
    "prediction"
  )
  
  outdir_pred_AOA <- file.path(
    out_base,
    season,
    "prediction_AOA"
  )
  
  outdir_AOA <- file.path(
    out_base,
    season,
    "AOA"
  )
  
  outdir_LPD <- file.path(
    out_base,
    season,
    "LPD"
  )
  
  walk(
    list(
      outdir_pred,
      outdir_pred_AOA,
      outdir_AOA,
      outdir_LPD
    ),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # names
  
  prefix <- paste0(
    sensor,
    "_",
    date_chr
  )
  
  out_pred <- file.path(
    outdir_pred,
    paste0(prefix, "_pred.tif")
  )
  
  out_pred_AOA <- file.path(
    outdir_pred_AOA,
    paste0(prefix, "_pred_AOA.tif")
  )
  
  out_AOA <- file.path(
    outdir_AOA,
    paste0(prefix, "_AOA.tif")
  )
  
  out_LPD <- file.path(
    outdir_LPD,
    paste0(prefix, "_LPD.tif")
  )
  
  # skip existing
  if (
    file.exists(out_pred) &&
    file.exists(out_pred_AOA) &&
    file.exists(out_AOA) &&
    file.exists(out_LPD)
  ) {
    
    message("Skipping existing: ", prefix)
    return(NULL)
  }
  
  # load stack
  
  pred_stack <- rast(
    stack_path
  )
  
  pred_stack <- mask(
    pred_stack,
    AOI,
    inverse = T
  )
  
  pred_stack <- mask(
    pred_stack,
    outline
  )
  # prediction
  
  pred <- terra::predict(
    pred_stack,
    model,
    na.rm = TRUE
  )
  
  names(pred) <- "prediction"
  
  # AOA/LPD
  
  aoa_res <- aoa(
    pred_stack,
    model,
    trainDI = tdi,
    LPD = TRUE
  )
  
  AOA <- aoa_res$AOA$AOA
  LPD <- aoa_res$LPD$LPD
  
  pred_AOA <- mask(pred, AOA, maskvalues = 0)
  
  # save
  
  writeRaster(
    pred,
    out_pred,
    overwrite = TRUE
  )
  
  writeRaster(
    pred_AOA,
    out_pred_AOA,
    overwrite = TRUE
  )
  
  writeRaster(
    AOA,
    out_AOA,
    overwrite = TRUE
  )
  
  writeRaster(
    LPD,
    out_LPD,
    overwrite = TRUE
  )
  
  rm(
    pred_stack,
    pred,
    pred_AOA,
    AOA,
    LPD,
    aoa_res
  )
  
  gc()
  
  message(
    "Wrote: ",
    prefix
  )
}


stack_paths <- list.files(
  "data/processed_data",
  pattern = "predictors_.*\\.tif$",
  recursive = TRUE,
  full.names = TRUE
)

walk(
  stack_paths[1],
  process_prediction_stack,
  model = model_k4_p5,
  tdi = tdi_k4_p5,
  AOI = AOI_mask,
  outline = Planet_outline 
)
