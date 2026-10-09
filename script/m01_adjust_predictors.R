#reclassify snow_dist to overcome AOA problem
library(purrr)
library(terra)
library(lubridate)
library(sf)
#test_stack <- rast("data/processed_data/predictor_stacks_AOI_mask_tif/pred_AOI_20211226.tif")
reclass_matrix <- matrix(
  c(0, 10,  1,
    10, 20,  2,
    20, 30, 3,
    30,40,  4,
    40,50,  5,
    50, 60, 6,
    60, 70, 7,
    70, 80, 8,
    80, 90, 9,
    90, 100, 10, 
    100, 110, 11,
    110, 120, 12,
    120, 130, 13, 
    130, 140, 14, 
    140, 150, 15, 
    150, 160, 16,
    160, 170, 17, 
    170, 180, 18, 
    180, 190, 19,
    190, 200, 20,
    200, 210, 21,
    210,220, 22,
    220, 230, 23,
    230, 240, 24,
    240, 250, 25,
    250,Inf,  26
  ),
  ncol = 3,
  byrow = TRUE
)

#test_stack$snow_dist_class <- classify(test_stack$snow_dist,
  #                                     reclass_matrix,
  #                                     include.lowest = T,
  #                                     right = F)

snow_class_function <- function(stack, matrix = reclass_matrix){
  new_layer <- classify(stack[[1]],                # take $snow_dist if stack contains more than one layer
                        matrix,
                        include.lowest = TRUE,
                        right = FALSE)
  
  names(new_layer) <- "snow_dist_class"
  stack <- c(stack, new_layer)
  
  return(stack)
}

pred_stack_path <- list.files("data/processed_data/predictor_stacks_AOI_mask_tif/", full.names = T)
pred_stack_names <- list.files("data/processed_data/predictor_stacks_AOI_mask_tif/", full.names = F)
pred_stack_names <- substr(pred_stack_names, start = 1, stop = 17)
pred_stack_list <- map(pred_stack_path, rast)
adj_stack_list <- map(pred_stack_list, snow_class_function)

#### classify snow distances for all time steps
season_path <- list.files("data/processed_data/snow_stacks/", full.names = T)
season_path <- season_path[-length(season_path)]

# single tif
process_snow_tif <- function(tif_path, season_dir, matrix = reclass_matrix) {
  
  stack <- rast(tif_path)
  
  comb_stack <- snow_class_function(
    stack,
    matrix = matrix
  )
  
  name <- basename(tif_path)
  date <- substr(name, start  = 11, stop = 18)
  
  outdir <- file.path(
    season_dir,
    "snow_class_stack"
  )
  
  dir.create(
    outdir,
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  outfile <- file.path(
    outdir,
    paste0(
      "snow_class_stack_",
      date,
      ".tif"
    )
  )
  
  writeRaster(
    comb_stack,
    outfile,
    overwrite = TRUE
  )
  
  rm(stack, comb_stack)
  gc()
  
  message("Wrote: ", date)
  
  invisible(NULL)
}


#seasons 
process_season <- function(season_dir, matrix = reclass_matrix) {
  
  tif_paths <- list.files(
    file.path(season_dir, "snow_dist_stack"),
    pattern = "\\.tif$",
    full.names = TRUE
  )
  
  walk(
    tif_paths,
    process_snow_tif,
    season_dir = season_dir,
    matrix = matrix
  )
  
  message("Finished season: ",
          basename(season_dir))
  
  invisible(NULL)
}

#try first season
walk(season_path[1], process_season)
# run all other seasons

walk(season_path[2:length(season_path)], process_season)
###########################

####

#add daily mean temperature from the Lake Fryxell LTER station
dates <- substr(pred_stack_names, start = 10, stop = 17)
dates <- ymd(dates)
met_data <- read.csv("data/raw_data/met_data_LTER/mcmlter-clim-frlm_daily-2026-02-28.csv", header = T)
met_data$date <- as.Date(met_data$date_time)

met_data_sub <- met_data[met_data$date %in% dates,]
air_data <- met_data_sub$airtemp_3m_degc
air_data_list <- as.list(air_data)

air_layer_function <- function(stack, met_data){
  template <- stack[[1]]
  values(template) <- met_data
  names(template) <- "air_temp"
  out <- c(stack, template)
return(out)
}

#
air_stack_list <- map2(adj_stack_list, air_data_list, air_layer_function)


#
names(adj_stack_list) <- pred_stack_names

outdir <- "data/processed_data/predictor_stacks_adj_AOI_tif/"
dir.create(outdir, showWarnings = FALSE)

for (nm in names(air_stack_list)) {
  
  r <- air_stack_list[[nm]]
  
  outpath <- file.path(outdir, paste0("adj_", nm, ".tif"))
  
  writeRaster(r, outpath, overwrite = TRUE)
}


#for model training: 
###additional test set
pred_stack_test_path <- list.files("data/processed_data/planet_test_dataset/test_predictor_stacks_tif/", full.names = T)
pred_stack_test_names <- list.files("data/processed_data/planet_test_dataset/test_predictor_stacks_tif/", full.names = F)
pred_stack_test_names <- substr(pred_stack_test_names, start = 14, stop = 21)
pred_stack_test_list <- map(pred_stack_test_path, rast)
adj_stack_test_list <- map(pred_stack_test_list, snow_class_function)

#add daily mean temperature from the Lake Fryxell LTER station
dates <- pred_stack_test_names
dates <- ymd(dates)
met_data_test <- read.csv("data/raw_data/met_data_LTER/mcmlter-clim-frlm_daily-2026-02-28.csv", header = T)
met_data_test$date <- as.Date(met_data_test$date_time)

met_data_test_sub <- met_data_test[met_data_test$date %in% dates,]
air_data_test <- met_data_test_sub$airtemp_3m_degc
air_data_test_list <- as.list(air_data_test)

air_stack_test_list <- map2(adj_stack_test_list, air_data_test_list, air_layer_function)

#mask to AOI
glaciers <- st_read("C:/Users/mini4/Desktop/Uni/Master/antarctica/data/raw_data/Antarctic_GIS/MDV_GLACIERS_STREAMS_POLYGONS/MDV_GLACIERS.shp")
lakes_ponds <- st_read("C:/Users/mini4/Desktop/Uni/Master/antarctica/data/raw_data/GIS layers/MCMLTER_GIS_Export_Layers/USGS_Map_1970/Lakes_and_Poonds_epsg32758.gpkg")
glaciers <- st_transform(glaciers, crs = st_crs(lakes_ponds))
AOI_mask <- st_union(glaciers, lakes_ponds)
plot(AOI_mask)
AOI_mask <- st_transform(AOI_mask, crs = crs(air_stack_test_list[[1]]))

#
names(air_stack_test_list) <- pred_stack_test_names
air_stack_test_list_AOI<- map(air_stack_test_list, ~ mask (.x, AOI_mask, inverse = T))

outdir <- "data/processed_data/predictor_stacks_test_adj_AOI_tif/"
dir.create(outdir, showWarnings = FALSE)

for (nm in names(air_stack_test_list_AOI)) {
  r <- air_stack_test_list_AOI[[nm]]
  outpath <- file.path(outdir, paste0("adj_test_", nm, ".tif"))
  writeRaster(r, outpath, overwrite = TRUE)
}

####################################
#for all time steps; add met data from Lake Fryxell LTER station and TWI
#load TWI and resample to 10m (Seninel resolution)
TWI <- rast("data/processed_data/DEM_predictors/TWI_1m_coast.tif")
AOI_mask <- st_transform(AOI_mask, crs(TWI))
TWI_c <- crop(TWI, AOI_mask)

TWI_3 <- rast("data/processed_data/DEM_predictors/TWI_3m_coast_.tif")
 
TWI_10_32758 <- project(
  TWI_c,
  "EPSG:32758",
  method = "bilinear"
)


target_10m <- rast(
  ext = ext(TWI_10_32758),
  resolution = 10,
  crs = crs(TWI_10_32758)
)

TWI_10_32758 <- resample(
  TWI_10_32758,
  target_10m,
  method = "bilinear"
)

TWI_10_32758 <- crop(TWI_10_32758, TWI_3)

writeRaster(TWI_10_32758, "data/processed_data/DEM_predictors/TWI_10m_32758_taylor.tif")
TWI_10 <- rast("data/processed_data/DEM_predictors/TWI_10m_32758_taylor.tif")
############################

rad_dir <- "D:/antarctica/data/processed_data/rsun_daily"

glob_files <- list.files(
  rad_dir,
  pattern = "^glob_rad_day_\\d+\\.tif$",
  full.names = TRUE
)

insol_files <- list.files(
  rad_dir,
  pattern = "^insol_day_\\d+\\.tif$",
  full.names = TRUE
)

glob_lookup <- setNames(
  glob_files,
  stringr::str_extract(basename(glob_files), "\\d+")
)

insol_lookup <- setNames(
  insol_files,
  stringr::str_extract(basename(insol_files), "\\d+")
)

process_predictor_season <- function(
    season_dir,
    met_data,
    twi_10,
    twi_3,
    glob_lookup,
    insol_lookup
) {
  
  season_name <- basename(season_dir)
  
  message("Processing season: ", season_name)
  
  # season wise folders (data and output)
  
  stack_dir <- file.path(
    season_dir,
    "snow_class_stack"
  )
  
  outdir <- file.path(
    season_dir,
    "pred_stack"
  )
  
  dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  # list season specific input files (snow_dist_class)
  
  stack_paths <- list.files(
    stack_dir,
    pattern = "\\.tif$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  #check whether its working
  if (length(stack_paths) == 0) {
    message("No tif files found in: ", stack_dir)
    return(NULL)
  }
  
  # adjust met data date format 
  
  met_data$date <- as.Date(met_data$date)
  
  # loop processing for every timestep per season
  
  walk(stack_paths, function(path) {
    
    file_name <- basename(path)
    
    date_chr <- stringr::str_extract(
      file_name,
      "\\d{8}"
    )
    
    outfile <- file.path(
      outdir,
      paste0(
        "predictors_",
        date_chr,
        ".tif"
      )
    )
    
    # skip existing files
    if (file.exists(outfile)) {
      
      message(
        "Skipping existing: ",
        basename(outfile)
      )
      
      return(NULL)
    }
    
    message(
      "Processing: ",
      file_name
    )
    
    #load snow_dist_class
    
    stack <- rast(path)
    
    if (!"snow_dist_class" %in% names(stack)) {
      message("Missing layer: ", file_name)
      return(NULL)
    }
    
    stack <- stack[["snow_dist_class"]]
    
    # use different TWI (3 or 10) depending on the sensor because of the different
    # spatial resolutions
    
    if (grepl("PS", season_name)) {
      
      twi <- twi_3
      
      if (!compareGeom(stack, twi, stopOnError = FALSE)) {
        
        message("Aligning geometry: ", file_name)
        
        stack <- project(
          stack,
          twi,
          method = "near"
        )
      }
    }
    
    if (grepl("S2", season_name)) {
      
      twi <- twi_10
      
      if (!compareGeom(stack, twi, stopOnError = FALSE)) {
        
        message("Aligning geometry: ", file_name)
        
        stack <- project(stack, twi, method = "near")
      }
    }
    
    # identifying day of year (doy) for the insolation layers 
    
    date_chr <- stringr::str_extract(file_name, "\\d{8}")
    
    if (is.na(date_chr)) {
      message("No date found: ", file_name)
      return(NULL)
    }
    
    date_obj <- lubridate::ymd(date_chr)
    
    doy <- as.integer(format(date_obj, "%j"))
    if (doy == 366) doy <- 365 #for leap years
    
    # prepare air temperature layer from Lake Fryxell LTER met station (daily mean)
    
    air_temp <- met_data |>
      dplyr::filter(date == date_obj) |>
      dplyr::pull(airtemp_3m_degc)
    
  
    if (length(air_temp) == 0 || is.na(air_temp[1])) {
      message("No met data for ", date_chr)
      return(NULL)
    }
    
    air_temp <- as.numeric(air_temp[1])
    
    
    # add the corresponding global radiation and daily insolation layers
    
    glob_path <- unname(
      glob_lookup[
        as.character(doy)
      ]
    )
    
    insol_path <- unname(
      insol_lookup[
        as.character(doy)
      ]
    )
    
    # safety check
    if (
      length(glob_path) == 0 ||
      is.na(glob_path) ||
      length(insol_path) == 0 ||
      is.na(insol_path)
    ) {
      
      message(
        "Missing radiation for DOY ",
        doy,
        " (",
        file_name,
        ")"
      )
      
      return(NULL)
    }
    
    glob_rad <- rast(
      glob_path
    )
    
    insol <- rast(
      insol_path
    )
    
    glob_rad <- project(
      glob_rad,
      twi,
      method = "bilinear"
    )
    
    insol <- project(
      insol,
      twi,
      method = "bilinear"
    )
    
    names(glob_rad) <- "daily_glob_rad"
    names(insol) <- "daily_insol"
    
    # add air temp layer
    
    air_layer <- stack
    values(air_layer) <- air_temp
    names(air_layer) <- "air_temp"
    
    # final layer names
    
    names(stack) <- "snow_dist_class"
    names(twi) <- "TWI"
    
    # combine to stack
    
    out_stack <- c(
      stack,
      twi,
      air_layer,
      glob_rad,
      insol
    )
    
    # output name
    
    outfile <- file.path(
      outdir,
      paste0("predictors_", date_chr, ".tif")
    )
    
    writeRaster(
      out_stack,
      outfile,
      overwrite = TRUE
    )
    
    rm(stack, twi, air_layer, glob_rad, insol, out_stack)
    gc() 
    
    message("Wrote: ", outfile)
  })
  
  message("Finished season: ", season_name)
  tmpFiles(remove = TRUE)
  gc()
}

#set tempdir
dir.create(
  "D:/terra_temp",
  showWarnings = FALSE
)

terraOptions(
  tempdir = "D:/terra_temp",
  memfrac = 0.8,
  progress = 10
)


# check for missing radiation layers 

missing_doy <- purrr::map_dfr(
  season_path,
  function(season_dir) {
    
    stack_dir <- file.path(
      season_dir,
      "snow_class_stack"
    )
    
    stack_paths <- list.files(
      stack_dir,
      pattern = "\\.tif$",
      full.names = TRUE,
      ignore.case = TRUE
    )
    
    tibble::tibble(
      season = basename(season_dir),
      file_name = basename(stack_paths)
    ) |>
      dplyr::mutate(
        date_chr = stringr::str_extract(
          file_name,
          "\\d{8}"
        ),
        
        date = lubridate::ymd(
          date_chr
        ),
        
        doy = as.integer(
          format(date, "%j")
        ),
        
        glob_exists = as.character(doy) %in% names(glob_lookup),
        
        insol_exists = as.character(doy) %in% names(insol_lookup)
      ) |>
      dplyr::filter(
        !glob_exists |
          !insol_exists
      )
  }
)

missing_doy

# run
# PS seasons
walk(
  season_path[1:3],
  process_predictor_season,
  met_data = met_data,
  twi_10 = TWI_10,
  twi_3 = TWI_3,
  glob_lookup = glob_lookup,
  insol_lookup = insol_lookup
)

# S2 seasons 
walk(
  season_path[4:length(season_path)],
  process_predictor_season,
  met_data = met_data,
  twi_10 = TWI_10,
  twi_3 = TWI_3,
  glob_lookup = glob_lookup,
  insol_lookup = insol_lookup
)
