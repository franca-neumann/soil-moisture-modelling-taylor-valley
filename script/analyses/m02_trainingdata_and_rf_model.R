library(CAST)
library(caret)
library(terra)
library(MLmetrics)
library(fastshap)
library(ggbeeswarm)
library(nestedcv)
library(sf)
library(dplyr)
library(lubridate)
library(tidyr)
library(purrr)
library(Cairo)

VWC_daily_sf <- st_read("C:/Users/mini4/Desktop/Uni/Master/antarctica/data/processed_data/training_data/planet_daily_mean.gpkg")

excluded_sensor <- st_read("C:/Users/mini4/Desktop/Uni/Master/antarctica/data/processed_data/training_data/final_combined_trainig_data.gpkg")
setdiff(excluded_sensor$VWC_daily,VWC_daily_sf$VWC_daily) #this is the correct data, exposed sensor data is already excluded

#extract training data for new stacks, directly as data frame  
VWC_daily_sf$date <- as.Date(VWC_daily_sf$date)
dates <- unique(VWC_daily_sf$date)

result_list <- lapply(dates, function(d) {
 
  pts <- VWC_daily_sf %>% dplyr::filter(date == d)
  
  d_str <- format(d, "%Y%m%d")
  
  r <- rast(paste0("D:/antarctica/data/processed_data/predictor_stacks_adj_AOI_tif/adj_pred_AOI_", d_str, ".tif"))

  pts <- st_transform(pts, crs(r))
  pts <- st_as_sf(pts)
  
  vals <- terra::extract(r, vect(pts))
  vals <- vals[, -1]
  dplyr::bind_cols(st_drop_geometry(pts), vals)
})

final_df <- bind_rows(result_list)
saveRDS(final_df,"D:/antarctica/data/processed_data/training_data/traindat_VWC_mean.RDS")

# again as option with geometry
result_list <- lapply(dates, function(d) {
  pts <- dplyr::filter(VWC_daily_sf, date == d)
  d_str <- format(d, "%Y%m%d")
  r <- rast(paste0("D:/antarctica/data/processed_data/predictor_stacks_adj_AOI_tif/adj_pred_AOI_", d_str, ".tif"))
  pts <- st_transform(pts, crs(r))
  vals <- terra::extract(r, pts)
  vals <- vals[, -1]  #remove ID
  cbind(pts, vals)
})

final_sf <- do.call(rbind, result_list)
st_write(final_sf,"D:/antarctica/data/processed_data/training_data/traindat_VWC_mean.gpkg", delete_dsn = TRUE)

###################
traindat <- final_df
traindat <- readRDS("D:/antarctica/data/processed_data/training_data/traindat_VWC_mean.RDS")
traindat <- traindat[!is.na(traindat$VWC_daily), ]

traindat_sf <- st_read("D:/antarctica/data/processed_data/training_data/traindat_VWC_mean.gpkg")
traindat_sf <- traindat_sf[!is.na(traindat_sf$VWC_daily), ]
#rf Modelle
set.seed(123)
###########
trainids <- CreateSpacetimeFolds(traindat, 
                                 spacevar = "site_id",
                                 k = 4,
                                 seed = 123)

ctrl <- trainControl(method = "cv",
                     index = trainids$index,
                     indexOut = trainids$indexOut,
                     summaryFunction = defaultSummary,
                     savePredictions = T)
selected_vars <- c("NDWI", "daily_insol", "snow_dist_class", 
                "daily_glob_rad", "TWI", "air_temp", "BI_imean", "BI_omean")

model_ffs_daily <- ffs(traindat[,selected_vars],
                            traindat$VWC_daily,
                            importance = TRUE,
                            metric = "RMSE",
                            tuneLength = 2,
                            globalval = T,
                            ntree=50,
                            method = "rf",
                            trControl = ctrl,
                            na.action = "na.omit")
model_ffs_daily
global_validation(model_ffs_daily)
plot(varImp(model_ffs_daily, scale = F))

predictors <- c( "daily_insol", "snow_dist_class", 
                "daily_glob_rad", "TWI", "air_temp")

model_daily<- train(traindat[, predictors],
                    traindat$VWC_daily,
                    method= "rf",
                    importance= TRUE,
                    ntree=500,
                    nodesize = 5,
                    tuneGrid = data.frame("mtry" = c(2:5)),
                    trControl=ctrl)


model_daily
global_validation(model_daily)
plot(varImp(model_daily, scale = F))
varImp(model_daily, scale = F)
saveRDS(model_daily, "D:/antarctica/data/processed_data/model_k4_p5/model_k4_p5_mtry2.RDS")
model_daily <- readRDS("D:/antarctica/data/processed_data/model_k4_p5/model_k4_p5_mtry2.RDS")
#ranger option
#model_daily<- train(traindat[, predictors],
                  #  traindat$VWC_daily,
                  #  method= "ranger",
                  #  importance= "impurity",
                  #  seed = 42,
                   # num.trees=500,
                  #  tuneGrid = expand.grid(
                  #    mtry = 2:5,
                  #    splitrule = "variance",
                  #    min.node.size = 2:5
                  #  ),
                  #  trControl=ctrl)
####no real improvement; dependent on the seed if better or not, inconsistent 
#results for min.node.size, therefore using the original "rf" model instead of ranger

#for exploration only 
par(mfrow = c(1,1))
boxplot(model_daily$resample$Rsquared,
        main = "R² from model$resample$Rsquared")


df_metrics <- model_daily$resample %>%
  dplyr::select(MAE, RMSE) %>%
  pivot_longer(cols = everything(), names_to = "metric", values_to = "value")

ggplot(df_metrics, aes(x = metric, y = value)) +
  geom_boxplot() +
  theme_bw(base_size = 14) +
  labs(title = "Final Model performance metrics",
       subtitle = "Retrieved from (model$resample$metric)",
       x = "Metric",
       y = "")

#Shapley-values

## Define prediction wrapper function
predict_function <- function(object, newdata) {
  predict(object, newdata = newdata)
}

## Compute Shapley values

shap_values <- fastshap::explain(
  object = model_daily,
  X = traindat[, predictors],
  pred_wrapper = predict_function,
  nsim = 100,               # Number of permutations
  adjust = TRUE            
)

### Beeswarm plotting
plot_shap_beeswarm(shap_values, traindat[,predictors], size = 1)

### save Beeswarm plot as png

CairoPNG(
  "plots/shap_beeswarm_daily_model.png",
  width = 6,
  height = 5,
  units = "in",
  dpi = 600
)

plot_shap_beeswarm(shap_values, traindat[, predictors], size = 1)

dev.off()
######################################################################
#AOA
all_stacks_AOI_files <- list.files("D:/antarctica/data/processed_data/predictor_stacks_adj_AOI_tif/", full.names = T)
all_stacks_AOI <- map(all_stacks_AOI_files, rast)
plot(all_stacks_AOI[[1]])

dates <- unique(traindat$date)
dates <- as.numeric(gsub("-", "", dates))

names(all_stacks_AOI) <- dates 
all_stacks_AOI[[1]]
names(all_stacks_AOI)

tdi_daily = trainDI(model_daily, verbose = FALSE)

#try out
#aoa_daily <- aoa(all_stacks_AOI[[5]], model_daily, trainDI = tdi_daily, LPD = F)
#aoa_daily$AOA
#aoa_list <- map(all_stacks_AOI, ~ aoa(.x, model_daily, trainDI = tdi_daily, LPD = F))
#names(aoa_list)

#save aoa's

#for (nm in names(aoa_list)) {
  
#  r <- aoa_list[[nm]]
  
#  outpath <- file.path(outdir, paste0("aoa", nm, ".tif"))
  
#  writeRaster(r$AOA, outpath, overwrite = TRUE)
# }


#LPD, AOA stepwise
outdir_LPD <- "D:/antarctica/data/processed_data/model_k4_p5/LPD_tif"
dir.create(outdir_LPD, showWarnings = T)

outdir_AOA <- "D:/antarctica/data/processed_data/model_k4_p5/AOA_tif"
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


#error profile is static; just once for the model
errorPro_daily <- errorProfiles(model_daily, trainDI = tdi_daily, variable = "DI")
errorPro_daily
plot(errorPro_daily)



#mask all_stack_AOI to indiv aoa
aoa_path <- list.files("D:/antarctica/data/processed_data/model_k4_p5/AOA_tif", full.names = T)
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

saveRDS(dist_list,"D:/antarctica/data/processed_data/model_k4_p5/dist_list.RDS")
#save all feature space distance plots- not working automatically 

out_dir_fd <- "D:/antarctica/data/processed_data/model_k4_p5/feature_dist_plots/"
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
outdir <- "D:/antarctica/data/processed_data/model_k4_p5/pred_in_AOI/"
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
outdir <- "D:/antarctica/data/processed_data/model_k4_p5/pred_in_AOI_and_AOA/"
dir.create(outdir, showWarnings = T)

for (nm in names(pred_aoa_list)) {
  
  r <- pred_aoa_list[[nm]]
  
  outpath <- file.path(outdir, paste0("pred_in_AOI_AOA", nm, ".tif"))
  
  writeRaster(r, outpath, overwrite = TRUE)
}

####################
library(terra)
library(magick)
library(viridisLite)


r_list <- pred_aoa_list   
dates  <- names(r_list)
dates <- ymd(dates)
dates <- format(dates, "%d.%m.%Y")  
dates
# one legend/scale for all layers
ncol <- 200
cols <- viridis(ncol)

mins <- map_dbl(r_list, ~ {
  gmin <- terra::global(.x, "min", na.rm = TRUE)
  min(unlist(gmin))
})

maxs <- map_dbl(r_list, ~ {
  gmax <- terra::global(.x, "max", na.rm = TRUE)
  max(unlist(gmax))
})

all_min <- min(mins)
all_max <- max(maxs)

outdir <- "D:/antarctica/data/processed_data/model_k4_p5/gif_frames_pred_in_aoa_aoi/"
dir.create(outdir, showWarnings = FALSE)

frames <- list()

#PNG

for (i in seq_along(r_list)) {
  
  r <- r_list[[i]]
  date_label <- dates[i]
  
  outfile <- file.path(outdir, sprintf("frame_%03d.png", i))
  
  png(outfile, width=900, height=800)
  par(mar=c(4,4,5,6))
  
  plot(r,
       col = cols,
       main = date_label,
       range = c(all_min, all_max),
       axes = TRUE,
       cex.main = 2)
  
  dev.off()
  
  frames[[i]] <- magick::image_read(outfile)
}

#gif
gif <- image_animate(image_join(frames), fps = 2)

image_write(gif, "D:/antarctica/data/processed_data/model_k4_p5/pred_in_aoa_aoi_animationfps.gif")

##################
par(mfrow = c(1,1))
values_pred_d <- values(pred_aoa_d)

df_plot_d <- bind_rows(
  tibble(group = "VWC traindat", value = traindat$VWC_daily),
  tibble(group = "prediction",  value = values_pred_d)
)

ggplot(df_plot_d, aes(x = group, y = value)) +
  geom_boxplot() +
  xlab("") +
  ylab("") +
  theme_minimal()


quantile(values_pred_d, na.rm = T)
quantile(traindat$VWC_daily)

##################################################################
#additional scenes were only P4 is available as additional test set
VWC_daily_test_sf <- st_read("data/processed_data/planet_test_dataset/planet_test_daily_mean_VWC_only.gpkg")

# again as option with geometry
result_list <- lapply(dates, function(d) {
  pts <- dplyr::filter(VWC_daily_test_sf, date == d)
  d_str <- format(d, "%Y%m%d")
  r <- rast(paste0("D:/antarctica/data/processed_data/predictor_stacks_test_adj_AOI_tif/","adj_test_", d_str, ".tif"))
  pts <- st_transform(pts, crs(r))
  vals <- terra::extract(r, pts)
  vals <- vals[, -1]  #remove ID
  cbind(pts, vals)
})

final_sf <- do.call(rbind, result_list)
st_write(final_sf,"D:/antarctica/data/processed_data/training_data/test_traindat_VWC_mean.gpkg", delete_dsn = TRUE)

#read test data and calculate test metrics
testdata_sf <- st_read("D:/antarctica/data/processed_data/training_data/test_traindat_VWC_mean.gpkg")
testdata <- st_drop_geometry(testdata_sf)
testdata <- testdata[!is.na(testdata$VWC_daily), ]
testdata <- testdata[complete.cases(testdata[, predictors]), ]
test_prediction <- predict(model_daily, newdata = testdata)
test_prediction_df <- data.frame(
  pred = as.vector(test_prediction)
)

testdata$predicted <- as.vector(test_prediction)
names(testdata)

#MAE, RMSE, R2 
MAE_test <- MAE(testdata$predicted, testdata$VWC_daily)
RMSE_test <- RMSE(testdata$predicted, testdata$VWC_daily)
R2_test <- R2_Score(testdata$predicted, testdata$VWC_daily)

test_metrics <- data.frame(
  RMSE_test = RMSE_test,
  Rsquared_test = R2_test,
  MAE_test = MAE_test
)

test_metrics

#AOA 
all_test_path  <- list.files("D:/antarctica/data/processed_data/predictor_stacks_test_adj_AOI_tif/", full.names = T)
all_test_stack_AOI <- map(all_test_path, rast)

dates <- unique(testdata$date)
dates <- as.numeric(gsub("-", "", dates))

names(all_test_stack_AOI) <- pred_stack_test_names# from script m01.
all_test_stack_AOI[[1]]
names(all_test_stack_AOI)

tdi_daily = trainDI(model_daily, verbose = FALSE)

aoa_list <- map(all_test_stack_AOI, ~ aoa(.x, model_daily, trainDI = tdi_daily, LPD = F))
names(aoa_list)

#save aoa's
outdir <- "D:/antarctica/data/processed_data/test_period_AOA_tif/"
#dir.create(outdir, showWarnings = T)

for (nm in names(aoa_list)) {
  
  r <- aoa_list[[nm]]
  
  outpath <- file.path(outdir, paste0("aoa_", nm, ".tif"))
  
  writeRaster(r$AOA, outpath, overwrite = TRUE)
}

#predictions
prediction_test_list <- map(all_test_stack_AOI, ~predict(.x, model_daily, na.rm = T))
#save predictions
outdir <- "D:/antarctica/data/processed_data/test_period_predictions_AOI_tif/"
#dir.create(outdir, showWarnings = T)

for (nm in names(prediction_test_list)) {
  
  r <- prediction_test_list[[nm]]
  
  outpath <- file.path(outdir, paste0("pred_in_AOI_", nm, ".tif"))
  
  writeRaster(r, outpath, overwrite = TRUE)
}


#restrict to AOA
pred_test_aoa_list <- map2(prediction_test_list, aoa_list, ~ {
  x <- .x
  x[.y$AOA == 0] <- NA
  x
})

names(pred_test_aoa_list)
#save predictions in AOI and AOA
outdir <- "D:/antarctica/data/processed_data/test_period_predictions_AOI_AOA_tif/"
dir.create(outdir, showWarnings = T)

for (nm in names(pred_test_aoa_list)) {
  
  r <- pred_test_aoa_list[[nm]]
  
  outpath <- file.path(outdir, paste0("pred_in_AOI_AOA", nm, ".tif"))
  
  writeRaster(r, outpath, overwrite = TRUE)
}



######try ffs with developer version of CAST (try what happens if more predictors were added)
#devtools::install_github("HannaMeyer/CAST")
#library(CAST)
#?CAST::ffs

#model_ffs_eS_daily <- ffs(traindat[,selected_vars],
#                       traindat$VWC_daily,
#                       importance = TRUE,
#                       metric = "RMSE",
#                       tuneLength = 2,
#                       globalval = T,
 #                      ntree=50,
#                       method = "rf",
#                       earlyStopping = F,
#                       trControl = ctrl,
#                       na.action = "na.omit")
#model_ffs_eS_daily
#global_validation(model_ffs_eS_daily)
#plot(varImp(model_ffs_eS_daily, scale = F))
# makes no difference for the end result! Same predictors were selected


