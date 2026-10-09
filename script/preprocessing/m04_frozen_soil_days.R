library(lubridate)

met_data <- read.csv("D:/antarctica/data/raw_data/met_data_LTER/mcmlter-clim-frlm_daily-2026-02-28.csv", header = T)
met_data$date <- as.Date(met_data$date_time)

above0deg_days <- unique(
  met_data$date[met_data$soiltemp1_5cm_degc > 0]
)

above0deg_days

met_20_21 <- met_data[
  met_data$date >= as.Date("2020-10-01") &
    met_data$date <= as.Date("2021-03-31"),
]

met_20_21$soiltemp1_10cm_degc
tage_pos <- unique(
  met_20_21$date[
    met_20_21$soiltemp1_10cm_degc > 0 
  ]
)
tage_pos

tage_pos_5 <- unique(
  met_20_21$date[
    met_20_21$soiltemp1_5cm_degc > 0 
  ]
)
tage_pos_5 #thats the correct approach, because thats similar to the temperature measurement position of the TMS 

met_19_20 <- met_data[
  met_data$date >= as.Date("2019-10-01") &
    met_data$date <= as.Date("2020-03-31"),
]


tage_pos_5 <- unique(
  met_19_20$date[
    met_19_20$soiltemp1_5cm_degc > 0 
  ]
)
tage_pos_5 

met_23_24 <- met_data[
  met_data$date >= as.Date("2023-10-01") &
    met_data$date <= as.Date("2024-03-31"),
]


tage_pos_5 <- unique(
  met_23_24$date[
    met_23_24$soiltemp1_5cm_degc > 0 
  ]
)
tage_pos_5 

met_18_19 <- met_data[
  met_data$date >= as.Date("2018-10-01") &
    met_data$date <= as.Date("2019-03-31"),
]


tage_pos_5 <- unique(
  met_18_19$date[
    met_18_19$soiltemp1_5cm_degc > 0 
  ]
)

tage_pos_5 

met_17_18 <- met_data[
  met_data$date >= as.Date("2017-10-01") &
    met_data$date <= as.Date("2018-03-31"),
]


tage_pos_5 <- unique(
  met_17_18$date[
    met_17_18$soiltemp1_5cm_degc > 0 
  ]
)

tage_pos_5 

met_16_17 <- met_data[
  met_data$date >= as.Date("2016-10-01") &
    met_data$date <= as.Date("2017-03-31"),
]


tage_pos_5 <- unique(
  met_16_17$date[
    met_16_17$soiltemp1_5cm_degc > 0 
  ]
)

tage_pos_5 

met_15_16 <- met_data[
  met_data$date >= as.Date("2015-10-01") &
    met_data$date <= as.Date("2016-03-31"),
]

tage_pos_5 <- unique(
  met_15_16$date[
    met_15_16$soiltemp1_5cm_degc > 0 
  ]
)

tage_pos_5 

met_21_22 <- met_data[
  met_data$date >= as.Date("2021-10-01") &
    met_data$date <= as.Date("2022-03-31"),
]

tage_pos_5 <- unique(
  met_21_22$date[
    met_21_22$soiltemp1_5cm_degc > 0 
  ]
)

tage_pos_5 

met_22_23 <- met_data[
  met_data$date >= as.Date("2022-10-01") &
    met_data$date <= as.Date("2023-03-31"),
]

tage_pos_5 <- unique(
  met_22_23$date[
    met_22_23$soiltemp1_5cm_degc > 0 
  ]
)

tage_pos_5 

met_25_26 <- met_data[
  met_data$date >= as.Date("2025-10-01") &
    met_data$date <= as.Date("2026-03-31"),
]

tage_pos_5 <- unique(
  met_25_26$date[
    met_25_26$soiltemp1_5cm_degc > 0 
  ]
)

tage_pos_5 

met_24_25 <- met_data[
  met_data$date >= as.Date("2024-10-01") &
    met_data$date <= as.Date("2025-03-31"),
]

tage_pos_5 <- unique(
  met_24_25$date[
    met_24_25$soiltemp1_5cm_degc > 0 
  ]
)

tage_pos_5 

met_23_24 <- met_data[
  met_data$date >= as.Date("2023-10-01") &
    met_data$date <= as.Date("2024-03-31"),
]

tage_pos_5 <- unique(
  met_23_24$date[
    met_23_24$soiltemp1_5cm_degc > 0 
  ]
)

tage_pos_5 
#Foehn
foehn_marte <- read.csv("data/raw_data/data_and_info_Eva_Bendix/foehn_dataset_Marte/fryxell_foehn_results.csv")
head(foehn_marte)
s_w_foehn <- foehn_marte[foehn_marte$foehn_Speirs == 1 & foehn_marte$foehn_Wiesenekker == 1,]
s_w_foehn$date <- date(s_w_foehn$date_time)
unique(s_w_foehn$date)

s_w_foehn <- s_w_foehn[
  s_w_foehn$date >= as.Date("2015-12-18") 
]

#check whether the days for we have available multispectral data are above 0 deg

multispec_info <- read.csv("D:/antarctica/data/raw_data/multispec_data_availability.csv", header = T, sep = ";")

all(as.Date(multispec_info$date) %in% as.Date(above0deg_days))
setdiff(as.character(multispec_info$date), as.character(above0deg_days))

multispec_info$date <- as.Date(multispec_info$date)
str(multispec_info$date)

library(dplyr)

multispec_info <- multispec_info %>%
  left_join(
    met_data %>% select(date, soiltemp1_5cm_degc),
    by = "date"
  )

multispec_info$soiltemp1_5cm_degc

multispec_info[multispec_info$soiltemp1_5cm_degc < 0,]$date

multispec_sub <- multispec_info[multispec_info$soiltemp1_5cm_degc >= 0 | is.na(multispec_info$soiltemp1_5cm_degc),] 

multispec_sub$doy <- yday(multispec_sub$date) 
unique(multispec_sub$doy)

#identify the doys that have not been calculated yet
PS_doy_name <- list.files("C:/Users/mini4/Desktop/Uni/Master/antarctica/data/processed_data/rsun_daily/", pattern = "^glob_rad.*\\.tif$", full.names = F)
PS_doy <- as.numeric(gsub(".*_day_([0-9]+)\\.tif$", "\\1", PS_doy_name))

doy_rsun <- setdiff(as.numeric(multispec_sub$doy), PS_doy)

setdiff( PS_doy,as.numeric(multispec_sub$doy))
