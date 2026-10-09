library(sf)
library(terra)
library(MLmetrics)
library(lubridate)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(rstatix)
library(patchwork)
library(DHARMa)
library(emmeans)
library(effects)
library(car)
library(writexl)
library(glmmTMB)

# load and prepare data

data_p2 <- readRDS("data/processed_data/TMS_VWC/all_data_p2_long_filtered.RDS")
data_p4 <- readRDS("data/processed_data/TMS_VWC/all_data_p4_long_filtered.RDS")

setdiff(names(data_p2), names(data_p4))
data_p2$season_day <- NULL
data_p2$season_year <- NULL

VWC_full <- rbind(data_p2,data_p4)
VWC_full$date <- as.Date(VWC_full$datetime)
VWC_full <- VWC_full[VWC_full$locality_id != "94204670", ] #considered as wrong because exposed (see research project for more information)
#filter to period where I have data (and not only NA for the moisture)
VWC_full <- VWC_full %>%
  filter(month(datetime) >= 11 | month(datetime) <= 2)

VWC_full_daily_mean <- VWC_full %>%
  group_by(site_id, date) %>%
  summarise(
    VWC_daily = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

daily_df <- VWC_full_daily_mean %>%
  mutate(
    year = year(date),
    month = month(date),
    season_year =
      ifelse(month >= 11,
             year + 1,
             year)
  ) %>%
  filter(
    season_year == 2022,
    month %in% c(11,12,1,2,3)
  )

gt_season_mean <- daily_df %>%
  group_by(site_id) %>%
  summarise(
    gt_mean =
      mean(VWC_daily,
           na.rm = TRUE),
    .groups = "drop"
  )

# available dates
all_dates <- daily_df %>%
  group_by(date) %>%
  summarise(
    n_obs =
      sum(!is.na(VWC_daily)),
    .groups = "drop"
  ) %>%
  filter(n_obs > 0) %>%
  pull(date) %>%
  sort()

cat(
  "Available dates:",
  length(all_dates),
  "\n"
)

# sensor-wise error function

calc_error_sensorwise <- function(
    sampled_df,
    gt_table
){
  
  sampled_df %>%
    group_by(site_id) %>%
    summarise(
      pred_mean =
        mean(
          VWC_daily,
          na.rm = TRUE
        ),
      .groups = "drop"
    ) %>%
    left_join(
      gt_table,
      by = "site_id"
    ) %>%
    mutate(
      
      MAE =
        abs(
          pred_mean -
            gt_mean
        ),
      
      rel_error =
        100 *
        MAE /
        abs(gt_mean)
    )
}

# 3.4 temporal metrics
calc_temporal <- function(dates){
  
  dates <- sort(
    unique(
      as.Date(dates)
    )
  )
  
  n_dates <- length(dates)
  
  if(n_dates == 1){
    
    return(
      tibble(
        span_days = 0,
        mean_gap = 0,
        max_gap = 0,
        n_months =
          length(
            unique(
              month(dates)
            )
          )
      )
    )
  }
  
  diffs <-
    as.numeric(
      diff(dates)
    )
  
  tibble(
    span_days =
      as.numeric(
        max(dates) -
          min(dates)
      ),
    
    mean_gap =
      mean(diffs),
    
    max_gap =
      max(diffs),
    
    n_months =
      length(
        unique(
          month(dates)
        )
      )
  )
}

# 3.5 sampling function
generate_random <- function(
    k,
    n = 500
){
  
  replicate(
    n,
    sort(
      sample(
        all_dates,
        k
      )
    ),
    simplify = FALSE
  ) %>%
    unique()
}



# 3.6 samples

results_sensor <- list()
counter <- 1

for(k in 1:25){
  
  cat(
    "\nRunning k =",
    k,
    "\n"
  )
  
  sampling_sets <- list(
    
    random =
      generate_random(
        k,
        500
      )
  )
  
  for(
    design in
    names(
      sampling_sets
    )
  ){
    
    combos <-
      sampling_sets[[design]]
    
    if(
      length(
        combos
      ) == 0
    ){
      next
    }
    
    for(i in seq_along(combos)){
      
      selected_dates <-
        combos[[i]]
      
      sampled <-
        daily_df %>%
        filter(
          date %in%
            selected_dates
        )
      
      perf <-
        calc_error_sensorwise(
          sampled,
          gt_season_mean
        )
      
      temporal <-
        calc_temporal(
          selected_dates
        )
      
      perf <- perf %>%
        mutate(
          design =
            design,
          k = k,
          rep = i,
          
          span_days =
            temporal$span_days,
          
          mean_gap =
            temporal$mean_gap,
          
          max_gap =
            temporal$max_gap,
          
          n_months =
            temporal$n_months
        )
      
      results_sensor[[counter]] <- perf
      
      counter <-
        counter + 1
    }
  }
}

analysis_random <-
  bind_rows(
    results_sensor
  )

analysis_random$k_factor <-
  factor(
    analysis_random$k
  )

analysis_random$n_month_factor <-
  factor(
    analysis_random$n_months
  )

# assign classes for span days to overcome missing values (non occuring spans)
analysis_random$span_class <- cut(
  analysis_random$span_days,
  breaks = c(0,10,20,30,40,50,60,Inf)
)
#model
mod_random_rel_error <-
  glmmTMB(
    rel_error ~
      span_class +
      #mean_gap +
      #max_gap +
      #n_month_factor +
      k_factor +
      (1 | site_id),
    
    data =
      analysis_random,
    
    family =
      Gamma(link = "log")
  )

saveRDS(mod_random_rel_error,"data/processed_data/model_timesteps_mean.RDS")
mod_random_rel_error <- readRDS("data/processed_data/model_timesteps_mean.RDS")
simresrdrel <- simulateResiduals(mod_random_rel_error)

#plot(simresme)
testUniformity(simresrdrel)
testDispersion(simresrdrel)

car::Anova(
  mod_random_rel_error,
  type = "II"
)

emm_rel_rd <- emmeans(
  mod_random_rel_error,
  ~ k_factor,
  type = "response"
)

summary(emm_rel_rd)



pairs(
  emm_rel_rd,
  adjust = "tukey",
  type = "response"
)

consec <- contrast(
  emm_rel_rd,
  method = "consec",
  adjust = "sidak"
)

consec

plot_df_rel <- as.data.frame(emm_rel_rd) %>%
  as_tibble() %>%
  mutate(
    k = as.numeric(as.character(k_factor))
  ) %>%
  filter(!is.na(response))



#publication plot
windowsFonts(
  Calibri = windowsFont("Calibri")
)

p_k_rel_err <- ggplot(plot_df_rel, aes(x = k, y = response)) +
  
  
  geom_vline(
    xintercept = 8.2,
    linetype = "22",
    linewidth = 0.5,
    colour = "grey40"
  ) +
 
  geom_line(
    linewidth = 0.8,
    colour = "black"
  ) +
  
  geom_point(
    size = 2.2,
    colour = "black"
  ) +
  
  scale_x_continuous(
    breaks = seq(2, 24, by = 2)
  ) +
  
  scale_y_continuous(
    limits = c(0, 4),
    breaks = seq(1, 4, by = 1),
    labels = c("1.0", "2.0", "3.0", "4.0"),
    expand = c(0, 0)
  ) +
  
  labs(
    x = "Number of included observations",
    y = "Relative error [%]"
  ) +
  
  theme_classic(base_family = "Calibri") +
  
  theme(
    
    text = element_text(family = "Calibri"),
    
    axis.title = element_text(
      size = 12,
      colour = "black"
    ),
    
    axis.text = element_text(
      size = 11,
      colour = "black"
    ),
    
    axis.line = element_line(
      linewidth = 0.5,
      colour = "black"
    ),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    )
    
  )

p_k_rel_err

ggsave(
  filename = "plots/relative_error_included_timesteps.png",
  plot = p_k_rel_err,
  width = 160,
  height = 100,
  units = "mm",
  dpi = 600,
  bg = "white"
)

### for span_days
emm_rel_span_rd <- emmeans(
  mod_random_rel_error,
  ~ span_class,
  type = "response"
)

summary(emm_rel_span_rd)



pairs(
  emm_rel_span_rd,
  adjust = "tukey",
  type = "response"
)

consec_span <- contrast(
  emm_rel_span_rd,
  method = "consec",
  adjust = "sidak"
)

consec_span

plot_df_rel_span <- as.data.frame(emm_rel_span_rd) %>%
  as_tibble() %>%
  mutate(
    k = as.numeric((span_class))
  ) %>%
  filter(!is.na(response))


p_k_rel_err_span <- ggplot(
  plot_df_rel_span,
  aes(x = span_class, y = response)
) +
  
  geom_vline(
    xintercept = 3.5,
    linetype = "22",
    linewidth = 0.5,
    colour = "grey40"
  ) +
  
  geom_errorbar(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    width = 0.15,
    linewidth = 0.5,
    colour = "black"
  ) +
  
  geom_point(
    size = 2.8,
    colour = "black"
  ) +
  
  scale_x_discrete(
    labels = c(
      "1–10",
      "11–20",
      "21–30",
      "31–40",
      "41–50",
      "51–60",
      ">60"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 4),
    breaks = seq(1, 4, by = 1),
    labels = c("1.0", "2.0", "3.0", "4.0"),
    expand = c(0, 0)
  ) +
  labs(
    x = "Temporal span [days]",
    y = "Relative error [%]"
  ) +
  
  theme_classic(base_family = "Calibri") +
  
  theme(
    
    text = element_text(family = "Calibri"),
    
    axis.title = element_text(
      size = 12,
      colour = "black"
    ),
    
    axis.text = element_text(
      size = 11,
      colour = "black"
    ),
    
    axis.line = element_line(
      linewidth = 0.5,
      colour = "black"
    ),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    )
  )

p_k_rel_err_span

ggsave(
  filename = "plots/relative_error_spandays.png",
  plot = p_k_rel_err_span,
  width = 160,
  height = 100,
  units = "mm",
  dpi = 600,
  bg = "white"
)

#evaluation of sd

rel_error_sd_site_k <- analysis_random %>%
  group_by(site_id, k) %>%
  summarise(
    sd_rel_error = sd(rel_error, na.rm = TRUE),
    mean_rel_error = mean(rel_error, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
rel_error_sd_site_k



ggplot(
  rel_error_sd_site_k,
  aes(
    x = k,
    y = sd_rel_error
  )
) +
  geom_line() +
  geom_point(size = 1.5) +
  facet_wrap(~site_id) +
  labs(
    x = "Number of included observations",
    y = "SD of relative error [%]"
  ) +
  theme_classic()


rel_error_sd_k <- rel_error_sd_site_k %>%
  group_by(k) %>%
  summarise(
    mean_sd = mean(sd_rel_error, na.rm = TRUE),
    sd_sd = sd(sd_rel_error, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    lower = mean_sd - sd_sd,
    upper = mean_sd + sd_sd
  )

p_sd_rel_err <- ggplot(
  rel_error_sd_k,
  aes(
    x = k,
    y = mean_sd
  )
) +
  
  geom_vline(
    xintercept = 8.2,
    linetype = "22",
    linewidth = 0.5,
    colour = "grey40"
  ) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.15
  ) +
  geom_line(
    linewidth = 0.8,
    colour = "black"
  ) +
  
  geom_point(
    size = 2.2,
    colour = "black"
  ) + 
  
  scale_x_continuous(
    breaks = seq(2, 24, by = 2)
  ) +
  
  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.1),
    expand = c(0, 0)
  ) +
  
  labs(
    x = "Number of included observations",
    y = "SD of relative error [%]"
  ) +
  
  theme_classic(
    base_family = "Calibri"
  ) +
  
  theme(
    
    text = element_text(
      family = "Calibri"
    ),
    
    axis.title = element_text(
      size = 12,
      colour = "black"
    ),
    
    axis.text = element_text(
      size = 11,
      colour = "black"
    ),
    
    axis.line = element_line(
      linewidth = 0.5,
      colour = "black"
    ),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.5
    )
  )

p_sd_rel_err

# combined 2
combined_plot <-
  p_k_rel_err +
  p_k_rel_err_span +
  p_sd_rel_err+
  plot_layout(
    ncol = 3,
    widths = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      text = element_text(family = "Calibri"),
      plot.tag = element_text(
        family = "Calibri",
        face = "bold",
        size = 13
      )
    )
  ) &
  theme(
    text = element_text(family = "Calibri")
  )
combined_plot

ggsave(
  filename = "plots/relative_error_combined_plot_incl_sd.png",
  plot = combined_plot,
  width = 280,
  height = 80,
  units = "mm",
  dpi = 600,
  bg = "white"
)
