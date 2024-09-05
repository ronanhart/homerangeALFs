# X===========================================X
# X--------------- HR Analysis ---------------X
# X------------ Model Validation -------------X
# X===========================================X

library(tidyverse)
library(lme4)
library(MuMIn)
library(lmerTest)
library(wCorr)

# ---------------X
# ---- SET UP ----
# ---------------X
data_dir <- "prepped_data/"
result_dir <- "results/"

# ... Load & Prep the data ----
dat_used_avail_wide <- read.csv(paste0(in_dir, "used_avail_wide.csv"))
hr_info <- read.csv(paste0(in_dir, "hr_info_full_si.csv"))
subset_info <- read.csv(paste0(in_dir, "subset_info.csv"))

# Join into one data set and select only the columns that are necessary for modeling
dat <- left_join(dat_used_avail_wide, hr_info, by = "hr_id") %>%
  left_join(subset_info, by = c("hr_id", "species", "sex", "season"))
head(dat)
(cols <- colnames(dat))
(avail_logsr_cols <- cols[grepl("_sc", cols) & !grepl("used", cols)])

dat <- dat %>%
  select(subset, hr_id, animal_id, species, is_M, season, n_days, n.pts_X_n.days_norm, 
         log_area, log_log_shape, all_of(avail_logsr_cols))
head(dat)
colnames(dat)

# Subset the data based on testing flags (from 02a_data_subsetting)
dat_test <- filter(dat, subset == "Test")
summary(dat_test)


# ...Load the models ----
all_mods <- readRDS(paste0(out_dir, "all_models.rds"))
View(all_mods)

# -------------------------X
# ---- MODEL VALIDATION ----
# -------------------------X
# Function to pull the observed & predicted values and weights 
# for each species-season-model
obs_pred_func <- function(model_name, data, spp, ssn){
  if(!spp %in% c("Mule Deer", "Pronghorn") | 
     !ssn %in% c("Winter", "Summer")){
    stop("`spp` must be either `Mule Deer` or `Pronghorn`
       `ssn` must be either `Winter` or `Summer`")
  }
  # skip if the season is summer and the model name is logSR_snd
  # (there is no snow in the summer so we did not fit a logSR snow depth model)
  if(!(ssn == "Summer" & model_name == "logSR_snd")){
    
    # Filter the data to the given species and season
    data <- filter(data, species == spp & season == ssn)
    
    # Pull the observed column name
    cols <- colnames(data)
    resp_col <- str_replace(model_name, "logSR", "logSR_sc")
    if(!resp_col %in% cols){stop(model_name, ": check the response name")}
    
    # pull the given species-season-model
    model <- all_mods[[spp]][[ssn]][[model_name]]
    
    # pull the predicted and observed values (and the weights)
    pred <- predict(model, data, re.form = NA)
    obs <- data[, resp_col]
    wgt <- data$n.pts_X_n.days_norm
    
    # return a data frame of the species-season-model
    # and the model's observed & predicted values and weights
    df <- data.frame(species = spp, season = ssn, response = model_name, 
                     obs = obs, pred = pred, wgt = wgt)
    row.names(df) <- NULL
    return(df)
  }
}

# Function to 
#  1. calculate the Pearson's correlation coefficient for each model
#  2. Fit a linear model to the observed and predicted values
#      a. pull the beta coefficients, std. error, and p-value
#      b. pull the r-squared and adjusted r-squared values
valid_func <- function(x){
  corr <- weightedCorr(x = x$obs, y = x$pred, weights = x$wgt, method = "Pearson")
  lm <- lm(obs ~ pred, data = x)
  rsq <- summary(lm)$r.squared
  adj_rsq <- summary(lm)$adj.r.squared
  coef <- summary(lm)$coefficients %>%
    as.data.frame() %>%
    mutate(corr = corr, rsq = rsq, adj_rsq = adj_rsq)
  coef$param <- row.names(coef)
  coef <- select(coef, param, estimate = Estimate, se = `Std. Error`, 
                 p_val = `Pr(>|t|)`, corr, rsq, adj_rsq) %>%
    mutate(param = ifelse(param == "(Intercept)", "intercept", param))
  row.names(coef) <- NULL
  return(coef)
}

# Loop through each species-season-model combination 
# to get a data frame of all models' observed & predicted values
# For TESTING data
obs_pred <- lapply(c("Mule Deer", "Pronghorn"), function(spp){
  lapply(c("Winter", "Summer"), function(ssn){
    lapply(names(response_vars), obs_pred_func, 
           data = dat_test, spp = spp, ssn = ssn) %>%
      bind_rows()
  }) %>% bind_rows()
}) %>% bind_rows() %>%
  # turn the species and season columns into factors
  mutate(species = factor(species),
         season = factor(season, levels = c("Winter", "Summer")),
         # add a column for human-readable model names
         response_name = case_when(
           response == "log_area" ~ "log(Area)",
           response == "log_log_shape" ~ "log(log(Shape))",
           response == "logSR_rd.pvd" ~ "Paved Roads",
           response == "logSR_rd.unpvd" ~ "Unpaved Roads",
           response == "logSR_fence" ~ "Fences",
           response == "logSR_elev" ~ "Elevation",
           response == "logSR_rough" ~ "Roughness",
           response == "logSR_forage" ~ "Forage",
           response == "logSR_cvr.shrub" ~ "Shrub Cover",
           response == "logSR_cvr.tree" ~ "Tree Cover",
           response == "logSR_snd" ~ "Snow Depth"),
         response_name = factor(
           response_name, levels = c("log(Area)", "log(log(Shape))", "Paved Roads",
                                     "Unpaved Roads", "Fences", "Elevation", 
                                     "Roughness", "Forage", "Shrub Cover", 
                                     "Tree Cover", "Snow Depth"))
  ) %>%
  relocate(response_name, .after = response)
head(obs_pred)

# Put the weights into bins
(wgt_bins <- BAMMtools::getJenksBreaks(var = obs_pred$wgt, k = 4) %>% 
    signif(1)) # round to a number with 2 significant figures
obs_pred <- obs_pred %>%
  mutate(wgt_bin = case_when(
    wgt >= min(wgt) & 
      wgt < wgt_bins[2] ~ "Low",
    wgt >= wgt_bins[2] & 
      wgt < wgt_bins[3] ~ "Medium",
    wgt >= wgt_bins[3] & 
      wgt <= max(wgt) ~ "High"),
    wgt_bin = factor(wgt_bin, levels = c("Low", "Medium", "High")))

# nest by species-season-model and calculate model fit
model_fit <- obs_pred %>%
  nest(data = -c(species, season, response, response_name)) %>%
  mutate(model_fit = map(data, valid_func)) %>%
  select(-data) %>%
  unnest(cols = model_fit)
head(model_fit)
summary(model_fit)

# ... PLOT ----
model_fit %>%
  select(species, season, response_name, param, estimate) %>%
  pivot_wider(names_from = param, values_from = estimate) %>%
  distinct() %>% as.data.frame()
lapply(unique(model_fit$response_name), function(resp){
  model_fit_resp <- model_fit %>%
    filter(response_name == resp) %>%
    pivot_wider(names_from = param, 
                values_from = c(estimate, se, p_val)) 
  obs_pred_resp <- filter(obs_pred, response_name == resp)
  min <- min(c(obs_pred_resp$obs, obs_pred_resp$pred))
  max <- max(c(obs_pred_resp$obs, obs_pred_resp$pred))
  
  ggplot(model_fit_resp) +
    facet_grid(species~season) +
    geom_abline(aes(slope = estimate_pred, intercept = estimate_intercept,
                    col = season), linewidth = 1) +
    scale_alpha_manual(values = c(0.15, 0.4, 0.75), name = "Weights") +
    geom_abline(slope = 1, linetype = 2) +
    geom_text(aes(x = min, y = max, 
                  label = paste("Intercept:", round(estimate_intercept, 3))),
              hjust = 0) +
    geom_text(aes(x = min, y = max, 
                  label = paste("Slope:", round(estimate_pred, 3))),
              hjust = 0, vjust = 2) +
    geom_text(aes(x = min, y = max,
                  label = paste("Pearson's Correlation:", round(corr, 3))),
              hjust = 0, vjust = 3.5) +
    geom_text(aes(x = min, y = max,
                  label = paste("R-squared:", round(rsq, 3))),
              hjust = 0, vjust = 5) +
    theme_bw() +
    labs(title = resp, y = "Observed", x = "Predicted") +
    theme(legend.position = "none") +
    xlim(c(min, max)) + ylim(c(min, max))
})