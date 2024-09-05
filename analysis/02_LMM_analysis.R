# X===========================================X
# X--------------- HR Analysis ---------------X
# X----------- Linear Mixed Models -----------X
# X===========================================X

# Response variables:
#   - log(Area)
#   - log(log(Shape Index))
#   - logSR ALFs (paved roads, unpaved roads, fences)
#   - logSR Environment (elevation, roughness, forage, shrub cover,
#                       tree cover, snow depth)

# Model Structure:
#   - log(Area) ~ 
#       available LFs + available environmental covariates + sex + 
#       sex:road_paved + sex:road_unpaved + sex:fence + n_days + (1|animal_id)
#   - log(log(shape)) ~
#       available LFs + available rough + available snow depth + sex + 
#       sex:road_paved + sex:road_unpaved + sex:fence + n_days + (1|animal_id)
#   - logSR_j ~ 
#       available_j (ie only include the availability of the focal
#       response variable) + logSR other predictors + sex + 
#       sex:road_paved + sex:road_unpaved + sex:fence + (1|animal_id)

# Weights: N points x N days (normalized) 

# Separate models for species and season

#   - mule deer x winter 
#   - mule deer x summer 
#   - pronghorn x winter 
#   - pronghorn x summer 

library(tidyverse)
library(lme4)
library(MuMIn)
library(lmerTest)
library(wCorr)


# ---------------X
# ---- SET UP ----
# ---------------X
in_dir <- "prepped_data/"

# Load the data 
dat_used_avail_wide <- read.csv(paste0(in_dir, "used_avail_wide.csv"))
hr_info <- read.csv(paste0(in_dir, "hr_info_full_si.csv"))
subset_info <- read.csv(paste0(in_dir, "subset_info.csv"))

# Join into one dataset and select only the columns that are necessary for modelling
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

# Subset the data based on training/testing flags (from 02a_data_subsetting)
dat_train <- filter(dat, subset == "Train")
summary(dat_train)

dat_test <- filter(dat, subset == "Test")
summary(dat_test)


# ----------------------X
# ---- MODEL SET UP  ----
# ----------------------X

# Set up model structure
# These will go in the formula() function
# availability and logSR predictors
avail <- "avail_sc_"
logSR <- "logSR_sc_"
env <- c("elev", "rough", "forage", "cvr.shrub", "cvr.tree", "snd", "pdsi")
lf <- c("rd.pvd", "rd.unpvd", "fence")

# function to remove the given attribute_j from the list of logSR
rm_j <- function(j, env){
  x <- paste0(logSR, c(lf, env))
  return(x[!grepl(j, x)])
}

# interact sex indicator variables with LF variables ('int' means 'interaction')
ind <- "is_M"
(int.lf <- paste(ind, paste0(avail, lf), sep = ":"))

# random effect
rand <- "(1|animal_id)"

# response variables
response_vars <- paste(c("log_area", "log_log_shape", 
                         paste0(logSR, lf), paste0(logSR, env)), "~ ")
names(response_vars) <- c("log_area", "log_log_shape", 
                          paste0("logSR_", lf), paste0("logSR_", env))

# remove pdsi from a response variable (it's only in there as a predictor)
response_vars <- response_vars[!grepl("pdsi", response_vars)]
response_vars


# Function that will run the given model with the given the species-season combo
model_func <- function(model_name, data, spp, ssn){
  if(!spp %in% c("Mule Deer", "Pronghorn") | 
     !ssn %in% c("Winter", "Summer")){
    stop("`spp` must be either `Mule Deer` or `Pronghorn`
       `ssn` must be either `Winter` or `Summer`")
  }
  
  # Filter the data to the given species and season
  data <- filter(data, species == spp & season == ssn)
  
  # If it's summer, remove snow depth from a list of environmental predictors
  if(ssn == "Summer") env <- env[!env %in% "snd"] 
  
  # ... a. log(Area) ----
  if(model_name == "log_area"){
    formula.area <- formula(
      paste(response_vars["log_area"],
            paste(c(paste0(avail, c(lf, env)), ind, int.lf, "n_days", rand), 
                  collapse = " + "),
            collapse = " "))
    
    tryCatch({
      mod.area <- lmer(formula.area,
                       weights = n.pts_X_n.days_norm, REML = T, na.action = na.omit,
                       data = data)
      return(mod.area)
    }, warning = function(w){
      cat("log(area)", spp, ssn, ":", conditionMessage(w), "\n")
    })
  }
  
  # ... b. log(log(Shape)) ----
  if(model_name == "log_log_shape"){
    # if it's summer, remove snow depth from the list of predictors
    shape_avail <- if(ssn == "Winter"){
      c("rough", "snd")
    } else {
      "rough"
    }
    formula.shape <- formula(
      paste(response_vars["log_log_shape"],
            paste(c(paste0(avail, c(lf, shape_avail)), ind, int.lf,
                    "n_days", rand), # include N days
                  collapse = " + "),
            collapse = " "))
    
    tryCatch({
      mod.shape <- lmer(formula.shape,
                        weights = n.pts_X_n.days_norm, REML = T, na.action = na.omit,
                        data = data)
      return(mod.shape)
    }, warning = function(w){
      cat("log(log(shape index))", spp, ssn, ":", conditionMessage(w), "\n")
    })
  }
  
  # ... c. logSR Paved Road ----
  if(model_name == "logSR_rd.pvd"){
    formula.road_pvd <- formula(
      paste(response_vars["logSR_rd.pvd"],
            paste(c(paste0(avail, "rd.pvd"), rm_j("rd.pvd", env), ind,
                    # only have interactions with the focal ALF
                    int.lf[grepl("rd.pvd", int.lf)], rand),
                  collapse = " + "),
            collapse = " "))
    
    tryCatch({
      mod.road_pvd <- lmer(formula.road_pvd,
                           weights = n.pts_X_n.days_norm, REML = T, na.action = na.omit,
                           data = data)
      return(mod.road_pvd)
    }, warning = function(w){
      cat("logSR paved roads", spp, ssn, ":", conditionMessage(w), "\n")
    })
  }
  
  # ... d. logSR Unpaved Road ----
  if(model_name == "logSR_rd.unpvd"){
    formula.road_unpvd <- formula(
      paste(response_vars["logSR_rd.unpvd"],
            paste(c(paste0(avail, "rd.unpvd"), rm_j("rd.unpvd", env), ind,
                    int.lf[grepl("rd.unpvd", int.lf)], rand),
                  collapse = " + "),
            collapse = " "))
    
    tryCatch({
      mod.road_unpvd <- lmer(formula.road_unpvd,
                             weights = n.pts_X_n.days_norm, REML = T, na.action = na.omit,
                             data = data)
      return(mod.road_unpvd)
    }, warning = function(w){
      cat("logSR unpaved roads", spp, ssn, ":", conditionMessage(w), "\n")
    })
  }
  
  # ... e. logSR Fence ----
  if(model_name == "logSR_fence"){
    formula.fence <- formula(
      paste(response_vars["logSR_fence"], 
            paste(c(paste0(avail, "fence"), rm_j("fence", env), ind, 
                    int.lf[grepl("fence", int.lf)], rand), 
                  collapse = " + "), 
            collapse = " "))
    
    tryCatch({
      mod.fence <- lmer(formula.fence,
                        weights = n.pts_X_n.days_norm, REML = T, na.action = na.omit,
                        data = data)
      return(mod.fence)
    }, warning = function(w){
      cat("logSR fence", spp, ssn, ":", conditionMessage(w), "\n")
    })
  }
  
  # ... f. logSR Elevation ----
  if(model_name == "logSR_elev"){
    formula.elev <- formula(
      paste(response_vars["logSR_elev"],
            paste(c(paste0(avail, c("elev", lf)), rm_j("elev", env), ind,
                    int.lf, rand),
                  collapse = " + "),
            collapse = " "))
    
    tryCatch({
      mod.elev <- lmer(formula.elev,
                       weights = n.pts_X_n.days_norm, REML = T, na.action = na.omit,
                       data = data)
      return(mod.elev)
    }, warning = function(w){
      cat("elevation", spp, ssn, ":", conditionMessage(w), "\n")
    })
  }
  
  # ... g. logSR Roughness ----
  if(model_name == "logSR_rough"){
    formula.rough <- formula(
      paste(response_vars["logSR_rough"],
            paste(c(paste0(avail, c("rough", lf)), rm_j("rough", env), ind,
                    int.lf, rand),
                  collapse = " + "),
            collapse = " "))
    
    tryCatch({
      mod.rough <- lmer(formula.rough,
                        weights = n.pts_X_n.days_norm, REML = T, na.action = na.omit,
                        data = data)
      return(mod.rough)
    }, warning = function(w){
      cat("roughness", spp, ssn, ":", conditionMessage(w), "\n")
    })
  }
  
  # ... h. logSR Forage ----
  if(model_name == "logSR_forage"){
    formula.forage <- formula(
      paste(response_vars["logSR_forage"],
            paste(c(paste0(avail, "forage"), paste0(avail, lf),
                    rm_j("forage", env), ind, int.lf, rand),
                  collapse = " + "),
            collapse = " "))
    
    tryCatch({
      mod.forage <- lmer(formula.forage,
                         weights = n.pts_X_n.days_norm, REML = T, na.action = na.omit,
                         data = data)
      return(mod.forage)
    }, warning = function(w){
      cat("forage", spp, ssn, ":", conditionMessage(w), "\n")
    })
  }
  
  # ... i. logSR Shrub Cover ----
  if(model_name == "logSR_cvr.shrub"){
    formula.shrub <- formula(
      paste(response_vars["logSR_cvr.shrub"],
            paste(c(paste0(avail, "cvr.shrub"), paste0(avail, lf),
                    rm_j("cvr.shrub", env), ind, int.lf, rand),
                  collapse = " + "),
            collapse = " "))
    
    tryCatch({
      mod.shrub <- lmer(formula.shrub,
                        weights = n.pts_X_n.days_norm, REML = T, na.action = na.omit,
                        data = data)
      return(mod.shrub)
    }, warning = function(w){
      cat("shrub cover", spp, ssn, ":", conditionMessage(w), "\n")
    })
  }
  
  # ... j. logSR Tree Cover ----
  if(model_name == "logSR_cvr.tree"){
    formula.tree <- formula(
      paste(response_vars["logSR_cvr.tree"],
            paste(c(paste0(avail, "cvr.tree"), paste0(avail, lf),
                    rm_j("cvr.tree", env), ind, int.lf, rand),
                  collapse = " + "),
            collapse = " "))
    
    tryCatch({
      mod.tree <- lmer(formula.tree,
                       weights = n.pts_X_n.days_norm, REML = T, na.action = na.omit,
                       data = data)
      return(mod.tree)
    }, warning = function(w){
      cat("tree cover", spp, ssn, ":", conditionMessage(w), "\n")
    })
  }
  
  # ... k. logSR Snow Depth ----
  if(model_name == "logSR_snd" & ssn == "Winter"){
    formula.snd <- formula(
      paste(response_vars["logSR_snd"],
            paste(c(paste0(avail, "snd"), paste0(avail, lf), 
                    rm_j("snd", env), ind, int.lf, rand),
                  collapse = " + "),
            collapse = " "))
    
    tryCatch({
      mod.snd <- lmer(formula.snd,
                      weights = n.pts_X_n.days_norm, REML = T, na.action = na.omit,
                      data = data)
      return(mod.snd)
    }, warning = function(w){
      cat("snow depth", spp, ssn, ":", conditionMessage(w), "\n")
    })
  }
}


# --------------X
# ---- MODEL ----
# --------------X

# ... Mule Deer ----
# ... ... Winter ----
md_win_mods <- lapply(names(response_vars), model_func, 
                      data = dat_train, spp = "Mule Deer", ssn = "Winter")
names(md_win_mods) <- names(response_vars)
md_win_mods
lapply(md_win_mods, summary)

# ... ... Summer ----
md_sum_mods <- lapply(names(response_vars), model_func, 
                      data = dat_train, spp = "Mule Deer", ssn = "Summer")
names(md_sum_mods) <- names(response_vars)
md_sum_mods
lapply(md_sum_mods, summary)


# ... Pronghorn ----
# ... ... Winter ----
pr_win_mods <- lapply(names(response_vars), model_func, 
                      data = dat_train, spp = "Pronghorn", ssn = "Winter")
names(pr_win_mods) <- names(response_vars)

# the tree cover model failed to converge because of the random effects
# re-run but as a linear model instead of LMM
pr_win_mods[["logSR_cvr.tree"]] 
pr_win_tree <- lm(
  formula(paste(response_vars["logSR_cvr.tree"],
                paste(c(paste0(avail, "cvr.tree"), paste0(avail, lf),
                        rm_j("cvr.tree", env), ind, int.lf),
                      collapse = " + "),
                collapse = " ")),
  weights = n.pts_X_n.days_norm, na.action = na.omit, 
  data = filter(dat_train, species == "Pronghorn" & season == "Winter")
)
pr_win_mods[["logSR_cvr.tree"]] <- pr_win_tree
pr_win_mods
lapply(pr_win_mods, summary)

# ... ... Summer ----
pr_sum_mods <- lapply(names(response_vars), model_func, 
                      data = dat_train, spp = "Pronghorn", ssn = "Summer")
names(pr_sum_mods) <- names(response_vars)
pr_sum_mods
lapply(pr_sum_mods, summary)


# ... Bind into one list of lists ----
md_mods <- list(md_win_mods, md_sum_mods)
names(md_mods) <- c("Winter", "Summer")

pr_mods <- list(pr_win_mods, pr_sum_mods)
names(pr_mods) <- c("Winter", "Summer")

all_mods <- list(md_mods, pr_mods)
names(all_mods) <- c("Mule Deer", "Pronghorn")
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
obs_pred <- lapply(c("Mule Deer", "Pronghorn"), function(spp){
  lapply(c("Winter", "Summer"), function(ssn){
    lapply(names(response_vars), obs_pred_func, 
           data = dat_train, spp = spp, ssn = ssn) %>%
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
    geom_text(aes(x = min, y = max,
                  label = paste("Adj. r-squared:", round(adj_rsq, 3))),
              hjust = 0, vjust = 6.5) +
    theme_bw() +
    labs(title = resp) +
    xlim(c(min, max)) + ylim(c(min, max))
})


