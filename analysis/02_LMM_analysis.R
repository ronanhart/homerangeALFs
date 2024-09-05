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


# ---------------X
# ---- SET UP ----
# ---------------X
in_dir <- "prepped_data/"
out_dir <- "results/"
if(!dir.exists(out_dir)) dir.create(out_dir)

# Load the data 
dat_used_avail_wide <- read.csv(paste0(in_dir, "used_avail_wide.csv"))
hr_info <- read.csv(paste0(in_dir, "hr_info_full_si.csv"))
subset_info <- read.csv(paste0(in_dir, "subset_info.csv"))

# Join into one dataset and select only the columns that are necessary for modeling
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

# Subset the data based on training flags (from 02a_data_subsetting)
dat_train <- filter(dat, subset == "Train")
summary(dat_train)


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

# Save
saveRDS(all_mods, paste0(out_dir, "all_models.rds"))






