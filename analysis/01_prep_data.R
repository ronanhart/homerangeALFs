# X================================================X
# X------------------ HR Analysis -----------------X
# X---------- Prepare Data for Analysis -----------X
# X================================================X

library(tidyverse)

# for calculating perimeter and shape index of the HRs
library(sf)
library(terra)

# ---------------X
# ---- SET UP ----
# ---------------X
in_dir <- "raw_data/"
out_dir <- "prepped_data/"
if(!dir.exists(out_dir)) dir.create(out_dir)

# Load the used/available data frame
dat_used_avail <- read.csv(paste0(in_dir, "used_avail_long.csv"))
head(dat_used_avail)
summary(dat_used_avail)
table(dat_used_avail$attribute)

# Load the HR info data frame
hr_info <- read.csv(paste0(in_dir, "hr_info.csv"))
head(hr_info)
summary(hr_info)

# Load the HR sf
hr_sf <- st_read(paste0(in_dir, "HRs"))
hr_sf


# --------------------------------------X
# ---- PREPARE USED/AVAIL DATA FRAME ----
# --------------------------------------X
#   a. Normalize used/avail values from 0-1
#   b. Scale and center available values
#   c. Widen the data


# ... a. Normalize used/avail values from 0-1 ----
#    based on their natural range
# First, keep the original used/avail values
dat_used_avail <- dat_used_avail %>%
  # "og" = original
  rename(used_og = used) %>%
  rename(avail_og = avail) %>%
  # divide used/avail (don't take log)
  mutate(used.avail_og = used_og/avail_og)

# Find the absolute min and max across used and all avail for all attributes
abs_min_max <- dat_used_avail %>%
  group_by(attribute) %>%
  summarize(used_min = min(used_og, na.rm = T),
            avail_min = min(avail_og, na.rm = T),
            used_max = max(used_og, na.rm = T),
            avail_max = max(avail_og, na.rm = T)) %>%
  group_by(attribute) %>%
  mutate(abs_min = min(used_min, avail_min),
         abs_max = max(used_max, avail_max)) %>%
  select(-c(used_min, avail_min, 
            used_max, avail_max)) %>%
  # cover (shrub/tree), LFs, forage, and PDSI already have ranges
  # (so make those the ranges instead of the min/max)
  # elevation, roughness, and snow depth don't have ranges 
  # (so the min/max are their ranges)
  mutate(
    oldmin = case_when(
      attribute == "pdsi" ~ -15,
      grepl("cvr", attribute) ~ 0,
      # snow, and roughness don't have an upper bound but they have a lower bound
      attribute == "snd" | attribute == "rough" ~ 0,
      attribute == "forage" ~ -100,
      grepl("rd", attribute) | grepl("fence", attribute) ~ 0,
      TRUE ~ floor(abs_min)),
    oldmax = case_when(
      grepl("cvr", attribute) ~ 100,
      attribute == "pdsi" ~ 15,
      attribute == "forage" ~ 100,
      grepl("road", attribute) | grepl("fence", attribute) ~ 1,
      TRUE ~ ceiling(abs_max))) %>%
  select(-c(abs_min, abs_max))
abs_min_max

# OldRange = (OldMax - OldMin)  
# NewRange = (NewMax - NewMin)  
# NewValue = (((OldValue - OldMin) * NewRange) / OldRange) + NewMin
dat_used_avail <- dat_used_avail %>%
  left_join(abs_min_max, by = "attribute") %>%
  mutate(used_new = (((used_og - oldmin)) / (oldmax - oldmin)),
         avail_new = (((avail_og - oldmin)) / (oldmax - oldmin)),
         # add 1 (to avoid errors with log)
         used_new = 1 + used_new,
         avail_new = 1 + avail_new,
         # calculate a new log(used/avail)
         used.avail_new = (used_new/avail_new),
         logSR_new = log(used.avail_new)) %>%
  select(-c(oldmin, oldmax))

# ... b. Scale and center available values ----
# Used values are not going in the models 
# The available values going in the model (and thus need to be scaled) will ONLY
#   be the original available values (not the ones with any manipulations)
#   UNLESS the original values are NA, then use the nonzero values
# Create a data.frame with means and SDs for all continuous covariates
scale_df <- dat_used_avail  %>%
  group_by(attribute) %>%
  summarize(
    mean_avail = mean(avail_og, na.rm = T),
    sd_avail = sd(avail_og, na.rm = T),
    mean_used = mean(used_og, na.rm = T),
    sd_used = sd(used_og, na.rm = T),
    mean_logSR = mean(logSR_new, na.rm = T),
    sd_logSR = sd(logSR_new, na.rm = T)
  )
scale_df


# Now a function that takes 'scale_df' and scales and centers all variables
# that appear in 'scale_df'.
scale_func <- function(x, scale_df){
  lapply(1:nrow(scale_df), function(i){
    attr <- scale_df$attribute[i]
    m_used <- scale_df$mean_used[i]
    s_used <- scale_df$sd_used[i]
    m_avail <- scale_df$mean_avail[i]
    s_avail <- scale_df$sd_avail[i]
    m_log.u.a_buffer <- scale_df$mean_logSR[i]
    s_log.u.a_buffer <- scale_df$sd_logSR[i]
    
    x_attr <- filter(x, attribute == attr)
    
    x_scale <- x_attr %>%
      mutate(used_sc = (used_og - m_used)/s_used,
             avail_sc = (avail_og - m_avail)/s_avail,
             logSR_sc = 
               (logSR_new - m_log.u.a_buffer)/s_log.u.a_buffer)
    return(x_scale)
  }) %>%
    bind_rows() %>%
    return()
}

# Put the used/avail df in the scaled function to get columns of scaled &
#   centered values
# (_sc means "scaled & centered))
dat_used_avail <- scale_func(dat_used_avail, scale_df)

# ... c. Widen the data ----
# make the used (og), available (og and sc), used/avail,
#   and log(used/avail) for every attribute their own column
# (putting it through an lapply so I don't have to reorder all the columns)
dat_used_avail <- split(dat_used_avail, dat_used_avail$attribute)
dat_used_avail_wide <- lapply(dat_used_avail, function(x){
  k <- unique(x$attribute)
  
  x_wide <- x %>%
    group_by(hr_id) %>%
    mutate(row = row_number()) %>%
    pivot_wider(id_cols = hr_id, names_from = attribute, 
                values_from = 
                  c(used_og, avail_og, used_new, avail_new, 
                    used.avail_new, logSR_new,
                    used_sc, avail_sc, logSR_sc))
  
  # make sure HR IDs are arranged numerically
  arrange(x_wide, hr_id) %>%
    return()
}) %>%
  bind_cols() 

# binding the columns together added all of the HR IDs too, so there are repeats
# remove these
dat_used_avail_wide <- rename(dat_used_avail_wide, hr_id = hr_id...1)
dat_used_avail_wide <- dat_used_avail_wide %>%
  select(-which(grepl("hr_id...", colnames(dat_used_avail_wide))))

# make sure there aren't any repeat HRs
nrow(dat_used_avail_wide) == length(unique(dat_used_avail_wide$hr_id))

head(dat_used_avail_wide)
summary(dat_used_avail_wide)


# ---------------------------X
# ---- HR INFO DATA FRAME ----
# ---------------------------X
# Make sure there are only HRs with more than 1 day and at least 20 points
filter(hr_info, n_days <= 1 | n_pts < 20)

# make sure there aren't any repeat HRs
nrow(hr_info) == length(unique(hr_info$hr_id))

# make sure all the HR IDs in the used/avail data frame are in the HR info table
all(hr_info$hr_id %in% dat_used_avail_wide$hr_id)

# Create a 'is_M' column 
# to create a binary indicator for if the individual is a male or not
hr_info <- hr_info %>%
  mutate(is_M = ifelse(sex == "M", 1, 0))
summary(hr_info$is_M)
table(hr_info$is_M)

# ... Calculate perimeter and shape index ----
# SHAPE INDEX:
#   Area of a Circle (A) = pi*r^2
#   r = sqrt(A/pi)
#   Perimeter of a Circle = 2*pi*r = 2*pi*sqrt(A/pi)
#   Perimeter of a Polygon/Perimeter of a Cirle = 
#       Shape Index of Compactness & Elongation
#   Shape Index = P/(2*pi*sqrt(A/pi))

#   if P_polygon = P_circle, then SI = 1, 
#       meaning perfect compactness and no elongation
#   if P_polygon > P_circle, then SI > 1, 
#     meaning the polygon's shape is more elongated and less compact 
#     compared to a circle of equal area
#   When SI = 1, the shape is a circle, so SI ~ 1 would be approximately a circle
#   As SI > 1, the shape becomes more elongated & less compact

# function to find perimeter and shape index for a given HR
perim_SI <- function(x){
  perim_m <- x %>%
    filter(level == max(x$level)) %>%
    terra::vect() %>%
    terra::perim()
  area_m <- max(x$area %>% max())
  
  # calculate the shape index
  shape_index <- perim_m/(2*pi*sqrt(area_m/pi))
  
  p_a <- data.frame(perim_m, area_m, shape_index)
  return(p_a)
}

# Nest the HR sf by hr_id
hr_nest <- hr_sf %>%
  nest(data = -hr_id)
head(hr_nest)

# loop through each row in the nested data to apply the perim_SI function
# (takes ~1.5 min)
hr_info <- hr_nest %>%
  mutate(si = map(data, perim_SI)) %>%
  unnest(cols = si) %>%
  select(-data) %>%
  left_join(hr_info, ., by = "hr_id") %>%
  distinct()

# Make log Area & log(log(Shape)) columns
hr_info <- hr_info %>%
  mutate(log_area = log(area_m),
         log_shape = log(shape_index),
         log_log_shape = log(log_shape)) 

# Multiply n pts x n days and normalize
hr_info <- hr_info %>%
  mutate(n.pts_X_n.days = n_pts * n_days,
         n.pts_X_n.days_norm = (n_pts * n_days) / sum(n_pts * n_days)) 

# Relocate columns 
hr_info <- hr_info %>%
  relocate(is_M, .after = sex) %>%
  relocate(log_area, .after = area_m) %>%
  relocate(c(log_shape, log_log_shape), .after = shape_index) %>%
  relocate(c(n.pts_X_n.days, n.pts_X_n.days_norm), .after = n_days) 
head(hr_info)
summary(hr_info)


# -----------------------------------------X
# ---- SUBSET DATA FOR TESTING/TRAINING ----
# -----------------------------------------X
# For the LMMs, we want to try keeping out 10% of the data as a testing set
#   (10% for each species-season-sex datafold)
# Check how much data would be left (training) and how much you'd have in testing set
# If there are at least 5 in each species-season-sex, then use this split

# nest by species-sex-season
dat_nest <- hr_info %>%
  select(c(hr_id, species, sex, season)) %>%
  nest(cols = hr_id)
dat_nest

# check how many rows of data there would be if we removed 10% of each datafold
lapply(dat_nest$cols, function(x){
  return(round(nrow(x) * 0.1))
}) %>% unlist()
# Since there are at least 5 rows for each split, we'll use this 

subset_info <- lapply(1:nrow(dat_nest), function(i){
  x <- unnest(dat_nest[i,], cols = cols)
  
  # Pull a random set of 10% of the data of this datafold
  set.seed(1234)
  subset_i <- sample(1:nrow(x), nrow(x) * 0.1)
  
  # Make a column called 'subset' that flags if that row should be 
  # testing or training
  x$subset <- NA
  x$subset[subset_i] <- "Test"
  x$subset[-subset_i] <- "Train"
  
  return(x)
}) %>%
  bind_rows() %>%
  arrange(hr_id)

# Check how many rows are in each set
table(subset_info$subset)
table(subset_info$subset)["Test"]/table(subset_info$subset)["Train"]

# Check how many rows per species-season
subset_info %>%
  group_by(species, season, subset) %>%
  tally()


# -------------X
# ---- SAVE ----
# -------------X
write.csv(dat_used_avail_wide, paste0(out_dir, "used_avail_wide.csv"), row.names = F)
write.csv(hr_info, paste0(out_dir, "hr_info_full_si.csv"), row.names = F)
write.csv(subset_info, paste0(out_dir, "subset_info.csv"), row.names = F)
write.csv(scale_df, paste0(out_dir, "scale_df.csv"), row.names = F)
