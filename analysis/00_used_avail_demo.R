# X================================================X
# X------------------ HR Analysis -----------------X
# X---- Demonstrate Used/Available Calculation ----X
# X================================================X

# NOTE: This is just a demonstration to show how the used/available framework was done
# For that reason, I have generated a simulated environmental raster stack
# The output of this script is not saved, as this is just a demonstration


library(terra)
library(sf)
library(tidyverse)

# ---------------X
# ---- SET UP ---- 
# ---------------X
# Load the home ranges
hrs <- st_read("raw_data/HRs")
head(hrs)

# Load the home range info dataset
hr_info <- read.csv("raw_data/hr_info.csv")
head(hr_info)

# ---------------------------------------------------X
# ---- CREATE A SIMULATED STACK OF FORAGE RASTERS ---- 
# ---------------------------------------------------X
# For each year-season combination, 
#   create a simulated NDVI x herbaceous cover (forage) raster
#   and stack together

# First, create a template raster 
# set the extent to the extnd of the home ranges and the cell resolution to 30m
template <- rast(crs = "epsg:32612", extent = ext(st_bbox(hrs)), resolution = 30)
template
(ncell <- ncell(template))

# Loop through every year and season
(years <- seq(min(hr_info$year), max(hr_info$year)))
(seasons <- unique(hr_info$month))

system.time({
  forage <- rast(lapply(years, function(year){
    r <- rast(lapply(seasons, function(ssn){
      cat(year, ssn, "\n")
      set.seed(year * ssn)
      
      # create a vector of ndvi values from around -1 to around 1
      ndvi <- seq(rnorm(1, -1, 0.5), rnorm(1, 1, 0.5), 
                  length.out = ncell)
      
      # add some random noise
      rand_ndvi <- rnorm(ncell, mean = 0, sd = abs(rnorm(1, 0, 0.1)))
      ndvi <- ndvi + rand_ndvi
      
      # re-set the bounds to -1 to 1
      ndvi <- ifelse(ndvi < -1, -1, 
                     ifelse(ndvi > 1, 1, ndvi))
      
      # create a vector of herbaceous cover values from around 0 to around 100
      herb_cover <- seq(rnorm(1, 0, 1), 
                        rnorm(1, 100, 1), length.out = ncell)
      
      # add some random noise
      rand_herb <- rnorm(ncell, mean = 0, sd = abs(rnorm(1)))
      herb_cover <- herb_cover + rand_herb
      
      # re-set the bounds to 0 to 100
      
      herb_cover <- ifelse(herb_cover < 0, 0, 
                           ifelse(herb_cover > 100, 100, herb_cover))
      
      # multiply together to create a vector of forage values
      forage <- ndvi * herb_cover
      
      # set these forage values as this year-season's raster values
      r <- template
      names(r) <- paste(year, ssn, sep = "_")
      values(r) <- forage
      
      return(r)
    }))
    return(r)
  }))
}) # took my system ~30 min
forage
plot(forage)
gc()


# -----------------------------------X
# ---- PERFORM USED-AVAIL ANALYIS ---- 
# -----------------------------------X
## Set up functions ----
# Funciton to rasterize the home range
rasterize_hr <- function(hr){
  # Add a new column called "intensity"
  # this is so the smaller isopleths are attributed with a higher intensity of use
  hr <- hr %>%
    arrange(level) %>%
    mutate(intensity = seq(n(), 1))
  
  # crop the template raster to the extent of the HR
  temp <- crop(template, hr)
  temp[] <- 1 # fill with values, the values don't matter
  
  # convert HR to raster
  hr_rast <- rasterize(x = hr, 
                       y = temp, 
                       # use the col "intensity" as the raster's values
                       field = hr$intensity,
                       fun = "count")
  
  # Check that the output raster isn't empty 
  # (this could happen )
  if(!all(is.na(values(hr_rast)))){
    # function to normalize data so that it sums to the max value
    # default max value is 1
    norm <- function(x, sum, max = 1){
      n_x <- (x / sum) * max
      return(n_x)
    }
    
    total <- sum(values(hr_rast), na.rm = T)
    rast_out <- norm(hr_rast, total)
    
    # check that all values sum to the max value of the UD
    if(!all.equal(sum(values(rast_out), na.rm = T), 1)){
      stop("HR raster not normalized correctly")
    }
    return(rast_out)
  }
}

# Function to generate a buffer (availability domain) around the home range
# (default buffer width = 10km = 100km^2)
generate_avail <- function(hr, buff_width = 10000){
  # find the centroid of the 95% isopleth
  center <- hr %>%
    filter(level == max(hr$level)) %>%
    st_centroid()
  
  # (half) width of buffer (width from the center) in meters
  # (side length is 2x)
  # *---------------*
  # \       \       \
  # \     x \       \
  # \       \   x   \
  # \       O-------\ 2x
  # \               \
  # \               \
  # \      2x       \
  # *---------------*
  # buf_width = 2x, so need to divide by 2
  # create blank raster of same crs and resolution as the HR,  
  # extent is the center +/- the given buffer width divided by 2)
  
  # first create an sf buffer (circular) and then find the extent (square)
  buff <- st_buffer(center, (buff_width/2))
  buff_out <- ext(buff)
  
  # return the extent of the buffer
  return(buff_out)
}

# Function to calculate used and avail of the given environmental/attribute (attr) raster
# In this demo, we are using the given year-season 'forage' raster
extract_used_avail <- function(attr_stack, used_rast, avail_rast){
  # Crop attribute raster stack to weighted used space and unweighted avail space
  attr_crop_used <- crop(attr_stack, used_rast)
  attr_crop_avail <- crop(attr_stack, avail_rast)
  
  # 1) Calculate weighted average use
    # multiply cropped environmental raster with the used space
    #   (this is because the values of the used space raster are weighted by 
    #   intensity of use, so the attribute values that fall within areas that were 
    #   used more intensely will be weighted higher)
  attr_mult_used <- attr_crop_used * used_rast
  
    # sum the values together
    # because the HR raster cells were weighted by intensity of use
    # and then normalized to sum to 1,
    # this sum will result in a weighted average
  wgt_avg_used <- sum(values(attr_mult_used), na.rm = T)
  
  # 2) Calculate average availability
    # we simply need to calculate the mean value of the attribute's values 
    # within the available space
  avg_avail <- mean(values(attr_crop_avail), na.rm = T)
  
  # return the used and available values
  used_avail <- list(wgt_avg_used, avg_avail)
  names(used_avail) <- c("used", "available")
  return(used_avail)
}

# Loop through each year-season
system.time({
  used_avail_df <- lapply(years, function(yr){
    lapply(seasons, function(ssn){
      cat(yr, ssn, "\n")
      
      # Pull the corresponding environmental raster
      r <- forage[[paste(yr, ssn, sep = "_")]]
      
      # Filter home ranges to those that were in this year and season
      hr_info_yr_ssn <- filter(hr_info, year == yr & month == ssn)
      hrs_yr_ssn <- filter(hrs, hr_id %in% unique(hr_info_yr_ssn$hr_id))
      
      # split by HR ID
      hrs_yr_ssn <- split(hrs_yr_ssn, hrs_yr_ssn$hr_id)
      
      # Loop through each home range
      lapply(hrs_yr_ssn, function(hr){
        
        # rasterize the HR and normalize so all cell values sum to 1
        hr_rast <- rasterize_hr(hr)
        
        # find the extent of the availability domain (100-km2 square around HR center)
        avail <- generate_avail(hr)
        
        # Calculate the weighted average of use (weighted by intensity of use)
        used_avail <- extract_used_avail(attr_stack = r, 
                                         used_rast = hr_rast, 
                                         avail_rast = avail)
        
        # fill use/avail table
        used_avail_df <- list(unique(hr$hr_id), used_avail[[1]], used_avail[[2]]) %>% 
          as.data.frame(optional = T, col.names = c("hr_id", "used", "available"))
        rownames(used_avail_df) <- NULL
        
        # Return
        return(used_avail_df)
      }) %>% bind_rows()
    }) %>% bind_rows()
  }) %>% bind_rows()
}) # ~15 min

head(used_avail_df)
used_avail_df %>%
  mutate(sr = used/available,
         logSR = log(sr)) %>%
  summary()
