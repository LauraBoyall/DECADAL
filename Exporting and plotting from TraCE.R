library(ggplot2) # used for plotting
library(ncdf4) # used for reading and extracting data from NetCDF model file
library(raster) # helpful package for plotting maps
library(dplyr) # helpful tool for filtering data
library(maps) # to create the map outlines
library(sf) # essential for spatial data handing
library(patchwork)
setwd("~/Mirror/DECADAL/Ashs models/")


annual_data <- ("~/Mirror/TraCE files/trace.01-36.22000BP.cam2.TS.22000BP_decavg_400BCE.nc")
JJA_data <- ("trace.01-36.22000BP.cam2.TS.22000BP_decavgJJA_400BCE.nc")
DJF_data <- ("trace.01-36.22000BP.cam2.TS.22000BP_decavgDJF_400BCE.nc")


# I have made a loop so it is much easier but feel free to run outside if easier!
extract_regional_timeseries <- function(varname, model_file) {
  nc <- nc_open(model_file)
  
  # Extract coords and time
  lat <- ncvar_get(nc, "lat")
  lon <- ncvar_get(nc, "lon")
  Time <- ncvar_get(nc, "time") * -1000  # Convert from ka BP to BP
  
  # Define regional bounds
  lon_idx <- which(lon >= -0 | lon <= 30)
  lat_idx <- which(lat >= 40 & lat <= 60)
  
  # Extract the variable
  var_data <- ncvar_get(nc, varname)
  
  # Mean over region
  regional_mean <- apply(var_data[lon_idx, lat_idx, ], 3, mean, na.rm = TRUE)
  
  # Convert to °C if it's temperature
  if (varname == "TS") {
    regional_mean <- regional_mean - 273.15
  }
  
  # Close file
  nc_close(nc)
  
  # Return as data frame
  data.frame(AgeBP = Time, regional_mean = regional_mean)
}

Annual<- extract_regional_timeseries("TS", annual_data) # TS is the temperature name of temp in the file
JJA<- extract_regional_timeseries("TS", JJA_data)
DJF<- extract_regional_timeseries("TS", DJF_data)


###############################################################################
################################ Plotting #####################################
###############################################################################
# now that we have regional_temp we can plot anything we want!

########### Creating Anomalies #####

calc_anomaly <- function(df) {
  baseline <- df %>%
    filter(AgeBP >= 0 & AgeBP <= 30) %>%
    summarise(ref_mean = mean(regional_mean, na.rm = TRUE)) %>%
    pull(ref_mean)
  
  df %>%
    mutate(anomaly = regional_mean - baseline)
}

# Apply to each time series
Annual_Holocene <- Annual %>%
  filter(AgeBP < 11700) %>%
  calc_anomaly()

DJF_Holocene <- DJF %>%
  filter(AgeBP < 11700) %>%
  calc_anomaly()

JJA_Holocene <- JJA %>%
  filter(AgeBP < 11700) %>%
  calc_anomaly()

ggplot() + 
  geom_line(data = Annual_Holocene, aes(x = AgeBP, y = anomaly), color = "black") +
  geom_line(data = DJF_Holocene, aes(x = AgeBP, y = anomaly), color = "navy") +
  geom_line(data = JJA_Holocene, aes(x = AgeBP, y = anomaly), color = "orange") +
  labs(
    title = "TraCE Holocene Temperature Anomalies ",
    x = "Age (years BP)",
    y = expression("Temperature Anomaly ("*degree*C*")")
  ) +
  ggpubr::theme_pubr()



### A plot zooming into a time window - just for JJA for now, but obvs can do for any

ggplot()+ 
  geom_line(data = JJA_Holocene, aes(x=AgeBP, y=anomaly), col = "orange")+ # Change to annuual or DJF
  labs(
    title = "JJA TraCE 21ka temperatures", # change title to match 
    x = "Age (years BP)",
    y = expression("Temperature "(degree * C))
  ) +
  scale_x_continuous(limits=c(5000,6000))+ # changing between 5-6 ka
  ggpubr::theme_pubr()


##### plotting together 
ggplot() + 
  geom_line(data = Annual_Holocene, aes(x = AgeBP, y = anomaly), color = "black") +
  geom_line(data = DJF_Holocene, aes(x = AgeBP, y = anomaly), color = "navy") +
  geom_line(data = JJA_Holocene, aes(x = AgeBP, y = anomaly), color = "orange") +
  labs(
    title = "TraCE Holocene Temperature Anomalies",
    x = "Age (years BP)",
    y = expression("Temperature Anomaly ("*degree*C*")")
  ) +
  scale_y_continuous(limits=c(-2.5,1.5))+
  scale_x_continuous(limits=c(5000,6000))+ # changing between 5-6 ka
  ggpubr::theme_pubr()



######### saving the data 
#### exporting the data 
write.csv(JJA_Holocene, "JJA Holocene TraCE.csv")
write.csv(DJF_Holocene, "DJF Holocene TraCE.csv")
write.csv(Annual_Holocene, "Annual Holocene TraCE.csv") 

###############################################################################
########################## Equilibrium plots / maps ###########################
###############################################################################
# for the equilibrium maps / time slices we are going to be using the same 
# model files but we load in again for ease
rm(list = ls()) # just clearing the environment

annual_data <- ("~/Mirror/TraCE files/trace.01-36.22000BP.cam2.TS.22000BP_decavg_400BCE.nc")
JJA_data <- ("trace.01-36.22000BP.cam2.TS.22000BP_decavgJJA_400BCE.nc")
DJF_data <- ("trace.01-36.22000BP.cam2.TS.22000BP_decavgDJF_400BCE.nc")
world_map <- map_data("world") # this is for the basemap




################## Plotting ##################################
###### ANNUAL MEAN TEMPERATURE <-- I have just repeated the same code, but highlight where I change to the different seasons
# Setup 
nc_file <- annual_data # can also be Annual_data or DJF_data
varname <- "TS"

# Defining the time windows for the plots
time_windows <- list(
  c(5000, 5500),
  c(5500, 6000),
  c(6000, 6500),
  c(6500, 7000)
)

# FUNCTION TO EXTRACT THE DATA AND CREATE ANOMALIES
anomaly_timeslice <- function(varname, nc_file, start_year, end_year) {
  nc <- nc_open(nc_file)
  Time <- ncvar_get(nc, "time") * -1000
  lat <- ncvar_get(nc, "lat")
  lon <- ncvar_get(nc, "lon")
  lon <- ifelse(lon > 180, lon - 360, lon)
  
  time_slice <- which(Time >= start_year & Time <= end_year)
  ref_slice  <- which(Time >= -30 & Time <= 0)
  
  var_data <- ncvar_get(nc, varname)
  if (varname == "TS") {
    var_data <- var_data - 273.15
  }
  
  slice_mean <- apply(var_data[, , time_slice], c(1, 2), mean, na.rm = TRUE)
  ref_mean   <- apply(var_data[, , ref_slice], c(1, 2), mean, na.rm = TRUE)
  anomaly <- slice_mean - ref_mean
  
  nc_close(nc)
  
  expand.grid(lon = lon, lat = lat) %>%
    mutate(value = as.vector(anomaly)) %>%
    arrange(lat, lon)
}

# Color scale
limits <- c(-2, 2) #<----- this is important here, some plots you will want the colour scale to be -3,+3, others more or less, annual is maybe +/- 2
col_palette <- colorRampPalette(c("navy", "white", "darkred"))(100)

# Store plots
plot_list <- list()

for (i in seq_along(time_windows)) {
  start <- time_windows[[i]][1]
  end <- time_windows[[i]][2]
  
  ts_data <- anomaly_timeslice(varname, nc_file, start, end)
  
  # Interpolate for smoothing
  interp_result <- with(ts_data, interp(
    x = lon, y = lat, z = value,
    xo = seq(min(lon), max(lon), length = 200),
    yo = seq(min(lat), max(lat), length = 200),
    duplicate = "mean"
  ))
  
  interp_df <- expand.grid(
    lon = interp_result$x,
    lat = interp_result$y
  ) %>%
    mutate(value = as.vector(interp_result$z))
  
  # Create plot without legend
  p <- ggplot(interp_df, aes(x = lon, y = lat, fill = value)) +
    geom_raster(interpolate = TRUE) +
    geom_path(data = world_map, aes(x = long, y = lat, group = group),
              color = "black", size = 0.2, inherit.aes = FALSE) +
    scale_fill_gradientn(
      colours = col_palette,
      limits = limits,
      name = "Temp Anomaly (°C)",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = 15,
        barheight = 0.5,
        ticks.colour = "black"
      )
    ) +
    coord_fixed(xlim = c(-25, 40), ylim = c(35, 70)) +
    labs(
      title = paste0(start, "–", end, " BP")
    ) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 12, hjust = 0.5),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  plot_list[[i]] <- p
}

# Combine with patchwork and shared legend
 Annual_plots <- (plot_list[[1]] + plot_list[[2]]) /
  (plot_list[[3]] + plot_list[[4]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
 
 # just adding the title etc
 Annual_plots + plot_annotation(
   title = "Annual Mean Temperature Anomalies",
   theme = theme(
     plot.title = element_text(size = 16, face = "bold", hjust = 0.5,
                               margin = margin(b = 10))
   )
 )

###############################
##############################
 ############ JJA

 nc_file <- JJA_data # can also be Annual_data or DJF_data
 varname <- "TS"
 
 
 # Defining the time windows for the plots
 time_windows <- list(
   c(5000, 5500),
   c(5500, 6000),
   c(6000, 6500),
   c(6500, 7000)
 )
 
 # FUNCTION TO EXTRACT THE DATA AND CREATE ANOMALIES
 anomaly_timeslice <- function(varname, nc_file, start_year, end_year) {
   nc <- nc_open(nc_file)
   Time <- ncvar_get(nc, "time") * -1000
   lat <- ncvar_get(nc, "lat")
   lon <- ncvar_get(nc, "lon")
   lon <- ifelse(lon > 180, lon - 360, lon)
   
   time_slice <- which(Time >= start_year & Time <= end_year)
   ref_slice  <- which(Time >= -30 & Time <= 0)
   
   var_data <- ncvar_get(nc, varname)
   if (varname == "TS") {
     var_data <- var_data - 273.15
   }
   
   slice_mean <- apply(var_data[, , time_slice], c(1, 2), mean, na.rm = TRUE)
   ref_mean   <- apply(var_data[, , ref_slice], c(1, 2), mean, na.rm = TRUE)
   anomaly <- slice_mean - ref_mean
   
   nc_close(nc)
   
   expand.grid(lon = lon, lat = lat) %>%
     mutate(value = as.vector(anomaly)) %>%
     arrange(lat, lon)
 }
 
 # Color scale
 limits <- c(-3, 3) #<----- this is important here, some plots you will want the colour scale to be -3,+3, others more or less, annual is maybe +/- 2
 col_palette <- colorRampPalette(c("navy", "white", "darkred"))(100)
 
 # Store plots
 plot_list <- list()
 
 for (i in seq_along(time_windows)) {
   start <- time_windows[[i]][1]
   end <- time_windows[[i]][2]
   
   ts_data <- anomaly_timeslice(varname, nc_file, start, end)
   
   # Interpolate for smoothing
   interp_result <- with(ts_data, interp(
     x = lon, y = lat, z = value,
     xo = seq(min(lon), max(lon), length = 200),
     yo = seq(min(lat), max(lat), length = 200),
     duplicate = "mean"
   ))
   
   interp_df <- expand.grid(
     lon = interp_result$x,
     lat = interp_result$y
   ) %>%
     mutate(value = as.vector(interp_result$z))
   
   # Create plot without legend
   p <- ggplot(interp_df, aes(x = lon, y = lat, fill = value)) +
     geom_raster(interpolate = TRUE) +
     geom_path(data = world_map, aes(x = long, y = lat, group = group),
               color = "black", size = 0.2, inherit.aes = FALSE) +
     scale_fill_gradientn(
       colours = col_palette,
       limits = limits,
       name = "Temp Anomaly (°C)",
       guide = guide_colorbar(
         title.position = "top",
         title.hjust = 0.5,
         barwidth = 15,
         barheight = 0.5,
         ticks.colour = "black"
       )
     ) +
     coord_fixed(xlim = c(-25, 40), ylim = c(35, 70)) +
     labs(
       title = paste0(start, "–", end, " BP")
     ) +
     theme_void() +
     theme(
       legend.position = "none",
       plot.title = element_text(size = 12, hjust = 0.5),
       plot.margin = margin(5, 5, 5, 5)
     )
   
   plot_list[[i]] <- p
 }
 
 # Combine with patchwork and shared legend
JJA_plots <- (plot_list[[1]] + plot_list[[2]]) /
   (plot_list[[3]] + plot_list[[4]]) +
   plot_layout(guides = "collect") &
   theme(legend.position = "bottom")

# just adding the title etc
JJA_plots + plot_annotation(
  title = "JJA Temperature Anomalies",
  theme = theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5,
                              margin = margin(b = 10))
  )
)
 
 #####################################
#####################################
##################################### 
#DJF plots

nc_file <- DJF_data # can also be Annual_data or DJF_data
varname <- "TS"

# Load world map
world_map <- map_data("world")

# Defining the time windows for the plots
time_windows <- list(
  c(5000, 5500),
  c(5500, 6000),
  c(6000, 6500),
  c(6500, 7000)
)

# FUNCTION TO EXTRACT THE DATA AND CREATE ANOMALIES
anomaly_timeslice <- function(varname, nc_file, start_year, end_year) {
  nc <- nc_open(nc_file)
  Time <- ncvar_get(nc, "time") * -1000
  lat <- ncvar_get(nc, "lat")
  lon <- ncvar_get(nc, "lon")
  lon <- ifelse(lon > 180, lon - 360, lon)
  
  time_slice <- which(Time >= start_year & Time <= end_year)
  ref_slice  <- which(Time >= -30 & Time <= 0)
  
  var_data <- ncvar_get(nc, varname)
  if (varname == "TS") {
    var_data <- var_data - 273.15
  }
  
  slice_mean <- apply(var_data[, , time_slice], c(1, 2), mean, na.rm = TRUE)
  ref_mean   <- apply(var_data[, , ref_slice], c(1, 2), mean, na.rm = TRUE)
  anomaly <- slice_mean - ref_mean
  
  nc_close(nc)
  
  expand.grid(lon = lon, lat = lat) %>%
    mutate(value = as.vector(anomaly)) %>%
    arrange(lat, lon)
}

# THIS IS IMPORTANT, I HAVE CHANGED THE LIMITS E.G. LEDGEND WILL GO FROM -4 TO +4 BUT NOT IN JJA / ANNUAL E.G. 
limits <- c(-4, 4) #<----- this is important here, some plots you will want the colour scale to be -3,+3, others more or less, annual is maybe +/- 2
col_palette <- colorRampPalette(c("navy", "white", "darkred"))(100)

# making a list of each of the timeslice maps
plot_list <- list()

for (i in seq_along(time_windows)) {
  start <- time_windows[[i]][1]
  end <- time_windows[[i]][2]
  
  ts_data <- anomaly_timeslice(varname, nc_file, start, end)
  
  # Interpolate for smoothing
  interp_result <- with(ts_data, interp(
    x = lon, y = lat, z = value,
    xo = seq(min(lon), max(lon), length = 200),
    yo = seq(min(lat), max(lat), length = 200),
    duplicate = "mean"
  ))
  
  interp_df <- expand.grid(
    lon = interp_result$x,
    lat = interp_result$y
  ) %>%
    mutate(value = as.vector(interp_result$z))
  
  # Create plot without legend
  p <- ggplot(interp_df, aes(x = lon, y = lat, fill = value)) +
    geom_raster(interpolate = TRUE) +
    geom_path(data = world_map, aes(x = long, y = lat, group = group),
              color = "black", size = 0.2, inherit.aes = FALSE) +
    scale_fill_gradientn(
      colours = col_palette,
      limits = limits,
      name = "Temp Anomaly (°C)",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = 15,
        barheight = 0.5,
        ticks.colour = "black"
      )
    ) +
    coord_fixed(xlim = c(-25, 40), ylim = c(35, 70)) +
    labs(
      title = paste0(start, "–", end, " BP")
    ) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 12, hjust = 0.5),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  plot_list[[i]] <- p
}

# plotting each of them together
DJF_plots <- (plot_list[[1]] + plot_list[[2]]) /
  (plot_list[[3]] + plot_list[[4]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

 DJF_plots + plot_annotation(
  title = "DJF Temperature Anomalies",
  theme = theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5,
                              margin = margin(b = 10))
  )
)

