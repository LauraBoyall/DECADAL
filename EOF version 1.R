#### A code to download and plot the ERA5 data and eventually create and compare the EOFs
######## Testing creating an EOF
Diss <- read.csv("~/Mirror/Bayesian/1700 new reconstruction/Diss Mere reconstruction.csv")
Naut <- read.csv("~/Mirror/Bayesian/Naut reconstruction final.csv")
TFS <- read.csv("~/Mirror/Bayesian/Bayesian Project/TieferSee old.csv")
Cimera <- read.csv("~/Mirror/Bayesian/Bayesian Project/Cimera.csv")
ZAB <- read.csv("~/Mirror/Bayesian/Bayesian Project/ZAB old.csv")


library(dplyr)
library(ggplot2)
library(maps)
library(gstat)
library(sp)
library(raster)

library(rasterVis)


Diss <- Diss %>%
  filter(yAD >= 0 & yAD<=  1930)
Diss$yAD <- round(Diss$yAD)

Naut <- Naut %>%
  filter(Age_AD >= 0 & Age_AD<=  1930)
Naut$Age_AD <- round(Naut$Age_AD)

TFS$AgeAD <- 1950 - TFS$Age..cal.BP.

TFS <- TFS %>%
  filter(AgeAD >= 0 & AgeAD<=  1930)
TFS$AgeAD <- round(TFS$AgeAD)


Cimera <- Cimera %>%
  filter(Age_AD >= 0 & Age_AD<=  1930)
Cimera$Age_AD <- round(Cimera$Age_AD)


ZAB <- ZAB %>%
  filter(Age..CE. >= 0 & Age..CE.<=  1930)



X <- cbind(Naut$temp_med, Diss$temp_med, TFS$MAT, Cimera$temp_med, ZAB$MAT)

C <- cov(X)
eof <- prcomp(X, center = FALSE, scale. = FALSE)

lake_coords$EOF2 <- eof$rotation[, 2]
lake_coords <- data.frame(
  lon = c(24.24, 1.6, 12.31, 4.18, 21.9836), # naut = 1, TFS = 2, Cimera = 4, ZAB = 5
  lat = c(61.48, 52.2, 53.35, 40.15, 54.1318), # Diss = 2
  EOF1 = eof$rotation[,1]
)


europe_map <- map_data("world", region = c("UK", "France", "Germany", "Poland", "Italy", "Spain", "Norway", "Sweden", "Finland", "Austria", "Switzerland", "Czech Republic", "Slovakia", "Hungary"))

ggplot() +
  geom_polygon(data = europe_map, aes(x = long, y = lat, group = group), fill = "grey90", color = "grey40") +
  geom_point(data = lake_coords, aes(x = lon, y = lat, color = EOF1), size = 6) +
  scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  #coord_fixed(1.3, xlim = c(0, 30), ylim = c(40, 60)) +  # adjust as needed
  theme_minimal() +
  labs(title = "Spatial EOF1 Loadings from Lake Temperature Records",
       x = "Longitude", y = "Latitude", color = "EOF1")




#### interpolating 

coordinates(lake_coords) <- ~lon+lat
proj4string(lake_coords) <- CRS("+proj=longlat +datum=WGS84")

# europe grid
lon_range <- seq(-10, 35, by = 0.5)
lat_range <- seq(35, 70, by = 0.5)
grid <- expand.grid(lon = lon_range, lat = lat_range)
coordinates(grid) <- ~lon+lat
gridded(grid) <- TRUE
proj4string(grid) <- CRS("+proj=longlat +datum=WGS84")

idw_result <- idw(formula = EOF2 ~ 1, locations = lake_coords, newdata = grid, idp = 2.0)
idw_result <- idw(formula = EOF1 ~ 1, locations = lake_coords, newdata = grid, idp = 2.0)
r <- raster(idw_result)

r_df <- as.data.frame(r, xy = TRUE)
names(r_df)[3] <- "EOF1"
names(r_df)[3] <- "EOF2"

# Europe map
europe_map <- map_data("world", region = c("UK", "France", "Germany", "Poland", "Italy", "Spain", "Norway", "Sweden", "Finland"))

ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = EOF2)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  geom_polygon(data = europe_map, aes(x = long, y = lat, group = group), fill = NA, color = "grey40") +
  coord_fixed(xlim = c(-10, 35), ylim = c(35, 70)) +
  theme_minimal() +
  labs(title = "Interpolated EOF1 Surface", x = "Longitude", y = "Latitude", fill = "EOF1")


##### now the ERA5 reanalysis


#Load the libraries

library("ecmwfr")
library(terra)
library(sf)
library(raster)


# set the key 
wf_set_key(key = "3e7b54fd-4971-4c7e-9cab-0488a600b3ee")


# Define years, months, and days
Year <- 1991:2020
Month <- sprintf("%02d", 1:12)  # Format months as "01" to "12"


# Define the request for ERA5 monthly averaged data
request <- list(
  dataset_short_name = "reanalysis-era5-land-monthly-means",
  product_type = "reanalysis",
  variable = "2m_temperature",
  year = Year,
  month = Month,
  time = "00:00",
  format = "netcdf",               
  area = c(90, -10, 30, 30),  #Europe, north, west, south and east
  target = "era5-temp2m.nc"         
)

# call:

file_1 <- wf_request(
  request  = request,  # the request
  transfer = TRUE,     # download the file
  path     = "."       # store data in current working directory
)

#Open the NC file

tp_raster_stack <- brick(unzip(file_1), varname=t2m)


nlayers(tp_raster_stack)
plot(tp_raster_stack )

#Visualize the mean layer

mean_layer <- mean(tp_raster_stack , na.rm = TRUE)



# Plot the mean raster layer 
plot(mean_layer, main = "ERA-5 Reanalysis Data (total_precipitation 2020-21)")


# Add world map and shapefile layers
maps::map("world", add = TRUE)


lake_coords_df <- data.frame(
  lon = c(24.24, 1.6, 12.31, 4.18, 21.9836),
  lat = c(61.48, 52.2, 53.35, 40.15, 54.1318)
)

coordinates(lake_coords_df) <- ~lon+lat
proj4string(lake_coords_df) <- CRS("+proj=longlat +datum=WGS84")

# Extract ERA5 time series at each site
reanalysis_matrix <- raster::extract(tp_raster_stack, lake_coords_df)  # rows = sites, columns = time

# Find columns without NA standard deviation
valid_cols <- which(!is.na(apply(X_reanalysis, 2, sd)))

# Filter the reanalysis matrix
X_reanalysis_clean <- X_reanalysis[, valid_cols]
X_std_reanalysis <- scale(X_reanalysis_clean)

# Run EOF
eof_reanalysis <- prcomp(X_std_reanalysis, center = FALSE, scale. = FALSE)
EOF1_re <- eof_reanalysis$rotation[, 1]
EOF2_re <- eof_reanalysis$rotation[, 2]

valid_coords <- data.frame(
  lon = c(24.24, 12.31, 21.9836)[valid_cols],  # match to valid columns
  lat = c(61.48, 53.35, 54.1318)[valid_cols],
  EOF1 = EOF1_re
)


# Convert raster to dataframe
r_df_era5 <- as.data.frame(mean_layer, xy = TRUE)
names(r_df_era5)[3] <- "ERA5_temp"

###### plot mapping the surface air temperature 
ggplot() +
  geom_raster(data = r_df_era5, aes(x = x, y = y, fill = ERA5_temp)) +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = mean(r_df_era5$ERA5_temp, na.rm = TRUE)
  ) +
  geom_polygon(data = europe_map, aes(x = long, y = lat, group = group),
               fill = NA, color = "grey40") +
  coord_fixed(xlim = c(-10, 35), ylim = c(35, 70)) +
  theme_minimal() +
  labs(
    title = "ERA5 Mean 2m Temperature (1991–2020)",
    x = "Longitude", y = "Latitude", fill = "Temp (K)"
  )




