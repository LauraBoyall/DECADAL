
#
#### Libraries ####
library(magrittr)
library(tidyverse)
library(geoChronR)
library(lubridate)
library(progress)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)



setwd("~/Documents/Documents – Laura’s MacBook Pro/Hydroclimate/")
load("~/Downloads/Pages2kTemperature2_1_2.RData")
# Remove extraneous objects
rm(D, TS)

shapes = c("GlacierIce" = 20, "LakeSediment" = 17, "Speleothem" = 15, "MolluskShells" = 15,
           "MarineSediment" = 17, "Wood" = 18, "TerrestrialSediement" = 17)


world <- ne_countries(scale = "medium", returnclass = "sf")


library(purrr)

temp_ts <- keep(sTS, function(x) {
  isTRUE(x[["climateInterpretation1_variable"]] == "T") &&
    !is.null(x[["paleoData_values"]]) &&
    !is.null(x[["year"]])
})



# Convert to a dataframe or list of dataframes
temp_data <- map(temp_ts, function(x) {
  tibble(
    year = x$year,
    temperature = x$paleoData_values,
    site = x$geo_siteName %||% NA,
    archive = x$paleoData_archiveType %||% NA,
    proxy = x$paleoData_proxy %||% NA,
    tsid = x$paleoData_TSid %||% NA
  )
})

# Combine into a single dataframe
df_temp <- bind_rows(temp_data, .id = "record_id")



site_meta <- map_dfr(temp_ts, function(x) {
  tibble(
    site = x$paleoData_TSid %||% NA,
    lat = x$geo_latitude %||% NA,
    lon = x$geo_longitude %||% NA,
    archive = x$paleoData_archiveType %||% NA
  )
})

# Plot on a world map
world_map <- map_data("world")
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "gray90", color = "gray70") +
  geom_point(data = site_meta, aes(x = lon, y = lat, color = archive), size = 2) +
  theme_minimal() +
  labs(title = "Temperature Proxy Site Locations", x = "Longitude", y = "Latitude")



is_high_resolution <- function(x) {
  yrs <- x$year
  if (length(yrs) < 2 || any(is.na(yrs))) return(FALSE)
  median(diff(sort(yrs))) <= 5
}

temp_ts_highres <- keep(temp_ts, is_high_resolution)



library(tidyr)
library(ggplot2)

df_temp_stack <- map_dfr(temp_ts_highres, function(x) {
  tibble(
    year = x$year,
    value = x$paleoData_values,
    site = x$paleoData_TSid %||% NA
  )
}, .id = "record_id")

# Optional: normalize or offset (z-score)
df_temp_stack <- df_temp_stack %>%
  group_by(site) %>%
  mutate(value_z = scale(value)[, 1]) %>%
  ungroup()

df_temp_stack<- df_temp_stack%>%
  filter(year >=0)
  

# Plot as stacked timeseries (z-scored)
ggplot(df_temp_stack, aes(x = year, y = value_z, group = site)) +
  geom_line(alpha = 0.6) +
  theme_minimal() +
  labs(title = "Stacked Temperature Time Series (Z-scored)",
       x = "Year", y = "Z-scored Temperature Anomaly")

###### Annual data
is_annual <- function(x) {
  yrs <- x$year
  if (length(yrs) < 2 || any(is.na(yrs))) return(FALSE)
  all(diff(sort(yrs)) == 1)
}

temp_ts_annual <- keep(temp_ts_highres, is_annual)

excluded_archives <- c("MarineSediment") #, "Speleothem")

temp_ts_filtered <- keep(temp_ts_annual, function(x) {
  !(x$paleoData_archiveType %||% "") %in% excluded_archives
})


europe_ts <- keep(temp_ts_filtered, function(x) {
  lat <- x$geo_latitude %||% NA
  lon <- x$geo_longitude %||% NA
  !is.na(lat) && !is.na(lon) && lat >= 35 && lat <= 70 && lon >= -10 && lon <= 40
})



###### EOF
eof_df <- map_dfr(europe_ts, function(x) {
  interp_vals <- approx(x$year, x$paleoData_values, xout = year_grid)$y
  tibble(year = year_grid, value = interp_vals, site = x$paleoData_TSid %||% NA)
})

eof_matrix <- eof_df %>%
  pivot_wider(names_from = site, values_from = value) %>%
  column_to_rownames("year") %>%
  as.matrix()

# Remove columns (sites) with too many NAs
eof_matrix <- eof_matrix[, colMeans(is.na(eof_matrix)) < 0.2]
eof_matrix <- scale(eof_matrix, center = TRUE, scale = TRUE)
eof_matrix <- na.approx(eof_matrix, na.rm = FALSE)  # Linear interpolate NAs



library(zoo)

library(zoo)
library(dplyr)
library(tidyr)

# 1. Rebuild the wide matrix from long data (eof_df from earlier step)
eof_matrix <- eof_df %>%
  pivot_wider(names_from = site, values_from = value) %>%
  arrange(year)  # very important!

# 2. Save years as a separate vector
eof_years <- eof_matrix$year

# 3. Drop 'year' column to get numeric matrix
eof_values <- eof_matrix %>%
  select(-year) %>%
  as.matrix()

# 4. Interpolate and clean missing values
eof_filled <- apply(eof_values, 2, function(col) {
  col <- na.approx(col, na.rm = FALSE)
  col <- na.locf(col, na.rm = FALSE)
  col <- na.locf(col, fromLast = TRUE, na.rm = FALSE)
  col[is.na(col)] <- mean(col, na.rm = TRUE)
  col
})

# 5. Make sure everything is finite
eof_filled <- eof_filled[rowSums(is.finite(eof_filled)) == ncol(eof_filled), ]
eof_years_clean <- eof_years[rowSums(is.finite(eof_filled)) == ncol(eof_filled)]

# 6. Run PCA
eof_res <- prcomp(eof_filled, center = TRUE, scale. = TRUE)

# 7. Plot PC1
plot(eof_years_clean, eof_res$x[, 1], type = "l",
     main = "First EOF", xlab = "Year", ylab = "PC1")


########
# Get site names in same order as in EOF matrix
site_names <- colnames(eof_filled)

# Match them to metadata (from earlier)
site_meta_map <- map_dfr(temp_ts_filtered, function(x) {
  tibble(
    site = x$paleoData_TSid %||% NA,
    lat = x$geo_latitude %||% NA,
    lon = x$geo_longitude %||% NA
  )
})

# Join loadings to coordinates
loading_df <- tibble(
  site = site_names,
  loading = eof_res$rotation[, 1]
) %>%
  left_join(site_meta_map, by = "site")

library(ggplot2)
library(maps)

world_map <- map_data("world")

ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70") +
  geom_point(data = loading_df, aes(x = lon, y = lat, color = loading),
             size = 3) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  coord_quickmap(xlim = c(-10, 40), ylim = c(35, 70)) +  # Europe bounds
  theme_minimal() +
  labs(title = "EOF1 Spatial Loadings temp2k", x = "Longitude", y = "Latitude",
       color = "Loading")



################## Now with SCUBIDO
######## Testing creating an EOF
Diss <- read.csv("~/Mirror/Bayesian/1700 new reconstruction/Diss Mere reconstruction.csv")
Naut <- read.csv("~/Mirror/Bayesian/Naut reconstruction final.csv")
TFS <- read.csv("TieferSee old.csv")
Cimera <- read.csv("Cimera.csv")
ZAB <- read.csv("ZAB old.csv")


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



# Make a list of standardized lake datasets
lake_list <- list(
  Diss = Diss %>% rename(year = yAD, value = temp_med),
  Naut = Naut %>% rename(year = Age_AD, value = temp_med),
  TFS = TFS %>% rename(year = AgeAD, value = MAT),
  Cimera = Cimera %>% rename(year = Age_AD, value = temp_med),
  ZAB = ZAB %>% rename(year = Age..CE., value = MAT)
)

# Add site names
lake_list <- imap(lake_list, ~ mutate(.x, site = .y))

lake_long <- bind_rows(lake_list)

lake_wide <- lake_long %>%
  group_by(year, site) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = site, values_from = value) %>%
  arrange(year)


lake_years <- lake_wide$year
lake_matrix <- lake_wide %>% select(-year) %>% as.matrix()

lake_matrix <- lake_wide %>%
  select(-year) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

lake_matrix_filled <- apply(lake_matrix, 2, function(col) {
  col <- zoo::na.approx(col, na.rm = FALSE)
  col <- zoo::na.locf(col, na.rm = FALSE)
  col <- zoo::na.locf(col, fromLast = TRUE, na.rm = FALSE)
  col[is.na(col)] <- mean(col, na.rm = TRUE)
  col
})

lake_matrix_filled <- lake_matrix_filled[complete.cases(lake_matrix_filled), ]
lake_years_clean <- lake_years[complete.cases(lake_matrix_filled)]

lake_eof <- prcomp(lake_matrix_filled, center = TRUE, scale. = TRUE)

plot(lake_years_clean, lake_eof$x[, 1], type = "l",
     main = "EOF1: Lakes Only", xlab = "Year", ylab = "PC1")


site_names <- colnames(lake_matrix_filled)


lake_coords <- data.frame(
  site = site_names,
  lon = c(-5.18, 24.24, 12.31,  21.9836, 1.6),   # Naut, Diss, TFS, Cimera, ZAB
  lat = c( 40.15, 61.48,  53.35, 54.1318, 52.2), # Cimera, Naut, TFS, ZAB, Diss
  EOF1 = lake_eof$rotation[, 1]
)


ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70") +
  geom_point(data = lake_coords, aes(x = lon, y = lat, color = EOF1), size = 4) +
  geom_text(data = lake_coords, aes(x = lon, y = lat, label = site), vjust = -1.2, size = 3) +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red") +
  coord_quickmap(xlim = c(-10, 40), ylim = c(35, 70)) +
  theme_minimal() +
  labs(title = "EOF1 Spatial Loadings (Lake Sites)", x = "Longitude", y = "Latitude", color = "Loading")



##### comparing : 
loading_df <- loading_df %>%
  mutate(dataset = "Pages2k")

lake_coords_labeled <- lake_coords %>%
  select(site, lon, lat, EOF1) %>%
  rename(loading = EOF1) %>%
  mutate(dataset = "Lake")

combined_loadings <- bind_rows(loading_df, lake_coords_labeled)
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70") +
  geom_point(data = combined_loadings,
             aes(x = lon, y = lat, color = loading, shape = dataset),
             size = 3) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_shape_manual(values = c("Pages2k" = 16, "Lake" = 15)) +  # circle and square
  coord_quickmap(xlim = c(-10, 40), ylim = c(35, 70)) +
  theme_minimal()+
  labs(title = "EOF1 Loadings: Pages2k vs Lake Records",
       x = "Longitude", y = "Latitude",
       color = "Loading", shape = "Dataset")




###### interpolation: 

library(akima)   # For interpolation
library(ggplot2)

# Only interpolate the Pages2k data (not lake sites)
interp_data <- loading_df %>%
  filter(!is.na(lat), !is.na(lon), !is.na(loading))

interp_grid <- with(interp_data, interp(
  x = lon,
  y = lat,
  z = loading,
  duplicate = "mean",
  nx = 100, ny = 100, extrap = TRUE
))

# Convert to data frame for ggplot
interp_df <- expand.grid(lon = interp_grid$x, lat = interp_grid$y)
interp_df$loading <- as.vector(interp_grid$z)

ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70") +
  
  # Interpolated background
  library(akima)

# Combine and filter to ensure all values needed for interpolation are present
interp_all <- combined_loadings %>%
  filter(!is.na(lat), !is.na(lon), !is.na(loading))

# Perform 2D interpolation over the region
interp_grid_all <- with(interp_all, interp(
  x = lon,
  y = lat,
  z = loading,
  duplicate = "mean",
  nx = 100, ny = 100, 
  extrap = TRUE 
))

# Convert interpolation result to dataframe for ggplot
interp_df_all <- expand.grid(lon = interp_grid_all$x, lat = interp_grid_all$y)
interp_df_all$loading <- as.vector(interp_grid_all$z)

ggplot() +
  # World map background
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70") +
  
  # Interpolated background (smoothed color)
  geom_raster(data = interp_df_all, aes(x = lon, y = lat, fill = loading), alpha = 0.7) +
  
  # Points from both datasets
  geom_point(data = combined_loadings,
             aes(x = lon, y = lat, color = loading, shape = dataset),
             size = 3) +
  
  # Site labels
  geom_text(data = combined_loadings,
            aes(x = lon, y = lat, label = site),
            size = 2.5, vjust = -0.5, hjust = -0.1, check_overlap = TRUE) +
  
  # Color and shape styling
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "Interpolated\nEOF1") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        name = "Point\nEOF1") +
  scale_shape_manual(values = c("Pages2k" = 16, "Lake" = 15)) +  # circle vs square
  
  coord_quickmap(xlim = c(-10, 40), ylim = c(35, 70)) +
  theme_minimal() +
  labs(title = "Interpolated EOF1 Field (Pages2k + Lakes)",
       x = "Longitude", y = "Latitude",
       shape = "Dataset")


##### trying to fill the space: 
library(sp)
library(gstat)

# Create SpatialPointsDataFrame
coordinates(interp_all) <- ~ lon + lat
proj4string(interp_all) <- CRS("+proj=longlat +datum=WGS84")

# Define a full grid over the map extent
grid <- expand.grid(
  lon = seq(-25, 45, by = 0.5),
  lat = seq(30, 75, by = 0.5)
)
coordinates(grid) <- ~ lon + lat
gridded(grid) <- TRUE
proj4string(grid) <- CRS("+proj=longlat +datum=WGS84")


# Create gstat object
g_model <- gstat(formula = loading ~ 1, data = interp_all)

# Automatic variogram (or define manually)
vario <- variogram(loading ~ 1, interp_all)
fit <- fit.variogram(vario, model = vgm("Sph"))

# Update model with fit
g_model$variogram.model <- fit

# Predict across the grid
kriged <- predict(g_model, grid)

krig_df <- as.data.frame(kriged)
names(krig_df)[1:3] <- c("lon", "lat", "loading")


ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70") +
  geom_raster(data = krig_df, aes(x = lon, y = lat, fill = loading), alpha = 0.7) +
  geom_point(data = combined_loadings, aes(x = lon, y = lat, color = loading, shape = dataset), size = 3) +
  geom_text(data = combined_loadings, aes(x = lon, y = lat, label = site), size = 2.5, vjust = -0.5, hjust = -0.1, check_overlap = TRUE) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "EOF1 (Interpolated)") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "EOF1 (Point)") +
  scale_shape_manual(values = c("Pages2k" = 16, "Lake" = 15)) +
  coord_quickmap(xlim = c(-25, 45), ylim = c(30, 75)) +
  theme_minimal() +
  labs(title = "EOF1 Field Kriging Interpolation",
       x = "Longitude", y = "Latitude", shape = "Dataset")



###### timeseries
pages_ts <- tibble(
  year = eof_years_clean,
  PC1 = eof_res$x[, 1],
  dataset = "Pages2k"
)

lake_ts <- tibble(
  year = lake_years_clean,
  PC1 = lake_eof$x[, 1],
  dataset = "Lake"
)


combined_ts <- bind_rows(pages_ts, lake_ts)

# Set factor levels so "Lake" is drawn last (on top)
combined_ts$dataset <- factor(combined_ts$dataset, levels = c("Pages2k", "Lake"))

ggplot() +
  geom_line(data = pages_ts, 
            aes(x = year, y = rollmean(PC1, 5, na.pad = TRUE), color = "Pages2k"), 
            size = 0.6) +
  geom_line(data = lake_ts, 
            aes(x = year, y = PC1, color = "Lake"), 
            size = 0.6) +
  scale_color_manual(values = c("Pages2k" = "gray60", "Lake" = "blue")) +
  theme_minimal() +
  labs(title = "EOF1 Time Series",
       x = "Year", y = "PC1 Score", color = "Dataset")





######### EOF2
eof2_df <- tibble(
  site = colnames(eof_filled),
  EOF2 = eof_res$rotation[, 2]
)

lake2_df <- tibble(
  site = colnames(lake_matrix),
  EOF2 = lake_eof$rotation[, 2]
)

combined_loadings <- combined_loadings %>%
  left_join(bind_rows(eof2_df, lake2_df), by = "site")


interp_all2 <- combined_loadings %>%
  filter(!is.na(EOF2)) %>%
  select(lon, lat, loading = EOF2)

coordinates(interp_all2) <- ~ lon + lat
proj4string(interp_all2) <- CRS("+proj=longlat +datum=WGS84")

g_model2 <- gstat(formula = loading ~ 1, data = interp_all2)
vario2 <- variogram(loading ~ 1, interp_all2)
fit2 <- fit.variogram(vario2, model = vgm("Sph"))
g_model2$variogram.model <- fit2

kriged2 <- predict(g_model2, grid)

krig_df2 <- as.data.frame(kriged2)
names(krig_df2)[1:3] <- c("lon", "lat", "loading")

ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70") +
  geom_raster(data = krig_df2, aes(x = lon, y = lat, fill = loading), alpha = 0.7) +
  geom_point(data = combined_loadings, aes(x = lon, y = lat, color = EOF2, shape = dataset), size = 3) +
  geom_text(data = combined_loadings, aes(x = lon, y = lat, label = site), size = 2.5, vjust = -0.5, hjust = -0.1, check_overlap = TRUE) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "EOF2 (Interpolated)") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "EOF2 (Point)") +
  scale_shape_manual(values = c("Pages2k" = 16, "Lake" = 15)) +
  coord_quickmap(xlim = c(-25, 45), ylim = c(30, 75)) +
  theme_minimal() +
  labs(title = "EOF2 Field Kriging Interpolation",
       x = "Longitude", y = "Latitude", shape = "Dataset")




#### comparison
# EOF1 map
map_eof1 <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70") +
  geom_raster(data = krig_df, aes(x = lon, y = lat, fill = loading), alpha = 0.7) +
  geom_point(data = combined_loadings, aes(x = lon, y = lat, color = loading, shape = dataset), size = 3) +
  geom_text(data = combined_loadings, aes(x = lon, y = lat, label = site), size = 2.5, vjust = -0.5, hjust = -0.1, check_overlap = TRUE) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "EOF1 (Interpolated)") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "EOF1 (Point)") +
  scale_shape_manual(values = c("Pages2k" = 16, "Lake" = 15)) +
  coord_quickmap(xlim = c(-25, 45), ylim = c(30, 75)) +
  theme_minimal() +
  labs(title = "EOF1 Field", x = "Longitude", y = "Latitude", shape = "Dataset")

# EOF2 map
map_eof2 <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70") +
  geom_raster(data = krig_df2, aes(x = lon, y = lat, fill = loading), alpha = 0.7) +
  geom_point(data = combined_loadings, aes(x = lon, y = lat, color = EOF2, shape = dataset), size = 3) +
  geom_text(data = combined_loadings, aes(x = lon, y = lat, label = site), size = 2.5, vjust = -0.5, hjust = -0.1, check_overlap = TRUE) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "EOF2 (Interpolated)") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "EOF2 (Point)") +
  scale_shape_manual(values = c("Pages2k" = 16, "Lake" = 15)) +
  coord_quickmap(xlim = c(-25, 45), ylim = c(30, 75)) +
  theme_minimal() +
  labs(title = "EOF2 Field", x = "Longitude", y = "Latitude", shape = "Dataset")

map_eof1 <- map_eof1 + theme(legend.position = "none")
map_eof2 <- map_eof2 + theme(legend.position = "none")

library(patchwork)
(map_eof1 + map_eof2) +
  plot_layout(guides = "collect") 
