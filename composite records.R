Diss <- read.csv("~/Mirror/Bayesian/initial Diss Mere reconstruction.csv")
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
library(tidyr)
library(rasterVis)

lake_coords <- data.frame(
  site_id = c("Naut", "Diss", "TFS", "Cimera", "ZAB"),
  lon = c(24.24, 1.6, 12.31, 4.18, 21.9836),
  lat = c(61.48, 52.2, 53.35, 40.15, 54.1318)
)

diss_tidy <- Diss %>%
  filter(yBP >= 0, yBP <= 1930) %>%
  mutate(year = round(yBP),
         site_id = "Diss",
         temp_anomaly = temp_med) %>%
  dplyr::select(site_id, year, temp_anomaly)  

naut_tidy <- Naut %>%
  filter(Age_AD >= 0, Age_AD <= 1930) %>%
  mutate(year = round(Age_AD),
         site_id = "Naut",
         temp_anomaly = temp_med) %>%  
  dplyr::select(site_id, year, temp_anomaly) 

TFS <- TFS %>%
  mutate(AgeAD = round(1950 - Age..cal.BP.)) %>%
  filter(AgeAD >= 0, AgeAD <= 1930)

tfs_tidy <- TFS %>%
  mutate(year = AgeAD,
         site_id = "TFS",
         temp_anomaly = MAT) %>%  
  dplyr::select(site_id, year, temp_anomaly) 

cimera_tidy <- Cimera %>%
  filter(Age_AD >= 0, Age_AD <= 1930) %>%
  mutate(year = round(Age_AD),
         site_id = "Cimera",
         temp_anomaly = temp_med) %>%  # Replace TEMP with actual temp column name
  dplyr::select(site_id, year, temp_anomaly) 

zab_tidy <- ZAB %>%
  filter(Age..CE. >= 0, Age..CE. <= 1930) %>%
  mutate(year = round(Age..CE.),
         site_id = "ZAB",
         temp_anomaly = MAT) %>%  # Replace TEMP with actual temp column name
  dplyr::select(site_id, year, temp_anomaly) 




all_tidy <- bind_rows(
  diss_tidy,
  naut_tidy,
 # tfs_tidy,
  cimera_tidy,
  zab_tidy
)

# Calculate the mean temperature anomaly for each year
mean_anomaly <- all_tidy %>%
  group_by(year) %>%
  summarise(mean_temp_anomaly = mean(temp_anomaly, na.rm = TRUE),
            n_sites = n()) %>%
  ungroup()

ggplot()+
  geom_line(data = diss_tidy, aes(x=year, y=temp_anomaly), colour = "grey")+
  geom_line(data = naut_tidy, aes(x=year, y=temp_anomaly), colour = "grey")+
  geom_line(data = cimera_tidy, aes(x=year, y=temp_anomaly), colour = "grey")+
  geom_line(data = zab_tidy, aes(x=year, y=temp_anomaly), colour = "grey")+
  geom_line(data = mean_anomaly, aes(x=year, y=rollmean(mean_temp_anomaly, 10, na.pad =T)), colour = "navy")+
  ggpubr::theme_pubr()
  
  


plot(mean_anomaly$year, mean_anomaly$mean_temp_anomaly, t='l')
