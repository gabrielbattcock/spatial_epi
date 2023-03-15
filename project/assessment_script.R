# Assessment script

#read in all the libraries
library(geoR)
library(wesanderson) 
source("r/map_maker.R")
source('3135_packages.R')
tsource("../spatial_custom_functions.R")
pal <- wes_palette("Zissou1", 100, type = "continuous")
library(tidyr)

#-data exploration-------------------------------------------------------------------
# 
# load in the shape files
country <- st_read("vector/gadm41_LBR_0.shp")
region <-st_read("vector/gadm41_LBR_1.shp")
district <- st_read("vector/gadm41_LBR_2.shp")
x <- st_read("vector/gadm41_LBR_3.shp")

head(region)

liberia <- read.csv("data/Liberia_STH.csv")
head(liberia)
#plot hists of the prevelance 
par( mfrow= c(1,3) )
hist(liberia$HK_examined)
hist(liberia$Asc_examined)
hist(liberia$TT_examined)
# data needs to be logged to be normal
liberia$HK_logit <- log(liberia$HK_positive+0.5/(liberia$HK_examined-liberia$HK_positive+0.5))
liberia$Asc_logit <- log(liberia$Asc_positive+0.5/(liberia$Asc_examined-liberia$Asc_positive+0.5))
liberia$TT_logit <- log(liberia$TT_positive+0.5/(liberia$TT_examined-liberia$TT_positive+0.5))

par( mfrow= c(1,3) )
hist(liberia$HK_logit)
hist(liberia$Asc_logit)
hist(log(liberia$TT_logit))

# create a map of the area for each disease
liberia$TT_PREV_BIN <- cut2(liberia$TT_prevalence, cuts = c(0.2, 0.4, 0.6))
liberia$HK_PREV_BIN <- cut2(liberia$HK_prevalence, cuts = c(0.2, 0.4, 0.6))
liberia$Asc_PREV_BIN <- cut2(liberia$Asc_prevalence, cuts = c(0.2, 0.4, 0.6))

#create a map of the area
qmplot(x = Longitude, y = Latitude,
       data = liberia, source = "stamen", maptype = "toner-lite") +
  geom_point(aes(fill = TT_PREV_BIN), shape = 21, alpha = 0.8) +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 3))

qmplot(x = Longitude, y = Latitude,
       data = liberia, source = "stamen", maptype = "toner-lite") +
  geom_point(aes(fill = HK_PREV_BIN), shape = 21, alpha = 0.8) +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 4))
qmplot(x = Longitude, y = Latitude,
       data = liberia, source = "stamen", maptype = "toner-lite") +
  geom_point(aes(fill = Asc_PREV_BIN), shape = 21, alpha = 0.8) +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 4))

####
# transform the data into a spatial object

# Transform to a sf object
liberia_sp <- liberia %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_transform(crs = "+init=epsg:3857 +units=m")
head(liberia)

par(mfrow=c(1,1))
plot(liberia_sp$geometry)


plot(country$geometry, color = "white")
plot(liberia_sp["HK_logit"], add =T) 
plot(country$geometry, color = "white")
plot(liberia_sp["Asc_logit"], add = T)
plot(country$geometry, color = "white")
plot(liberia_sp["TT_logit"], add = T)


# Continuous processes #################################################################
# 

liberia_sp <- st_transform(liberia_sp, crs = "+init=epsg:3857 +units=m")
head(liberia_sp)

## plot a semivariance for hookworm 
ggvario(coords = st_coordinates(liberia_sp), data = liberia_sp$HK_logit)
ggvario(coords = st_coordinates(liberia_sp), data = liberia_sp$TT_logit)
ggvario(coords = st_coordinates(liberia_sp), data = liberia_sp$Asc_logit)

#### reset the boundaries
country <- country %>% 
st_set_crs(4326) %>%
  st_transform(crs = "+init=epsg:3857 +units=m")

region <- region %>% 
  st_set_crs(4326) %>% 
  st_transform(crs = "+init=epsg:3857 +units=m")

district <- district %>% 
  st_set_crs(4326) %>% 
  st_transform(crs = "+init=epsg:3857 +units=m")

### make a prediction for HK based on the semivariogram


pred_locations <- st_make_grid(country, cellsize = 5000, what = "centers")
pred_locations <- st_intersection(pred_locations, country)#resize the prediction locations to just the area
plot(country$geometry, col = "white")
plot(pred_locations, cex = .02, col = "red", pch = 19, add = T) 


# hook worm 
variofit <- autoKrige(formula = HK_logit ~ 1,
                      input_data = as(liberia_sp, "Spatial"),
                      new_data = as(pred_locations, "Spatial"),
                        model = c("Sph", "Exp", "Gau")) # types of models to try

HK_pred_prevalence <- plogis(variofit$krige_output$var1.pred)

HK_prevalence <- rasterFromXYZ(xyz = data.frame(st_coordinates(pred_locations), HK_pred_prevalence),
                            crs = crs(liberia_sp))
plot(HK_prevalence)
plot(country$geometry, add = T)
plot(region$geometry, add = T)

writeRaster(HK_prevalence, "outputs/HK_pred.tif",
            format = 'GTiff', overwrite = TRUE)

# TT

variofit <- autoKrige(formula = TT_logit ~ 1,
                      input_data = as(liberia_sp, "Spatial"),
                      new_data = as(pred_locations, "Spatial"),
                      model = c("Sph", "Exp", "Gau")) 

TT_pred_prevalence <- plogis(variofit$krige_output$var1.pred)

TT_prevalence <- rasterFromXYZ(xyz = data.frame(st_coordinates(pred_locations), TT_pred_prevalence),
                               crs = crs(liberia_sp))
plot(TT_prevalence)
plot(country$geometry, add = T)
plot(region$geometry, add = T)

writeRaster(TT_prevalence, "outputs/TT_pred.tif",
            format = 'GTiff', overwrite = TRUE)

# Asc

variofit <- autoKrige(formula = Asc_logit ~ 1,
                      input_data = as(liberia_sp, "Spatial"),
                      new_data = as(pred_locations, "Spatial"),
                      model = c("Sph", "Exp", "Gau")) 

Asc_pred_prevalence <- plogis(variofit$krige_output$var1.pred)

Asc_prevalence <- rasterFromXYZ(xyz = data.frame(st_coordinates(pred_locations), Asc_pred_prevalence),
                               crs = crs(liberia_sp))
plot(Asc_prevalence)
plot(country$geometry, add = T)
plot(region$geometry, add = T)

writeRaster(Asc_prevalence, "outputs/Asc_pred.tif",
            format = 'GTiff', overwrite = TRUE)


#-Discrete ---------------------------------------------------------------------
# 

# head(liberia)

# region <- region %>% rename(ADMIN1 = NAME_1)

region_data <- region %>% st_join(liberia_sp, join = st_intersects) %>%
  group_by(ADMIN1) %>%
  summarise(HK_positive = sum(HK_positive),
            TT_positive = sum(TT_positive),
            Asc_positive = sum(Asc_positive),
            HK_examined = sum(HK_examined),
            TT_examined = sum(TT_examined),
            Asc_examined = sum(Asc_examined),
            HK_prev = sum(HK_positive)/sum(HK_examined),
            TT_prev = sum(TT_positive)/sum(TT_examined),
            Asc_prev = sum(Asc_positive)/sum(Asc_examined)) 

  # select(ADMIN1, HK_prev, TT_prev, Asc_prev, geom_mean)
st_crs(liberia_sp)

district_data <-  district %>% st_join(liberia_sp, join = st_intersects) %>%
  group_by(ADMIN2) %>%
  summarise(HK_positive = sum(HK_positive),
            TT_positive = sum(TT_positive),
            Asc_positive = sum(Asc_positive),
            HK_examined = sum(HK_examined),
            TT_examined = sum(TT_examined),
            Asc_examined = sum(Asc_examined),
            HK_prev = sum(HK_positive)/sum(HK_examined),
            TT_prev = sum(TT_positive)/sum(TT_examined),
            Asc_prev = sum(Asc_positive)/sum(Asc_examined)) %>% 
   na.omit()



#   
tm_shape(region_data) +
  tm_polygons(col = "HK_prev",
              title = "Hook worm prevalence by region",
              style = "quantile",
              palette = "-Spectral",
              alpha = 0.7,
              lwd = 0.2,
              n = 4,
              border.col = 1,
              legend.hist = TRUE) +
  tm_layout(legend.outside = TRUE)

tm_shape(district_data) +
  tm_polygons(col = "HK_prev",
              title = "Hook worm prevalence by region",
              style = "quantile",
              palette = "-Spectral",
              alpha = 0.7,
              lwd = 0.2,
              n = 4,
              border.col = 1,
              legend.hist = TRUE) +
  tm_layout(legend.outside = TRUE)


# region -----------------------------------------------------------------------

region_neighbour <- poly2nb(region_data)
region_centroids <- region_data %>%
  summarise(geom_mean = st_centroid(geometry))

# region_centroids <- region_centroids %>% 
#   st_set_crs(4326) %>% 

plot(region_data$geometry, lwd = 0.2)
plot(region_centroids,
     col = 'blue', cex = 0.5, lwd = 1, add = TRUE)

# Update the neighbour pattern to a row standardised weights matrix.
weights_queen <- nb2listw(region_neighbour, style = "W", zero.policy = T)
# compare neighbour pattern defined by queen contiguity with...
weights_queen$neighbours

#### Global moran's I
region_moran <- moran.test(region_data$HK_prev, listw = weights_queen, zero.policy = T)
region_moran

# Plot a Moran's Scatter plot using scaled data
region_data$scaled_HKprev <- scale(region_data$HK_prev)
# create a lag vector from the neighbour list and the scaled  values
region_data$lag_scaled_HKprev <- lag.listw(weights_queen,
                                         region_data$scaled_HKprev, zero.policy = T)
# plot the output
ggplot(data = region_data, aes(x = scaled_HKprev, y = lag_scaled_HKprev)) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linetype = "dashed") +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("Raw HK prevalence (scaled)") +
  ylab("HK prevalence lag vector") +
  theme_bw()

# Simulate the ASMR UK data with a random distribution with monte carlo
mc_moran <- moran.mc(region_data$HK_prev, listw = weights_queen,
                     nsim = 9999, zero.policy = T)
# Plot a simple histogram of the data
hist(mc_moran$res, breaks = 50)


#####
# extract the residuals from the MC simulation output above
mc_moran_df <- as.data.frame(mc_moran$res)
# change the column name to something easier to write

colnames(mc_moran_df) <- "mcresiduals"
# plot MC residuals and compare against pur Moran's I estimate
ggplot(mc_moran_df, aes(x = mcresiduals)) +
  geom_histogram(colour = "gray1", fill = "seagreen",
                 alpha = 0.7, binwidth = 0.01) +
  geom_vline(xintercept = region_moran$estimate[1],
             colour = "red", linetype = "dashed") +
  theme_bw()
mc_moran

# Calculate the local Moran statistic for each region using queen contiguity
local_moran <- localmoran(region_data$HK_prev,
                          weights_queen, # our weights object
                          zero.policy = TRUE,
                          alternative = "two.sided")
# Replace the column names to make it clearer
colnames(local_moran) <- c("local_I", "expected_I",
                           "variance_I", "z_statistic", "p_value")
# Summarise the first few rows so you can see what has happened
head(local_moran)

# Join Local Moran's I information back on to the spatial data for plotting
local_moran <- cbind(region_data, local_moran) %>%
  dplyr::select(scaled_HKprev, lag_scaled_HKprev, p_value, geometry)
# plot the data to look for hot and cold spot clusters based on p_values
tm_shape(local_moran) +
  tm_polygons(col = "p_value", title = "Local Moran",
              style = "fixed", breaks = c(0, 0.2, 0.4, 0.6, 0.8,1),
              lwd = 0.2, border.col = 1, palette = "-YlOrRd")

# Join our local moran outputs to the main data and select columns for analysis
local_map <- region_data %>%
  cbind(., local_moran) %>% #
  dplyr::select(scaled_HKprev, lag_scaled_HKprev, p_value, geometry)
# Run the map_maker function
local_map <- map_maker(local_map)
# Plot the number of each cluster type
table(local_map$cluster)

# Plot the data to look for hot and cold spot clusters based on p_values
tm_shape(local_map) +
  tm_polygons(col = "cluster",
              title = "HK Local Moran",
              style = "fixed",
              palette = "Spectral",
              lwd = 0.2,
              border.col = 1)+
  tm_layout(legend.outside = TRUE)
st_write(local_map, "outputs/HK_region.shp")

# District ---------------------------------------------------------------------

district_neighbour <- poly2nb(district_data)
district_centroids <- district_data %>%
  summarise(geom_mean = st_centroid(geometry))

plot(district_data$geometry, lwd = 0.2)
plot(district_neighbour,district_centroids$geom_mean, 
     col = 'blue', cex = 0.5, lwd = 1, add = TRUE)

# Update the neighbour pattern to a row standardised weights matrix.
weights_queen <- nb2listw(district_neighbour, style = "W", zero.policy = T)
# compare neighbour pattern defined by queen contiguity with...
weights_queen$neighbours

#### Global moran's I
district_moran <- moran.test(district_data$HK_prev, listw = weights_queen, zero.policy = T)
district_moran

# Plot a Moran's Scatter plot using scaled data
district_data$scaled_HKprev <- scale(district_data$HK_prev)
# create a lag vector from the neighbour list and the scaled  values
district_data$lag_scaled_HKprev <- lag.listw(weights_queen,
                                             district_data$scaled_HKprev, zero.policy = T)
# plot the output
ggplot(data = district_data, aes(x = scaled_HKprev, y = lag_scaled_HKprev)) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linetype = "dashed") +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("Raw HK prevalence (scaled)") +
  ylab("HK prevalence lag vector") +
  theme_bw()

# Simulate the ASMR UK data with a random distribution with monte carlo
mc_moran <- moran.mc(district_data$HK_prev, listw = weights_queen,
                     nsim = 9999, zero.policy = T)
# Plot a simple histogram of the data
hist(mc_moran$res, breaks = 50)


#####
# extract the residuals from the MC simulation output above
mc_moran_df <- as.data.frame(mc_moran$res)
# change the column name to something easier to write

colnames(mc_moran_df) <- "mcresiduals"
# plot MC residuals and compare against pur Moran's I estimate
ggplot(mc_moran_df, aes(x = mcresiduals)) +
  geom_histogram(colour = "gray1", fill = "seagreen",
                 alpha = 0.7, binwidth = 0.01) +
  geom_vline(xintercept = district_moran$estimate[1],
             colour = "red", linetype = "dashed") +
  theme_bw()
mc_moran

# Calculate the local Moran statistic for each region using queen contiguity
local_moran <- localmoran(district_data$HK_prev,
                          weights_queen, # our weights object
                          zero.policy = TRUE,
                          alternative = "two.sided")
# Replace the column names to make it clearer
colnames(local_moran) <- c("local_I", "expected_I",
                           "variance_I", "z_statistic", "p_value")
# Summarise the first few rows so you can see what has happened
head(local_moran)

# Join Local Moran's I information back on to the spatial data for plotting
local_moran <- cbind(district_data, local_moran) %>%
  dplyr::select(scaled_HKprev, lag_scaled_HKprev, p_value, geometry)
# plot the data to look for hot and cold spot clusters based on p_values
tm_shape(local_moran) +
  tm_polygons(col = "p_value", title = "Local Moran",
              style = "fixed", breaks = c(0, 0.2, 0.4, 0.6, 0.8,1),
              lwd = 0.2, border.col = 1, palette = "-YlOrRd")

# Join our local moran outputs to the main data and select columns for analysis
local_map <- district_data %>%
  cbind(., local_moran) %>% #
  dplyr::select(scaled_HKprev, lag_scaled_HKprev, p_value, geometry)
# Run the map_maker function
local_map <- map_maker(local_map)
# Plot the number of each cluster type
table(local_map$cluster)

# Plot the data to look for hot and cold spot clusters based on p_values
tm_shape(local_map) +
  tm_polygons(col = "cluster",
              title = "HK Local Moran",
              style = "fixed",
              palette = "Spectral",
              lwd = 0.2,
              border.col = 1)+
  tm_layout(legend.outside = TRUE)

st_write(local_map, "outputs/HK_district.shp")
