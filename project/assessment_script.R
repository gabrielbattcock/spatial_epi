# Assessment script
#candidate number 221352

#read in all the libraries
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("r/map_maker.R")
source('3135_packages.R')
source("../spatial_custom_functions.R") ## this is saved locally for me...
library(wesanderson)### you may need to download this, will also change your life
library(tidyr)

#-data exploration-------------------------------------------------------------------
# 
# load in the shape files
country <- st_read("vector/gadm41_LBR_0.shp")
region <-st_read("vector/gadm41_LBR_1.shp")
district <- st_read("vector/gadm41_LBR_2.shp")
level3 <- st_read("vector/gadm41_LBR_3.shp")

head(region)

liberia <- read.csv("data/Liberia_STH.csv")
head(liberia)



# number of pupils per school

summary(liberia$HK_examined)
summary(liberia$TT_examined)
summary(liberia$Asc_examined)

sum(liberia$HK_examined)
sum(liberia$TT_examined)
sum(liberia$Asc_examined)

# liberia %>% count(distinct(ADMIN2))

#
liberia %>% group_by(ADMIN2) %>%
  count(ADMIN2) %>% summary

sum(liberia$HK_examined)
sum(liberia$TT_examined)
sum(liberia$Asc_examined)

par(mfrow = c(1,3))
boxplot(liberia$HK_prevalence, ylim = c(0,1), main = "Hookworms" , ylab = "Prevalence")
boxplot(liberia$Asc_prevalence, ylim = c(0,1), main = "Roundworm", ylab = "Prevalence")
boxplot(liberia$TT_prevalence, ylim = c(0,1), main = "Whipworm", ylab = "Prevalence")

#descriptive stats

#plot hists of the prevelance 
par( mfrow= c(1,3) )
hist(liberia$HK_prevalence, xlab = "Prevalence", main = "Hookworms")
hist(liberia$Asc_prevalence,xlab = "Prevalence",main = "Roundworm")
hist(liberia$TT_prevalence,xlab = "Prevalence",main = "Whipworm")
# data needs to be logged to be normal
liberia$HK_logit <- log((liberia$HK_positive+0.5)/(liberia$HK_examined-liberia$HK_positive+0.5))
liberia$Asc_logit <- log((liberia$Asc_positive+0.5)/(liberia$Asc_examined-liberia$Asc_positive+0.5))
liberia$TT_logit <- log((liberia$TT_positive+0.5)/(liberia$TT_examined-liberia$TT_positive+0.5))

par( mfrow= c(1,3))
hist(liberia$HK_logit, bins = 30,xlab = "logit(prevalence)", main = "Hookworms")
hist(liberia$Asc_logit, bins = 30,xlab = "logit(prevalence)", main = "Roundworm")
hist(liberia$TT_logit, bins = 30, xlab = "logit(prevalence)", main = "Whipworm")

# create a map of the area for each disease
liberia$TT_PREV_BIN <- cut2(liberia$TT_prevalence, cuts = c(0.2, 0.4, 0.6))
liberia$HK_PREV_BIN <- cut2(liberia$HK_prevalence, cuts = c(0.2, 0.4, 0.6))
liberia$Asc_PREV_BIN <- cut2(liberia$Asc_prevalence, cuts = c(0.2, 0.4, 0.6))

#create a map of the area
qmplot(x = Longitude, y = Latitude,
       data = liberia, source = "stamen", maptype = "toner-lite") +
  geom_point(aes(fill = TT_PREV_BIN), shape = 21, alpha = 0.8) +
  scale_fill_manual(values = wes_palette("GrandBudapest1", n = 4))

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

variofitHK <- autofitVariogram(formula = HK_logit ~ 1,
                               input_data = as(liberia_sp, "Spatial"),
                               model = c("Sph", "Exp", "Gau"))
plot(variofitHK)

variofitTT <- autofitVariogram(formula = TT_logit ~ 1,
                               input_data = as(liberia_sp, "Spatial"),
                               model = c("Sph", "Exp", "Gau"))
plot(variofitTT)

variofitAsc <- autofitVariogram(formula = Asc_logit ~ 1,
                               input_data = as(liberia_sp, "Spatial"),
                               model = c("Sph", "Exp", "Gau"))
plot(variofitAsc)



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

level3 <- level3 %>% 
  st_set_crs(4326) %>% 
  st_transform(crs = "+init=epsg:3857 +units=m")


### make a prediction for HK based on the semivariogram


pred_locations <- st_make_grid(country, cellsize = 5000, what = "centers")
pred_locations <- st_intersection(pred_locations, country)#resize the prediction locations to just the area
plot(country$geometry, col = "white")
plot(pred_locations, cex = .2, col = "red", pch = 19, add = T) 


# hook worm 

#ordinary krigging
variofit <- autoKrige(formula = HK_logit ~ 1,
                      input_data = as(liberia_sp, "Spatial"),
                      new_data = as(pred_locations, "Spatial"),
                        model = c("Sph", "Exp", "Gau")) # types of models to try
plot(variofit)
HK_pred_prevalence <- plogis(variofit$krige_output$var1.pred)

variofit$krige_output$var1.var

HK_prevalence <- rasterFromXYZ(xyz = data.frame(st_coordinates(pred_locations), HK_pred_prevalence),
                            crs = crs(liberia_sp))
plot(HK_prevalence)
plot(country$geometry, add = T)
plot(district$geometry, add = T)

writeRaster(HK_prevalence, "outputs/HK_pred.tif",
            format = 'GTiff', overwrite = TRUE)

# TT

variofit <- autoKrige(formula = TT_logit ~ 1,
                      input_data = as(liberia_sp, "Spatial"),
                      new_data = as(pred_locations, "Spatial"),
                      model = c("Sph", "Exp", "Gau")) 
plot(variofit)
TT_pred_prevalence <- plogis(variofit$krige_output$var1.pred)

TT_prevalence <- rasterFromXYZ(xyz = data.frame(st_coordinates(pred_locations), TT_pred_prevalence),
                               crs = crs(liberia_sp))
plot(TT_prevalence)
plot(country$geometry, add = T)
plot(district$geometry, add = T)

writeRaster(TT_prevalence, "outputs/TT_pred.tif",
            format = 'GTiff', overwrite = TRUE)

# Asc

variofit <- autoKrige(formula = Asc_logit ~ 1,
                      input_data = as(liberia_sp, "Spatial"),
                      new_data = as(pred_locations, "Spatial"),
                      model = c("Sph", "Exp", "Gau")) 

Asc_pred_prevalence <- plogis(variofit$krige_output$var1.pred)
asc_prev_points <-  data.frame(variofit$krige_output) %>% 
  st_as_sf(coords = c("coords.x1","coords.x2"), crs = 3857) 

tm_shape(asc_prev_points) +
  tm_dots(col = "var1.var",
              title = "Hook worm prevalence by region",
              style = "quantile",
              palette = "-Spectral",
              alpha = 1,
              size = .2,
              n = 4,
              border.col = 1,
              legend.hist = TRUE) +
  tm_layout(legend.outside = TRUE)


Asc_prevalence <- rasterFromXYZ(xyz = data.frame(st_coordinates(pred_locations), Asc_pred_prevalence),
                               crs = crs(liberia_sp))
plot(Asc_prevalence)
plot(country$geometry, add = T)
plot(district$geometry, add = T)

writeRaster(Asc_prevalence, "outputs/Asc_pred.tif",
            format = 'GTiff', overwrite = TRUE)


#-Discrete ---------------------------------------------------------------------
# 
# doing spatial joins


region_data <- region %>% st_join(liberia_sp) %>%
  group_by(NAME_1) %>%
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
st_write(region_data, "outputs/region.shp", append = F)

district_data <-  district %>% st_join(liberia_sp, join = st_intersects) %>%
  group_by(NAME_2) %>%
  summarise(HK_positive = sum(HK_positive),
            TT_positive = sum(TT_positive),
            Asc_positive = sum(Asc_positive),
            HK_examined = sum(HK_examined),
            TT_examined = sum(TT_examined),
            Asc_examined = sum(Asc_examined),
            HK_prev = sum(HK_positive)/sum(HK_examined),
            TT_prev = sum(TT_positive)/sum(TT_examined),
            Asc_prev = sum(Asc_positive)/sum(Asc_examined)) %>% 
  na.omit()#valid for the analysis otherwise moran's I doesn't work

st_write(district_data, "outputs/district.shp", append = F)

level3_data <- level3 %>%  st_join(liberia_sp, join = st_intersects) %>%
  group_by(NAME_3) %>%
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

# st_write(local_map_HK, "outputs/HK_district.shp", append = F)

summary(region_data)
#   
tm_shape(region_data) +
  tm_polygons(col = "Asc_prev",
              title = "Hook worm prevalence by region",
              style = "quantile",
              palette = "-Spectral",
              alpha = 0.8,
              lwd = 2,
              n = 4,
              border.col = 1,
              legend.hist = TRUE) +
  tm_layout(legend.outside = TRUE)

tm_shape(district_data) +
  tm_polygons(col = "Asc_prev",
              title = "Hook worm prevalence by region",
              style = "quantile",
              breaks = 6,
              palette = "Spectral",
              alpha = 1,
              lwd = 0.2,
              n = 4,
              border.col = 1,
              legend.hist = TRUE) +
  tm_layout(legend.outside = TRUE)

tm_shape(level3_data) +
  tm_polygons(col = "HK_prev",
              title = "Hook worm prevalence by region",
              style = "quantile",
              breaks = 4,
              palette = "Spectral",
              alpha = 1,
              lwd = 0.2,
              n = 4,
              border.col = 1,
              legend.hist = TRUE) +
  tm_layout(legend.outside = TRUE)


# region -----------------------------------------------------------------------

region_neighbour <- poly2nb(region_data)
region_centroids <- region_data %>%
  summarise(geom_mean = st_centroid(geometry))

#check the distribution
ggplot(region_data, aes(x = HK_prev)) +
  geom_histogram(colour = "gray1",
                 fill = "seagreen",
                 alpha = 1,
                 binwidth = 0.02) +
  theme_bw() 

ggplot(region_data, aes(x = TT_prev)) +
  geom_histogram(colour = "gray1",
                 fill = "seagreen",
                 alpha = 1,
                 binwidth = 0.02) +
  theme_bw() 

ggplot(region_data, aes(x = Asc_prev)) +
  geom_histogram(colour = "gray1",
                 fill = "seagreen",
                 alpha = 1,
                 binwidth = 0.02) +
  theme_bw() 

#use the logit scale
region_data$logit_HK <- log((region_data$HK_positive+0.5)/(region_data$HK_examined-region_data$HK_positive+0.5))
region_data$logit_Asc <- log((region_data$Asc_positive+0.5)/(region_data$Asc_examined-region_data$Asc_positive+0.5))
region_data$logit_TT <- log((region_data$TT_positive+0.5)/(region_data$TT_examined-region_data$TT_positive+0.5))
## check to see the prevalence after logit
ggplot(region_data, aes(x = logit_HK)) +
  geom_histogram(colour = "gray1",
                 fill = "seagreen",
                 alpha = 1,
                 binwidth = 1) +
  theme_bw() 
ggplot(region_data, aes(x = logit_Asc)) +
  geom_histogram(colour = "gray1",
                 fill = "seagreen",
                 alpha = 1,
                 binwidth = 1) +
  theme_bw() 
ggplot(region_data, aes(x = logit_TT)) +
  geom_histogram(colour = "gray1",
                 fill = "seagreen",
                 alpha = 1,
                 binwidth = 1) +
  theme_bw() 


#plot the neighbourgood

plot(region_data$geometry, lwd = 0.2)
plot(region_centroids,
     col = 'blue', cex = 0.5, lwd = 1, add = TRUE)

# Update the neighbour pattern to a row standardised weights matrix.
weights_queen_region <- nb2listw(region_neighbour, style = "W", zero.policy = T)
# compare neighbour pattern defined by queen contiguity with...
weights_queen_region$neighbours

#### Global moran's I
region_moran_HK <- moran.test(region_data$logit_HK,  listw = weights_queen_region, zero.policy = T)
region_moran_HK
region_moran_TT <- moran.test(region_data$logit_TT,  listw = weights_queen_region, zero.policy = T)
region_moran_TT
region_moran_Asc <- moran.test(region_data$logit_Asc , listw = weights_queen_region, zero.policy = T)
region_moran_Asc

# Plot a Moran's Scatter plot using scaled data
region_data$scaled_HKprev <- scale(region_data$logit_HK)
region_data$scaled_TTprev <- scale(region_data$logit_TT)
region_data$scaled_Ascprev <- scale(region_data$logit_Asc)
# create a lag vector from the neighbour list and the scaled  values
region_data$lag_scaled_HKprev <- lag.listw(weights_queen_region,
                                         region_data$scaled_HKprev, zero.policy = T)
region_data$lag_scaled_TTprev <- lag.listw(weights_queen_region,
                                           region_data$scaled_TTprev, zero.policy = T)
region_data$lag_scaled_Ascprev <- lag.listw(weights_queen_region,
                                           region_data$scaled_Ascprev, zero.policy = T)
# plot the output
ggplot(data = region_data, aes(x = scaled_HKprev, y = lag_scaled_HKprev)) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linetype = "dashed") +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("Raw HK prevalence (scaled)") +
  ylab("HK prevalence lag vector") +
  theme_bw()

ggplot(data = region_data, aes(x = scaled_TTprev, y = lag_scaled_TTprev)) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linetype = "dashed") +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("Raw TT prevalence (scaled)") +
  ylab("TT prevalence lag vector") +
  theme_bw()

ggplot(data = region_data, aes(x = scaled_Ascprev, y = lag_scaled_Ascprev)) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linetype = "dashed") +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("Raw Asc prevalence (scaled)") +
  ylab("Asc prevalence lag vector") +
  theme_bw()

# Simulate the data with a random distribution with monte carlo
mc_moran_reg_HK <- moran.mc(region_data$logit_HK, listw = weights_queen_region,
                     nsim = 9999, zero.policy = T)
mc_moran_reg_TT <- moran.mc(region_data$logit_TT, listw = weights_queen_region,
                     nsim = 9999, zero.policy = T)
mc_moran_reg_Asc <- moran.mc(region_data$logit_Asc, listw = weights_queen_region,
                     nsim = 9999, zero.policy = T)

# local
# extract the residuals from the MC simulation output above
mc_moran_reg_HK_df <- as.data.frame(mc_moran_reg_HK$res)
mc_moran_reg_TT_df <- as.data.frame(mc_moran_reg_TT$res)
mc_moran_reg_Asc_df <- as.data.frame(mc_moran_reg_Asc$res)

# change the column name to something easier to write

colnames(mc_moran_reg_HK_df) <- "mcresiduals"
colnames(mc_moran_reg_TT_df) <- "mcresiduals"
colnames(mc_moran_reg_Asc_df) <- "mcresiduals"
# plot MC residuals and compare against pur Moran's I estimate
ggplot(mc_moran_reg_HK_df, aes(x = mcresiduals)) +
  geom_histogram(colour = "gray1", fill = "seagreen",
                 alpha = 0.7, binwidth = 0.01) +
  geom_vline(xintercept = region_moran_HK$estimate[1],
             colour = "red", linetype = "dashed") +
  theme_bw()
mc_moran_reg_HK_df

ggplot(mc_moran_reg_TT_df, aes(x = mcresiduals)) +
  geom_histogram(colour = "gray1", fill = "seagreen",
                 alpha = 0.7, binwidth = 0.01) +
  geom_vline(xintercept = region_moran_TT$estimate[1],
             colour = "red", linetype = "dashed") +
  theme_bw()
mc_moran_reg_TT_df

ggplot(mc_moran_reg_Asc_df, aes(x = mcresiduals)) +
  geom_histogram(colour = "gray1", fill = "seagreen",
                 alpha = 0.7, binwidth = 0.01) +
  geom_vline(xintercept = region_moran_Asc$estimate[1],
             colour = "red", linetype = "dashed") +
  theme_bw()
mc_moran_reg_Asc_df

#there is not enough evidence against the null hypothesis here

# Calculate the local Moran statistic for each region using queen contiguity
local_moran_reg_HK <- localmoran(region_data$logit_HK,
                          weights_queen_region, 
                          zero.policy = TRUE,
                          alternative = "two.sided")
local_moran_reg_TT <- localmoran(region_data$logit_TT,
                                 weights_queen_region, 
                          zero.policy = TRUE,
                          alternative = "two.sided")
local_moran_reg_Asc <- localmoran(region_data$logit_Asc,
                                  weights_queen_region,
                          zero.policy = TRUE,
                          alternative = "two.sided")
# Replace the column names to make it clearer
colnames(local_moran_reg_HK) <- c("local_I", "expected_I",
                           "variance_I", "z_statistic", "p_value")
colnames(local_moran_reg_TT) <- c("local_I", "expected_I",
                           "variance_I", "z_statistic", "p_value")
colnames(local_moran_reg_Asc) <- c("local_I", "expected_I",
                           "variance_I", "z_statistic", "p_value")
# Summarise the first few rows so you can see what has happened
head(local_moran_reg_HK)
head(local_moran_reg_TT)
head(local_moran_reg_Asc)
# Join Local Moran's I information back on to the spatial data for plotting
local_moran_reg_HK <- cbind(region_data, local_moran_reg_HK) %>%
  dplyr::select(scaled_HKprev, lag_scaled_HKprev, p_value, geometry)
local_moran_reg_TT <- cbind(region_data, local_moran_reg_TT) %>%
  dplyr::select(scaled_TTprev, lag_scaled_TTprev, p_value, geometry)
local_moran_reg_Asc <- cbind(region_data, local_moran_reg_Asc) %>%
  dplyr::select(scaled_Ascprev, lag_scaled_Ascprev, p_value, geometry)
# plot the data to look for hot and cold spot clusters based on p_values
tm_shape(local_moran_reg_HK) +
  tm_polygons(col = "p_value", title = "Local Moran HK" ,
              style = "fixed", breaks = c(0, 0.2, 0.4, 0.6, 0.8,1),
              lwd = 0.2, border.col = 1, palette = "-YlOrRd")
tm_shape(local_moran_reg_TT) +
  tm_polygons(col = "p_value", title = "Local Moran TT",
              style = "fixed", breaks = c(0, 0.2, 0.4, 0.6, 0.8,1),
              lwd = 0.2, border.col = 1, palette = "-YlOrRd")
tm_shape(local_moran_reg_Asc) +
  tm_polygons(col = "p_value", title = "Local Moran Asc",
              style = "fixed", breaks = c(0, 0.2, 0.4, 0.6, 0.8,1),
              lwd = 0.2, border.col = 1, palette = "-YlOrRd")

# Join our local moran outputs to the main data and select columns for analysis
local_map_reg_HK <- region_data %>%
  cbind(., local_moran_reg_HK) %>% #
  dplyr::select(scaled_HKprev, lag_scaled_HKprev, p_value, geometry)

local_map_reg_TT <- region_data %>%
  cbind(., local_moran_reg_TT) %>% #
  dplyr::select(scaled_TTprev, lag_scaled_TTprev, p_value, geometry)

local_map_reg_Asc <- region_data %>%
  cbind(., local_moran_reg_Asc) %>% #
  dplyr::select(scaled_Ascprev, lag_scaled_Ascprev, p_value, geometry)
# Run the map_maker function
local_map_reg_HK <- map_maker(local_map_reg_HK)
local_map_reg_TT <- map_maker(local_map_reg_TT)
local_map_reg_Asc <- map_maker(local_map_reg_Asc)
# Plot the number of each cluster type
table(local_map_reg_HK$cluster)
table(local_map_reg_TT$cluster)
table(local_map_reg_Asc$cluster)
# Plot the data to look for hot and cold spot clusters based on p_values
tm_shape(local_map_reg_HK) +
  tm_polygons(col = "cluster",
              title = "HK Local Moran",
              style = "fixed",
              palette = "Spectral",
              lwd = 0.2,
              border.col = 1)+
  tm_layout(legend.outside = TRUE)
st_write(local_map_reg_HK, "outputs/HK_region.shp", append = FALSE)

tm_shape(local_map_reg_TT) +
  tm_polygons(col = "cluster",
              title = "TT Local Moran",
              style = "fixed",
              palette = "Spectral",
              lwd = 0.2,
              border.col = 1)+
  tm_layout(legend.outside = TRUE)
st_write(local_map_reg_TT, "outputs/TT_region.shp",append = FALSE)

tm_shape(local_map_reg_Asc) +
  tm_polygons(col = "cluster",
              title = "Asc Local Moran",
              style = "fixed",
              palette = "Spectral",
              lwd = 0.2,
              border.col = 1)+
  tm_layout(legend.outside = TRUE)
st_write(local_map_reg_Asc, "outputs/Asc_region.shp",append = FALSE)

# District ---------------------------------------------------------------------

district_neighbour <- poly2nb(district_data)
district_centroids <- district_data %>%
  summarise(geom_mean = st_centroid(geometry))
## check to see the prevalence 
ggplot(district_data, aes(x = HK_prev)) +
  geom_histogram(colour = "gray1",
                 fill = "seagreen",
                 alpha = 1,
                 binwidth = 0.02) +
  theme_bw() 

ggplot(district_data, aes(x = TT_prev)) +
  geom_histogram(colour = "gray1",
                 fill = "seagreen",
                 alpha = 1,
                 binwidth = 0.02) +
  theme_bw() 

ggplot(district_data, aes(x = Asc_prev)) +
  geom_histogram(colour = "gray1",
                 fill = "seagreen",
                 alpha = 1,
                 binwidth = 0.02) +
  theme_bw() 

#use the logit scale
district_data$logit_HK <- log((district_data$HK_positive+0.5)/(district_data$HK_examined-district_data$HK_positive+0.5))
district_data$logit_Asc <- log((district_data$Asc_positive+0.5)/(district_data$Asc_examined-district_data$Asc_positive+0.5))
district_data$logit_TT <- log((district_data$TT_positive+0.5)/(district_data$TT_examined-district_data$TT_positive+0.5))
## check to see the prevalence after logit
ggplot(district_data, aes(x = logit_HK)) +
  geom_histogram(colour = "gray1",
                 fill = "seagreen",
                 alpha = 1,
                 binwidth = 1) +
  theme_bw() 
ggplot(district_data, aes(x = logit_Asc)) +
  geom_histogram(colour = "gray1",
                 fill = "seagreen",
                 alpha = 1,
                 binwidth = 1) +
  theme_bw() 
ggplot(district_data, aes(x = logit_TT)) +
  geom_histogram(colour = "gray1",
                 fill = "seagreen",
                 alpha = 1,
                 binwidth = 1) +
  theme_bw() 



plot(district_data$geometry, lwd = 0.2)
plot(district_neighbour,district_centroids$geom_mean, 
     col = 'blue', cex = 0.5, lwd = 1, add = TRUE)

# Update the neighbour pattern to a row standardised weights matrix.
weights_queen <- nb2listw(district_neighbour, style = "W", zero.policy = T)
# compare neighbour pattern defined by queen contiguity with...
weights_queen$neighbours

#### Global moran's I
district_moran_HK <- moran.test(district_data$logit_HK, listw = weights_queen, zero.policy = T)
district_moran_HK

district_moran_Asc <- moran.test(district_data$logit_Asc, listw = weights_queen, zero.policy = T)
district_moran_Asc

district_moran_TT <- moran.test(district_data$logit_TT, listw = weights_queen, zero.policy = T)
district_moran_TT

# Plot a Moran's Scatter plot using scaled data
district_data$scaled_HKprev <- scale(district_data$logit_HK)
district_data$scaled_TTprev <- scale(district_data$logit_TT)
district_data$scaled_Ascprev <- scale(district_data$logit_Asc)

# create a lag vector from the neighbour list and the scaled  values
district_data$lag_scaled_HKprev <- lag.listw(weights_queen,
                                             district_data$scaled_HKprev, zero.policy = T)
district_data$lag_scaled_TTprev <- lag.listw(weights_queen,
                                             district_data$scaled_TTprev, zero.policy = T)
district_data$lag_scaled_Ascprev <- lag.listw(weights_queen,
                                             district_data$scaled_Ascprev, zero.policy = T)

# plot the output
ggplot(data = district_data, aes(x = scaled_HKprev, y = lag_scaled_HKprev)) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linetype = "dashed") +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("Raw HK prevalence (scaled)") +
  ylab("HK prevalence lag vector") +
  theme_bw()

ggplot(data = district_data, aes(x = scaled_TTprev, y = lag_scaled_TTprev)) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linetype = "dashed") +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("Raw TT prevalence (scaled)") +
  ylab("TT prevalence lag vector") +
  theme_bw()

ggplot(data = district_data, aes(x = scaled_Ascprev, y = lag_scaled_Ascprev)) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linetype = "dashed") +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("Raw Asc prevalence (scaled)") +
  ylab("Asc prevalence lag vector") +
  theme_bw()

# Simulate data with a random distribution with monte carlo
mc_moran_HK <- moran.mc(district_data$logit_HK, listw = weights_queen,
                     nsim = 9999, zero.policy = T)
mc_moran_TT <- moran.mc(district_data$logit_TT, listw = weights_queen,
                        nsim = 9999, zero.policy = T)
mc_moran_Asc <- moran.mc(district_data$logit_Asc, listw = weights_queen,
                        nsim = 9999, zero.policy = T)

##
# extract the residuals from the MC simulation output above
mc_moran_df_HK <- as.data.frame(mc_moran_HK$res)
mc_moran_df_Asc <- as.data.frame(mc_moran_Asc$res)
mc_moran_df_TT <- as.data.frame(mc_moran_TT$res)

# change the column name to something easier to write
colnames(mc_moran_df_HK) <- "mcresiduals"
colnames(mc_moran_df_Asc) <- "mcresiduals"
colnames(mc_moran_df_TT) <- "mcresiduals"

# plot MC residuals and compare against pur Moran's I estimate
ggplot(mc_moran_df_HK, aes(x = mcresiduals)) +
  geom_histogram(colour = "gray1", fill = "seagreen",
                 alpha = 0.7, binwidth = 0.01) +
  geom_vline(xintercept = district_moran_HK$estimate[1],
             colour = "red", linetype = "dashed") +
  theme_bw()
mc_moran_df_HK

ggplot(mc_moran_df_Asc, aes(x = mcresiduals)) +
  geom_histogram(colour = "gray1", fill = "seagreen",
                 alpha = 0.7, binwidth = 0.01) +
  geom_vline(xintercept = district_moran_Asc$estimate[1],
             colour = "red", linetype = "dashed") +
  theme_bw()
mc_moran_df_Asc

ggplot(mc_moran_df_TT, aes(x = mcresiduals)) +
  geom_histogram(colour = "gray1", fill = "seagreen",
                 alpha = 0.7, binwidth = 0.01) +
  geom_vline(xintercept = district_moran_TT$estimate[1],
             colour = "red", linetype = "dashed") +
  theme_bw()
mc_moran_df_TT
# Calculate the local Moran statistic for each region using queen contiguity
local_moran_HK <- localmoran(district_data$logit_HK,
                          weights_queen, # our weights object
                          zero.policy = TRUE,
                          alternative = "two.sided")

local_moran_TT <- localmoran(district_data$logit_TT,
                             weights_queen, # our weights object
                             zero.policy = TRUE,
                             alternative = "two.sided")

local_moran_Asc <- localmoran(district_data$logit_Asc,
                             weights_queen, # our weights object
                             zero.policy = TRUE,
                             alternative = "two.sided")
# Replace the column names to make it clearer
colnames(local_moran_HK) <- c("local_I", "expected_I",
                           "variance_I", "z_statistic", "p_value")
colnames(local_moran_Asc) <- c("local_I", "expected_I",
                              "variance_I", "z_statistic", "p_value")
colnames(local_moran_TT) <- c("local_I", "expected_I",
                              "variance_I", "z_statistic", "p_value")
# Summarise the first few rows so you can see what has happened
head(local_moran_HK)
head(local_moran_Asc)
head(local_moran_TT)
# Join Local Moran's I information back on to the spatial data for plotting

local_moran_HK <- cbind(district_data, local_moran_HK) %>%
  dplyr::select(scaled_HKprev, lag_scaled_HKprev, p_value, geometry)
local_moran_TT <- cbind(district_data, local_moran_TT) %>%
  dplyr::select(scaled_TTprev, lag_scaled_TTprev, p_value, geometry)
local_moran_Asc <- cbind(district_data, local_moran_Asc) %>%
  dplyr::select(scaled_Ascprev, lag_scaled_Ascprev, p_value, geometry)

# plot the data to look for hot and cold spot clusters based on p_values

tm_shape(local_moran_HK) +
  tm_polygons(col = "p_value", title = "Local Moran - Hookworm",
              style = "fixed", breaks = c(0, 0.2, 0.5,1),
              lwd = 0.2, border.col = 1, palette = "-YlOrRd")

tm_shape(local_moran_TT) +
  tm_polygons(col = "p_value", title = "Local Moran - Whipworm",
              style = "fixed", breaks = c(0, 0.2, 0.5,1),
              lwd = 0.2, border.col = 1, palette = "-YlOrRd")
tm_shape(local_moran_Asc) +
  tm_polygons(col = "p_value", title = "Local Moran - Roundworm",
              style = "fixed", breaks = c(0, 0.2, 0.5, 1),
              lwd = 0.2, border.col = 1, palette = "-YlOrRd")

# Join our local moran outputs to the main data and select columns for analysis
local_map_HK <- district_data %>%
  cbind(., local_moran_HK) %>% #
  dplyr::select(scaled_HKprev, lag_scaled_HKprev, p_value, geometry)

local_map_TT <- district_data %>%
  cbind(., local_moran_TT) %>% #
  dplyr::select(scaled_TTprev, lag_scaled_TTprev, p_value, geometry)

local_map_Asc <- district_data %>%
  cbind(., local_moran_Asc) %>% #
  dplyr::select(scaled_Ascprev, lag_scaled_Ascprev, p_value, geometry)
# Run the map_maker function
local_map_HK <- map_maker(local_map_HK)
local_map_TT <- map_maker(local_map_TT)
local_map_Asc <- map_maker(local_map_Asc)
# Plot the number of each cluster type
table(data.frame(local_map_HK$cluster))
table(local_map_TT$cluster)
table(local_map_Asc$cluster)

# Plot the data to look for hot and cold spot clusters based on p_values
tm_shape(local_map_HK) +
  tm_polygons(col = "cluster",
              title = "Hookworm",
              style = "fixed",
              palette = "Spectral",
              lwd = 0.2,
              border.col = 1)+
  tm_layout(legend.outside = TRUE,
            legend.text.size = 1.2,
            legend.title.size = 1.5)

tm_shape(local_map_Asc) +
  tm_polygons(col = "cluster",
              title = "Roundworm",
              style = "fixed",
              palette = "Spectral",
              lwd = 0.2,
              border.col = 1)+
  tm_layout(legend.outside = TRUE,
            legend.text.size = 1.2,
            legend.title.size = 1.5)

tm_shape(local_map_TT) +
  tm_polygons(col = "cluster",
              title = "Whipworm",
              style = "fixed",
              palette = "Spectral",
              lwd = 0.2,
              border.col = 1)+
  tm_layout(legend.outside = TRUE,
            legend.text.size = 1.2,
            legend.title.size = 1.5)


local_map_HK$cluster <- as.character(local_map_HK$cluster)
local_map_TT$cluster <- as.character(local_map_TT$cluster)
local_map_Asc$cluster <- as.character(local_map_Asc$cluster)



#get names of districts with high-high
names_HK <- local_map_HK  %>% st_centroid %>% 
  st_transform(crs = "+init=epsg:3857 +units=m") %>%
  st_join(district_data) %>% 
  filter(cluster == "high-high") %>% 
  group_by(NAME_2) %>% 
  slice(1)%>% 
  dplyr::select(NAME_2, cluster, HK_prev)



names_Asc <- local_map_Asc  %>%  st_centroid %>% 
  st_transform(crs = "+init=epsg:3857 +units=m") %>%
  st_join(district_data) %>% 
  filter(cluster == "high-high") %>% 
  group_by(NAME_2) %>% 
  slice(1) %>% 
  dplyr::select(NAME_2, cluster, Asc_prev)


names_TT <- local_map_TT %>% st_centroid %>% 
  st_transform(crs = "+init=epsg:3857 +units=m") %>%
  st_join(district_data) %>% 
  filter(cluster == "high-high") %>% 
  group_by(NAME_2) %>% 
  slice(1)%>% 
  dplyr::select(NAME_2, cluster, TT_prev)



st_write(names_Asc, "outputs/names_Asc.shp", append = F)
st_write(names_HK, "outputs/names_HK.shp", append = F)
st_write(names_TT, "outputs/names_TT.shp", append = F)

st_write(local_map_HK, "outputs/HK_district.shp", append = F)
st_write(local_map_TT, "outputs/TT_district.shp", append = F)
st_write(local_map_Asc, "outputs/Asc_district.shp", append = F)

#--------------------------- extension with different boundaries ----


level3_neighbour <- poly2nb(level3_data)
level3_centroids <- level3_data %>%
  summarise(geom_mean = st_centroid(geometry))
## check to see the prevalence 
ggplot(level3_data, aes(x = HK_prev)) +
  geom_histogram(colour = "gray1",
                 fill = "seagreen",
                 alpha = 1,
                 binwidth = 0.02) +
  theme_bw() 

ggplot(level3_data, aes(x = TT_prev)) +
  geom_histogram(colour = "gray1",
                 fill = "seagreen",
                 alpha = 1,
                 binwidth = 0.02) +
  theme_bw() 

ggplot(level3_data, aes(x = Asc_prev)) +
  geom_histogram(colour = "gray1",
                 fill = "seagreen",
                 alpha = 1,
                 binwidth = 0.02) +
  theme_bw() 

#use the logit scale
level3_data$logit_HK <- log((level3_data$HK_positive+0.5)/(level3_data$HK_examined-level3_data$HK_positive+0.5))
level3_data$logit_Asc <- log((level3_data$Asc_positive+0.5)/(level3_data$Asc_examined-level3_data$Asc_positive+0.5))
level3_data$logit_TT <- log((level3_data$TT_positive+0.5)/(level3_data$TT_examined-level3_data$TT_positive+0.5))
## check to see the prevalence after logit
ggplot(level3_data, aes(x = logit_HK)) +
  geom_histogram(colour = "gray1",
                 fill = "seagreen",
                 alpha = 1,
                 binwidth = 1) +
  theme_bw() 
ggplot(level3_data, aes(x = logit_Asc)) +
  geom_histogram(colour = "gray1",
                 fill = "seagreen",
                 alpha = 1,
                 binwidth = 1) +
  theme_bw() 
ggplot(level3_data, aes(x = logit_TT)) +
  geom_histogram(colour = "gray1",
                 fill = "seagreen",
                 alpha = 1,
                 binwidth = 1) +
  theme_bw() 



plot(level3_data$geometry, lwd = 0.2)
plot(level3_neighbour,level3_centroids$geom_mean, 
     col = 'blue', cex = 0.5, lwd = 1, add = TRUE)

# Update the neighbour pattern to a row standardised weights matrix.
weights_queen_l3 <- nb2listw(level3_neighbour, style = "W", zero.policy = T)
# compare neighbour pattern defined by queen contiguity with...
weights_queen_l3$neighbours

#### Global moran's I
level3_moran_HK <- moran.test(level3_data$logit_HK, listw = weights_queen_l3, zero.policy = T)
level3_moran_HK

level3_moran_Asc <- moran.test(level3_data$logit_Asc, listw = weights_queen_l3, zero.policy = T)
level3_moran_Asc

level3_moran_TT <- moran.test(level3_data$logit_TT, listw = weights_queen_l3, zero.policy = T)
level3_moran_TT

# Plot a Moran's Scatter plot using scaled data
level3_data$scaled_HKprev <- scale(level3_data$logit_HK)
level3_data$scaled_TTprev <- scale(level3_data$logit_TT)
level3_data$scaled_Ascprev <- scale(level3_data$logit_Asc)

# create a lag vector from the neighbour list and the scaled  values
level3_data$lag_scaled_HKprev <- lag.listw(weights_queen_l3,
                                             level3_data$scaled_HKprev, zero.policy = T)
level3_data$lag_scaled_TTprev <- lag.listw(weights_queen_l3,
                                           level3_data$scaled_TTprev, zero.policy = T)
level3_data$lag_scaled_Ascprev <- lag.listw(weights_queen_l3,
                                            level3_data$scaled_Ascprev, zero.policy = T)

# plot the output
ggplot(data = level3_data, aes(x = scaled_HKprev, y = lag_scaled_HKprev)) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linetype = "dashed") +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("Raw HK prevalence (scaled)") +
  ylab("HK prevalence lag vector") +
  theme_bw()

ggplot(data = level3_data, aes(x = scaled_TTprev, y = lag_scaled_TTprev)) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linetype = "dashed") +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("Raw TT prevalence (scaled)") +
  ylab("TT prevalence lag vector") +
  theme_bw()

ggplot(data = level3_data, aes(x = scaled_Ascprev, y = lag_scaled_Ascprev)) +
  geom_smooth(method = "lm", se = FALSE, colour = "red", linetype = "dashed") +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  xlab("Raw Asc prevalence (scaled)") +
  ylab("Asc prevalence lag vector") +
  theme_bw()

# Simulate data with a random distribution with monte carlo
mc_moran_HK_l3 <- moran.mc(level3_data$logit_HK, listw = weights_queen_l3,
                        nsim = 9999, zero.policy = T)
mc_moran_TT_l3 <- moran.mc(level3_data$logit_TT, listw = weights_queen_l3,
                        nsim = 9999, zero.policy = T)
mc_moran_Asc_l3<- moran.mc(level3_data$logit_Asc, listw = weights_queen_l3,
                         nsim = 9999, zero.policy = T)

##
# extract the residuals from the MC simulation output above
mc_moran_df_HK_l3 <- as.data.frame(mc_moran_HK_l3$res)
mc_moran_df_Asc_l3 <- as.data.frame(mc_moran_Asc_l3$res)
mc_moran_df_TT_l3 <- as.data.frame(mc_moran_TT_l3$res)

# change the column name to something easier to write
colnames(mc_moran_df_HK_l3) <- "mcresiduals"
colnames(mc_moran_df_Asc_l3) <- "mcresiduals"
colnames(mc_moran_df_TT_l3) <- "mcresiduals"

# plot MC residuals and compare against pur Moran's I estimate
ggplot(mc_moran_df_HK_l3, aes(x = mcresiduals)) +
  geom_histogram(colour = "gray1", fill = "seagreen",
                 alpha = 0.7, binwidth = 0.01) +
  geom_vline(xintercept = level3_moran_HK$estimate[1],
             colour = "red", linetype = "dashed") +
  theme_bw()
mc_moran_df_HK_l3

ggplot(mc_moran_df_Asc_l3, aes(x = mcresiduals)) +
  geom_histogram(colour = "gray1", fill = "seagreen",
                 alpha = 0.7, binwidth = 0.01) +
  geom_vline(xintercept = level3_moran_Asc$estimate[1],
             colour = "red", linetype = "dashed") +
  theme_bw()
mc_moran_df_Asc_l3

ggplot(mc_moran_df_TT_l3, aes(x = mcresiduals)) +
  geom_histogram(colour = "gray1", fill = "seagreen",
                 alpha = 0.7, binwidth = 0.01) +
  geom_vline(xintercept = level3_moran_TT$estimate[1],
             colour = "red", linetype = "dashed") +
  theme_bw()
mc_moran_df_TT_l3
# Calculate the local Moran statistic for each region using queen contiguity
local_moran_HK_l3 <- localmoran(level3_data$logit_HK,
                             weights_queen_l3, # our weights object
                             zero.policy = TRUE,
                             alternative = "two.sided")

local_moran_TT_l3 <- localmoran(level3_data$logit_TT,
                                weights_queen_l3, # our weights object
                             zero.policy = TRUE,
                             alternative = "two.sided")

local_moran_Asc_l3 <- localmoran(level3_data$logit_Asc,
                                weights_queen_l3, # our weights object
                              zero.policy = TRUE,
                              alternative = "two.sided")
# Replace the column names to make it clearer
colnames(local_moran_HK_l3) <- c("local_I", "expected_I",
                              "variance_I", "z_statistic", "p_value")
colnames(local_moran_TT_l3) <- c("local_I", "expected_I",
                               "variance_I", "z_statistic", "p_value")
colnames(local_moran_Asc_l3) <- c("local_I", "expected_I",
                              "variance_I", "z_statistic", "p_value")
# Summarise the first few rows so you can see what has happened
head(local_moran_HK_l3)
head(local_moran_TT_l3)
head(local_moran_Asc_l3)
# Join Local Moran's I information back on to the spatial data for plotting

local_moran_HK_l3 <- cbind(level3_data, local_moran_HK_l3) %>%
  dplyr::select(scaled_HKprev, lag_scaled_HKprev, p_value, geometry)
local_moran_TT_l3 <- cbind(level3_data, local_moran_TT_l3) %>%
  dplyr::select(scaled_TTprev, lag_scaled_TTprev, p_value, geometry)
local_moran_Asc_l3 <- cbind(level3_data, local_moran_Asc_l3) %>%
  dplyr::select(scaled_Ascprev, lag_scaled_Ascprev, p_value, geometry)

# plot the data to look for hot and cold spot clusters based on p_values

tm_shape(local_moran_HK_l3) +
  tm_polygons(col = "p_value", title = "Local Moran - Hookworm",
              style = "fixed", breaks = c(0, 0.05, 0.5,1),
              lwd = 0.2, border.col = 1, palette = "-YlOrRd")

tm_shape(local_moran_TT_l3) +
  tm_polygons(col = "p_value", title = "Local Moran - Whipworm",
              style = "fixed", breaks = c(0, 0.05, 0.5,1),
              lwd = 0.2, border.col = 1, palette = "-YlOrRd")
tm_shape(local_moran_Asc_l3) +
  tm_polygons(col = "p_value", title = "Local Moran - Roundworm",
              style = "fixed", breaks = c(0, 0.05, 0.5, 1),
              lwd = 0.2, border.col = 1, palette = "-YlOrRd")

# Join our local moran outputs to the main data and select columns for analysis
local_map_HK_l3 <- level3_data %>%
  cbind(., local_moran_HK_l3) %>% #
  dplyr::select(scaled_HKprev, lag_scaled_HKprev, p_value, geometry)

local_map_TT_l3 <- level3_data %>%
  cbind(., local_moran_TT_l3) %>% #
  dplyr::select(scaled_TTprev, lag_scaled_TTprev, p_value, geometry)

local_map_Asc_l3 <- level3_data %>%
  cbind(., local_moran_Asc_l3) %>% #
  dplyr::select(scaled_Ascprev, lag_scaled_Ascprev, p_value, geometry)
# Run the map_maker function
local_map_HK_l3 <- map_maker(local_map_HK_l3)
local_map_TT_l3 <- map_maker(local_map_TT_l3)
local_map_Asc_l3 <- map_maker(local_map_Asc_l3)
# Plot the number of each cluster type
table(local_map_HK_l3$cluster)
table(local_map_TT_l3$cluster)
table(local_map_Asc_l3$cluster)

# Plot the data to look for hot and cold spot clusters based on p_values
tm_shape(local_map_HK_l3) +
  tm_polygons(col = "cluster",
              title = "Hookworm",
              style = "fixed",
              palette = "Spectral",
              lwd = 0.2,
              border.col = 1)+
  tm_layout(legend.outside = TRUE,
            legend.text.size = 1.2,
            legend.title.size = 1.5)

tm_shape(local_map_Asc_l3) +
  tm_polygons(col = "cluster",
              title = "Roundworm",
              style = "fixed",
              palette = "Spectral",
              lwd = 0.2,
              border.col = 1)+
  tm_layout(legend.outside = TRUE,
            legend.text.size = 1.2,
            legend.title.size = 1.5)

tm_shape(local_map_TT_l3) +
  tm_polygons(col = "cluster",
              title = "Whipworm",
              style = "fixed",
              palette = "Spectral",
              lwd = 0.2,
              border.col = 1)+
  tm_layout(legend.outside = TRUE,
            legend.text.size = 1.2,
            legend.title.size = 1.5)

local_map_HK_l3$cluster <- as.character(local_map_HK_l3$cluster)
local_map_TT_l3$cluster <- as.character(local_map_TT_l3$cluster)
local_map_Asc_l3$cluster <- as.character(local_map_Asc_l3$cluster)


st_write(local_map_HK_l3, "outputs/HK_level3.shp", append = F)
st_write(local_map_TT_l3, "outputs/TT_level3.shp", append = F)
st_write(local_map_Asc_l3, "outputs/Asc_level3.shp", append = F)



