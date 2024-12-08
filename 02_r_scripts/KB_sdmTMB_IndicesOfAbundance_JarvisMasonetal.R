# JarvisMason et al.
# Generate standardized index of larval abundance for Kelp Bass
# 2024


# load libraries ----------------------------------------------------------

library(dplyr)
library(sdmTMB)
library(sdmTMBextra)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthhires)
library(MASS)
library(terra)
library(sf)
library(marmap)
library(future)
future::plan(multisession)
library("png")                      # Load png package
library(patchwork)                  # Load patchwork package

pcla <- readPNG("03_plots/kb.png", native = TRUE)



## function to run sdmTMB models using different model specs

check_mods <- function(sp, sptmp, y, covar, time, mesh, dataset, xt, tv, reml) {
  fit <- sdmTMB(reformulate(
    termlabels = covar,
    response = y
  ),
  mesh = mesh,
  time = time,
  spatial = sp, # can be on / off. On estimates a shared spatial field that is the same across years
  spatiotemporal = sptmp, # should be iid / off
  time_varying = tv,
  extra_time = xt,
  family = tweedie(), # could also be one of the delta_ families
  data = dataset,
  reml = reml,
  silent = FALSE
  )
  san <- sanity(fit)
  AIC <- AIC(fit)
  tidy <- tidy(fit,"ran_pars",conf.int = TRUE)
  
  print(list(fit,san,AIC,tidy))
}

## function to create pub quality figs

theme_Publication <- function(base_size=20, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(0.5),vjust = 0), #hjust = 0.5),
            # text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border =  element_blank(),#element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major.y = element_line(colour="#f0f0f0"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            legend.text = element_text(size = rel(0.5)),
            plot.margin=unit(c(5,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}



# 01 load data ------------------------------------------------

load(file = "01_tidy_data/SCB_paralabrax.Rdata") 

names(df.intrp)
kb <- df.intrp
kb$moon.facnum <- as.numeric(kb$moon.fac)
kb$scale_day <- scale(kb$day)
kb$scale_dist2shore <- scale(kb$dist_km_calc)
kb$scale_moon.num <- scale(kb$moon.num)



## set seed------------

set.seed(022824)

# 02 make a mesh -------------------------------------------------------------

# make a mesh -- which is approximating the spatial surface
# note that the data frame 'bsb' already contains lat and lon in UTMs -- which is what
# we want to use

mesh = make_mesh(kb, xy_cols = c("Lon.km","Lat.km"), cutoff = 8, type = "kmeans")

# smaller cutoff values lead to more complexity -- and we can check the number of knots,
# or random effects being estimated with the following code
mesh$mesh$n
plot(mesh)

# account for any overlap of mesh with land

map_data <- rnaturalearth::ne_countries(
  scale = "large",
  returnclass = "sf", country = "United States of America")

# Crop the polygon for plotting and efficiency
st_bbox(map_data) # find the rough coordinates

na_coast <- suppressWarnings(suppressMessages(
  sf::st_crop(map_data,
              c(xmin = -122, ymin = 30, xmax = -116, ymax = 36))))

na_coast_proj <- sf::st_transform(na_coast, crs = 32611)

# Add on the barrier mesh component:
# Note that the original mesh was created using km coordinates,
# thus, the proj_scaling needs to be x1000 to convert to meters

bspde <- add_barrier_mesh(
  mesh, na_coast_proj, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)

# In the above, the grey dots are the centre of triangles that are in the
# ocean. The red crosses are centres of triangles that are over land. The
# spatial range will be assumed to be 0.1 (`range_fraction`) over land compared
# to over water.


# 03 fit models with spatiotemporal and spatial random effects on, but having different fixed effects ----------------------------------------
# NOTE: use reml = FALSE (here, AIC based on fixed effects, random effects must stay the same across model comparisons)


# test the effect of distance to shore, call it m1a
# specify some function arguments

names(kb)

# extra time argument (no data these years)
xt <- c(1967,1968,1970,1971,1973,1974,1976,1977,1979:1983,1985:1988,1990:1991,
        1993:1994,1996,1999,2008,2010:2011)


covar <- c("s(year)") # make sure your variables are named correctly; scale covariates where necessary
m1a <- check_mods("on", "iid", "KB", covar, "year", bspde, kb, xt,tv=NULL,reml=FALSE)

# now try adding in moon phase, call it m1b

names(kb)
covar <- c("s(year)","s(moon.facnum,k=3)")
m1b <- check_mods("on", "iid", "KB", covar, "year", bspde, kb, xt,tv=NULL,reml=FALSE)

covar <- c("s(year)","scale_moon.num") 
m1b2 <- check_mods("on", "iid", "KB", covar, "year", bspde, kb, xt,tv=NULL,reml=FALSE)

covar <- c("s(year)","s(scale_moon.num)")
m1b3 <- check_mods("on", "iid", "KB", covar, "year", bspde, kb, xt,tv=NULL,reml=FALSE)

covar <- c("s(year)","s(scale_moon.num,bs='cc',k=4)")
m1b4 <- check_mods("on", "iid", "KB", covar, "year", bspde, kb, xt,tv=NULL,reml=FALSE)



# now try adding in hour, call it m1c

covar <- c("s(year)","s(hour, bs ='cc',k = 24)")   # use single quotes to specify basis function
covar <- c("s(year)","s(hour, bs =\"cc\",k = 24)") # another way that works
m1c <- check_mods("on", "iid", "KB", covar, "year", bspde, kb, xt,tv=NULL,reml=FALSE)

# now try adding in day_night by net type, call it m1d

covar <- c("s(year)","s(day_night, by = net_type, bs = 'fs')")
m1d <- check_mods("on", "iid", "KB", covar, "year", bspde, kb, xt,tv=NULL,reml=FALSE)

covar <- c("s(year)","s(net_type, bs = 'fs')")
m1d2 <- check_mods("on", "iid", "KB", covar, "year", bspde, kb, xt,tv=formula(~1),reml=FALSE)

covar <- c("s(year)","net_type")
m1d3 <- check_mods("on", "iid", "KB", covar, "year", bspde, kb, xt,tv=formula(~1),reml=FALSE)

# now try adding in day, call it m1e

covar <- c("s(year)","scale_day")
m1e <- check_mods("on", "iid", "KB", covar, "year", bspde, kb, xt,tv=NULL,reml=FALSE)

covar <- c("s(year)", "s(scale_day)")
m1e2 <- check_mods("on", "iid", "KB", covar, "year", bspde, kb, xt,tv=NULL,reml=FALSE)

covar <- c("s(year)", "scale_day", "s(scale_moon.num,bs='cc',k=4)")
m2e1 <- check_mods("on", "iid", "KB", covar, "year", bspde, kb, xt,tv=NULL,reml=FALSE)

# now try adding in distance to shore, call it m1f

covar <- c("s(year)","scale_dist2shore")
m1f <- check_mods("on", "iid", "KB", covar, "year", bspde, kb, xt,tv=NULL,reml=FALSE)

covar <- c("s(year)","scale_dist2shore","scale_moon.num")
m1f3 <- check_mods("on", "iid", "KB", covar, "year", bspde, kb, xt,tv=NULL,reml=FALSE)

covar <- c("s(year)", "scale_dist2shore", "scale_day", "s(scale_moon.num,bs='cc',k=4)")
m1f2 <- check_mods("on", "iid", "KB", covar, "year", bspde, kb, xt,tv=NULL,reml=FALSE)

covar <- c("s(year)", "scale_dist2shore", "scale_day")
m2f1 <- check_mods("on", "iid", "KB", covar, "year", bspde, kb, xt,tv=NULL,reml=FALSE)

covar <- c("s(year)", "scale_dist2shore", "scale_day", "moon.num")
m2f2 <- check_mods("on", "iid", "KB", covar, "year", bspde, kb, xt,tv=NULL,reml=FALSE)

# 04 check model fit----------------------------------------------

# sanity checks

c(
  m1a[2],   #1  
  m1b[2],   #2  
  m1b2[2],  #3                                  
  m1b3[2],  #4  
  m1b4[2],  #5  
  m1c[2],   #6  
  m1d[2],   #7  
  m1d2[2],  #8  
  m1d3[2],  #9  
  m1e[2],   #10                               
  m1e2[2],  #11 
  m2e1[2],  #12                  
  m1f[2],   #13                        
  m1f2[2],  #14 
  m1f3[2],  #15 
  m2f1[2],  #16         
  m2f2[2]   #17 
)
# mods 1,3,5,10,12,14 passed sanity checks

# check parsimony
# AIC
as.data.frame(tibble(
  model = c(
  # c("s(year),scale_moon.num"),
  c("s(year),s(scale_moon.num,bs='cc',k=4)"), # between this and the one above, this moon model is slightly better
  c("s(year),scale_day"),
  c("s(year),scale_day,s(scale_moon.num,bs='cc',k=4)"),
  c("s(year),scale_dist2shore,scale_day,s(scale_moon.num,bs='cc',k=4)")
),
 AIC = c(
   # unlist(m1b2[3]),   
   unlist(m1b4[3]),   
   unlist(m1e[3]),                                      
   unlist(m2e1[3]),   
   unlist(m1f2[3])
 ))) %>% 
  arrange(AIC)

# m1f2 and m1e most parsimonious

# check if random effects structure with spatial off fits better 
# NOTE: use reml = TRUE (here, AIC based on random effects)

covar <- c("s(year)", "scale_dist2shore", "scale_day", "s(scale_moon.num,bs='cc',k=4)")
m1f2.spOff <- check_mods("off", "iid", y="BSB", covar, "year", bspde, kb, xt,tv=NULL,reml=TRUE)
covar <- c("s(year)", "scale_day")
m1e.spOff <- check_mods("off", "iid", y="BSB", covar, "year", bspde, kb, xt,tv=NULL,reml=TRUE)


## cross validation --------------------------------------------------------

## function to run models for cross validation (includes args for clust & k_folds) 
check_cv_mods <- function(sp, sptmp, y, covar, time, mesh, dataset, tv, clust,k_folds) {
  fit <- sdmTMB_cv(reformulate(
    termlabels = covar,
    response = y
  ),
  mesh = mesh,
  time = time,
  spatial = sp, 
  spatiotemporal = sptmp, 
  time_varying = tv,
  family = tweedie(), 
  data = dataset,
  fold_ids = clust,
  k_folds = k_folds,
  parallel = TRUE
  )
}

set.seed(40) # must run set.seed before assigning clust to ensure reproducibility
clust <- sample(1:8, size = nrow(kb), replace = T)
k_folds <- length(unique(clust))

# note that m1f2.cv does not complete
# covar <- c("s(year)", "scale_dist2shore", "scale_day", "s(scale_moon.num,bs='cc',k=4)")
# m1f2.cv <- check_cv_mods("on", "iid", "KB", covar, "year", bspde, kb, tv=NULL,clust,k_folds)

covar <- c("s(year)","scale_day")
m1e.cv <- check_cv_mods("on", "iid", "KB", covar, "year", bspde, kb, tv=NULL,clust,k_folds)

covar <- covar <- c("s(year)","s(scale_moon.num,bs='cc',k=4)")
m1b4.cv <- check_cv_mods("on", "iid", "KB", covar, "year", bspde, kb, tv=NULL,clust,k_folds)

covar <- c("s(year)","scale_day", "s(scale_moon.num,bs='cc',k=4)")    
m2e1.cv <- check_cv_mods("on", "iid", "KB", covar, "year", bspde, kb, tv=NULL,clust,k_folds)

# 8-fold cv - did the models converge?
c(#m2f2.cv$converged,  # no, model run couldn't complete
  m1e.cv$converged,   # "s(year)","scale_day"
  m1b4.cv$converged,  # "s(year)","s(scale_moon.num,bs='cc',k=4)" 
  m2e1.cv$converged)  # no 

# compare models
as.data.frame(tibble(
  model = c(
    c("s(year),scale_day"),
    c("s(year),s(scale_moon.num,bs='cc',k=4")
  ),
  AIC = c(
    unlist(m1e[3]),  
    unlist(m1b4[3])
  ),
  sum_LL = c(
    m1e.cv$sum_loglik,   
    m1b4.cv$sum_loglik
  )
)) %>% 
  arrange(AIC)


# 05 save best predictive mod as determined by cross-validation----------------------------------

# best mod
fit <- sdmTMB(KB ~ s(year) + scale_day,
              mesh = bspde,
              time = "year", 
              spatial = "on", 
              spatiotemporal = "iid", 
              extra_time = c(1967,1968,1970,1971,1973,1974,1976,1977,1979:1980,1982, 1985:1988,1990:1991,
                             1993:1994,1996,1999,2008,2010:2011),
              family = tweedie(), # could also be one of the delta_ families
              data = kb,
              reml = FALSE,
              silent = FALSE)


tidy(fit,"fixed",conf.int = TRUE)


# 06 create prediction grid -----------------------------

# create a polygon based off the mesh boundaries
new_poly <- st_convex_hull(st_union(bspde$mesh_sf))
plot(new_poly)
st_crs(new_poly)

# remove land from polygon
poly <- st_difference(new_poly, na_coast_proj)
st_crs(poly) 
plot(poly)
st_bbox(poly)


# create a grid within the polygon using the raster package
resolution <- 2000 # poly is in meters; 2000 is pretty fine: 2000m x 2000m = 4 km2
r <- raster::raster(as(poly, "Spatial"), resolution = resolution)
rr <- raster::rasterize(as(poly, "Spatial"), r, getCover = TRUE)
plot(rr)

# convert to data frame and calculate grid cell area for each row
grid <- as.data.frame(raster::rasterToPoints(rr))
grid$area <- grid$layer * resolution * resolution # meters squared
grid <- dplyr::filter(grid, area > 0) |> 
  dplyr::select(-layer)

grid$Lon.km <- grid$x / 1000
grid$Lat.km <- grid$y / 1000
grid$area.km <- grid$area / 1e06


# plot to see a display of grid cell areas coded by color
ggplot(grid, aes(Lon.km, Lat.km, colour = area.km)) +
  geom_tile(width = 10, height = 10, fill = NA) + # in meters
  scale_colour_viridis_c(direction = -1) +
  geom_point(size = 0.5) +
  coord_fixed()


# convert to a shapefile using the appropriate crs (we want UTMs) - this is so
# we can fill in the corresponding missing distance data using a distance to shore raster (could do the same for sea floor depth)
grid_sdf <- st_as_sf(grid, coords = c("x", "y"), crs = "epsg:32611")


map_data <- rnaturalearth::ne_coastline(
  scale = "large",
  returnclass = "sf")#, country = "United States of America")
# Crop the polygon for plotting and efficiency:
st_bbox(map_data) # find the rough coordinates

na_coast <- suppressWarnings(suppressMessages(
  st_crop(map_data,
          c(xmin = -122, ymin = 32, xmax = -116, ymax = 35.5))
)) %>% 
  # slice(1:(n() - 9)) # exclude islands
  filter(scalerank == 0) # exclude islands

na_coast_proj <- sf::st_transform(na_coast, crs = 32611)

# check it
ggplot() +
  geom_sf(data = na_coast_proj) 

# we need to set the Coordinate Reference System to be the same
st_crs(grid_sdf) <- st_crs(na_coast_proj)

# map coast
# Plot our coast and survey data:
ggplot() +
  geom_sf(data = na_coast_proj) +
  geom_sf(data = grid_sdf, size = 0.5) +
  coord_sf() +
  theme_minimal()

points_sf <- sample_n(grid_sdf, nrow(grid_sdf))


# get only the lines that are the coast (no islands),
# and combine into one multilinestring
polygon_sf <- 
  na_coast_proj 

# very fast!
dist <- st_distance(points_sf, 
                    polygon_sf)

# convert the distance to numeric
points_sf_dist <- 
  points_sf %>% 
  mutate(dist = as.numeric(dist) / 1000)

# get lines to show distances
dist_lines <- st_nearest_points(points_sf_dist, polygon_sf)

# visual check that distance-to-coast is 
# computed as expected
ggplot() +
  geom_sf(data = points_sf_dist,
          aes(colour = dist)) +
  geom_sf(data = polygon_sf)  +
  geom_sf(data = dist_lines,
          aes(colour = points_sf_dist$dist)) +
  scale_color_viridis_c() +
  coord_sf() +
  theme_minimal() 


df <- left_join(grid,points_sf_dist)

# now expand our data frame to include all time variables and any covariates in our model (year, Julian day,
# lat/lon in kilometers) - this is a necessary step for prediction

# replication factor  
yrs <- sort(as.vector(unique(kb$year))) 
n <- length(yrs)

# replicate our dataframe for year in our model, call our df, 'new_df'
new_df <- do.call("rbind", replicate( 
  n, df, simplify = FALSE)) 

# add year for designating time
new_df$year_f <-  as.factor(rep(yrs, each = nrow(df))) # if a factor in model
new_df$year <-  rep(yrs, each = nrow(df))


# add day column (choose a day to predict on) 
new_df$day <- 196         # July 15th


# grab mean and stdev of day column in kb 
# these will be used to scale the prediction grid distances

mu.yday <- mean(kb$day)
sd.yday <- sd(kb$day)
new_df <- new_df %>%
  mutate(scale_day = (day - mu.yday) / sd.yday) %>% # use same name as covariate name in model
  filter(area.km>3.9999) # we want 4 km sq grid cells, or one could create an area vector to predict on all grid cells, but in this case, the vec was too large for making predictions 


# 07 calculate index of abundance----------------------------------------------------

# make predictions
pred.index = predict(fit, newdata = new_df, return_tmb_object = TRUE)

# now make index
index = get_index(pred.index, bias_correct = TRUE, area = 4) # area_vec doesn't work??

save(kb,df,new_df,fit,pred.index,index, file="01_tidy_data/output_SIA_KB.RData")

## index with missing years (skip this section to move on to plot index above)---------------

# replication factor  
yrs <- seq(1963,2016,1) 
n <- length(yrs)

# replicate our dataframe for year in our model, call our df, 'new_df'
new_df <- do.call("rbind", replicate( 
  n, df, simplify = FALSE)) 

# add year for designating time
new_df$year_f <-  as.factor(rep(yrs, each = nrow(df))) # if a factor in model
new_df$year <-  rep(yrs, each = nrow(df))

# add day column 
new_df$day <- 196


# grab mean and stdev of day column in kb 
# these will be used to scale the prediction grid distances
mu.yday <- mean(kb$day)
sd.yday <- sd(kb$day)

new_df <- new_df %>% 
  dplyr::filter(dist > 0) %>% 
  mutate(scale_day = (day - mu.yday) / sd.yday) %>% # use same name as covariate name in model
  filter(area.km>3.9999)

# make predictions
pred.index = predict(fit, newdata = new_df, return_tmb_object = TRUE)

# now make index
index = get_index(pred.index, bias_correct = TRUE, area = 4) 

save(kb,df,new_df,fit,pred.index,index, file="01_tidy_data/output_SIA_KB_with_missing_yrs.RData")


# 08 Figure 2 plot temporal trends ---------------------------

load(file="01_tidy_data/output_SIA_KB.RData")
load(file="01_tidy_data/output_SIA_KB_with_missing_yrs.RData")


# on linear scale (no transformation)
p4 <- ggplot(index, aes(year, est)) + 
  geom_ribbon(aes(ymin=lwr,ymax=upr), fill = "gray",alpha=0.6) +
  # geom_pointrange(aes(ymin=lwr,ymax=upr))+
  geom_line(linewidth = 0.5, color = "grey30")+
  xlim(1960,2020)+
  ylim(c(0,1000000))+
  xlab("Year") + ylab("Abundance")+
  theme_Publication()  +
  labs(tag = 'C')
p4

# natural log
p5 <- ggplot(index, aes(year, log_est)) +
  geom_ribbon(aes(ymin=log_est-2*se,ymax=log_est+2*se), fill = "gray",alpha=0.6) +
  # geom_pointrange(aes(ymin=log_est-2*se,ymax=log_est+2*se)) +
  geom_line()+
  xlab("Year") + ylab("Ln abundance")+
  theme_Publication() +
  ylim(4,16)+
  xlim(1960,2020)+
  labs(tag = 'D')


p6 <-  inset_element(p = pcla,
                     left = 0.65,
                     bottom = 0.75,
                     right = 1,
                     top = 1,
                     on_top = TRUE,
                     align_to = 'full') +
  labs(tag = '')

KB <- p4 + p5 + plot_annotation(title = 'Kelp Bass',
                                theme = theme(plot.title = element_text(size = 24, hjust = 0.5))) # Draw combined plot

KB + p6



# 09 Figure 4 predicted density maps -------------------------------------------------------

load(file="01_tidy_data/output_SIA_KB.RData")

# predict on new data
set.seed(011823)

pred.spatial = predict(fit, newdata = new_df, return_tmb_object = FALSE)


# plot pretty maps

library(rnaturalearth)
library(rnaturalearthhires)

map_data <- rnaturalearth::ne_countries(
  scale = "large",
  returnclass = "sf", country = "United States of America")
# Crop the polygon for plotting and efficiency:
st_bbox(map_data) # find the rough coordinates

na_coast <- suppressWarnings(suppressMessages(
  st_crop(map_data,
          c(xmin = -122, ymin = 30, xmax = -116, ymax = 35.5))
))

na_coast_proj <- sf::st_transform(na_coast, crs = 32611)

ggplot(na_coast_proj) + geom_sf()

sf::st_boundary(na_coast_proj)

par(mfrow =c(1,1))

p6 <-   inset_element(p = pcla,
                      left = 0.85,
                      bottom = 0.8,
                      right = 1,
                      top = 1,
                      align_to = 'full')+
  labs(tag = '')

pred.spatial <- pred.spatial %>% 
  filter(year != 1983) # only two stations surveyed this year

## Supplemental Figure S10----------------------
# fixed effects and random effects
ggplot() + 
  geom_tile(data = pred.spatial, aes(x = Lon.km*1000, y = Lat.km*1000, fill = exp(est), color = exp(est))) +
  scale_fill_viridis_c( trans = "sqrt",
                        # trim extreme high values to make spatial variation more visible
                        na.value = "yellow", limits = c(0, quantile(exp(pred.spatial$est), 0.995))) +
  scale_color_viridis_c( trans = "sqrt",
                         # trim extreme high values to make spatial variation more visible
                         na.value = "yellow", limits = c(0, quantile(exp(pred.spatial$est), 0.995)),
                         guide = 'none') +
  geom_sf(data=na_coast_proj,fill = "antique white") +
  # geom_sf(data = mex_coast_proj,fill = "antique white")+
  coord_sf(xlim = c(164390.8,526929),ylim = c(3599377,3859781))+
  scale_y_continuous(breaks = seq(32.5,35.5,1)) +
  scale_x_continuous(breaks = seq(-120.5, -117, by = 4)) +
  facet_wrap(~year) +
  theme_Publication() +
  labs(fill = "Ln(abun)") +
  labs(x = "Longitude", y = "Latitude")+
  ggtitle("Kelp Bass",
          subtitle = "Prediction (fixed effects + random effects)")+
  theme(strip.background =element_rect(fill="white"))+
  p6

# just spatiotemporal random effects
ggplot() + 
  geom_raster(data = pred.spatial, aes(x = Lon.km*1000, y = Lat.km*1000, fill = epsilon_st)) +
  scale_fill_gradient2() +
  geom_sf(data=na_coast_proj,fill = "antique white") +
  coord_sf(xlim = c(164390.8,526929),ylim = c(3599377,3859781))+
  scale_y_continuous(breaks = seq(32.5,34.5,1)) +
  scale_x_continuous(breaks = seq(-120, -117, by = 4)) +
  facet_wrap(~year) +
  labs(x = "Longitude", y = "Latitude")+
  theme_Publication() +
  ggtitle("Kelp Bass",
          subtitle = "Prediction (spatiotemporal random effects)")+
  theme(strip.background =element_rect(fill="white"))+
  p6

# 10 Figure 5 center of gravity -------------------------------------------------------

load(file="01_tidy_data/output_SIA_KB.RData")

cog <- get_cog(pred.index, level = 0.75, format = "wide", area = 4)

save(cog, file="01_tidy_data/output_COG_KB.RData")

#---
load(file="01_tidy_data/output_COG_KB.RData")

pcog <- cog %>% 
  mutate(decade = case_when(year > 1969 & year < 1980 ~ "1970s",
                            year < 1970 ~ "1960s",
                            year > 1979 & year < 1990 ~ "1980s",
                            year > 1989 & year < 2000 ~ "1990s",
                            year > 1999 & year < 2010 ~ "2000s",
                            year > 2009 & year < 2020 ~ "2010s")) %>% 
  ggplot(aes(est_x, est_y)) +
  geom_pointrange(aes(xmin = lwr_x, xmax = upr_x),color="gray",linetype=2) +
  geom_pointrange(aes(ymin = lwr_y, ymax = upr_y),color="gray",linetype=2) +
  theme_classic()+
  scale_fill_viridis_c()+
  xlab("Longitude (km)")+
  ylab("Latitude (km)")

pcog

map_data <- rnaturalearth::ne_countries(
  scale = "large",
  returnclass = "sf", country = "United States of America")

# Crop the polygon for plotting and efficiency:
st_bbox(map_data) # find the rough coordinates

na_coast <- suppressWarnings(suppressMessages(
  sf::st_crop(map_data,
              c(xmin = -122, ymin = 30, xmax = -116, ymax = 35))
))

na_coast_proj <- sf::st_transform(na_coast, crs = 32611)


# Pt. Dume shore station sst data
# citation: Carter, Melissa L.; Flick, Reinhard E.; Terrill, Eric; Beckhaus, Elena C.; Martin, Kayla; Fey, Connie L.; Walker, Patricia W.; Largier, John L.; McGowan, John A. (2022). Shore Stations Program, Point Dume - Zuma Beach. In Shore Stations Program Data Archive: Current and Historical Coastal Ocean Temperature and Salinity Measurements from California Stations. UC San Diego Library Digital Collections. https://doi.org/10.6075/J00001XZ

df <- read_csv("00_raw_data/Zuma_TEMP_1956-2021.csv", 
               skip = 47)
df$SST <- as.numeric(df$SURF_TEMP_C)

# get July mean temp off Pt. Dume by year (we need to fill in missing CalCOFI temp data in 1981)
sstdf <- df %>% 
  dplyr::select(YEAR,MONTH,SST) %>% 
  drop_na("SST") %>% 
  dplyr::group_by(YEAR,MONTH) %>% 
  dplyr::filter(MONTH == 7) %>% 
  dplyr::rename(year = YEAR) %>% 
  dplyr::summarise(sst = mean(SST)) # July mean by year

cogwsst <- left_join(cog,sstdf)

# grab CalCOFI station-specific temperatures and join with cog data

sstdf.2 <- kb %>% 
  dplyr::select(year,temp.10m) %>% 
  dplyr::group_by(year) %>% 
  dplyr::summarise(meansst = mean(temp.10m, na.rm=TRUE))

cogwsst <- left_join(cog,sstdf.2)
cogwsst[9,10] <- 19.40968 # fill in missing CalCOFI temp data in 1981 from Pt. Dume


# assign annual COG estimates to corresponding decade

cogwdec <- cogwsst %>% 
  mutate(decade = case_when(year > 1969 & year < 1980 ~ "1970s",
                            year < 1970 ~ "1960s",
                            year > 1979 & year < 1990 ~ "1980s",
                            year > 1989 & year < 2000 ~ "1990s",
                            year > 1999 & year < 2010 ~ "2000s",
                            year > 2009 & year < 2020 ~ "2010s")) %>% 
  mutate(est_x_m = est_x * 1000,
         est_y_m = est_y * 1000,
         lwr_x_m = lwr_x * 1000,
         upr_x_m = upr_x * 1000,
         lwr_y_m = lwr_y * 1000,
         upr_y_m = upr_y * 1000) %>% 
  filter(year != 1983)
cogwdec
# Farthest northCOG northing: 1981, 3759.756 km; Farthest south COG northing: 1965 3661.771 km


# cool way to plot contours using density quantiles (from https://stackoverflow.com/questions/75598144/interpretation-of-2d-density-estimate-charts)

density_quantiles <- function(x, y, quantiles) {
  dens <- MASS::kde2d(x, y, n = 500)
  df   <- cbind(expand.grid(x = dens$x, y = dens$y), z = c(dens$z))
  r    <- terra::rast(df)
  ind  <- sapply(seq_along(x), function(i) cellFromXY(r, cbind(x[i], y[i])))
  ind  <- ind[order(-r[ind][[1]])]
  vals <- r[ind][[1]]
  ret  <- approx(seq_along(ind)/length(ind), vals, xout = quantiles)$y
  replace(ret, is.na(ret), max(r[]))
}


quantiles <- c(0, 0.5, 0.75, 0.9)

theme_Publication <- function(base_size=14, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(0.5),vjust = 0), #hjust = 0.5),
            # text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border =  element_blank(),#element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major.y = element_line(colour="#f0f0f0"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic",size = rel(0.75)),
            legend.text = element_text(size = rel(0.75)),
            plot.margin=unit(c(5,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


# plot annual COG estimates and kernel density quantiles (COG estimates color coded by decade)
p6 <-  inset_element(p = pcla,
                     left = 0.025,
                     bottom = 0.025,
                     right = 0.25,
                     top = 0.25,
                     # on_top = TRUE,
                     align_to = 'panel') +
  labs(tag = '')

e2 <- ggplot(na_coast_proj) +
  theme(panel.background = element_rect(fill = "aliceblue")) +
  geom_pointrange(data = cogwdec,aes(x=est_x_m, y=est_y_m,xmin = lwr_x_m, xmax = upr_x_m),alpha = 0) + # making these transparent because when we simply turn them off, the kde function cuts off the quartiles
  geom_pointrange(data = cogwdec,aes(x=est_x_m, y=est_y_m,ymin = lwr_y_m, ymax = upr_y_m),alpha = 0) +
  geom_density2d_filled(data = cogwdec, aes(x=est_x_m, y=est_y_m,fill = after_stat(level)),
                        alpha=0.5,
                        contour_var = "density",
                        breaks = density_quantiles(cogwdec$est_x_m, cogwdec$est_y_m, quantiles)) +
  scale_fill_viridis_d('Quantiles', 
                       labels = scales::percent(quantiles[-1]),
                       direction = -1,
                       option="G")+
  scale_alpha(guide = 'none')+
  guides(color = 'none',
         fill = "none")+
  geom_point(data = cogwdec,aes(x=est_x_m, y=est_y_m,color=year),size = 3,shape=20)+
  scale_color_viridis_b('Decade',
                        option = "viridis",
                        direction = 1)+
  labs(x = "Longitude", y = "Latitude")+
  labs(subtitle = "Kelp Bass")+
  theme_Publication()+
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    plot.background = element_blank(),
    panel.border = element_blank())+
  coord_sf(xlim = c(164390.8,526929),ylim = c(3599377,3859781))+
  scale_y_continuous(breaks = seq(32.5,34.5,1)) +
  geom_sf(fill = "antique white") +
  labs(tag = 'C') +
  p6 
e2

# scatter plot of annual COG estimates by mean July SST (estimates color coded by decade)
f<- ggplot(cogwdec, aes(x=meansst, y=est_y,color=year)) + 
  geom_point(size=1.5)+
  geom_smooth(method="lm", se=TRUE, col="red", size=0.5,linetype = 1,fill="gray39")  +
  scale_color_viridis_b('Decade',
                        option = "viridis",
                        direction = 1)+
  guides(color ="none")+
  scale_x_continuous(limits = c(15.5,19.7,by=2))+
  ylim(3650,3800) +
  xlab("Mean SST")+
  ylab("Northing (km)")+
  theme_classic()+
  theme(legend.position="right")+
  labs(tag = 'D')
f

mod <- lm(data=cogwdec, est_y ~ meansst)
summary(mod)

# plot Figure 5C,D
g2 <- e2+
  inset_element(f, left = 0.6, bottom = 0.55, right = 1.1, top = 1.2) 

g2

# get median COG for 1960s and 2010s
cogwdec %>% 
  group_by(decade) %>% 
  filter(decade == "1960s" | decade == "2010s") %>% 
  summarise(meanNorth = median(est_y))

library(ggpubr)
ggplot(data=subset(cogwdec, decade %in% "1960s" | decade %in% "2010s"),aes(decade,est_y))+
  geom_boxplot(outlier.shape = NA, draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_jitter(width=0.1,color="darkgray")+
  theme_Publication()+
  xlab("Species")+
  ylab("Ln abundance")+
  #  Add p-value
  stat_compare_means(method = "wilcox", size=8)

compare_means(est_y ~ decade, subset(cogwdec, decade %in% "1960s" | decade %in% "2010s"), method = "wilcox.test")
wilcox.test(est_y ~ decade, subset(cogwdec, decade %in% "1960s" | decade %in% "2010s"))


# this is plotting mean July sst by year

z <- ggplot(cogwdec,aes(year,meansst,color=year))+
  geom_point(size=4)+
  scale_color_viridis_b('Decade',
                        option = "viridis",
                        direction = 1)+
  geom_smooth(method="lm", se=TRUE, col="red", size=1,linetype =4,fill="darkgray")
z

zmod <- lm(data=cogwdec, year ~ meansst)
summary(zmod)
