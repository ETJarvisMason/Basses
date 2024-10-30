# Jarvis Mason, ET
# 25 OCT 2024


# load libraries ----------------------------------------------------------

library(dplyr)
library(sdmTMB)
library(sdmTMBextra)
library(ggplot2)
library(ggpubr)
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

pneb <- readPNG("03_plots/bsb.png", native = TRUE)



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
  
  print(list(fit,which(san==FALSE),AIC,tidy))
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

# additional formatting

names(df.intrp)
bsb <- df.intrp
bsb$moon.facnum <- as.numeric(bsb$moon.fac)
bsb$scale_day <- scale(bsb$day)
bsb$day2 <- bsb$day^2
bsb$scale_dist2shore <- scale(bsb$dist_km_calc)
bsb$scale_moon.num <- scale(bsb$moon.num)



## set seed------------

set.seed(022824)

# 02 make a mesh -------------------------------------------------------------

# make a mesh -- which is approximating the spatial surface
# note that the data frame 'bsb' already contains lat and lon in UTMs -- which is what
# we want to use

mesh = make_mesh(bsb, xy_cols = c("Lon.km","Lat.km"), cutoff = 8, type = "kmeans")

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



# 03 fit models having different fixed effects ----------------------------------------
# NOTE: use reml = FALSE (here, AIC based on fixed effects, random effects must stay the same across model comparisons)
# spatiotemporal and Tespatial random effects on


# test the effect of distance to shore, call it m1a
# specify some function arguments

names(bsb)

# extra time argument (no data these years)
xt <- c(1967, 1968, 1970, 1971, 1973, 1974, 1976, 1977,
        1979, 1980, 1982, 1985, 1986, 1987, 1988, 1990, 1991, 1993, 1994,
        1996, 1999, 2008, 2010, 2011)

covar <- c("0","scale_dist2shore") # make sure your variables are named correctly; scale covariates where necessary
m1a <- check_mods("on", "iid", "BSB", covar, "year", bspde, bsb, xt,tv=formula(~1),reml=FALSE)


# now try adding in moon phase

names(bsb)
covar <- c("0","scale_dist2shore","s(moon.facnum, k=4)")
m1b <- check_mods("on", "iid", "BSB", covar, "year", bspde, bsb, xt,tv=formula(~1),reml=FALSE)

covar <- c("0","scale_dist2shore","scale_moon.num")
m1b2 <- check_mods("on", "iid", "BSB", covar, "year", bspde, bsb, xt,tv=formula(~1),reml=FALSE)

covar <- c("0","scale_dist2shore","s(scale_moon.num)")
m1b3 <- check_mods("on", "iid", "BSB", covar, "year", bspde, bsb, xt,tv=formula(~1),reml=FALSE)

covar <- c("0","scale_dist2shore","s(scale_moon.num,bs='cc',k=4)")
m1b4 <- check_mods("on", "iid", "BSB", covar, "year", bspde, bsb, xt,tv=formula(~1),reml=FALSE)

# now try adding in hour, call it m1c

covar <- c("0","scale_dist2shore","s(hour, bs ='cc',k = 24)")   # use single quotes to specify basis function
covar <- c("0","scale_dist2shore","s(hour, bs =\"cc\",k = 24)") # another way that works
m1c <- check_mods("on", "iid", "BSB", covar, "year", bspde, bsb, xt,tv=formula(~1),reml=FALSE)

# now try adding in day_night by net type, call it m1d

covar <- c("0","scale_dist2shore","s(day_night, by = net_type, bs = 'fs')")
m1d <- check_mods("on", "iid", "BSB", covar, "year", bspde, bsb, xt,tv=formula(~1),reml=FALSE)

covar <- c("0","scale_dist2shore","s(net_type, bs = 'fs')")
m1d2 <- check_mods("on", "iid", "BSB", covar, "year", bspde, bsb, xt,tv=formula(~1),reml=FALSE)

covar <- c("0","scale_dist2shore","net_type")
m1d3 <- check_mods("on", "iid", "BSB", covar, "year", bspde, bsb, xt,tv=formula(~1),reml=FALSE)

# now try adding in Julian day, call it m1e

covar <- c("0","scale_dist2shore","scale_day")
m1e <- check_mods("on", "iid", "BSB", covar, "year", bspde, bsb, xt,tv=formula(~1),reml=FALSE)

covar <- c("0","scale_dist2shore","s(scale_day, k=3)")
m1f <- check_mods("on", "iid", "BSB", covar, "year", bspde, bsb, xt,tv=formula(~1),reml=FALSE)

covar <- c("0","s(scale_dist2shore)","day","day2")
m1g <- check_mods("on", "iid", "BSB", covar, "year", bspde, bsb, xt,tv=formula(~1),reml=FALSE)

covar <- c("0","s(scale_dist2shore)","scale_day","s(moon.num,bs='cc',k=4)")
m1h <- check_mods("on", "iid", "BSB", covar, "year", bspde, bsb, xt,tv=formula(~1),reml=FALSE)

covar <- c("0","s(scale_dist2shore)","scale_day","moon.num")
m1i <- check_mods("on", "iid", "BSB", covar, "year", bspde, bsb, xt,tv=formula(~1),reml=FALSE)

covar <- c("0","s(scale_day)")
m1j <- check_mods("on", "iid", "BSB", covar, "year", bspde, bsb, xt,tv=formula(~1),reml=FALSE)


# 04 check model fit ------------------------------------------------------

# sanity checks
c(
  m1a[2], #1  
  m1b[2], #2  
  m1b2[2],#3  
  m1b3[2],#4  
  m1b4[2],#5  
  m1c[2], #6
  m1d[2], #7
  m1d2[2],#8  
  m1d3[2],#9
  m1e[2], #10  
  m1f[2], #11 
  m1g[2], #12 
  m1h[2], #13
  m1i[2], #14
  m1j[2]) #15 

# mods 1-5,8,10 passed sanity checks


# check parsimony
# AIC
as.data.frame(tibble(
  model = c(
    c("0,scale_dist2shore"),
    c("0,scale_dist2shore,s(moon.facnum, k=4)"),
    c("0,scale_dist2shore,s(moon.num, k=4)"),
    c("0,scale_dist2shore,s(moon.num)"),
    c("0,scale_dist2shore,s(moon.num,bs='cc',k=4)"), # best moon model
    c("0,scale_dist2shore,s(net_type, bs = 'fs')"),
    c("0,scale_dist2shore,scale_day")
  ),
  AIC = c(
    unlist(m1a[3]),  
    unlist(m1b[3]), 
    unlist(m1b2[3]),
    unlist(m1b3[3]), 
    unlist(m1b4[3]),
    unlist(m1d2[3]),
    unlist(m1e[3])  
  ))) %>% 
  arrange(AIC)

# m1a and m1b most parsimonious

# check if random effects structure with spatial off fits better 
# NOTE: use reml = TRUE (here, AIC based on random effects)

covar <- c("0","scale_dist2shore")
m1a.spOff <- check_mods("off", "iid", y="BSB", covar, "year", bspde, bsb, xt,tv=formula(~1),reml=TRUE)
covar <- c("0","scale_dist2shore","s(moon.facnum, k=4)")
m1b.spOff <- check_mods("off", "iid", y="BSB", covar, "year", bspde, bsb, xt,tv=formula(~1),reml=TRUE)



## cross validation --------------------------------------------------------
set.seed(30) # must run set.seed before assigning clust to ensure reproducibility
clust <- sample(1:8, size = nrow(bsb), replace = T)
k_folds <- length(unique(clust))

## function to run models for cross validation (includes args for clust & k_folds) ---------------
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

covar <- c("0","scale_dist2shore")
m1.cv <- check_cv_mods("on", "iid", "BSB", covar, "year", bspde, bsb, tv=formula(~1),clust,k_folds)

covar <- c("0","scale_dist2shore","scale_day")
m2.cv <- check_cv_mods("on", "iid", "BSB", covar, "year", bspde, bsb, tv=formula(~1),clust,k_folds)

covar <- c("0","scale_dist2shore","s(moon.num,bs='cc',k=4)")
m3.cv <- check_cv_mods("on", "iid", "BSB", covar, "year", bspde, bsb, tv=formula(~1),clust,k_folds)

covar <- c("0","scale_dist2shore","s(net_type, bs = 'fs')")
m4.cv <- check_cv_mods("on", "iid", "BSB", covar, "year", bspde, bsb, tv=formula(~1),clust,k_folds)


# convergence across folds?
c(m1.cv$converged,
  m2.cv$converged,
  m3.cv$converged,
  m4.cv$converged) # no convergence

# compare models
as.data.frame(tibble(
  model = c(
    c("0,scale_dist2shore"),
    c("0,scale_dist2shore,scale_day"),
    c("0,scale_dist2shore,s(moon.num,bs='cc',k=4)"),
    c("0,scale_dist2shore,s(net_type, bs = 'fs')")
  ),
  AIC = c(
    unlist(m1a[3]),  
    unlist(m1e[3]), 
    unlist(m1b[3]),
    unlist(m1d2[3])
    ),
  sum_LL = c(
    m1.cv$sum_loglik,   
    m2.cv$sum_loglik, 
    m3.cv$sum_loglik, 
    m4.cv$sum_loglik
    )
  )) %>% 
  arrange(AIC)



# 05 save best predictive mod as determined by cross-validation----------------------------------

fit <- sdmTMB(BSB ~ 0 + scale_dist2shore,
              mesh = bspde,
              time = "year", # name of time variable
              spatial = "on", # can be on / off. On estimates a shared spatial field that is the same across years
              spatiotemporal = "iid", # should be iid / off
              time_varying = ~ 1,
              extra_time = c(1967, 1968, 1970, 1971, 1973, 1974, 1976, 1977,
                             1979, 1980, 1982, 1985, 1986, 1987, 1988, 1990, 1991, 1993, 1994,
                             1996, 1999, 2008, 2010, 2011),
              family = tweedie(), # could also be one of the delta_ families
              data = bsb,
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
# we can fill in the corresponding missing depth data using a bathymetry raster
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

# check it for no islands
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

# now expand our data frame to include all time variables and any covariates in our model (year, distance, moon phase, and
# lat/lon in kilometers) - this is a necessary step for prediction

# replication factor  
yrs <- sort(unique(bsb$year))
n <- length(yrs)

# replicate our dataframe for year in our model, call our df, 'new_df'
new_df <- do.call("rbind", replicate( 
  n, df, simplify = FALSE)) 

# add year for designating time
new_df$year_f <-  as.factor(rep(yrs, each = nrow(df))) # if a factor in model
new_df$year <-  rep(yrs, each = nrow(df))


# grab mean and stdev of distance column in bsb 
# these will be used to scale the prediction grid distances
mu.dist <- mean(bsb$dist_km_calc)     
sd.dist <- sd(bsb$dist_km_calc)      

new_df <- new_df %>% 
  dplyr::filter(dist > 0) %>% 
  mutate(scale_dist2shore = (dist - mu.dist) / sd.dist) %>% # use same name as covariate name in model
  filter(area.km>3.9999)



# 07 calculate index of abundance----------------------------------------------------

# make predictions
pred.index = predict(fit, newdata = new_df, return_tmb_object = TRUE)

# get index
index = get_index(pred.index, bias_correct = TRUE, area = 4)

save(bsb,df,new_df,fit,pred.index,index, file="01_tidy_data/output_SIA_BSB.RData")

## index with missing years (skip to plot index above)-----------------------------------------------------


# replication factor  
yrs <- seq(1963,2016,1) 
n <- length(yrs)

# replicate our dataframe for year in our model, call our df, 'new_df'
new_df <- do.call("rbind", replicate( 
  n, df, simplify = FALSE)) 

# add year for designating time
new_df$year_f <-  as.factor(rep(yrs, each = nrow(df))) # if a factor in model
new_df$year <-  rep(yrs, each = nrow(df))


# grab mean and stdev of distance column in bsb 
# these will be used to scale the prediction grid distances
mu.dist <- mean(bsb$dist_km_calc)     
sd.dist <- sd(bsb$dist_km_calc)      

new_df <- new_df %>% 
  dplyr::filter(dist > 0) %>% 
  mutate(scale_dist2shore = (dist - mu.dist) / sd.dist) %>% # use same name as covariate name in model
  filter(area.km>3.9999)

# make predictions
pred.index = predict(fit, newdata = new_df, return_tmb_object = TRUE)

# get index
index = get_index(pred.index, bias_correct = TRUE, area = 4)

save(bsb,df,new_df,fit,pred.index,index, file="01_tidy_data/output_SIA_BSB_with_missing_yrs.RData")

# 08 Figure 2 temporal trends ---------------------------

load(file="01_tidy_data/output_SIA_BSB.RData")

# load instead, to plot index with missing years # (Fig. S10)
load(file="01_tidy_data/output_SIA_BSB_with_missing_yrs.RData")


# on linear scale (no transformation)
p1 <- ggplot(index, aes(year, est)) + 
  geom_ribbon(aes(ymin=lwr,ymax=upr), fill = "gray",alpha=0.6) +
  # geom_pointrange(aes(ymin=lwr,ymax=upr))+
  geom_line(linewidth = 0.5, color = "grey30")+
  xlim(1960,2020)+
  ylim(c(0,1000000))+
  # ggtitle("Barred Sand Bass")+
  xlab("Year") + ylab("Abundance")+
  theme_Publication()  +
  labs(tag = 'A')


p1 


# natural log
p2 <- ggplot(index, aes(year, log_est)) +
  geom_ribbon(aes(ymin=log_est-2*se,ymax=log_est+2*se), fill ="gray",alpha=0.6) +
  # geom_pointrange(aes(ymin=log_est-2*se,ymax=log_est+2*se)) +
  geom_line()+
  # ggtitle("Barred Sand Bass")+
  xlab("Year") + ylab("Ln abundance")+
  theme_Publication() +
  ylim(4,16)+
  xlim(1960,2020)+
  labs(tag = 'B')

p3 <- inset_element(p = pneb,
                    left = 0.65,
                    bottom = 0.75,
                    right = 1,
                    top = 1,
                    align_to = 'full') +
  labs(tag = '')
BSB <- p1 + p2 + plot_annotation(title = 'Barred Sand Bass',
                                 theme = theme(plot.title = element_text(size = 24, hjust = 0.5))) # Draw combined plot

BSB + p3

## Supplemental Figure S8----------------------

load(file="01_tidy_data/output_SIA_KB.RData")
index.kb <- index %>% 
  mutate(sp = "Kelp Bass")

load(file="01_tidy_data/output_SIA_BSB.RData")
index.bsb <- index%>% 
  mutate(sp = "Barred Sand Bass")

index.both <- full_join(index.kb,index.bsb)

library(ggpubr)


ggplot(data=index.both,aes(sp,log_est))+
  geom_boxplot(outlier.shape = NA, draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_jitter(width=0.1,color="darkgray")+
  stat_summary(fun = "mean",
               geom = "point",
               color = "red",
               size = 3)+
  theme_Publication()+
  xlab("Species")+
  ylab("Ln abundance")+
  #  Add p-value
  stat_compare_means(method = "anova", size=8)

# Compute the analysis of variance
res.aov <- aov(log_est ~ sp, data = index.both)
# Summary of the analysis
summary(res.aov)

index.both %>% 
  group_by(sp) %>% 
  summarize(mean=mean(log_est),
            sd=sd(log_est))

# 09 Figure 3 predicted density maps -------------------------------------------------------

load(file="01_tidy_data/output_SIA_BSB.RData")

# predict on new data
set.seed(011823)

pred.spatial = predict(fit, newdata = new_df, return_tmb_object = FALSE)


# plot pretty maps

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

p3 <-   inset_element(p = pneb,
                      left = 0.85,
                      bottom = 0.8,
                      right = 1,
                      top = 1,
                      align_to = 'full')+
  labs(tag = '')

pred.spatial <- pred.spatial %>% 
  filter(year != 1983) # only 2 stations surveyed this year

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
  coord_sf(xlim = c(164390.8,526929),ylim = c(3599377,3859781))+
  scale_y_continuous(breaks = seq(32.5,35.5,1)) +
  scale_x_continuous(breaks = seq(-120.5, -117, by = 4)) +
  facet_wrap(~year) +
  theme_Publication() +
  labs(fill = "Ln(abun)") +
  labs(x = "Longitude", y = "Latitude")+
  ggtitle("Barred Sand Bass",
          subtitle = "Prediction (fixed effects + random effects)")+
  theme(strip.background =element_rect(fill="white"))+
  p3

## Supplemental Figure S9----------------------
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
  ggtitle("Barred Sand Bass",
          subtitle = "Prediction (spatiotemporal random effects)")+
  theme(strip.background =element_rect(fill="white"))+
  p3


# 10 Figure 5 center of gravity -------------------------------------------------------

# load data and predictions for index
load(file="01_tidy_data/output_SIA_BSB.RData")

cog <- get_cog(pred.index, level = 0.75, format = "wide", area = 4)

save(cog, file="01_tidy_data/output_COG_BSB.RData")

load(file="01_tidy_data/output_COG_BSB.RData")

map_data <- rnaturalearth::ne_countries(
  scale = "large",
  returnclass = "sf", country = "United States of America")

# Crop the polygon for plotting efficiency:
st_bbox(map_data) # first find the rough coordinates

na_coast <- suppressWarnings(suppressMessages(
  sf::st_crop(map_data,
              c(xmin = -122, ymin = 30, xmax = -116, ymax = 35))
))

na_coast_proj <- sf::st_transform(na_coast, crs = 32611)


# Pt. Dume shore station sst data
# cite:"Point Dume measurements collected by the Los Angeles County Lifeguard service. Data provided by the Shore Stations Program sponsored at Scripps Institution of Oceanography by California State Parks, Division of Boating and Waterways. Contact: shorestation@ucsd.edu."
# citation: Carter, Melissa L.; Flick, Reinhard E.; Terrill, Eric; Beckhaus, Elena C.; Martin, Kayla; Fey, Connie L.; Walker, Patricia W.; Largier, John L.; McGowan, John A. (2022). Shore Stations Program, Point Dume - Zuma Beach. In Shore Stations Program Data Archive: Current and Historical Coastal Ocean Temperature and Salinity Measurements from California Stations. UC San Diego Library Digital Collections. https://doi.org/10.6075/J00001XZ

df <- read_csv("00_raw_data/Zuma_TEMP_1956-2021.csv", 
               skip = 47)
df$SST <- as.numeric(df$SURF_TEMP_C)

# note that some months may have missing temp days, grab
# monthly means for summer months and then calculate an
# annual summer mean
pdsst <- df %>% 
  dplyr::select(YEAR,MONTH,SST) %>% 
  drop_na("SST") %>% 
  dplyr::group_by(YEAR,MONTH) %>% 
  dplyr::filter(MONTH == 7) %>% 
  dplyr::summarise(mean = mean(SST)) # monthly means by year

sstdf <- pdsst %>% 
  dplyr::group_by(YEAR) %>% 
  dplyr::rename(year = YEAR) %>% 
  dplyr::summarise(sst = mean(mean)) # annual mean of monthly means

cogwsst <- left_join(cog,sstdf)

sstdf.2 <- bsb %>% 
  dplyr::select(year,temp.10m) %>% 
  dplyr::group_by(year) %>% 
  dplyr::summarise(meansst = mean(temp.10m, na.rm=TRUE))

cogwsst <- left_join(cog,sstdf.2)
cogwsst[9,10] <- 19.40968 # fill in temp data from Pt. Dume


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
         upr_y_m = upr_y * 1000)%>% 
  filter(year != 1983)
cogwdec

pneb <- readPNG("03_plots/bsb.png", native = TRUE)

p7 <-  inset_element(p = pneb,
                     left = 0.025,
                     bottom = 0.025,
                     right = 0.25,
                     top = 0.25,
                     # on_top = TRUE,
                     align_to = 'panel') +
  labs(tag = '')

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
 
a2 <- ggplot(na_coast_proj) +
  theme(panel.background = element_rect(fill = "aliceblue")) +
  geom_density2d_filled(data = cogwdec, aes(x=est_x_m, y=est_y_m,fill = after_stat(level)),
                        alpha=0.5,
                        contour_var = "density",
                        breaks = density_quantiles(cogwdec$est_x_m, cogwdec$est_y_m, quantiles)) +
  scale_fill_viridis_d('Quantiles', 
                       labels = scales::percent(quantiles[-1]),
                       direction = -1,
                       option="G")+
  scale_alpha(guide = 'none')+
  geom_point(data = cogwdec,aes(x=est_x_m, y=est_y_m,color=year),size = 3,shape=20)+
  scale_color_viridis_b('Decade',
                        option = "viridis",
                        direction = 1)+
  labs(x = "Longitude", y = "Latitude")+
  labs(subtitle = "Barred Sand Bass")+
  theme_Publication()+
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    plot.background = element_blank(),
    panel.border = element_blank())+
  coord_sf(xlim = c(164390.8,526929),ylim = c(3599377,3859781))+
  scale_y_continuous(breaks = seq(32.5,34.5,1)) +
  geom_sf(fill = "antique white") +
  labs(tag = 'A') +
  p7 
a2



# Farthest north COG northing: 2006, 3753.554 km; Farthest south COG northing: 2012 3676.620 km 
b<- ggplot(cogwdec, aes(x=meansst, y=est_y,color=year)) + 
  geom_point(size=1.5)+
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
  labs(tag = 'B')
b

c2 <- a2+b + plot_annotation(title="Barred Sand Bass",
                             theme = theme(plot.title = element_text(size = 18, hjust = 0.5)))
c2 <- a2+
  inset_element(b, left = 0.6, bottom = 0.55, right = 1.1, top = 1.2) 

c2

# get median COG for 1960s and 2010s
cogwdec %>% 
  group_by(decade) %>% 
  filter(decade == "1960s" | decade == "2010s") %>% 
  summarise(meanNorth = median(est_y))

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

# load KB index for plotting with BSB
load(file="01_tidy_data/output_COG_KB.RData")
pcla <- readPNG("03_plots/kb.png", native = TRUE)

p6 <-  inset_element(p = pcla,
                     left = 0.025,
                     bottom = 0.025,
                     right = 0.25,
                     top = 0.25,
                     # on_top = TRUE,
                     align_to = 'panel') +
  labs(tag = '')

cogwsst <- left_join(cog,sstdf.2) # stdf.2 from section above
cogwsst[9,10] <- 19.40968 # fill in temp data from Pt. Dume


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


e2 <- ggplot(na_coast_proj) +
  theme(panel.background = element_rect(fill = "aliceblue")) +
  geom_pointrange(data = cogwdec,aes(x=est_x_m, y=est_y_m,xmin = lwr_x_m, xmax = upr_x_m),size = 0,color="gray",linetype=2,alpha = 0) +
  geom_pointrange(data = cogwdec,aes(x=est_x_m, y=est_y_m,ymin = lwr_y_m, ymax = upr_y_m),size = 0,color="gray",linetype=2,alpha = 0) +
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


# Farthest northCOG northing: 1981, 3759.756 km; Farthest south COG northing: 1965 3661.771 km (difference of 97.985 km)
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


g2 <- e2+
  inset_element(f, left = 0.6, bottom = 0.55, right = 1.1, top = 1.2) 

g2

(c2 / g2) + plot_layout(guides = "collect") + theme(legend.position = "bottom",
                                                    legend.direction = "horizontal")

