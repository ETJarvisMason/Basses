# ET Mason
# 29 OCT 2023

# load libraries ----------------------------------------------------------

library(INLA)
library(sdmTMB)
library(tidyverse)
library(ggplot2)
library(sf)
library(sp)
library(visreg)
library(marmap)
library(lubridate)
library(ggeffects)
library(sdmTMBextra)
library(rnaturalearth)
library(mgcv)
library(patchwork)
library(magick)
library(ggimage)
library(performance)

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


# Import fish drawings
library("png")                
pneb <- readPNG("03_plots/bsb.png", native = TRUE)
pneb.magick <- (file="03_plots/bsb.png")
pcla <- readPNG("03_plots/kb.png", native = TRUE)
pcla.magick <- (file="03_plots/kb.png")


# 01 load data --------------------------------------------------------

load(file="01_tidy_data/scaled_data.Rdata")

# additional formatting
head(bsb.scaled)
bsb.scaled$scale_day <- scale(bsb.scaled$day)

# plot the distribution of July cruise data
ggplot(bsb.scaled,aes(Lon.km,Lat.km, size = BSB, colour = log(BSB + 1)))+
  geom_point()+
  facet_wrap(~year)

# years with zero positives across stations (9 years: 1964, 66, 72, 78, 92, 95, 01, 02, 05)
yrs.no.BSB <- unlist(bsb.scaled %>%
                       dplyr::group_by(year) %>%
                       dplyr::summarise(density = sum(BSB)) %>%
                       dplyr::filter(density == 0) %>%
                       dplyr::select(year))



# Explore relationships among environmental covariates 

justvars <- bsb.scaled %>%
  ungroup() %>% 
  dplyr::select(temp.10m,npgo,oni,zoop,scale_day,dist_km_calc,bakun33,bakun30,south.ild,south.cuti,south.beuti,calpac)
source("http://www.sthda.com/upload/rquery_cormat.r")
library(corrplot)
par(mfrow=c(1,1))
rquery.cormat(justvars)

# From Hunsicker et al. 2022 (Tracking and forecasting community
# responses to climate perturbations in the
# California Current Ecosystem):

# isothermal layer depth (ILD):  The ILD is similar to mixed layer depth and defines the depth where temperature
# deviates by 0.5˚C from the surface value.

# CUTI and BEUTI are
# upwelling indices that quantify vertical transport and nitrate flux through the base of the
# mixed layer, respectively [44]. The data were annually averaged (July-June) from the coast to
# 100 km offshore, with the exception of CUTI and BEUTI, which capture coastal upwelling
# within 75 km of shore. In the alongshore direction, we calculated averages for two regions
# with a division at Point Conception, California, separating the southern portion of the CCE
# (31–34.5˚N) from the central region (34.5–40.5˚N, Fig 1). This is in response to the recognition
# of Point Conception as a major biogeographic boundary for the California Current System,
# with differing wind and current patterns north and south of that feature [45, 46]. The
# annual averages were taken from July to June to capture the influence of the El Niño–Southern
# Oscillation (ENSO), which peaks in winter and is the dominant mode of interannual variability
# influencing the California Current [47].

# "We found that nitrate flux through the base of the
# mixed layer (BEUTI) exhibited the strongest relationship with species and community-level
# responses." This was for the central coast. Ild was important off southern California.

# r : The table of correlation coefficients
# p : Table of p-values corresponding to the significance levels of the correlations
# sym : A representation of the correlation matrix in which coefficients are replaced by symbols according to the strength of the dependence. For more description, see this article: Visualize correlation matrix using symnum function
# In the generated graph, negative correlations are in blue and positive ones in red color.

# correlations above 0.4 are relatively strong,
# between 0.2 and 0.4 are moderate,
# below 0.2 are weak
rquery.cormat(justvars, type="flatten", graph=TRUE)

## set seed------------

set.seed(022824)

# 02 BSB Time Period 1 ------------------------------------------------------------------

# first drop nas

# NOTE:  only npgo and oni go back as far as CalCOFI data
# in Model 2 we will also include upwelling indices (bakun33, bakun30, or beuti/cuti) and isothermal layer depth for na drop
bsb.scaled2 <- bsb.scaled %>%
  drop_na(scale_day,npgo,zoop,temp.10m,dist_km_calc,oni)

# create new mesh based on new dataframe
map_data <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf", country = "United States of America")

na_coast <- suppressWarnings(suppressMessages(
  sf::st_crop(map_data,
              c(xmin = -122, ymin = 30, xmax = -116, ymax = 36))))

na_coast_proj <- sf::st_transform(na_coast, crs = 32611)

mesh2 = make_mesh(bsb.scaled2, xy_cols = c("Lon.km","Lat.km"), cutoff = 10)

bspde2 <- add_barrier_mesh(
  mesh2, na_coast_proj, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)

# specify some function arguments

names(bsb.scaled2)
covar <- c("scale_day","dist_km_calc", "temp.10m" , "npgo" , "zoop" , "oni") # make sure your variables are named correctly; scale covariates where necessary

fit1a.bsb <- sdmTMB(
  BSB ~ scale_day + dist_km_calc + temp.10m + npgo + zoop + oni,
  mesh = bspde2,
  time = "year", # name of time variable
  spatial = "off", # can be on / off. On estimates a shared spatial field that is the same across years
  spatiotemporal = "iid", # should be iid / off
  family = tweedie(), # could also be one of the delta_ families
  data = bsb.scaled2,
  reml = FALSE,
  silent = FALSE)

## check fit  ---------------------------------------------------------------

fit <- fit1a.bsb
sanity(fit)
AIC(fit)

check_collinearity(fit)

max(fit$gradients)
fit$sd_report
tidy(fit,"ran_pars",conf.int = TRUE)
tidy(fit,"fixed",conf.int = TRUE)

## k-fold CV - BSB Time Period 1--------------------------------

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

set.seed(56)
clust <- sample(1:8, size = nrow(bsb.scaled2), replace = T)
k_folds <- length(unique(clust))

covar <- c("scale_day","dist_km_calc","temp.10m","npgo","zoop","oni")
m1_cv.bsb <- check_cv_mods("off", "iid", "BSB", covar, "year", mesh=bspde2, dataset=bsb.scaled2, tv=NULL,clust,k_folds)

# results from model without cross-validation
m1_bsb <- check_mods("off", "iid", y="BSB", covar, "year", mesh=bspde2, dataset=bsb.scaled2, xt=NULL,tv=NULL,reml=FALSE) 


# now substitute zoop for calpac
dataset <- bsb.scaled %>% 
  drop_na(scale_day,npgo,calpac,temp.10m,dist_km_calc,oni) # cal pac
mesh6 <- make_mesh(dataset, xy_cols = c("Lon.km","Lat.km"), cutoff = 10)
mesh <- add_barrier_mesh(
  mesh6, na_coast_proj, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)
set.seed(42)
clust <- sample(1:8, size = nrow(dataset), replace = T)
k_folds <- length(unique(clust))
covar <- c("scale_day","dist_km_calc","temp.10m","npgo","calpac","oni")
m2_cv.bsb <- check_cv_mods("off", "iid", "BSB", covar, "year", mesh=mesh, dataset=dataset, tv=NULL,clust,k_folds)

# results from model without cross-validation
m2_bsb <- check_mods("off", "iid", y="BSB", covar, "year", mesh=mesh, dataset=dataset, xt=NULL,tv=NULL,reml=FALSE) 

# Compare models

as.data.frame(tibble(
  model = c(
    c("scale_day,zoop,dist_km_calc,temp.10m,npgo,oni"),
    c("scale_day,calpac,dist_km_calc,temp.10m,npgo,oni") # make sure your variables are named correctly; scale covariates where necessary
  ),
  AIC = c(
    unlist(m1_bsb[3]),   
    unlist(m2_bsb[3])
  ),
  sum_LL = c(
    m1_cv.bsb$sum_loglik,   
    m2_cv.bsb$sum_loglik
  ))) %>% 
  arrange(AIC)

# 03 BSB Time Period 2 ------------------------------------------------------------------


# let's incorporate other data available that relates to upwelling
# and mixed layer depth

# include bakun33/bakun30 and south.ild for na drop,
# this represents a shorter time series since south.ild is only available since 1980,
# but no CalCOFI temp data in 1981 or 1983 anyway

# call it bsb.scaled3
bsb.scaled3 <- bsb.scaled %>%
  drop_na(scale_day,bakun33,npgo,zoop,temp.10m,dist_km_calc,oni,south.ild)

# NOTE: no temp.10m from 198107 or 198307; also no bucket temps from those two cruises either
# thus, the earliest year is 1984

# create new mesh based on new dataframe
mesh3 = make_mesh(bsb.scaled3, xy_cols = c("Lon.km","Lat.km"), cutoff = 10) ## 8 is too low for Model 2

bspde3 <- add_barrier_mesh(
  mesh3, na_coast_proj, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)


covar <- c("scale_day","bakun33","npgo","zoop","temp.10m","dist_km_calc","oni","south.ild")
m3 <- check_mods("on", "iid", y="BSB", covar, "year", bspde3, bsb.scaled3, xt=NULL,tv=NULL,reml=TRUE)
m6 <- check_mods("off", "iid", y="BSB", covar, "year", bspde3, bsb.scaled3, xt=NULL,tv=NULL,reml=TRUE) 


fit2a.bsb <-
  sdmTMB(
    BSB ~
      scale_day + dist_km_calc + temp.10m + npgo + zoop + oni + bakun33 + south.ild,
    mesh = bspde3,
    time = "year",
    spatial = "off",
    spatiotemporal = "iid",
    family = tweedie(),
    data = bsb.scaled3,
    reml = FALSE,
    silent = FALSE
  )


## check fit ---------------------------------------------------------------
fit <- fit2a.bsb
fit
sanity(fit)
AIC(fit)

check_collinearity(fit)

tidy.coef <- tidy(fit,"fixed",conf.int = TRUE)

tidy.coef.noint <- tidy.coef %>% 
  filter(term != "(Intercept)")


tidy.coef2 <- tidy(fit1a.bsb,"fixed",conf.int = TRUE)

tidy.coef.noint2 <- tidy.coef2 %>% 
  filter(term != "(Intercept)")

tidy <- rbind(tidy.coef.noint,tidy.coef.noint2)

tidy$model <- "Time Period 1"
tidy[7:14,6] <- "Time Period 2"

cols <- c( 'Time Period 1'="#D95F02", 'Time Period 2'="#7570B3")
shapes <- c('Time Period 1'= 16, 'Time Period 2'= 17)
ONI <- expression(paste("Ocean Ni",tilde(n),"o Index"))

theme_Publication <- function(base_size=24, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(0.7),vjust = 0), #hjust = 0.5),
            # text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border =  element_blank(),#element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(0.7)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size = rel(0.6)), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major.y = element_line(colour="#f0f0f0"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.1, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            legend.text = element_text(size = rel(0.7)),
            plot.margin=unit(c(5,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


p5 <- ggplot(tidy,aes(estimate,term))+
  geom_point(aes(shape=model,size=2,color=model),alpha=0.5,
             position = position_dodge(width=0.5))+
  geom_errorbarh(aes(xmin=conf.low,xmax=conf.high,y=term,
                     shape=model,
                     color=model,
                     linewidth=0.025),
                 height=0.0,
                 alpha=0.5,
                 position = position_dodge(width=0.75))+
  scale_y_discrete(limits=c('scale_day','dist_km_calc','npgo','oni','temp.10m','zoop','bakun33','south.ild'),
                   labels=c('day of year','distance to mainland','NPGO',ONI,'surface temperature','zooplankton biomass','Bakun upwelling 33\u00b0N','isothermal layer depth'))+
  geom_vline(xintercept=0.0,linetype=2)+
  xlim(c(-6.2,2))+
  theme_Publication()+
  scale_color_manual(values = cols,
                     labels = c("Time Period 1\n(1963-2016)","Time Period 2\n(1984-2016)"))+
  scale_shape_manual(values = shapes,
                     labels = c("Time Period 1\n(1963-2016)","Time Period 2\n(1984-2016)"))+
  theme(legend.title = element_blank(),
        legend.key.spacing.y = unit(1, 'cm'))+
  guides(
    size = "none",
    color = "none",
    shape = "none",
    # color = guide_legend(override.aes = list(size = 5)),
    # shape = guide_legend(
    #   key.spacing.y = unit(1, 'cm')),
    linewidth = "none")+
  theme(legend.title = element_blank())+
  xlab("Coefficient estimate")+
  ylab("Covariate")+
  ggtitle("Barred Sand Bass")+
  labs(tag = 'A')
p5

p4 <- inset_element(p = pneb,
                    left = 0.75,
                    bottom = 0.85,
                    right = 1,
                    top = 1,
                    align_to = 'full')+
  labs(tag = '')

##k-fold CV - BSB Time Period 2-----------------------
# NOTE: no temp.10m from 198107 or 198307; also no bucket temps from those two cruises either
# thus, the earliest year is 1984


set.seed(42)
clust <- sample(1:8, size = nrow(bsb.scaled3), replace = T)
k_folds <- length(unique(clust))
mesh <- bspde3
dataset <- bsb.scaled3
covar <- c("scale_day","dist_km_calc","temp.10m","npgo","zoop","oni","bakun33","south.ild")
m1b_cv.bsb <- check_cv_mods("off", "iid", "BSB", covar, "year", mesh=mesh, dataset=dataset, tv=NULL,clust,k_folds)
m1b_bsb <- check_mods("off", "iid", y="BSB", covar, "year", mesh=mesh, dataset=dataset, xt=NULL,tv=NULL,reml=FALSE) 


# replace bakun33 w/upwelling south of Bight
dataset <- bsb.scaled %>% 
  drop_na(scale_day,npgo,zoop,temp.10m,dist_km_calc,oni,bakun30,south.ild) # upwelling south

# create new mesh based on new dataframe
mesh5 = make_mesh(dataset, xy_cols = c("Lon.km","Lat.km"), cutoff = 10) 
mesh <- add_barrier_mesh(
  mesh5, na_coast_proj, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)
set.seed(42)
clust <- sample(1:8, size = nrow(dataset), replace = T)
covar <- c("scale_day","dist_km_calc","temp.10m","npgo","zoop","oni","bakun30","south.ild")
m2b_cv.bsb <- check_cv_mods("off", "iid", "BSB", covar, "year", mesh=mesh, dataset=dataset, tv=NULL,clust,k_folds)
m2b_bsb <- check_mods("off", "iid", y="BSB", covar, "year", mesh=mesh, dataset=dataset, xt=NULL,tv=NULL,reml=FALSE) 


# replace bakun33 with beuti and cuti

# now, include south.beuti and south.ild for na drop,
dataset <- bsb.scaled %>%
  drop_na(scale_day,south.beuti,south.cuti,npgo,zoop,temp.10m,dist_km_calc,oni,south.ild)
# create new mesh based on new dataframe
mesh4 = make_mesh(dataset, xy_cols = c("Lon.km","Lat.km"), cutoff = 10)
mesh <- add_barrier_mesh(
  mesh4, na_coast_proj, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)

set.seed(42)
clust <- sample(1:8, size = nrow(dataset), replace = T)
covar <- c("scale_day","dist_km_calc","temp.10m","npgo","zoop","oni","south.beuti","south.cuti","south.ild")
m4b_cv.bsb <- check_cv_mods("off", "iid", "BSB", covar, "year", mesh=mesh, dataset=dataset, tv=NULL,clust,k_folds)

# results from model without cross-validation
m4b_bsb <- check_mods("off", "iid", y="BSB", covar, "year", mesh=mesh, dataset=dataset, xt=NULL,tv=NULL,reml=FALSE) 



# replace zoop with Calanus pacificus
# first drop nas

dataset<- bsb.scaled %>% 
  drop_na(scale_day,npgo,calpac,temp.10m,dist_km_calc,oni,bakun30,south.ild) # ratio.extent
mesh7 = make_mesh(dataset, xy_cols = c("Lon.km","Lat.km"), cutoff = 10)
mesh <- add_barrier_mesh(
  mesh7, na_coast_proj, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)

set.seed(42)
clust <- sample(1:8, size = nrow(dataset), replace = T)
covar <- c("scale_day","dist_km_calc","temp.10m","npgo","calpac","oni","bakun30","south.ild")
m5b_cv.bsb <- check_cv_mods("off", "iid", "BSB", covar, "year", mesh=mesh, dataset=dataset, tv=NULL,clust,k_folds)
m5b_bsb <- check_mods("off", "iid", y="BSB", covar, "year", mesh=mesh, dataset=dataset, xt=NULL,tv=NULL,reml=FALSE) 


# Compare models

as.data.frame(tibble(
  model = c(
    c("scale_day,zoop,dist_km_calc,temp.10m,npgo,oni,bakun33,south.ild"),
    c("scale_day,zoop,dist_km_calc,temp.10m,npgo,bakun30,oni"),
    c("scale_day,zoop,dist_km_calc,temp.10m,npgo,south.beuti,south.cuti,oni"),
    c("scale_day,calpac,dist_km_calc,temp.10m,npgo,bakun33,oni") # make sure your variables are named correctly; scale covariates where necessary
  ),
  AIC = c(
    unlist(m1b_bsb[3]),   
    unlist(m2b_bsb[3]),                                    
    unlist(m4b_bsb[3]),                                
    unlist(m5b_bsb[3])
  ),
  sum_LL = c(
    m1b_cv.bsb$sum_loglik,   
    m2b_cv.bsb$sum_loglik,   
    m4b_cv.bsb$sum_loglik,                                
    m5b_cv.bsb$sum_loglik
  ))) %>% 
  arrange(AIC)



# 04 KB Time Period 1 ------------------------------------------------------------------

kb <- bsb.scaled  # bsb.scaled contains KB data
kb$scale_day <- scale(kb$day)

# first drop nas

# NOTE:  only npgo and oni go back as far as CalCOFI data
# in Model 2 we will also include upwelling indices (bakun33, bakun30, or beuti/cuti) and isothermal layer depth for na drop
kb.scaled2 <- kb %>%
  drop_na(scale_day,npgo,zoop,temp.10m,dist_km_calc,oni)

# create new mesh based on new dataframe
map_data <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf", country = "United States of America")

na_coast <- suppressWarnings(suppressMessages(
  sf::st_crop(map_data,
              c(xmin = -122, ymin = 30, xmax = -116, ymax = 36))))

na_coast_proj <- sf::st_transform(na_coast, crs = 32611)

mesh2 = make_mesh(kb.scaled2, xy_cols = c("Lon.km","Lat.km"), cutoff = 10)

bspde2 <- add_barrier_mesh(
  mesh2, na_coast_proj, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)

# specify some function arguments

names(kb.scaled2)
covar <- c("day","dist_km_calc", "temp.10m" , "npgo" , "zoop" , "oni") # make sure your variables are named correctly; scale covariates where necessary

m2 <- check_mods("on", "iid", y="KB", covar, "year", bspde2, kb.scaled2, xt=NULL,tv=NULL,reml=TRUE) 
m3 <- check_mods("off", "iid", y="KB", covar, "year", bspde2, kb.scaled2, xt=NULL,tv=NULL,reml=TRUE)

fit1a.kb <- sdmTMB(KB ~
                  scale_day + dist_km_calc + temp.10m + npgo + zoop + oni,
                mesh = bspde2,
                time = "year", # name of time variable
                spatial = "off", # can be on / off. On estimates a shared spatial field that is the same across years
                spatiotemporal = "iid", # should be iid / off
                family = tweedie(), # could also be one of the delta_ families
                data = kb.scaled2,
                reml = FALSE,
                silent = FALSE)

## check fit  ---------------------------------------------------------------
fit <- fit1a.kb
sanity(fit)
AIC(fit)

library(performance)
check_collinearity(fit)
max(fit$gradients)
fit$sd_report
tidy(fit,"ran_pars",conf.int = TRUE)
tidy(fit,"fixed",conf.int = TRUE)


##k-fold CV - KB Time Period 1-----------------------

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


set.seed(42)
clust <- sample(1:8, size = nrow(kb.scaled2), replace = T)
k_folds <- length(unique(clust))

mesh <- bspde2
dataset <- kb.scaled2
covar <- c("scale_day","dist_km_calc","temp.10m","npgo","zoop","oni")
m1_cv.kb <- check_cv_mods("off", "iid", "KB", covar, "year", mesh=mesh, dataset=dataset, tv=NULL,clust,k_folds)
m1_kb <- check_mods("off", "iid", y="KB", covar, "year", mesh=mesh, dataset=dataset, xt=NULL,tv=NULL,reml=FALSE) 



# Same as Time Period 1 mod, but substitute zoop for calpac
dataset <- kb %>% 
  drop_na(scale_day,npgo,calpac,temp.10m,dist_km_calc,oni) # cal pac

mesh6 = make_mesh(dataset, xy_cols = c("Lon.km","Lat.km"), cutoff = 10)

mesh <- add_barrier_mesh(
  mesh6, na_coast_proj, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)

set.seed(42)
clust <- sample(1:8, size = nrow(dataset), replace = T)
k_folds <- length(unique(clust))
covar <- c("scale_day","dist_km_calc","temp.10m","npgo","calpac","oni")
m2_cv.kb <- check_cv_mods("off", "iid", "KB", covar, "year", mesh=mesh, dataset=dataset, tv=NULL,clust,k_folds)

# results from model without cross-validation
m2_kb <- check_mods("off", "iid", y="KB", covar, "year", mesh=mesh, dataset=dataset, xt=NULL,tv=NULL,reml=FALSE) 

as.data.frame(tibble(
  model = c(
    c("scale_day,zoop,dist_km_calc,temp.10m,npgo,oni"),
    c("scale_day,calpac,dist_km_calc,temp.10m,npgo,oni") # make sure your variables are named correctly; scale covariates where necessary
  ),
  AIC = c(
    unlist(m1_kb[3]),   
    unlist(m2_kb[3])
  ),
  sum_LL = c(
    m1_cv.kb$sum_loglik,   
    m2_cv.kb$sum_loglik
  ))) %>% 
  arrange(AIC)
# Compare log-likelihoods 
m1_cv.kb[["sum_loglik"]]# full mod with zoop: all folds converged
m2_cv.kb[["sum_loglik"]] # full mod with calpac: all folds converged

# 05 KB Time Period 2 ------------------------------------------------------------------

# let's incorporate other data available that relates to upwelling
# and mixed layer depth

# include bakun33/bakun30 and south.ild for na drop,
# this represents a shorter time series since south.ild is only available since 1980,
# but no CalCOFI temp data in 1981 or 1983 anyway

# call it kb.scaled3
kb.scaled3 <- kb %>%
  drop_na(scale_day,bakun33,npgo,zoop,temp.10m,dist_km_calc,oni,south.ild)

# NOTE: no temp.10m from 198107 or 198307; also no bucket temps from those two cruises either
# thus, the earliest year is 1984

# create new mesh based on new dataframe
mesh3 = make_mesh(kb.scaled3, xy_cols = c("Lon.km","Lat.km"), cutoff = 10) ## 8 is too low for Model 2

bspde3 <- add_barrier_mesh(
  mesh3, na_coast_proj, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)

fit2a.kb <-
  sdmTMB(
    KB ~  scale_day + dist_km_calc + temp.10m + npgo + zoop + oni + bakun33 + south.ild,
    mesh = bspde3,
    time = "year",
    # name of time variable
    spatial = "off",
    # can be on / off. On estimates a shared spatial field that is the same across years
    spatiotemporal = "iid",
    # should be iid / off
    family = tweedie(),
    # could also be one of the delta_ families
    data = kb.scaled3,
    reml = FALSE,
    silent = FALSE
  )


## check fit ---------------------------------------------------------------
fit <- fit2a.kb
fit
sanity(fit)
AIC(fit)

check_collinearity(fit)

tidy.coef <- tidy(fit1a.kb,"fixed",conf.int = TRUE)

tidy.coef.noint <- tidy.coef %>% 
  filter(term != "(Intercept)")


tidy.coef2 <- tidy(fit2a.kb,"fixed",conf.int = TRUE)

tidy.coef.noint2 <- tidy.coef2 %>% 
  filter(term != "(Intercept)")

tidy <- rbind(tidy.coef.noint,tidy.coef.noint2)

tidy$model <- "Time Period 1"
tidy[7:14,6] <- "Time Period 2"

cols <- c( 'Time Period 1'="#D95F02", 'Time Period 2'="#7570B3")
shapes <- c('Time Period 1'= 16, 'Time Period 2'= 17)
ONI <- expression(paste("Ocean Ni",tilde(n),"o Index"))


p9 <- ggplot(tidy,aes(estimate,term))+
  geom_point(aes(shape=model,size=2,color=model),alpha=0.5,
             position = position_dodge(width=0.5))+
  geom_errorbarh(aes(xmin=conf.low,xmax=conf.high,y=term,
                     shape=model,
                     color=model,
                     linewidth=0.025),
                 height=0.0,
                 alpha=0.5,
                 position = position_dodge(width=0.75))+
  scale_y_discrete(limits=c('scale_day','dist_km_calc','npgo','oni','temp.10m','zoop','bakun33','south.ild'),
                   # labels=c('day','shore_dist','npgo','oni','temp','zoop','upwelling','ild'))+
                   labels=c('day of year','distance to mainland','NPGO',ONI,'surface temperature','zooplankton biomass','Bakun upwelling 33\u00b0N','isothermal layer depth'))+
  geom_vline(xintercept=0.0,linetype=2)+
  xlim(c(-6,2))+
  theme_Publication()+
  scale_color_manual(values = cols,
                     labels = c("Time Period 1\n(1963-2016)","Time Period 2\n(1984-2016)"))+
  scale_shape_manual(values = shapes,
                     labels = c("Time Period 1\n(1963-2016)","Time Period 2\n(1984-2016)"))+
  theme(legend.title = element_blank(),
        legend.key.spacing.y = unit(1, 'cm'),
        legend.position = "bottom",
        legend.direction = "horizontal")+
  guides(
    size = "none",
    color = guide_legend(override.aes = list(size = 5)),
    # shape = guide_legend(
    #   key.spacing.y = unit(1, 'cm')),
    linewidth = "none")+
  theme(legend.title = element_blank())+
  xlab("Coefficient estimate")+
  ylab("Covariate")+
  ggtitle("Kelp Bass")+
  labs(tag = 'B')
p9


p8 <-   inset_element(p = pcla,
                      left = 0.75,
                      bottom = 0.85,
                      right = 1,
                      top = 1,
                      align_to = 'full')+
  labs(tag = '')

##k-fold CV - KB Time Period 2-----------------------
set.seed(42)
clust <- sample(1:8, size = nrow(kb.scaled3), replace = T)
k_folds <- length(unique(clust))

mesh <- bspde3
dataset <- kb.scaled3
covar <- c("scale_day","dist_km_calc","temp.10m","npgo","zoop","oni","bakun33","south.ild")
m1b_cv.kb <- check_cv_mods("off", "iid", "KB", covar, "year", mesh=mesh, dataset=dataset, tv=NULL,clust,k_folds)
m1b_kb <- check_mods("off", "iid", y="KB", covar, "year", mesh=mesh, dataset=dataset, xt=NULL,tv=NULL,reml=FALSE) 

# replace bakun33 with kelp

dataset <- kb %>% 
  drop_na(scale_day,npgo,zoop,temp.10m,dist_km_calc,oni,south.ild,area) # ratio.extent
mesh8 = make_mesh(dataset, xy_cols = c("Lon.km","Lat.km"), cutoff = 10)
mesh <- add_barrier_mesh(
  mesh8, na_coast_proj, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)

set.seed(42)
clust <- sample(1:8, size = nrow(dataset), replace = T)
k_folds <- length(unique(clust))
covar <- c("scale_day","zoop","dist_km_calc", "temp.10m" , "npgo" , "area" , "oni") # make sure your variables are named correctly; scale covariates where necessary
m2b_cv.kb <- check_cv_mods("off", "iid", "KB", covar, "year", mesh=mesh, dataset=dataset, tv=NULL,clust,k_folds)
m2b_kb <- check_mods("off", "iid", y="KB", covar, "year", mesh=mesh, dataset=dataset, xt=NULL,tv=NULL,reml=FALSE) 


# replace bakun33 with upwelling south of Bight

dataset <- kb %>% 
  drop_na(scale_day,npgo,zoop,temp.10m,dist_km_calc,oni,bakun30,south.ild) # upwelling south
mesh5 = make_mesh(dataset, xy_cols = c("Lon.km","Lat.km"), cutoff = 10)
mesh <- add_barrier_mesh(
  mesh5, na_coast_proj, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)

set.seed(42)
clust <- sample(1:8, size = nrow(dataset), replace = T)
k_folds <- length(unique(clust))
covar <- c("scale_day","zoop","dist_km_calc", "temp.10m" , "npgo" , "bakun30" , "oni") # make sure your variables are named correctly; scale covariates where necessary
m3b_cv.kb <- check_cv_mods("off", "iid", "KB", covar, "year", mesh=mesh, dataset=dataset, tv=NULL,clust,k_folds)
m3b_kb <- check_mods("off", "iid", y="KB", covar, "year", mesh=mesh, dataset=dataset, xt=NULL,tv=NULL,reml=FALSE) 


# now replace Bakun upwelling with south.beuti and south.ild
dataset <- kb %>%
  drop_na(scale_day,south.beuti,south.cuti,npgo,zoop,temp.10m,dist_km_calc,oni,south.ild)

# create new mesh based on new dataframe
mesh4 = make_mesh(dataset, xy_cols = c("Lon.km","Lat.km"), cutoff = 10)
mesh <- add_barrier_mesh(
  mesh4, na_coast_proj, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)

set.seed(42)
clust <- sample(1:8, size = nrow(dataset), replace = T)
k_folds <- length(unique(clust))
covar <- c("scale_day","zoop","dist_km_calc", "temp.10m" , "npgo" , "south.beuti", "south.cuti" , "oni") # make sure your variables are named correctly; scale covariates where necessary
m4b_cv.kb <- check_cv_mods("off", "iid", "KB", covar, "year", mesh=mesh, dataset=dataset, tv=NULL,clust,k_folds)
m4b_kb <- check_mods("off", "iid", y="KB", covar, "year", mesh=mesh, dataset=dataset, xt=NULL,tv=NULL,reml=FALSE) 


# keep bakun33 and add kelp extent
dataset<- kb %>% 
  drop_na(scale_day,npgo,zoop,temp.10m,dist_km_calc,oni,south.ild,bakun33,area) # ratio.extent
mesh10 = make_mesh(dataset, xy_cols = c("Lon.km","Lat.km"), cutoff = 10)
mesh <- add_barrier_mesh(
  mesh10, na_coast_proj, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)

set.seed(42)
clust <- sample(1:8, size = nrow(dataset), replace = T)
k_folds <- length(unique(clust))
covar <- c("scale_day","zoop","dist_km_calc", "temp.10m" , "npgo" , "bakun33", "area" , "oni") # make sure your variables are named correctly; scale covariates where necessary
m5b_cv.kb <- check_cv_mods("off", "iid", "KB", covar, "year", mesh=mesh, dataset=dataset, tv=NULL,clust,k_folds)
m5b_kb <- check_mods("off", "iid", y="KB", covar, "year", mesh=mesh, dataset=dataset, xt=NULL,tv=NULL,reml=FALSE) 


# substitute zoop for calpac

dataset <- kb %>% 
  drop_na(scale_day,npgo,calpac,temp.10m,dist_km_calc,oni,bakun33,south.ild) # ratio.extent
mesh7 = make_mesh(dataset, xy_cols = c("Lon.km","Lat.km"), cutoff = 10)
mesh <- add_barrier_mesh(
  mesh7, na_coast_proj, range_fraction = 0.1,
  proj_scaling = 1000, plot = TRUE
)

set.seed(42)
clust <- sample(1:8, size = nrow(dataset), replace = T)
k_folds <- length(unique(clust))
covar <- c("scale_day","calpac","dist_km_calc", "temp.10m" , "npgo" , "bakun33", "oni") # make sure your variables are named correctly; scale covariates where necessary
m6b_cv.kb <- check_cv_mods("off", "iid", "KB", covar, "year", mesh=mesh, dataset=dataset, tv=NULL,clust,k_folds)
m6b_kb <- check_mods("off", "iid", y="KB", covar, "year", mesh=mesh, dataset=dataset, xt=NULL,tv=NULL,reml=FALSE) 


# Compare models

as.data.frame(tibble(
  model = c(
    c("scale_day,zoop,dist_km_calc,temp.10m,npgo,oni,bakun33,south.ild"),
    c("scale_day,zoop,dist_km_calc,temp.10m,npgo,area,oni"),
    c("scale_day,zoop,dist_km_calc,temp.10m,npgo,bakun30,oni"),
    c("scale_day,zoop,dist_km_calc,temp.10m,npgo,south.beuti,south.cuti,oni"),
    c("scale_day,zoop,dist_km_calc,temp.10m,npgo,bakun33,area,oni"),
    c("scale_day,calpac,dist_km_calc,temp.10m,npgo,bakun33,oni") # make sure your variables are named correctly; scale covariates where necessary
  ),
  AIC = c(
    unlist(m1b_kb[3]),   
    unlist(m2b_kb[3]),   
    unlist(m3b_kb[3]),                                      
    unlist(m4b_kb[3]),                                
    unlist(m5b_kb[3]),
    unlist(m6b_kb[3])
    ),
  sum_LL = c(
    m1b_cv.kb$sum_loglik,   
    m2b_cv.kb$sum_loglik,   
    m3b_cv.kb$sum_loglik,                                      
    m4b_cv.kb$sum_loglik,                                
    m5b_cv.kb$sum_loglik,
    m6b_cv.kb$sum_loglik
  ))) %>% 
  arrange(AIC)

# 06 Figure 7-------------------------------------------------

(p5 + p4) / (p9 + p8)


# 07 Supplemental Figs 12-15 ----------------------------------------
theme_Publication <- function(base_size=18, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(0.7),vjust = 0), #hjust = 0.5),
            # text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border =  element_blank(),#element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(0.7)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size = rel(0.6)), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major.y = element_line(colour="#f0f0f0"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.1, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            legend.text = element_text(size=12),
            plot.margin=unit(c(5,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

# BSB Time Period 1
fit <- fit1a.bsb

par(mfrow=c(2,3))
plot_list <- visreg(fit,
                    line = list(col = "red"),
                    points = list(cex = 1, pch = 1),
                    gg=TRUE
)

library(patchwork)

h<-plot_list[[1]] + theme_Publication()
i<-plot_list[[2]] + theme_Publication()
j<-plot_list[[3]] + theme_Publication()
k<-plot_list[[4]] + theme_Publication()
l<-plot_list[[5]] + theme_Publication()
m<-plot_list[[6]] + theme_Publication()


BSB_TimePeriod1 <- h+i+j+k+l+m + plot_layout(widths = c(1, 1, 1),
                                             heights = c(1, 1, 1))
BSB_TimePeriod1

# create a df

df.points <- NULL
df.points$x <- as.vector(c(0))
df.points$y <- as.vector(c(-3.0))
df.points <- as.data.frame(df.points)
df.points$image <- pneb.magick


# create text/image plots
n <- ggplot(df.points, aes(x, y)) +
  geom_image(aes(image=image), size=1)+
  theme_void()+
  annotate("text", x=0, y=2.5, 
           fontface = "bold",
           # size = 12,
           label= "Barred Sand Bass\nTime Period 1\n(1963-2016)")+
  ylim(c(-5,5))

# pull it all together
h+plot_spacer()+ plot_spacer() + i+j+k+l+ m + n +plot_layout(widths = c(1, 1, 1),
                                                             heights = c(1, 1, 1))



# BSB Time Period 2
fit <- fit2a.bsb

test <- visreg(fit,
               line = list(col = "red"),
               points = list(cex = 1, pch = 1),
               gg=TRUE
)

library(patchwork)
a<-test[[1]] + theme_Publication()
b<-test[[2]] + theme_Publication()
c<-test[[3]] + theme_Publication()
d<-test[[4]] + theme_Publication()
e<-test[[5]] + theme_Publication()
f<-test[[6]] + theme_Publication()
g<-test[[7]] + theme_Publication()
h<-test[[8]] + theme_Publication()


BSB_TimePeriod2 <- a+b+c+d+e+f+g+h + plot_layout(ncol=3)
BSB_TimePeriod2

# create a df

df.points <- NULL
df.points$x <- as.vector(c(0))
df.points$y <- as.vector(c(-3.0))
df.points <- as.data.frame(df.points)
df.points$image <- pneb.magick


# create text/image plot
n <- ggplot(df.points, aes(x, y)) +
  geom_image(aes(image=image), size=1)+
  theme_void()+
  annotate("text", x=0, y=2.5, 
           fontface = "bold",
           # size = 12,
           label= "Barred Sand Bass\nTime Period 2\n(1984-2016)")+
  ylim(c(-5,5))

# pull it all together
BSB_TimePeriod2 + n + plot_layout(widths = c(1, 1, 1),
                                  heights = c(1, 1, 1))

# KB Time Period 1
fit <- fit1a.kb
plot_list <- visreg(fit,
                    line = list(col = "red"),
                    points = list(cex = 1, pch = 1),
                    gg=TRUE
)

h<-plot_list[[1]] + theme_Publication()
i<-plot_list[[2]] + theme_Publication()
j<-plot_list[[3]] + theme_Publication()
k<-plot_list[[4]] + theme_Publication()
l<-plot_list[[5]] + theme_Publication()
m<-plot_list[[6]] + theme_Publication()

# create a df

df.points <- NULL
df.points$x <- as.vector(c(0))
df.points$y <- as.vector(c(-3.0))
df.points <- as.data.frame(df.points)
df.points$image <- pcla.magick


# create text/image plot
n <- ggplot(df.points, aes(x, y)) +
  geom_image(aes(image=image), size=1)+
  theme_void()+
  annotate("text", x=0, y=2.5, 
           fontface = "bold",
           # size = 12,
           label= "Kelp Bass\nTime Period 1\n(1963-2016)")+
  ylim(c(-5,5))

# pull it all together                                                         
h+plot_spacer()+ plot_spacer() + i+j+k+l+ m + n +plot_layout(widths = c(1, 1, 1),
                                                             heights = c(1, 1, 1))

# KB Time Period 2
fit <- fit2a.kb
plot_list2 <- visreg(fit,
                     line = list(col = "red"),
                     points = list(cex = 1, pch = 1),
                     gg=TRUE
)

a<-plot_list2[[1]] + theme_Publication()
b<-plot_list2[[2]] + theme_Publication()
c<-plot_list2[[3]] + theme_Publication()
d<-plot_list2[[4]] + theme_Publication()
e<-plot_list2[[5]] + theme_Publication()
f<-plot_list2[[6]] + theme_Publication()
g<-plot_list2[[7]] + theme_Publication()
z<-plot_list2[[8]] + theme_Publication()


# create text/image plot
n <- ggplot(df.points, aes(x, y)) +
  geom_image(aes(image=image), size=1)+
  theme_void()+
  annotate("text", x=0, y=2.5, 
           fontface = "bold",
           # size = 12,
           label= "Kelp Bass\nTime Period 2\n(1984-2016)")+
  ylim(c(-5,5))

# pull it all together
a+b+c+d+e+f+g+z+ n + plot_layout(ncol=3,
                                 nrow=3,
                                 widths = c(1, 1, 1),
                                 heights = c(1, 1, 1))
