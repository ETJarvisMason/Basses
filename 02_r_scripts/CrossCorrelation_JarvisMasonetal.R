# Jarvis Mason, ET
# 29 OCT 2024



# load libraries
library(funtimes)
library(tidyverse)
library(sdmTMB)
library(sf)
library(marmap)
library(patchwork)

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
            legend.position = "top",
            legend.direction = "vertical",
            legend.key.size= unit(0.1, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(5,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

# Import fish drawings
library("png")                      # Load png package

pneb <- readPNG("03_plots/bsb.png", native = TRUE)
pcla <- readPNG("03_plots/kb.png", native = TRUE)

# 01 load data---------------------------------------------------
# load sdm_TMB indices for BSB
load(file="01_tidy_data/output_SIA_BSB_with_missing_yrs.RData")

# bring in scaled RecFIN estimates (total estimates)
est.scaled <- read.csv("01_tidy_data/est.scaled.csv")

# bring in CPFV harvest dat
bsbCPFV <- read.csv("00_raw_data/BSB_EstCPFVLands1947to2021.csv")
harv <- bsbCPFV %>% 
  filter(Year>1974) %>% 
  dplyr::rename(year = Year,
                bsan.c = BSB_CPFVlands) %>% 
  mutate(kbas.c = RockbassCPFVLands - bsan.c)

# scale data
index.scaled[,2] <- scale(index[,5])

index.scaled <- index.scaled[,1:2]         # we only want the first two cols
colnames(index.scaled)[2] <- "bsan.l"


# run cross correlation with RecFIN estimates (total catch [includes releases])
# run cross correlation with CPFV logbook data (total harvest)




# 02 Cross-correlation Analysis------------------
## Barred Sand Bass --------------------------------------------------------------------

# sand bass:     (bsan.c,bsan.l), (bsan.sum_catch,bsan.l)
ccm_dat.scaled.total <- full_join(index.scaled,est.scaled)

ccm_dat.scaled.harv <- full_join(index.scaled,harv)
ccm_dat.scaled.harv[,7:8]  <- scale(ccm_dat.scaled.harv[,7:8])

# assign column 1 to be the fishing data of interest (bsan.sum_catch or bsan.c)
col1 <- "bsan.sum_catch"
col2 <- "bsan.l"

# update dataframe of interest with appropriate col names
ccm_dat.scaled.total <- as.data.frame(ccm_dat.scaled.total) %>% 
  dplyr::select(year,!!as.symbol(col1),!!as.symbol(col2))


# total catch
# assign your df of interest as 'd'
d <- ccm_dat.scaled.total
attributes(d$bsan.l) <- NULL

# impute missing data
library(mice)
ccf.data <- mice(data=d,m=2,method="pmm",maxit=50,seed=500)

completedData <- complete(ccf.data,1) 

d <- completedData

# total estimates data available from 1980 (larval data go to 2016)
d <- as.data.frame(d) %>%
  filter(year > 1979 & year < 2017) 

# assign your fishing data type as 'd_source'
d_source <- "RecFIN"
sp <- "Barred Sand Bass"
pic <- pneb

###Figure 4A------------------------------------
d_source <- d_source  
sp <- sp
pic <- pic
d <- d
colnames(d) <- c("year","catch","larvae")
larvae <- d$larvae 
catch <- d$catch 


# par(mfrow=c(1,1))
ccf.boot <- ccf_boot(as.numeric(d$larvae),
                     as.numeric(d$catch),lag.max=10)

# total catch
ccf.plot.total <- ggplot(ccf.boot,aes(Lag,r_S))+  
  geom_ribbon(aes(x=Lag,ymin=lower_S,ymax=upper_S,alpha=0.3),fill="skyblue")+
  geom_point(aes(color = abs(r_S), size = 1))+
  # scale_color_viridis_b(option = "viridis", 
                        # limits = c(0,0.6),
                        # direction = -1)+
  scale_color_gradient(low = "yellow", high = "darkred",limits = c(0,0.55)) +
  # geom_segment(aes(x=Lag, xend=Lag, y=0, yend=r_S)) +
  # geom_vline(xintercept=0.0,color="red")+
  geom_hline(yintercept=0.0,linetype=2)+
  scale_x_continuous(limits = c(-10,0),breaks = seq(-10,0,2))+
  ylim(-0.4,0.6)+
  xlab("Lag")+
  ylab("Coefficient")+
  labs(subtitle = "Total catch")+
  theme_Publication()+
  # theme(legend.position = "") 
  labs(tag = 'A')+
  theme(legend.position = "")
ccf.plot.total

inset <-                  # Combine plot & image
  inset_element(p = pic,
                left = 0.65,
                bottom = 0.75,
                right = 1,
                top = 1,
                align_to = 'full')+
  labs(tag = '')

# harvest
col1 <- "bsan.c"
col2 <- "bsan.l"
ccm_dat.scaled.harv <- as.data.frame(ccm_dat.scaled.harv) %>% 
  dplyr::select(year,!!as.symbol(col1),!!as.symbol(col2))


d_source <- "CPFV"
d <- ccm_dat.scaled.harv
attributes(d$bsan.l) <- NULL

d <- as.data.frame(d) %>%
  filter(year > 1974 & year < 2017) %>% 
  dplyr::select(year,bsan.c,bsan.l)

###Figure 4B------------------------------------
d_source <- d_source  
sp <- sp
pic <- pic
d <- d
colnames(d) <- c("year","catch","larvae")
larvae <- d$larvae 
catch <- d$catch 

# par(mfrow=c(1,1))
ccf.boot <- ccf_boot(as.numeric(d$larvae),
                     as.numeric(d$catch),lag.max=10)

ccf.plot.harv <- ggplot(ccf.boot,aes(Lag,r_S))+
  geom_ribbon(aes(x=Lag,ymin=lower_S,ymax=upper_S,alpha=0.3),fill="skyblue")+
  geom_point(aes(color = abs(r_S), size = 1.75))+
  scale_color_gradient(low = "yellow", high = "darkred",limits = c(0,0.55)) +
  # geom_segment(aes(x=Lag, xend=Lag, y=0, yend=r_S)) +
  # geom_vline(xintercept=0.0,color="red")+
  geom_hline(yintercept=0.0,linetype=2)+
  scale_x_continuous(limits = c(-10,0),breaks = seq(-10,0,2))+
  ylim(-0.4,0.6)+
  xlab("Lag")+
  ylab("Coefficient")+
  labs(subtitle = "Harvest")+
  theme_Publication()+
  # theme(legend.position = "") 
  labs(tag = 'B')+
  theme(legend.position = "")
ccf.plot.harv

BSB <- ccf.plot.total + ccf.plot.harv + inset +
  plot_annotation(title = 'Barred Sand Bass',
                  theme = theme(plot.title = element_text(size = 18, hjust = 0.5))) # Draw combined plot
BSB      

## Kelp Bass --------------------------------------------------------------------

# load sdm_TMB indices for KB
load(file="01_tidy_data/output_SIA_KB_with_missing_yrs.RData")

# scale data
index.scaled[,2] <- scale(index[,5])

index.scaled <- index.scaled[,1:2]         # we only want the first two cols
colnames(index.scaled)[2] <- "kbas.l"
index.scaled$kbas.l <- as.vector(index.scaled$kbas.l)
str(index.scaled)

# kelp bass:     (kbas.c,kbas.l), (kbas.sum_catch,kbas.l)
ccm_dat.scaled.total <- full_join(index.scaled,est.scaled)

# assign column 1 to be the fishing data of interest (bsan.sum_catch or bsan.c)
col1 <- "kbas.sum_catch"
col2 <- "kbas.l"

# update dataframe of interest with appropriate col names
ccm_dat.scaled.total <- as.data.frame(ccm_dat.scaled.total) %>% 
  dplyr::select(year,!!as.symbol(col1),!!as.symbol(col2))



# total catch
# assign your df of interest as 'd'
d <- as.data.frame(ccm_dat.scaled.total)

attributes(d$kbas.l) <- NULL
d <- as.data.frame(d) %>%
  filter(year > 1979 & year < 2017)

library(mice)
ccf.data <- mice(data=d,m=2,method="pmm",maxit=50,seed=500)

completedData <- complete(ccf.data,1) 

d <- completedData

# assign your fishing data type as 'd_source'
d_source <- "RecFIN"

sp <- "Kelp Bass"
pic <- pcla


###Figure 4C------------------------------------
d_source <- d_source  
sp <- sp
pic <- pic
d <- d
colnames(d) <- c("year","catch","larvae")
larvae <- d$larvae 
catch <- d$catch 

# par(mfrow=c(1,1))
ccf.boot <- ccf_boot(as.numeric(d$larvae),
                     as.numeric(d$catch),lag.max=10)


ccf.plot.total <- ggplot(ccf.boot,aes(Lag,r_S))+  
  geom_ribbon(aes(x=Lag,ymin=lower_S,ymax=upper_S,alpha=0.3),fill="skyblue")+
  geom_point(aes(color = abs(r_S), size = 1))+
  scale_color_gradient(low = "yellow", high = "darkred",limits = c(0,0.55)) +
  # geom_segment(aes(x=Lag, xend=Lag, y=0, yend=r_S)) +
  # geom_vline(xintercept=0.0,color="red")+
  geom_hline(yintercept=0.0,linetype=2)+
  scale_x_continuous(limits = c(-10,0),breaks = seq(-10,0,2))+
  ylim(-0.4,0.6)+
  xlab("Lag")+
  ylab("Coefficient")+
  labs(subtitle = "Total catch")+
  theme_Publication()+
  # theme(legend.position = "") 
  labs(tag = 'C')+
  theme(legend.position = "")
ccf.plot.total

inset <-                  # Combine plot rsignrank()inset <-                  # Combine plot & image
  inset_element(p = pic,
                left = 0.65,
                bottom = 0.75,
                right = 1,
                top = 1,
                align_to = 'full')+
  labs(tag = '')

ccf.plot.total

###Figure 4D------------------------------------

# harvest

ccm_dat.scaled.harv <- full_join(index.scaled,harv)
col1 <- "kbas.c"
ccm_dat.scaled.harv <- as.data.frame(ccm_dat.scaled.harv) %>% 
  dplyr::select(year,!!as.symbol(col1),!!as.symbol(col2))

d <- ccm_dat.scaled.harv

d <- as.data.frame(d) %>%
  dplyr::filter(year > 1974 & year < 2017) %>% 
  dplyr::select(year,kbas.c,kbas.l)

# assign your fishing data type as 'd_source'

d_source <- "CPFV"

d_source <- d_source  
sp <- sp
pic <- pic
d <- d
colnames(d) <- c("year","catch","larvae")
larvae <- d$larvae 
catch <- d$catch 

# par(mfrow=c(1,1))
ccf.boot <- ccf_boot(as.numeric(d$larvae),
                     as.numeric(d$catch),lag.max=10)

ccf.plot.harv <- ggplot(ccf.boot,aes(Lag,r_S))+# ccf.plot.harv <- ggplot(ccf.boot,aes(Lag,r_S))+
  geom_ribbon(aes(x=Lag,ymin=lower_S,ymax=upper_S,alpha=0.3),fill="skyblue")+
  geom_point(aes(color = abs(r_S), size = 1))+
  scale_color_gradient(low = "yellow", high = "darkred",limits = c(0,0.55)) +
  # geom_segment(aes(x=Lag, xend=Lag, y=0, yend=r_S)) +
  # geom_vline(xintercept=0.0,color="red")+
  geom_hline(yintercept=0.0,linetype=2)+
  scale_x_continuous(limits = c(-10,0),breaks = seq(-10,0,2))+
  ylim(-0.4,0.6)+
  xlab("Lag")+
  ylab("Coefficient")+
  labs(subtitle = "Harvest")+
  theme_Publication()+
  # theme(legend.position = "") +
  labs(tag = 'D')+
  theme(legend.position = "")
ccf.plot.harv

KB <- ccf.plot.total + ccf.plot.harv + inset+
  plot_annotation(title = 'Kelp Bass',
                  theme = theme(plot.title = element_text(size = 18, hjust = 0.5))) # Draw combined plot
KB      


# 03 Figure 4------------------------------------
BSB / KB

# 04 Supp Fig. S3 catch trends -------------------------------------------------------
theme_Publication <- function(base_size=16, base_family="Helvetica") {
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
            panel.grid.major.y = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "top",
            legend.direction = "vertical",
            legend.key.size= unit(0.1, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(5,2,2,2),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}


# bring in un-scaled RecFIN estimates (total estimates)
estimates <- read.csv("01_tidy_data/est.NOTscaled.csv")

# bring in CPFV harvest dat
# Import original index data from back in March

bsbCPFV <- read.csv("00_raw_data/BSB_EstCPFVLands1947to2021.csv")
harv <- bsbCPFV %>% 
  filter(Year>1974) %>% 
  rename(year = Year,
         bsan.c = BSB_CPFVlands) %>% 
  mutate(kbas.c = RockbassCPFVLands - bsan.c)


# Barred Sand Bass
d.a <- estimates
d.a <- as.data.frame(d.a) %>%
  filter(year > 1979 & year < 2017) %>% 
  dplyr::select(year,bsan.sum_catch) %>% 
  mutate(bsan.est = bsan.sum_catch/1000) %>% 
  complete(year = 1980:2016, 
           fill = list(bsan.est = NA_real_))

mean.est.a <- mean(d.a$bsan.est,na.rm = TRUE)
std.est.a <- sd(d.a$bsan.est,na.rm = TRUE)

bsb.c <- ggplot()+
  geom_point(data=d.a,aes(x=year,y=bsan.est),shape=21)+
  geom_line(data=d.a,aes(x=year,y=bsan.est))+
  scale_x_continuous(limits = c(1970,2020),
                     breaks=seq(1970,2020,10))+
  scale_y_continuous("Thousands of fish",
                     sec.axis = sec_axis(~ (. - mean.est.a)/std.est.a, name = "Deviation"))+
  geom_hline(yintercept = mean.est.a,linetype="dotted")+
  theme_Publication()+
  xlab("Year")  +
  labs(subtitle = "Total catch")+
  theme(legend.position = "") +
  labs(tag = 'A')
# theme(legend.position = "") 
bsb.c

d.b <- harv
d.b <- as.data.frame(d.b) %>%
  filter(year < 2017)

mean.harv.b <- mean(d.b$bsan.c)
std.harv.b <- sd(d.b$bsan.c)

bsb.h <- ggplot()+
  geom_point(data=d.b,aes(x=year,y=bsan.c),shape=21)+
  geom_line(data=d.b,aes(x=year,y=bsan.c))+
  scale_y_continuous(limits = c(0,900),
                     breaks=seq(0,900,200),
                     sec.axis = sec_axis(~ (. - mean.harv.b)/std.harv.b, name = "Deviation"))+
  scale_x_continuous(limits = c(1970,2020),
                     breaks=seq(1970,2020,10))+
  geom_hline(yintercept = mean.harv.b,linetype="dotted")+
  theme_Publication()+
  ylab("Thousands of fish")+
  xlab("Year") +
  labs(subtitle = "Harvest")+
  theme(legend.position = "") +
  labs(tag = 'B')
# theme(legend.position = "") 


inset.bsb <-                  # Combine plot & image
  inset_element(p = pneb,
                left = 0.65,
                bottom = 0.75,
                right = 1,
                top = 1,
                align_to = 'full')+
  labs(tag = '')

BSB.catch <- bsb.c + bsb.h + plot_annotation(title = 'Barred Sand Bass',
                                             theme = theme(plot.title = element_text(size = 24, hjust = 0.5))) # Draw combined plot
BSB.catch + inset.bsb


# Kelp Bass
d.c <- estimates
d.c <- as.data.frame(d.c) %>%
  filter(year > 1979 & year < 2017) %>% 
  mutate(kbas.est = kbas.sum_catch/1000) %>% 
  dplyr::select(year,kbas.est)%>% 
  complete(year = 1980:2016, 
           fill = list(kbas.est = NA_real_))

mean.est.c <- mean(d.c$kbas.est,na.rm = TRUE)
std.est.c <- sd(d.c$kbas.est,na.rm = TRUE)


kb.c <- ggplot()+
  geom_point(data=d.c,aes(x=year,y=kbas.est),shape=21)+
  geom_line(data=d.c,aes(x=year,y=kbas.est))+
  scale_x_continuous(limits = c(1970,2020),
                     breaks=seq(1970,2020,10))+
  scale_y_continuous("Thousands of fish",
                     sec.axis = sec_axis(~ (. - mean.est.c)/std.est.c, name = "Deviation"))+
  geom_hline(yintercept = mean.est.c,linetype="dotted")+
  theme_Publication()+
  ylab("Deviation")+
  xlab("Year")  +
  labs(subtitle = "Total catch")+
  theme(legend.position = "") +
  labs(tag = 'C')
# theme(legend.position = "") 


d.d <- harv
d.d <- as.data.frame(d.d) %>%
  filter(year < 2017) 

mean.harv.d <- mean(d.d$kbas.c)
std.harv.d <- sd(d.d$kbas.c)


kb.h <- ggplot()+
  geom_point(data=d.d,aes(x=year,y=kbas.c),shape=21)+
  geom_line(data=d.d,aes(x=year,y=kbas.c))+
  scale_x_continuous(limits = c(1970,2020),
                     breaks=seq(1970,2020,10))+
  scale_y_continuous("Thousands of fish",
                     sec.axis = sec_axis(~ (. - mean.harv.d)/std.harv.d, name = "Deviation"))+
  geom_hline(yintercept = mean.harv.d,linetype="dotted")+
  theme_Publication()+
  ylab("Deviation")+
  xlab("Year") +
  labs(subtitle = "Harvest")+
  theme(legend.position = "") +
  labs(tag = 'D')



inset.kb <-                  # Combine plot & image
  inset_element(p = pcla,
                left = 0.65,
                bottom = 0.75,
                right = 1,
                top = 1,
                align_to = 'full')+
  labs(tag = '')

KB.catch <- kb.c + kb.h + plot_annotation(title = 'Kelp Bass',
                                          theme = theme(plot.title = element_text(size = 24, hjust = 0.5))) # Draw combined plot

(BSB.catch + inset.bsb) / (KB.catch + inset.kb)

ggsave(file="03_plots/Fig_S3_CatchTrends.png",
       dpi = 600,
       width=20,
       height=20,
       units="cm")
