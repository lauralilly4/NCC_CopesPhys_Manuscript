########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 1b: Plot nMDS scores - yearly timeseries
##  Laura E. Lilly
##  Updated: 10 Mar 2024
########################################
# From saved file ('Lilly_etal_NCC_CopePhys_S1_nMDS):
#   -> plot yearly timeseries of nMDS scores


library(tidyverse)
library(lubridate)
library(mgcv)
library(tidygam)
library(tidymv)


# Open datafile
# scrfl <- read.csv('NH05_Cope_biom_MDSscore_v4_CAM_RawDts.csv')
scrfl <- read.csv('NH05_CopeDens_log10_nMDSscr_Samp_k2.csv')
colnames(scrfl) <- c("Samp_Date","NMDS1","NMDS2")


# ### Step 0: Convert file days -> R Dates
scrs_tib <- as_tibble(scrfl) |>
  mutate(Samp_Date = as.Date(Samp_Date))


# ### Step 1: Create daily-resolution time structure
stdt <- as.Date(paste(min(year(scrs_tib$Samp_Date)),01,01,sep="-"),format="%Y-%m-%d")
endt <- as.Date(paste(max(year(scrs_tib$Samp_Date)),12,31,sep="-"),format="%Y-%m-%d")
mstrdt <- seq(stdt,endt,by="day")

# Remove Leap Days -> to make every timeseries the same length
lpdt <- format(as.Date(paste(2012,02,29,sep='-')),"%m-%d")
mstrmd <- format(mstrdt,"%m-%d")
lpidx <- which(mstrmd %in% lpdt)
mstrdt <- mstrdt[-lpidx]


# ### Step 2: Slot actual copepod data points into larger timeseries
# ***If a sample exists for a Leap Day, pause to deal with it
mstrdayall <- as_tibble(data.frame(mstrdt,matrix(nrow=length(mstrdt),ncol=ncol(scrs_tib)-1)))
sppidx <- vector(length=nrow(scrs_tib))

for(dn in 1:nrow(scrs_tib)){
  mstrvec = mstrdayall[,1, drop = TRUE]
  scrval = scrs_tib[dn,1, drop = TRUE]
  
  midx = which(mstrvec == scrval)
  if(length(midx) == 0){
    print(scrs_tib[dn,1])
    continue = FALSE
  } else {
    continue = TRUE
  }
  mstrdayall[midx,] = scrs_tib[dn,]
  sppidx[dn] = midx
}
names(mstrdayall) =  names(scrs_tib)



# ### Step 3: Reconfigure nMDS scores into yearly timeseries
mstrdts_tib <- as_tibble(mstrdayall) |>
  mutate(DOY = yday(Samp_Date),
         pltyr = year(Samp_Date),
         st_doy = yday(ymd(paste0(year(Samp_Date),"01_01"))),
         en_doy = yday(ymd(paste0(year(Samp_Date),"12_31")))
         )

yrsplt <- unique(mstrdts_tib$pltyr)


# ### Step 3b: Calculate various models: LM (timeseries), GAMs (yearly chunks)
lm1 <- lm(NMDS1 ~ Samp_Date, data = mstrdayall)

lm2 <- lm(NMDS2 ~ Samp_Date, data = mstrdayall)



# Calculate GAMs for NMDS1 & 2 -> for plotting
# NMDS1
gm1 <- gam(NMDS1 ~ s(DOY,bs='cc'), data = mstrdts_tib)
model_gm1 <- predict_gam(gm1)
# Rename columns of model so they can be recognized for plotting
colnames(model_gm1) <- c("DOY","fit","se.fit","lower_ci","upper_ci")

# NMDS2
gm2 <- gam(NMDS2 ~ s(DOY,bs='cc'), data = mstrdts_tib)
model_gm2 <- predict_gam(gm2)
colnames(model_gm2) <- c("DOY","fit","se.fit","lower_ci","upper_ci")






#############################
# #########  PLOTS  ######### 
# ### Step 4a: Plot entire 25-yr timeseries -> for trends
# 4.1: NMDS1
pltx1a <- ggplot(mstrdayall, aes(x = Samp_Date, y = NMDS1)) +
  
  geom_point(aes(x = Samp_Date, y = NMDS1), size = 0.7, color = 'grey40') +
  # geom_line(aes(x = Samp_Date, y = NMDS1), lwd = 0.5, linetype = "dotted") +
  geom_smooth(method = "lm") +

  annotate("text", label = "Adj. R^2 = 0.01", x = mstrdayall$Samp_Date[2000], 
           y = -0.8, size = 4, colour = "black") +
  annotate("text", label = "p < 0.01", x = mstrdayall$Samp_Date[2000], 
           y = -0.9, size = 4, colour = "black") +
  
  theme(axis.text = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "bottom", 
        axis.title = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.key=element_blank()) 

# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Plots_v4/PX1_nMDS1_LoadsTimeser.png", plot = pltx1a, width = 2000, height = 1000, units = 'px')

  
# 4.2: NMDS2
pltx2a <- ggplot(mstrdayall, aes(x = Samp_Date, y = NMDS2)) +
  
  geom_point(aes(x = Samp_Date, y = NMDS2), size = 0.7, color = 'grey40') +
  geom_smooth(method = "lm") +

  annotate("text", label = "Adj. R^2 = 0.03", x = mstrdayall$Samp_Date[8200], 
           y = 0.7, size = 4, colour = "black") +
  annotate("text", label = "p < 0.001", x = mstrdayall$Samp_Date[8200], 
           y = 0.6, size = 4, colour = "black") +
  
  theme(axis.text = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "bottom", 
        axis.title = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.key=element_blank()) 

# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Plots_v4/PX1_nMDS2_LoadsTimeser.png", plot = pltx2a, width = 2000, height = 1000, units = 'px')






# ### Step 4b: Plot all year-chunks of scores
## Color Scheme  -> four categories
symc <- c(15,16,17,18,12,3,4,8)
colel <- "orangered"
colwm <- "orange2"
colcd <- "skyblue"
colla <- "royalblue"
colnu <- "grey50"

# Set up color/symbol scheme based on year-classification (El Nino, La Nina
#           warm, cool, neutral)
#            1996,  1997,   1998,   1999,   2000,   2001,   2002,  
symsall <- c(symc[1],symc[1],symc[1],symc[1],symc[2],symc[1],symc[2],
            # 2003, 2004,   2005,   2006,   2007,   2008,   2009,
            symc[2],symc[2],symc[3],symc[3],symc[3],symc[3],symc[4],
            # 2010, 2011,   2012,   2013,   2014,   2015,   2016, 
            symc[2],symc[4],symc[5],symc[4],symc[5],symc[4],symc[3],
            # 2017, 2018,   2019,   2020
            symc[6],symc[6],symc[7],symc[7])
colsall <- c(colwm,colnu,colel,colla,colla,colcd,colcd,
            colwm,colnu,colwm,colcd,colnu,colla,colcd,
            colel,colla,colcd,colnu,colnu,colwm,colel,
            colnu,colcd,colnu,colcd)


# PLOT 1 - NMDS1
plt02_1 <- ggplot(data = model_gm1, aes(x = DOY, y = fit)) +  # For some reason, have to plot the GAM first...
  
  geom_segment(aes(x = 1, xend = 365, y = 0, yend = 0), color = "grey30", lwd = 0.3) + 
  geom_smooth_ci(linetype = 'dashed', lwd = 1) +
  geom_point(data = mstrdts_tib, aes(x = DOY, y = NMDS1, color = factor(pltyr), shape = factor(pltyr))) + 

  ylim(c(-1.1,1.1)) +
  scale_color_manual(name = 'Year',
                     labels = yrsplt,
                     values = colsall) + 
  scale_shape_manual(name = 'Year',
                     labels = yrsplt,
                     values = symsall) +
  scale_x_continuous(breaks = mstrdts_tib |>
                      group_by(month(Samp_Date)) |>
                       summarize(ndays = min(DOY)) |>
                       slice(c(1:12)) |>
                       pull(ndays),
                     labels = month.abb[c(1:12)]) +
  
  theme(axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(size = 14, colour = "black"),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "none", 
        legend.key = element_blank(),
        panel.background = element_blank(), 
        # panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  
  guides(color = guide_legend(ncol = 7),
         shape = guide_legend(ncol = 7))
# # Plot w legend
# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/P2_nMDS1_YrlyCyc.png", plot = plt02_1, width = 2000, height = 1600, units = 'px')
# # Plot w/o legend
# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/P2_nMDS1_YrlyCyc_noLgd.png", plot = plt02_1, width = 2000, height = 1200, units = 'px')


# PLOT 2 - NMDS2
plt02_2 <- ggplot(data = model_gm2, aes(x = DOY, y = fit)) +  # For some reason, have to plot the GAM first...
  
  geom_segment(aes(x = 1, xend = 365, y = 0, yend = 0), color = "grey30", lwd = 0.3) + 
  geom_smooth_ci(linetype = 'dashed', lwd = 1) +
  geom_point(data = mstrdts_tib, aes(x = DOY, y = NMDS2, color = factor(pltyr), shape = factor(pltyr))) + 
  
  ylim(c(-1.1,1.1)) +
  scale_color_manual(name = 'Year',
                     labels = yrsplt,
                     values = colsall) + 
  scale_shape_manual(name = 'Year',
                     labels = yrsplt,
                     values = symsall) +
  scale_x_continuous(breaks = mstrdts_tib |>
                       group_by(month(Samp_Date)) |>
                       summarize(ndays = min(DOY)) |>
                       slice(c(1:12)) |>
                       pull(ndays),
                     labels = month.abb[c(1:12)]) +
  
  theme(axis.title = element_text(size = 14, colour = "black"),
        axis.text = element_text(colour = "black", size = 12),
        legend.title = element_text(size = 14, colour = "black"),
        legend.text = element_text(size = 12, colour ="black"),
        # legend.position = "none",
        legend.position = "bottom",
        legend.key = element_blank(),
        panel.background = element_blank(), 
        # panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  guides(color = guide_legend(ncol = 7),
         shape = guide_legend(ncol = 7))

# # Plot w/ legend
# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/P2_nMDS2_YrlyCyc.png", plot = plt02_2, width = 2000, height = 1600, units = 'px')
# # Plot w/o legend
# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/P2_nMDS2_YrlyCyc_noLgd.png", plot = plt02_2, width = 2000, height = 1200, units = 'px')



