########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 2: Plot nMDS yearly scores
##  Laura E. Lilly
##  Updated: 18 May 2023
########################################
# From saved file ('Lilly_etal_NCC_CopePhys_S1_nMDS):
#   -> plot yearly timeseries of nMDS scores


library(tidyverse)
library(lubridate)
library(mgcv)
library(tidygam)
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Open datafile
# scrfl <- read.csv('NH05_Cope_biom_MDSscore_v4_CAM_RawDts.csv')
scrfl <- read.csv('NH05_CopeDens_log10_nMDSscr_Samp_k2.csv')
colnames(scrfl) <- c("Samp_Date","NMDS1","NMDS2")


# ### Step 0: Convert file days -> R Dates
scrs_tib <- as_tibble(scrfl) |>
  mutate(Samp_Date = as.Date(Samp_Date))

# # OLD code -> to convert dates to "Date" format
# scrdtfrm <- data.frame(scrfl[,3],scrfl[,1],scrfl[,2])
# scrdts <- as.Date(with(scrdtfrm,paste(scrfl[,3],scrfl[,1],scrfl[,2],sep="-")),format="%Y-%m-%d")
# scrtbl <- cbind(scrdts,scrfl[,4:6])


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

# # Delineate for NH05 vs. NH25
# if(yrsunq[1] == 1996){
#   symspal = symsall
#   colspal = colsall
# } else if(yrsunq[1] == 1998){
#   symspal = symsall[3:length(symsall)]
#   colspal = colsall[3:length(colsall)]
# }


# PLOT 1 - NMDS1
plt02_1 <- ggplot(data = model_gm1, aes(x = DOY, y = fit)) +  # For some reason, have to plot the GAM first...
  
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
  
  theme(axis.text = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "bottom", 
        axis.title = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.key=element_blank()) +
  guides(color = guide_legend(ncol = 7),
         shape = guide_legend(ncol = 7))

# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Plots_v4/P2_nMDS1_YrlyCyc.png", plot = plt02_1, width = 2000, height = 1600, units = 'px')



# PLOT 2 - NMDS2
plt02_2 <- ggplot(data = model_gm2, aes(x = DOY, y = fit)) +  # For some reason, have to plot the GAM first...
  
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
  
  theme(axis.text = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "bottom", 
        # legend.position = "none",
        axis.title = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.key=element_blank()) +
  guides(color = guide_legend(ncol = 7),
         shape = guide_legend(ncol = 7))

# # Plot w/ legend
# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Plots_v4/P2_nMDS2_YrlyCyc.png", plot = plt02_2, width = 2000, height = 1600, units = 'px')
# # Plot w/o legend
# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Plots_v4/P2_nMDS2_YrlyCyc_noLgd.png", plot = plt02_2, width = 2000, height = 1200, units = 'px')








  
  
################  OLD version -> Base R  ####################  
  
  
# # S3.1: Create Jan. 1 timeseries and find indices of 'mstrdayall' that match
# scrin <- readline(prompt = "Which nMDS dimension?")
# 
# jan1s <- seq(mstrdt[1],as.Date("2020-01-01"),by="year")
# idxj1 <- which(mstrdayall$Samp_Date %in% jan1s)
# idxend <- c(idxj1[2:length(idxj1)]-1,nrow(mstrdayall))
# 
# # S3.2: Get data.frame of year-long chunks for selected species
# yrsunq <- unique(year(mstrdayall$Samp_Date))
# dtswk <- format(mstrdayall$Samp_Date,"%m-%d")
# wksunq <- unique(dtswk)
# # Also get indices of Day 1 of each month -> for plotting purposes (below)
# mostunq <- unique(floor_date(mstrdayall$Samp_Date[1:365], unit = "month"))
# mostidx <- which(mstrdayall$Samp_Date %in% mostunq)
# mostdt <- format(as.Date(mostunq),"%m-%d")
# 
# 
# # For-loop to select each year-chunk
# scrchnks <- data.frame(matrix(nrow=length(yrsunq),ncol=length(wksunq)))
# for(jx in 1:length(idxj1)){
#   chnk = t(mstrdayall[idxj1[jx]:idxend[jx],as.numeric(scrin)+1])
#   scrchnks[jx,] = chnk
# }
# rownames(scrchnks) <- yrsunq
# colnames(scrchnks) <- wksunq




# # PLOT all yearly points
# dev.new()
# # plot(gam10,ylim=c(-1.2,1.2),xlab='',ylab='',xaxt='n',yaxt='n',bty='n',lwd=2,cex.axis = 1.7)
# plot(seq(1,ncol(scrchnks),1),scrchnks[1,],pch=symspal[1],col=colspal[1],ylim=c(-1.5,1.5),xlab='',ylab='',xaxt='n',yaxt='n',bty='n',cex.axis = 1.7)
# axis(1,labels=FALSE,tick=FALSE)
# for(pu in 2:nrow(scrchnks)){
#   points(seq(1,ncol(scrchnks),1),scrchnks[pu,],pch=symspal[pu],col=colspal[pu],cex=1)
# }
# lines(seq(1,ncol(scrchnks),1),rep(0,ncol(scrchnks)),col='grey50')
# axis(1,at=mostidx,labels=mostdt,cex.axis=1.7,las=2)
# axis(2,at=seq(-1,1,0.5),labels=c("-1","","0","","1"),cex.axis=1.7)
# legend("topright",legend=c(as.character(yrsunq)),col=colspal,pch=symspal,ncol=6)


