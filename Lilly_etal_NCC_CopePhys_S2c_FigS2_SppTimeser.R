###########    NOAA/OSU post-doc    ##########
###   Step 2b: Copepod spp. timeseries     ###
##############################################
# PURPOSE: 
#     - Plot long-term timeseries for copepod species and community
#     - Calculate long-term trends (or not)


library(dplyr)
library(tidyverse)
library(lubridate)

# ### Input copepod species file
copefl = read.csv(paste0('Biol_files_2025/NH05_CopeDens_log_subSpp_1996_2020_from_CAM.csv'))


# Select species
# PSEUDO,CALMAR,ACALON,CENABD,ACATON,CALPAC,CALSTY,CALTEN,CLASO,CORANG,CTNCAL,CALTENU,PARA,CLAARC
# sum Clauso = c(CLASO,CLAARC,CLAPER,CLAPAR,CLAPAU,CLALIV)
copespp = readline("Which species?   ")

sppcut <- copefl |>
  select(Mon,Day,Year,copespp)
  # select(Mon,Day,Year,CLASO,CLAARC,CLAPER,CLAPAR,CLAPAU,CLALIV)


sppdf <- sppcut |>
  mutate(Date = as.Date(paste(Year,Mon,Day, sep = "-"),format = "%Y-%m-%d"),
         across(4:ncol(sppcut), ~ 10^(.x-1)-0.1)) |>
  mutate(spp_sum = rowSums(across(4:ncol(sppcut))),
         spp_sum_log = log10(spp_sum+0.1)+1)


# Calculate linear model for plotting
lmspp <- lm(spp_sum_log ~ Date, data = sppdf)
summary(lmspp)



########## PLOT TIMESERIES ###########

# Cool spp. 
plts2_1 <- ggplot(data = sppdf, aes(x = Date, y = spp_sum_log)) + 
  # geom_line(color = "blue3", lwd = 0.5) + 
  geom_point(color = "blue3", size = 0.8) + 
  geom_smooth(method = "lm", se = FALSE, color = "grey50") +
  annotate("text", label = "Adj. R-sq = 0.006", x = sppdf$Date[600], 
           y = 0.6, size = 4, colour = "black") +
  annotate("text", label = "p < 0.05", x = sppdf$Date[600], 
           y = 0.3, size = 4, colour = "black") +
  
  scale_x_continuous(breaks = seq.Date(as.Date(paste(min(year(sppdf$Date)),01,01,sep = "-")),
                                       as.Date(paste(max(year(sppdf$Date)),01,01,sep = "-")),365),
                     labels = c("1996","","","","2000","","","",
                                "2004","","","","2008","","","",
                                "2012","","","","2016","","","","2020"),
                     limits = c(as.Date(paste(min(year(sppdf$Date)),01,01,sep = "-")),
                                as.Date(paste(max(year(sppdf$Date)),12,31,sep = "-"))),
                     name = "") +
  scale_y_continuous(breaks = c(0,1,2,3,4),
                     labels = c("0","","2","","4"),
                     limits = c(0,4.5),
                     name = "Log10(density)") + 
  
  theme(axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black", 
                                   angle = 90, vjust = 0.5, hjust = 1),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
  )

# ggsave(paste0("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/PS2a_",copespp,"_Timeser_v2.png"), plot = plts2_1, width = 2000, height = 1000, units = 'px')




# Warm spp. - warmer colors
plts2_2 <- ggplot(data = sppdf, aes(x = Date, y = spp_sum_log)) + 
  # geom_line(color = "orangered3", lwd = 0.5) + 
  geom_point(color = "orangered3", size = 0.8) + 
  geom_smooth(method = "lm", se = FALSE, color = "grey50") +
  annotate("text", label = "Adj. R-sq = 0.02", x = sppdf$Date[600], 
           y = 0.6, size = 4, colour = "black") +
  annotate("text", label = "p < 0.01", x = sppdf$Date[600], 
           y = 0.3, size = 4, colour = "black") +
  
  scale_x_continuous(breaks = seq.Date(as.Date(paste(min(year(sppdf$Date)),01,01,sep = "-")),
                                       as.Date(paste(max(year(sppdf$Date)),01,01,sep = "-")),365),
                     labels = c("1996","","","","2000","","","",
                                "2004","","","","2008","","","",
                                "2012","","","","2016","","","","2020"),
                     limits = c(as.Date(paste(min(year(sppdf$Date)),01,01,sep = "-")),
                                as.Date(paste(max(year(sppdf$Date)),12,31,sep = "-"))),
                     name = "") +
  scale_y_continuous(breaks = c(0,1,2,3,4),
                     labels = c("0","","2","","4"),
                     limits = c(0,4.5),
                     name = "Log10(density)") + 
  
  theme(axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black", 
                                   angle = 90, vjust = 0.5, hjust = 1),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
  )

# ggsave(paste0("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/PS2d_",copespp,"_Timeser_v2.png"), plot = plts2_2, width = 2000, height = 1000, units = 'px')

