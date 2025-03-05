###########    NOAA/OSU post-doc    ##########
###   Step 2b: Copepod spp. timeseries     ###
##############################################
# PURPOSE: 
#     - Plot long-term timeseries for copepod species and community
#     - Calculate long-term trends (or not)


library(tidyverse)

# ### Input copepod species file
# copefl = read.csv(paste0('../Datasets/Zoops_NHL_csv_v2_fromCAM/NH05_CopeDens_log_subSpp_1996_2020.csv'))
copefl = read.csv(paste0('NH05_CopeDens_log_subSpp_1996_2020.csv'))


# Select species
# PSEUDO,CALMAR,ACALON,CENABD,ACATON,CALPAC,CALSTY,CALTEN,CLASO,CORANG,CTNCAL,CALTENU,PARA,CLAARC
copespp = readline("Which species?   ")
sppcut = copefl[,copespp]

sppdts_df = data.frame(copefl$Mon,copefl$Day,copefl$Year)
sppdts = as.Date(with(sppdts_df,paste(copefl$Year,copefl$Mon,copefl$Day,sep="-")),format="%Y-%m-%d")

sppcut[sppcut == 0] = NA

sppdf_all <- data.frame(cbind(sppdts,sppcut)) |>
  mutate(Date = as.Date(sppdts)) |>
  select(Date, sppcut)


lmspp <- lm(sppcut ~ Date, data = sppdf_all)
summary(lmspp)



########## PLOT TIMESERIES ###########

# Cool spp. 
plts2_1 <- ggplot(data = sppdf_all, aes(x = Date, y = sppcut)) + 
  geom_line(color = "blue3", lwd = 0.5) + 
  geom_point(color = "blue3", size = 0.8) + 
  geom_smooth(method = "lm", se = FALSE, color = "grey50") +
  annotate("text", label = "Adj. R-sq = 0.005", x = sppdf_all$Date[600], 
           y = 0.6, size = 4, colour = "black") +
  annotate("text", label = "p < 0.10", x = sppdf_all$Date[600], 
           y = 0.3, size = 4, colour = "black") +
  
  scale_x_continuous(breaks = seq.Date(as.Date(paste(min(year(sppdf_all$Date)),01,01,sep = "-")),
                                       as.Date(paste(max(year(sppdf_all$Date)),01,01,sep = "-")),365),
                     labels = c("1996","","","","2000","","","",
                                "2004","","","","2008","","","",
                                "2012","","","","2016","","","","2020"),
                     limits = c(as.Date(paste(min(year(sppdf_all$Date)),01,01,sep = "-")),
                                as.Date(paste(max(year(sppdf_all$Date)),12,31,sep = "-"))),
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

# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/PS2a_Along_Timeser.png", plot = plts2_1, width = 2000, height = 1000, units = 'px')




# Warm spp. - warmer colors
plts2_2 <- ggplot(data = sppdf_all, aes(x = Date, y = sppcut)) + 
  geom_line(color = "orangered3", lwd = 0.5) + 
  geom_point(color = "orangered3", size = 0.8) + 
  geom_smooth(method = "lm", se = FALSE, color = "grey50") +
  annotate("text", label = "Adj. R-sq = 0.16", x = sppdf_all$Date[600], 
           y = 0.6, size = 4, colour = "black") +
  annotate("text", label = "p < 0.001", x = sppdf_all$Date[600], 
           y = 0.3, size = 4, colour = "black") +
  
  scale_x_continuous(breaks = seq.Date(as.Date(paste(min(year(sppdf_all$Date)),01,01,sep = "-")),
                                       as.Date(paste(max(year(sppdf_all$Date)),01,01,sep = "-")),365),
                     labels = c("1996","","","","2000","","","",
                                "2004","","","","2008","","","",
                                "2012","","","","2016","","","","2020"),
                     limits = c(as.Date(paste(min(year(sppdf_all$Date)),01,01,sep = "-")),
                                as.Date(paste(max(year(sppdf_all$Date)),12,31,sep = "-"))),
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

# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/PS2d_Claso_Timeser.png", plot = plts2_2, width = 2000, height = 1000, units = 'px')

