###########      NOAA/OSU post-doc     ##########
###   Step 12: BST date interannual trend    ###
################################################
# PURPOSE: 
#     - Evaluate whether the date of BST shows a trending shift across years
# NOTE: BST dates were pre-determined by Cheryl A. Morgan and given to me as such


library(tidyverse)
library(lubridate)


# Open file of BST dates by year (from CAM)
dtsfl <- read.csv('NH05_Cope_TransDts_from_CAM.csv')

dtsall <- dtsfl |>
  mutate(BST_dt = make_date(Year,Mo_spr,Dy_spr),
         BFT_dt = make_date(Year,Mo_fal,Dy_fal),
         BST_doy = yday(BST_dt),
         BFT_doy = yday(BFT_dt))



# ### Plot & add trendlines
# BST
lm12_1 <- lm(BST_doy ~ Year, data = dtsall)
summary(lm12_1)

plt12_1 <- ggplot(data = dtsall, aes(x = Year, y = BST_doy)) + 
  geom_point(aes(x = Year, y = BST_doy)) + 
  geom_smooth(method = 'lm', se = FALSE) + 
  annotate("text", x = 2010, y = 220, label = "Adj. R-sq = -0.01", size = 5) + 
  annotate("text", x = 2010, y = 210, label = "p > 0.05", size = 5) + 
  
  xlab("Year") + 
  ylab("DOY") + 
  
  theme(axis.title.y = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 16, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black", 
                                   # angle = 90, vjust = 0.5, hjust = 1
        ),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
  )
# ggsave(paste0("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/PX4a_BST_DOYtrend.png"), plot = plt12_1, width = 2400, height = 1600, units = 'px')



# BFT
lm12_2 <- lm(BFT_doy ~ Year, data = dtsall)
summary(lm12_2)

plt12_2 <- ggplot(data = dtsall, aes(x = Year, y = BFT_doy)) + 
  geom_point(aes(x = Year, y = BFT_doy)) + 
  geom_smooth(method = 'lm', se = FALSE) + 
  
  annotate("text", x = 2010, y = 220, label = "Adj. R-sq = -0.04", size = 5) + 
  annotate("text", x = 2010, y = 210, label = "p > 0.05", size = 5) + 
  
  xlab("Year") + 
  ylab("DOY") + 
  
  theme(axis.title.y = element_text(size = 16, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 16, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black", 
                                   # angle = 90, vjust = 0.5, hjust = 1
        ),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
  )
# ggsave(paste0("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/PX4b_BFT_DOYtrend.png"), plot = plt12_2, width = 2400, height = 1600, units = 'px')

