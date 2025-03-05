########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 5c - Plots for Supp. Fig. - Timeseries scatter of:
##      Year vs. BST DOY
##  Laura E. Lilly
##  Updated: 14 Dec 2024
########################################


library(tidyverse)


# # Input BST dates as vector
# bstdts <- as.Date(c(# "1996-07-01",
#                     "1997-05-01","1998-07-16",
#                     "1999-05-01","2000-04-01","2001-03-16",
#                     "2002-04-16","2003-06-01","2004-05-16",
#                     "2005-08-16","2006-05-16","2007-03-16",
#                     "2008-03-01","2009-03-01","2010-06-16",
#                     "2011-03-16","2012-05-01","2013-04-01",
#                     "2014-04-01","2017-07-01","2018-05-16",
#                     "2019-06-01","2020-04-01"),format="%Y-%m-%d")

bst_fl = read.csv("Biol_files_2025/NH05_Cope_BST_BFT_Dates_from_CAM_v2_Dec2024.csv")

bst_df <- bst_fl |>
  select(Year, Start.Date, End.Date) |>
  mutate(BST = as.Date(Start.Date, format = "%m/%d/%Y"),
         BFT = as.Date(End.Date, format = "%m/%d/%Y")) |>
  select(Year, BST, BFT) |>
  filter(Year %in% c(1997:2014,2017:2020))


### BST ###
bstdts <- bst_df$BST
bst_doy <- yday(bstdts)
bst_yr <- year(bstdts)
ninoyrs <- c("Neutral","Nino","Nina","Nina","Cool","Cool",
             "Warm","Neutral","Warm","Cool","Neutral","Nina","Cool",
             "Nino","Nina","Cool","Neutral","Neutral","Neutral","Cool","Neutral","Cool")

bst_df <- data.frame(cbind(bst_yr,bst_doy,ninoyrs)) |>
  mutate(bst_yr = as.numeric(bst_yr),
         bst_doy = as.numeric(bst_doy)) # |>
  # add_row(bst_yr = as.numeric(1996), bst_doy = NA, ninoyrs = "Warm", 
  #         .before = 1)


### BFT ###
bftdts <- bst_df$BFT
bft_doy <- yday(bftdts)
bft_yr <- year(bftdts)
ninoyrs <- c("Neutral","Nino","Nina","Nina","Cool","Cool",
             "Warm","Neutral","Warm","Cool","Neutral","Nina","Cool",
             "Nino","Nina","Cool","Neutral","Neutral","Neutral","Cool","Neutral","Cool")

bft_df <- data.frame(cbind(bft_yr,bft_doy,ninoyrs)) |>
  mutate(bft_yr = as.numeric(bft_yr),
         bft_doy = as.numeric(bft_doy))



# ### Plot: Year vs. BST DOY -> check for long-term trend
# ## Create a version without 1996 info

colel <- "orangered"
colwm <- "orange2"
colcd <- "skyblue2"
colla <- "royalblue"
colnu <- "grey50"
pltcolsing <- c(colcd,colnu,colla,colel,colwm)



#### Plot BST ####
# Calculate LM
lm_5c <- lm(bst_doy ~ bst_yr, data = bst_df)
summary(lm_5c)

# Plot data
plt01 <- ggplot(data = bst_df, aes(x = bst_yr, y = bst_doy)) + 
  
  geom_point(aes(x = bst_yr, y = bst_doy, color = factor(ninoyrs))) +
  geom_smooth(method = 'lm', se = FALSE, color = "black", lwd = 0.7) + 
  
  # labs(x = "Year", y  = "BST date (DOY)") + 
  
  scale_color_manual(name = "Year-type",
                     values = pltcolsing) + 
  scale_x_continuous(breaks = seq(1996,2020,1),
                     labels = c("1996","","","",
                                "2000","","","",
                                "2004","","","",
                                "2008","","","",
                                "2012","","","",
                                "2016","","","",
                                "2020"),
                     name = "") +
  scale_y_continuous(name = "") + 
  
  theme(axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black", 
                                   angle = 90, vjust = 0.5, hjust = 1),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        legend.position = 'none'
  )
            
# # Plot w/ legend
# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/PSX1_Year_v_BSTDOY_no1996_v2.png", plot = plt01, width = 2000, height = 1000, units = 'px')



#### Plot BFT ####
# Calculate LM
lm_5d <- lm(bft_doy ~ bft_yr, data = bft_df)
summary(lm_5d)

# Plot data
plt02 <- ggplot(data = bft_df, aes(x = bft_yr, y = bft_doy)) + 
  
  geom_point(aes(x = bft_yr, y = bft_doy, color = factor(ninoyrs))) +
  geom_smooth(method = 'lm', se = FALSE, color = "black", lwd = 0.7) + 
  
  scale_color_manual(name = "Year-type",
                     values = pltcolsing) + 
  scale_x_continuous(breaks = seq(1996,2020,1),
                     labels = c("1996","","","",
                                "2000","","","",
                                "2004","","","",
                                "2008","","","",
                                "2012","","","",
                                "2016","","","",
                                "2020"),
                     name = "") +
  scale_y_continuous(name = "") + 
  
  theme(axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black", 
                                   angle = 90, vjust = 0.5, hjust = 1),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        legend.position = 'none'
  )

# # Plot w/ legend
# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/PSX2_Year_v_BFTDOY_no1996_v2.png", plot = plt01, width = 2000, height = 1000, units = 'px')


