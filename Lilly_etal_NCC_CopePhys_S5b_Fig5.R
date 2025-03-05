########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 5b - Plots for Fig. 5 - Winter PSI vs. BST date,
#                 Winter PSI vs. Summer PSI
##  Laura E. Lilly
##  Updated: 10 Mar 2024
########################################
# Compare:
#   - Winter PSI vs. BST date
#   - Winter PSI vs. Summer PSI (BST mag)
# ### OLD - NOTE: Must run 'Lilly_etal_NCC_CopePhys_S5a_PSISeason.R' *prior* to this script

library(tidyverse)
library(ggrepel)


###################### LOAD FILES ########################
scrfl <- read.csv('NH05_CopeDens_log10_nMDSscr_Samp_k2.csv')
colnames(scrfl) <- c("Samp_Date","NMDS1","NMDS2")


# ### Step 0: Convert file days -> R Dates
scrs_tib <- as_tibble(scrfl) |>
  mutate(Date = as.Date(Samp_Date),
         Year = year(Date),
         DOY2 = yday(Date)) |>
  select(Date, Year, DOY2, NMDS1, NMDS2)


# Then input BST dates as vector
bstdts <- as.Date(c("1996-07-01","1997-05-01","1998-07-16",
                    "1999-05-01","2000-04-01","2001-03-16",
                    "2002-04-16","2003-06-01","2004-05-16",
                    "2005-08-16","2006-05-16","2007-03-16",
                    "2008-03-01","2009-03-01","2010-06-16",
                    "2011-03-16","2012-05-01","2013-04-01",
                    "2014-04-01","2017-07-01","2018-05-16",
                    "2019-06-01","2020-04-01"),format="%Y-%m-%d")

# Also create vector of BST dates WITHOUT 1996, since we don't want to
# include that date in these analyses
bstdts_cut <- as.Date(c("1997-05-01","1998-07-16",
                    "1999-05-01","2000-04-01","2001-03-16",
                    "2002-04-16","2003-06-01","2004-05-16",
                    "2005-08-16","2006-05-16","2007-03-16",
                    "2008-03-01","2009-03-01","2010-06-16",
                    "2011-03-16","2012-05-01","2013-04-01",
                    "2014-04-01","2017-07-01","2018-05-16",
                    "2019-06-01","2020-04-01"),format="%Y-%m-%d")

# Second, calculate Winter and Summer date vectors by subtracting/adding 6 weeks
dtinv <- 6 # Set the 'number of weeks' interval to calculate Winter and Summer 
# dates

bst_df <- data.frame(bstdts_cut) |>
  mutate(Date = as.Date(bstdts_cut),
         BST_win = Date-(dtinv*7),
         BST_sum = Date+((dtinv+2)*7),
         DOY1 = yday(Date),
         Year = year(Date)) |>
  select(Date,Year,DOY1,BST_win,BST_sum)


# ### Step 1: Merge nMDS and BST DFs -> Separate for *winter* and *summer*
# WINTER
bst_win <- bst_df |>
  mutate(DOY_winBST = yday(BST_win), 
         BST_Date = Date,
         Date = BST_win) |>
  select(BST_Date, Year, Date, DOY_winBST)
  
win_merge <- left_join(bst_win, scrs_tib, 
                       join_by(Year, closest(Date <= Date))) |>
  mutate(BST_DOY = yday(BST_Date)) |>
  select(BST_Date, Year, BST_DOY, Date.y, NMDS1, NMDS2)
colnames(win_merge)[c(4)] <- c("Date_WinMDS")


# SUMMER
bst_sum <- bst_df |>
  mutate(DOY_sumBST = yday(BST_sum), 
         BST_Date = Date,
         Date = BST_sum) |>
  select(BST_Date, Year, Date, DOY_sumBST)

sum_merge <- left_join(bst_sum, scrs_tib, 
                       join_by(Year, closest(Date <= Date))) |>
  mutate(BST_DOY = yday(BST_Date)) |>
  select(BST_Date, Year, BST_DOY, Date.y, NMDS1, NMDS2)
colnames(sum_merge)[c(4)] <- c("Date_SumMDS")



# # ### Step 2: Calculate avg PSI for each year within a season
# # ALT #1: Don't need to add placeholder values for 2015, 2016, because
# #   they are already removed from PSIs DF
# bstyrdy_plt = bstyrdy
# # # ALT #2: Add 'placeholder' values for 2015, 2016 --> if I'm using a
# # #     BST dates vector *without* 2015, 2016
# # bstyrdy_plt = append(bstyrdy,c(50,50),after=19) # Add values of '50' so they fall on x-axis
# 
# # WINTER (each year compared to all other years -> avg.)
# win_avg_psi = colMeans(win_psi)
# 
# # SUMMER
# sum_avg_psi = colMeans(sum_psi)
# 
# 



#######################
# ninoyrs <- c("Warm","Neutral","Nino","Nina","Nina","Cool","Cool",
#              "Warm","Neutral","Warm","Cool","Neutral","Nina","Cool",
#              "Nino","Nina","Cool","Neutral","Neutral","Neutral","Cool","Neutral","Cool")
# ## Again, create a versoin without 1996 info
ninoyrs <- c("Neutral","Nino","Nina","Nina","Cool","Cool",
             "Warm","Neutral","Warm","Cool","Neutral","Nina","Cool",
             "Nino","Nina","Cool","Neutral","Neutral","Neutral","Cool","Neutral","Cool")

win_lbl <- cbind(win_merge, ninoyrs)
colnames(win_lbl)[7] <- "Nino_type"

sum_lbl <- cbind(sum_merge, ninoyrs)
colnames(sum_lbl)[7] <- "Nino_type"


# ### Step 3: PLOTS
colel <- "orangered"
colwm <- "orange2"
colcd <- "skyblue2"
colla <- "royalblue"
colnu <- "grey50"

pltcolsing <- c(colcd,colnu,colla,colel,colwm)


# FIG. 5a:  Winter nMDS score vs. BST yearday
lm_5a <- lm(BST_DOY ~ NMDS1, data = win_lbl)

# Have to keep "old" plot name (07a) for organization purposes among
#   plots generated for figures
plt07a <- ggplot(win_lbl, aes(x = NMDS1, y = BST_DOY)) + 
  geom_point(aes(x = NMDS1, y = BST_DOY, color = factor(Nino_type))) + 
  geom_smooth(method = 'lm', se = FALSE, color = "black", lwd = 0.3) + 
  
  geom_text_repel(aes(label = Year), # color = factor(psi_df$Nino_lbls)), 
                size = 3.5, nudge_y = -2) +
  labs(x = "nMDS score", y  = "BST date (DOY)") + 
  
  scale_color_manual(name = "Year-type",
                     values = pltcolsing) + 
  
  theme(aspect.ratio = 1,
        axis.text = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "right", 
        axis.title = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "grey50", fill = NA, lwd = 0.5),
        legend.key=element_blank()) 

# # Plot w/ legend
# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/P7a_2_WinMDS1scr_v_BST_v2no1996.png", plot = plt07a, width = 1600, height = 1600, units = 'px')



# FIG. 7b: BST yearday vs. Summer nMDS1 score
lm_5b <- lm(NMDS1 ~ BST_DOY, data = sum_lbl)

plt07b <- ggplot(sum_lbl, aes(x = BST_DOY, y = NMDS1)) + 
  geom_point(aes(x = BST_DOY, y = NMDS1, color = factor(Nino_type))) + 
  geom_smooth(method = 'lm', se = FALSE, color = "black", lwd = 0.3) + 
  
  geom_text_repel(aes(label = Year),
                  size = 3.5) +
  labs(x  = "BST date (DOY)", y = "nMDS score") + 
  
  
  scale_color_manual(name = "Year-type",
                     values = pltcolsing) + 
  
  theme(aspect.ratio = 1,
        axis.text = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "right", 
        axis.title = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "grey50", fill = NA, linewidth = 0.5),
        legend.key=element_blank()) 

# # Plot w/ legend
# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/P7b_2_SumMDS1scr_v_BST.png", plot = plt07b, width = 1600, height = 1600, units = 'px')

