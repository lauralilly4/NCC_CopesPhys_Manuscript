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
# NOTE: Must run 'Lilly_etal_NCC_CopePhys_S5a_PSISeason.R' *prior* to this script


library(tidyverse)


# ### Step 1: Convert Dates -> Yearday
modys = c(31,28,31,30,31,30,31,31,30,31,30,31) # Number of days in each month -> to multiply by

# BST dates --> convert to yearday
bstyrdy = vector()
for(w in 1:length(bstdts)){
  mno = month(bstdts[w])-1 # Subtract 1 because you only want number of *whole* months prior
  if(mno > 0){
    mdsum = sum(modys[1:mno])
  } else if (mno == 0){
    mdsum = 0
  }
  dsum = mdsum+day(bstdts[w])
  bstyrdy = c(bstyrdy,dsum)
}
yrlbls = year(bstdts)


# ### Step 2: Calculate avg PSI for each year within a season
# ALT #1: Don't need to add placeholder values for 2015, 2016, because
#   they are already removed from PSIs DF
bstyrdy_plt = bstyrdy
# # ALT #2: Add 'placeholder' values for 2015, 2016 --> if I'm using a
# #     BST dates vector *without* 2015, 2016
# bstyrdy_plt = append(bstyrdy,c(50,50),after=19) # Add values of '50' so they fall on x-axis


# WINTER (each year compared to all other years -> avg.)
win_avg_psi = colMeans(win_psi)

# SUMMER
sum_avg_psi = colMeans(sum_psi)


ninoyrs <- c("Warm","Neutral","Nino","Nina","Nina","Cool","Cool",
             "Warm","Neutral","Warm","Cool","Neutral","Nina","Cool",
             "Nino","Nina","Cool","Neutral","Neutral","Neutral","Cool","Neutral","Cool")

# Combine PSI values and year-labels into DF
psi_df <- as_tibble(data.frame(yrlbls,as.numeric(bstyrdy_plt),win_avg_psi,sum_avg_psi,ninoyrs))
colnames(psi_df) <- c("Year","BST_yrdy","Win_PSI","Sum_PSI","Nino_lbl")



#######################
# ### Step 3: PLOTS
colel <- "orangered"
colwm <- "orange2"
colcd <- "skyblue2"
colla <- "royalblue"
colnu <- "grey50"

pltcolsing <- c(colcd,colnu,colla,colel,colwm)


# FIG. 5b_1:  Winter PSI avg. vs. BST yearday
plt05b_1 <- ggplot(psi_df, aes(x = Win_PSI, y = BST_yrdy)) + 
  geom_point(aes(x = Win_PSI, y = BST_yrdy, color = factor(Nino_lbl))) + 
  geom_smooth(method = 'lm', se = FALSE, color = "black", lwd = 0.3) + 
  
  geom_text_repel(aes(label = yrlbls), # color = factor(psi_df$Nino_lbls)), 
                size = 3.5, nudge_y = -2) +
  labs(x = "PSI", y  = "BST Yearday") + 
  
  
  scale_color_manual(name = "Year-type",
                     values = pltcolsing) + 
  
  theme(aspect.ratio = 1,
        axis.text = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "right", 
        axis.title = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "grey50", fill = NA, size = 0.5),
        legend.key=element_blank()) 

# # Plot w/ legend
# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Plots_v4/P5b_1_WinPSI_v_BST.png", plot = plt05b_1, width = 1600, height = 1600, units = 'px')




# FIG. 5b_2:  Winter PSI avg. vs. Summer PSI avg
plt05b_1 <- ggplot(psi_df, aes(x = Win_PSI, y = BST_yrdy)) + 
  geom_point(aes(x = Win_PSI, y = BST_yrdy, color = factor(Nino_lbl))) + 
  geom_smooth(method = 'lm', se = FALSE, color = "black", lwd = 0.3) + 
  
  geom_text_repel(aes(label = yrlbls), # color = factor(psi_df$Nino_lbls)), 
                  size = 3.5, nudge_y = -2) +
  labs(x = "PSI", y  = "BST Yearday") + 
  
  
  scale_color_manual(name = "Year-type",
                     values = pltcolsing) + 
  
  theme(aspect.ratio = 1,
        axis.text = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "right", 
        axis.title = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "grey50", fill = NA, size = 0.5),
        legend.key=element_blank()) 


# # Plot w/ legend
# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Plots_v4/P5b_1_WinPSI_v_BST.png", plot = plt05b_1, width = 1600, height = 1600, units = 'px')


