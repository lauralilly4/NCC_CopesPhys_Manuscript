###########      NOAA/OSU post-doc     ##########
###   Step 11a: Phenomix - species shifts     ###
##############################################
# PURPOSE: 
#     - Evaluate whether any species shows interannual shifts in peaks

# # Phenomix package downloaded via: 
# remotes::install_github("nwfsc-cb/phenomix",build_vignettes = TRUE)

# ### RUN AFTER - for *Species-level* analyses: 
# 'Lilly_etal_NCC_CopePhys_S2b_FigS2_SppTimeser.R'
# See below for 'nMDS scores' option


library(phenomix)
library(tidyverse)
library(ggplot2)


############ Data setup - for phenomix ############ 
### OPTION 1: Species density DF
sppdf_phen <- sppdf_all |>
 mutate(Year = year(Date),
        DOY = yday(Date)) |>
        # sppcut = 10^sppcut) |>
   na.omit()

# ### OPTION 2: nMDS score (population) - uncomment and run this section: 
# scrfl <- read.csv('NH05_CopeDens_log10_nMDSscr_Samp_k2.csv')
# colnames(scrfl) <- c("Date","NMDS1","NMDS2")
# 
# scrs_tib <- data.frame(scrfl) |>
#   mutate(Date = as.Date(Date)) |>
#   mutate(sppcut = NMDS2) |>  # SELECT: DIM 1 or 2
#   select(Date, sppcut)
# 
# sppdf_phen <- scrs_tib |>
#   mutate(Year = year(Date),
#          DOY = yday(Date),
#          sppcut = sppcut + 1) |>
#   na.omit()


####################################
# ### Set up DF for Fit Function ###
cov_dat <- data.frame(nyear = unique(sppdf_phen$Year))
cov_dat$nyear <- cov_dat$nyear - min(cov_dat$nyear)

spp_create <- create_data(data = sppdf_phen, 
                          variable = "sppcut", 
                          time = "Year", 
                          date = "DOY",
                          # asymmetric_model = TRUE, # For species-level
                          asymmetric_model = FALSE, # For nMDS
                          est_sigma_re = TRUE,
                          est_mu_re = TRUE,
                          mu = ~ nyear,
                          sigma = ~ nyear,
                          covar_data = cov_dat,
                          tail_model = "gaussian",
                          family = "gaussian"
                          )

# ### Run Fit Function
spp_fit <- fit(spp_create) # Create the fit

# Get names and values of fit predictions
fit_names <- names(spp_fit$sdreport$value)
pred_vals <- spp_fit$sdreport$value


####### Create Dataframe ######
# Create dataframe from 'create' object (which will be used for predictions and is for some reason a different length than original df)
create_df <- data.frame(spp_create$years,spp_create$x, spp_create$y)
colnames(create_df) <- c("Year","DOY","sppcut")

# Add '1995' to each year to get it to its proper value (right now: 1-25)
createdf_yrs <- create_df |>
  mutate(Year = Year + 1995)


# Combine original DF and 'Create' and see how they compare...
combo_df <- left_join(sppdf_phen,create_df)
# Get predicted values -> to combine with 'create' input values
createdf_yrs$pred <- as.numeric(spp_fit$sdreport$value[grep("pred", names(spp_fit$sdreport$value))])

# initvals <- spp_fit$init_vals


# ### Get Dates of 50% cumu. sum by Year ###
cumusumdf <- createdf_yrs |>
  group_by(Year) |>
  mutate(csum_spp = cumsum(sppcut),
         csum_pred = cumsum(pred))

# SPP DATA: Get 50% cumu. sum values
cumu50_spp <- cumusumdf |>
  select(Year,DOY,csum_spp) |>
  summarise(cumu50 = max(csum_spp)*0.5)

# SPP DATA: Then find DOY dates that match closest to 50% cumu
cumu50dts_spp <- data.frame(matrix(nrow = nrow(cumu50_spp), ncol = ncol(cumusumdf)))
for(s in 1:nrow(cumu50_spp)){
  mtchid = which.min(abs(cumusumdf$csum_spp[cumusumdf$Year == cumu50_spp$Year[s]]-cumu50_spp$cumu50[s]))
  stid = min(which(cumusumdf$Year == cumu50_spp$Year[s])) # Get first index of that year -> to add to subset index of matched cumu50
  idtot = mtchid+stid-1
  
  cumu50dts_spp[s,] <- cumusumdf[idtot,]
}
colnames(cumu50dts_spp) <- colnames(cumusumdf)



# PRED DATA: Get 50% cumu. sum values
cumu50_pred <- cumusumdf |>
  select(Year,DOY,csum_pred) |>
  summarise(cumu50 = max(csum_pred)*0.5)

# PRED DATA: Then find DOY dates that match closest to 50% cumu
cumu50dts_pred <- data.frame(matrix(nrow = nrow(cumu50_pred), ncol = ncol(cumusumdf)))
for(s in 1:nrow(cumu50_pred)){
  mtchid = which.min(abs(cumusumdf$csum_pred[cumusumdf$Year == cumu50_pred$Year[s]]-cumu50_pred$cumu50[s]))
  stid = min(which(cumusumdf$Year == cumu50_pred$Year[s])) # Get first index of that year -> to add to subset index of matched cumu50
  idtot = mtchid+stid-1
  
  cumu50dts_pred[s,] <- cumusumdf[idtot,]
}
colnames(cumu50dts_pred) <- colnames(cumusumdf)



# ######## PLOTS ########

# PLOT 0: Diagnostics
plt11_0 <- plot_diagnostics(spp_fit, type = "timing", logspace = FALSE)


# ### Plot 1: Species predictions
plt11_1 <- ggplot(data = createdf_yrs, aes(x = DOY, y = sppcut, group = factor(Year))) + 
  geom_point(aes(x = DOY, y = sppcut), color = "grey30", fill = "grey60", size = 0.7) + 
  geom_line(aes(x = DOY, y = pred), color = "magenta3") + 
  # geom_segment(data = cumu50dts_spp, aes(x = DOY, xend = DOY, y = 0, yend = 5), color = "black") + 
  # geom_segment(data = cumu50dts_pred, aes(x = DOY, xend = DOY, y = 0, yend = 5), color = "magenta3", linetype = "dashed") + 
  # scale_y_continuous(breaks = c(0,1,2,3,4,5),
  #                    labels = c("0","","2","","4",""),
  #                    limits = c(0,5), 
  #                    name = "Log(spp. density [no. m-2])") +
  
  
  geom_segment(data = cumu50dts_spp, aes(x = DOY, xend = DOY, y = 0, yend = 2), color = "black") +
  geom_segment(data = cumu50dts_pred, aes(x = DOY, xend = DOY, y = 0, yend = 2), color = "magenta3", linetype = "dashed") +
  ### For NMDS scores
  scale_y_continuous(breaks = c(0,0.5,1,1.5,2),
                     labels = c("0","","1","","2"),
                     limits = c(0,2),
                     name = "nMDS score") +
  
  
  theme(axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 9, colour = "black", 
                                   # angle = 90, vjust = 0.5, hjust = 1
                                   ),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
  ) + 
  
  facet_wrap(~ Year)
# ggsave(paste0("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/P5a_Yearly_RealvPred_",copespp,".png"), plot = plt11_1, width = 1600, height = 1200, units = 'px')



# PLOT 2: Peak population time by year
lm_spp <- lm(DOY ~ Year, data = cumu50dts_spp)
summary(lm_spp)
lm_pred <- lm(DOY ~ Year, data = cumu50dts_pred)
summary(lm_pred)

plt11_2 <- ggplot(data = cumu50dts_spp, aes(x = Year, y = DOY)) +
  
  geom_point(aes(x = Year, y = DOY), color = "grey40", fill = "grey70", size = 2) + 
  geom_smooth(data = cumu50dts_spp, method = "lm", se = TRUE, color = "black") +
  geom_point(data = cumu50dts_pred, aes(x = Year, y = DOY), color = "magenta4", fill = "magenta3", size = 2) + 
  geom_smooth(data = cumu50dts_pred, method = "lm", se = TRUE, color = "magenta4", linetype = "dashed") +
  annotate("text", x = 2016, y = 255, label = "Spp: Adj. R-sq. = 0.17") + 
  annotate("text", x = 2016, y = 250, label = "p-val < 0.05") + 
  annotate("text", x = 2016, y = 240, label = "Pred: Adj. R-sq. = 0.21") + 
  annotate("text", x = 2016, y = 235, label = "p-val < 0.05") + 

  ylim(c(140,260)) + 
  
  theme(axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black", 
                                   # angle = 90, vjust = 0.5, hjust = 1
        ),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
  )
# ggsave(paste0("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/P5b_DOYTrends_",copespp,".png"), plot = plt11_2, width = 1600, height = 1600, units = 'px')


# PLOT 3: Real vs. Modeled 50% cumu. 
# First, combine two dates DFs
peakdoy_cmb <- data.frame(cumu50dts_spp$Year, cumu50dts_spp$DOY, cumu50dts_pred$Year, cumu50dts_pred$DOY)
colnames(peakdoy_cmb) <- c("Year_spp", "DOY_spp", "Year_pred", "DOY_pred")
lm_peakdoy <- lm(DOY_pred ~ DOY_spp, peakdoy_cmb)


plt11_3 <- ggplot(data = peakdoy_cmb, aes(x = DOY_spp, y = DOY_pred)) + geom_point(aes(x = DOY_spp, y = DOY_pred)) + 
  geom_smooth(method = 'lm', se = FALSE) + 
  
  xlab("DOY of 50% cumu. population") + 
  ylab("DOY of 50% cumu. prediction") + 
  
  theme(axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black", 
                                   # angle = 90, vjust = 0.5, hjust = 1
        ),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
  )
# ggsave(paste0("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/P5c_DOYcomp_",copespp,".png"), plot = plt11_3, width = 1600, height = 1600, units = 'px')


