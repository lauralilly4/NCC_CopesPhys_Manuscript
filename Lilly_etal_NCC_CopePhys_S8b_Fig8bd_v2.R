########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 8a (formerly 6b), v2 - Get dates of Inst and Cumu ALF transitions
##  Plots for Fig. 8b,d - CUMULATIVE ALF vs. BST date & mag
##  Laura E. Lilly
##  Updated: 20 Jul 2024
## *SAME as 'S8a_Fig8ac', but for *cumulative* flow comparisons
########################################
# RUN AFTER: 
#   - 'Lilly_etal_NCC_CopePhys_S3a_SppProps.R'
#   - 'Lilly_etal_NCC_CopePhys_S5a_PSISeason.R'

# NOTE: Flow data start in 1997, NOT 1996 (copepod data start year)

library(tidyverse)
library(lubridate)
library(reshape2)
library(ggrepel)
library(zoo)


########## Data In ########## 
flwsfl <- read.csv('NH10_Flows_Inst_Cumu.csv')

# Dataset cleanup and date conversion
flws_data <- flwsfl |>
  mutate_all(~ifelse(is.nan(.), NA, .)) |>
  mutate(Date = as.Date(Date, format = "%d-%b-%Y"),
         Year = year(Date),
         .after = Date)

bst_doy <- yday(bstdts) # it only takes one line now, from former for-loop...


#############################
# ### DF #2: Cumulative flow
cumu_flw <- flws_data |>
  mutate(dt_mody = format(Date, format = "%m-%d"),
         DOY = yday(Date)) |>
  select(DOY,Year,Cumu_flow) |>
  mutate(cumu_sub = c(0,Cumu_flow[1:nrow(flws_data)-1]),
         cumu_diff = Cumu_flow-cumu_sub)

# Pivot DF to group values by year -> to get first value from each year where 
# slope changes to zero (then negative)
cumu_flw_yr <- cumu_flw |>
  select(DOY,Year,Cumu_flow) |>
  pivot_wider(names_from = Year, values_from = Cumu_flow)

# # Do the same for 'diffs' -> to get smallest value for each year
# cumu_diffs_yr <- cumu_flw |>
#   select(DOY,Year,cumu_diff) |>
#   pivot_wider(names_from = Year, values_from = cumu_diff)



# From Stack Overflow: https://stackoverflow.com/questions/41061140/how-to-calculate-the-average-slope-within-a-moving-window-in-r
#   Get rolling average slope across specified no. days (e.g., 14)
slopedys <- 14
cumu_slopes_yr <- data.frame(matrix(nrow = nrow(cumu_flw_yr)-(slopedys-1),
                                    ncol = ncol(cumu_flw_yr)))


######## LEFT OFF HERE: Make the 'Coef' function work for all-NA windows...
 
# Version 1 - returns the first error about can't run lm with NA, blah blah
# Coef <- function(Z) if (all(is.na(df$y))) NA else coef(lm(x ~ y, as.data.frame(Z))) 
Coef <- function(Z) ifelse (!all(is.na(df$y)), coef(lm(x ~ y, as.data.frame(Z), na.action = na.exclude)), NA) 


# # Test code with just a subset of data
# coef_sub <- df[251:(251+14),]
# ifelse (!all(is.na(coef_sub$y)), coef(lm(x ~ y, as.data.frame(coef_sub), na.action = na.exclude)), NA) 
# # Version 2 - returns a different error when I set the " == TRUE" part
# Coef <- function(Z) if ((all(is.na(df$y))) == TRUE) NA else coef(lm(x ~ y, as.data.frame(Z))) 


for(o in 2:ncol(cumu_flw_yr)){
  df = cumu_flw_yr[,c(1,o)]
  colnames(df) = c("x","y")
  if(sum(is.na(df)) > slopedys*3){
    slopesyr <- rep(NA,nrow(cumu_flw_yr)-(slopedys-1))
  } else {
    # slopesyr <- rollapplyr(zoo(df), slopedys, Coef, by.column = FALSE)
    slopesyr <- rollapplyr(df, slopedys, Coef, by.column = FALSE)
    
    # slopesyr <- rollapplyr(zoo(df), slopedys, Coef)
  }
  cumu_slopes_yr[,o] <- slopesyr
}
cumu_slopes_yr[,1] <- cumu_flw_yr[366-353+1:353,1]
colnames(cumu_slopes_yr) <- colnames(cumu_flw_yr)


# Now get minimum slope value for each year
cumu_slopes_long <- cumu_slopes_yr |>
  pivot_longer(!DOY, names_to = "Year", values_to = "Slope")

cumu_slopesmin <- cumu_slopes_long |>
  mutate(Year = as.numeric(Year)) |>
  group_by(Year) |>
  summarise(min = which.max(Slope < 0))


# Combine vectors into DF -> for plotting
# First, get avg Summer PSI
sum_avg_psi = colMeans(sum_psi)
# REMOVE from copepod-oriented datasets: 1996:1998,2004:2016 (because cumu_flow doesn't have
#     values for the first dates, and BST doesn't have dates for 2015, 2016)
sumpsi_sub <- sum_avg_psi[-c(1:3,9:19,23)]
bst_sub <- bst_doy[-c(1:3,9:19,23)]

# Remove 2015, 2016 from cumu_flow (because no BST dates)
cumu_sub <- cumu_slopesmin |>
  filter(!(Year %in% c(2015,2016)))

cumu_comps_df <- data.frame(cbind(cumu_sub$Year,cumu_sub$min,
                                  sumpsi_sub,bst_sub))
colnames(cumu_comps_df) <- c("Year","Cumu_DOY","PSI_DOY","BST_DOY")


########## Plots ##########
symc <- c(15,16,17,18,12,3,4,8)
colel <- "orangered"
colwm <- "orange2"
colcd <- "skyblue"
colla <- "royalblue"
colnu <- "grey50"

# Set up color/symbol scheme based on year-classification (El Nino, La Nina
#           warm, cool, neutral)
# #            1997,   1998,   1999,   2000,   2001,   2002,  
# symsall <- c(symc[1],symc[1],symc[1],symc[2],symc[1],symc[2],
#              # 2003, 2004,   2005,   2006,   2007,   2008,   2009,
#              symc[2],symc[2],symc[3],symc[3],symc[3],symc[3],symc[4],
#              # 2010, 2011,   2012,   2013,   2014,    
#              symc[2],symc[4],symc[5],symc[4],symc[5],
#              # 2017, 2018,   2019,   2020
#              symc[6],symc[6],symc[7],symc[7])
colsall <- c(colel,colla,colla,colcd,colcd,
             colwm,colnu,colwm,colcd,colnu,colla,colcd,
             colel,colla,colcd,colnu,colnu,
             colnu,colcd,colnu,colcd)


# ### Plot 20: Cumulative flow - year-chunks
plt20 <- ggplot(data = cumu_flw, aes(x = DOY, y = Cumu_flow)) + 
  
  geom_line(aes(color = factor(Year), group = Year)) + 
  geom_segment(aes(x = DOY[1], xend = DOY[length(DOY)],
                   y = 0, yend = 0), color = "black", lwd = 1) + 
  
  ylab("Flow (m/s)") + 
  xlab("Date") 
# scale_x_date(date_labels = "%b", date_breaks = "1 month")  


####
# ### Plot 21: moving *slopes* for all years
cumu_slopes_long <- cumu_slopes_yr |>
  pivot_longer(!DOY, names_to = "Year", values_to = "Slope")

plt21 <- ggplot(data = cumu_slopes_long, aes(x = DOY, y = Slope)) + 
  geom_line(aes(color = factor(Year), group = Year)) 



###
# ### Correlation: Cumu_Slope_Neg Yrdy vs. BST Yrdy
lm22 <- lm(BST_DOY ~ Cumu_DOY, data = cumu_comps_df)
summary(lm22)

cols2 <- c(colla,colla,colcd,colcd,colwm,
             colnu,colcd,colnu)

plt22 <- ggplot(data = cumu_comps_df,
                aes(x = Cumu_DOY, y = BST_DOY)) + 
  geom_point(aes(x = Cumu_DOY, y = BST_DOY, color = factor(Year)), size = 3) + 
  geom_smooth(method = 'lm', fill = NA) + 
  geom_text_repel(aes(label = Year)) + 
  annotate('text', x = 175, y = 210, label = "R adj. = -0.17") + 
  annotate('text', x = 175, y = 200, label = "p > 0.10") + 
  
  xlab("Yearday - negative Cumulative Slope") +
  ylab("BST Yearday") + 
  scale_color_manual(name = 'Year',
                     # labels = factor(Year),
                     values = cols2,
                     guide = 'none') +
  
  theme(axis.text = element_text(colour = "black", size = 12),
        # legend.text = element_text(size = 12, colour ="black"),
        # legend.position = "bottom", 
        axis.title = element_text(size = 14, colour = "black"),
        # legend.title = element_text(size = 14, colour = "black"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "grey50", fill = NA, size = 0.4),
        # legend.key=element_blank()
  ) 
# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/P8b_cumuFlw_v_BST.png", plot = plt22, width = 2000, height = 1600, units = 'px')



# ### Correlation: Cumu_flow_neg Yrdy vs. Summer PSI
lm23 <- lm(PSI_DOY ~ Cumu_DOY, data = cumu_comps_df)
summary(lm23)

plt23 <- ggplot(data = cumu_comps_df,
                aes(x = Cumu_DOY, y = PSI_DOY)) + 
  geom_point(aes(x = Cumu_DOY, y = PSI_DOY, color = factor(Year)), size = 3) + 
  geom_smooth(method = 'lm', fill = NA) + 
  geom_text_repel(aes(label = Year)) + 
  annotate('text', x = 174, y = 55, label = "R adj. = -0.07") + 
  annotate('text', x = 174, y = 52, label = "p > 0.10") + 
  
  xlab("Yearday - negative Cumulative Flow") +
  ylab("Summer avg. PSI") + 
  scale_color_manual(name = 'Year',
                     # labels = instyrs2,
                     values = cols2,
                     guide = 'none') +
  
  theme(axis.text = element_text(colour = "black", size = 12),
        # legend.text = element_text(size = 12, colour ="black"),
        # legend.position = "bottom", 
        axis.title = element_text(size = 14, colour = "black"),
        # legend.title = element_text(size = 14, colour = "black"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "grey50", fill = NA, size = 0.4),
        # legend.key=element_blank()
  ) 
# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/P8d_CumuFlw_v_SumPSI.png", plot = plt23, width = 2000, height = 1600, units = 'px')







