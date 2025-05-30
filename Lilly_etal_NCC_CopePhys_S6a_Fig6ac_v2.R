########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 8a (formerly 6b), v2 - Get dates of Inst and Cumu ALF transitions
##  Plots for Fig. 8a,c - INSTANTANEOUS ALF vs. BST date & mag
##  Laura E. Lilly
##  Updated: 20 Jul 2024
########################################
# RUN AFTER: 
#   - 'Lilly_etal_NCC_CopePhys_S3a_SppProps.R'
#   - 'Lilly_etal_NCC_CopePhys_S5a_PSISeason.R'

# NOTE: Flow data start in 1997, NOT 1996 (copepod data start year)

# ### RUN AFTER: 'Lilly_etal_NCC_CopePhys_S5b_Fig5.R' -> to get nMDS scores 
#     selected for summer dates

library(tidyverse)
library(lubridate)
library(reshape2)
library(ggrepel)



########## Data In ########## 
flwsfl <- read.csv('NH10_Flows_Inst_Cumu.csv')

# Flow data: Dataset cleanup and date conversion
flws_data <- flwsfl |>
    mutate_all(~ifelse(is.nan(.), NA, .)) |>
    mutate(Date = as.Date(Date, format = "%d-%b-%Y"),
      Year = year(Date),
         .after = Date)


########## Reshape DFs and get info #########
# ### DF #1: Instantaneous flow
# First, reconfigure dates
inst_flw <- flws_data |>
  mutate(dt_mody = format(Date, format = "%m-%d"),
         # Will have to remove Leap Days if I want to add a repeat Date seq for plotting
         # plt_dt = rep(seq.Date(as.Date("1999-01-01"), 
         #                  as.Date("1999-12-31"),
         #                  by = "day")
         ) |>
  select(dt_mody,Year,Inst_flow)


# Second, get date of first negative value for each year
# Pivot DF to group values by year -> to select one value from each year
inst_flw_yrly <- inst_flw |>
  pivot_wider(names_from = Year, values_from = Inst_flow)

inst_neg1 <- data.frame(matrix(nrow = 1, ncol = ncol(inst_flw_yrly)-1))

for(i in 2:ncol(inst_flw_yrly)){
  idx1 = which(inst_flw_yrly[,i] < 0)[1]
  dt1 = inst_flw_yrly[idx1,1]
  
  inst_neg1[1,i-1] <- dt1
}
colnames(inst_neg1) <- unique(inst_flw$Year)

bstyrdy <- yday(bstdts) # it only takes one line now, from former for-loop...

# # ### Convert Dates Yearday -> both BST and ALF_neg
# modys = c(31,28,31,30,31,30,31,31,30,31,30,31) # Number of days in each month -> to multiply by
# 
# # BST dates --> convert to yearday
# bstyrdy = vector()
# for(w in 1:length(bstdts)){
#   mno = month(bstdts[w])-1 # Subtract 1 because you only want number of *whole* months prior
#   if(mno > 0){
#     mdsum = sum(modys[1:mno])
#   } else if (mno == 0){
#     mdsum = 0
#   }
#   dsum = mdsum+day(bstdts[w])
#   bstyrdy = c(bstyrdy,dsum)
# }
# yrlbls = year(bstdts)


# Inst ALF dates --> convert to yearday
inst_matdts <- as.Date(as.matrix(inst_neg1), format = "%m-%d")
inst_negdts <- as.Date(paste(as.numeric(colnames(inst_neg1)),
                             month(inst_matdts),
                             day(inst_matdts), sep = "-"), format = "%Y-%m-%d")

instyrdy = vector()
for(w in 1:length(inst_negdts)){
  mno = month(inst_negdts[w])-1 # Subtract 1 because you only want number of *whole* months prior
  if(is.na(mno)){
    mdsum = NA
  } else if (mno > 0){
    mdsum = sum(modys[1:mno])
  } else if (mno == 0){
    mdsum = 0
  }
  dsum = mdsum+day(inst_negdts[w])
  instyrdy = c(instyrdy,dsum)
}
yrlbls <- year(bstdts)


# # Combine vectors into DF -> for plotting
# # First, get avg Summer PSI
# sum_avg_psi = colMeans(sum_psi)
# # REMOVE from *flow* datasets: 1997, 2015, 2016, 2021-24
# instyrdy2 <- instyrdy[-c(1,19,20,25:27)]
# instyrs2 <- as.numeric(colnames(inst_neg1))[-c(1,19,20,25:27)]
# REMOVE from *copepod* datasets: 1996,1997  - because flow doesn't have that year
# bstyrdy2 <- bstyrdy[-c(1:2)]
# sum_avg_psi2 <- sum_avg_psi[-c(1:2)]
# inst_comps_df <- data.frame(cbind(instyrs2,instyrdy2,bstyrdy2,sum_avg_psi2))


# Combine Instantaneous flow, BST date, and nMDS1 score
# REMOVE from Flow datasets: 1997, 2015, 2016, 2021-24
instyrdy2 <- instyrdy[-c(1,19,20,25:27)]
instyrs2 <- as.numeric(colnames(inst_neg1))[-c(1,19,20,25:27)]

# REMOVE from BST and nMDS datasets: 1996,1997  - because flow doesn't have those years
bstyrdy2 <- bstyrdy[-c(1:2)]
sum_nmds2 <- sum_merge[-c(1:2),]

# Now combined all of the above vectors
inst_comps_df <- data.frame(cbind(instyrs2,instyrdy2,bstyrdy2,sum_nmds2))


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


# ### Plot 10: Instantaneous flow - year-chunks
plt10 <- ggplot(data = inst_flw, aes(x = dt_mody, y = Inst_flow)) + 
  
  geom_line(aes(color = factor(Year), group = Year)) + 
  geom_segment(aes(x = dt_mody[1], xend = dt_mody[length(dt_mody)],
                   y = 0, yend = 0), color = "black", lwd = 1) + 
  
  ylab("Flow (m/s)") + 
  xlab("Date") 
  # scale_x_date(date_labels = "%b", date_breaks = "1 month")  


# ### Correlation: Inst_flow_neg Yrdy vs. BST Yrdy
lm11 <- lm(bstyrdy2 ~ instyrdy2, data = inst_comps_df)
summary(lm11)

plt11 <- ggplot(data = inst_comps_df,
                aes(x = instyrdy2, y = bstyrdy2)) + 
  geom_point(aes(x = instyrdy2, y = bstyrdy2, color = factor(instyrs2)), size = 3) + 
  geom_smooth(method = 'lm', fill = NA) + 
  geom_text_repel(aes(label = instyrs2)) + 
  annotate('text', x = 175, y = 210, label = "R adj. = -0.06") + 
  annotate('text', x = 175, y = 200, label = "p > 0.10") + 
  
  xlab("Yearday - Instantaneous Flow Transition") +
  ylab("BST Yearday") + 
  scale_color_manual(name = 'Year',
                     labels = instyrs2,
                     values = colsall,
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

# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/P8a_InstFlw_v_BST_v2.png", plot = plt11, width = 2000, height = 1600, units = 'px')



# ### Correlation: Inst_flow_neg Yrdy vs. Summer PSI
lm12 <- lm(NMDS1 ~ instyrdy2, data = inst_comps_df)
summary(lm12)

plt12 <- ggplot(data = inst_comps_df,
                aes(x = instyrdy2, y = NMDS1)) + 
  geom_point(aes(x = instyrdy2, y = NMDS1, color = factor(instyrs2)), size = 3) + 
  geom_smooth(method = 'lm', fill = NA) + 
  geom_text_repel(aes(label = instyrs2)) + 
  annotate('text', x = 174, y = 1, label = "R adj. = 0.07") + 
  annotate('text', x = 174, y = 0.9, label = "p > 0.05") + 
  
  xlab("Yearday - Instantaneous Flow Transition") +
  ylab("Summer nMDS1") + 
  scale_color_manual(name = 'Year',
                     labels = instyrs2,
                     values = colsall,
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

# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/P8c_2_InstFlw_v_SumNMDS.png", plot = plt12, width = 2000, height = 1600, units = 'px')











# ### Plot 2: Instantaneous flow - year-chunks
plt20 <- ggplot(data = cumu_flw, aes(x = dt_mody, y = Cumu_flow)) + 
  
  geom_line(aes(color = factor(Year), group = Year)) + 
  geom_segment(aes(x = dt_mody[1], xend = dt_mody[length(dt_mody)],
                   y = 0, yend = 0), color = "black", lwd = 1) + 
  
  ylab("Flow (m/s)") + 
  xlab("Date")


