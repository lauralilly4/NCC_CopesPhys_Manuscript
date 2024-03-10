########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 5a: Calculate Percent Similarity Index (PSI) between years - by season
##  Laura E. Lilly
##  Updated: 10 Mar 2024
########################################
# Percent Similarity Index (PSI) of interannual community variability
#   within a season
# 'Seasonal dates' are as follows:
#   - Biological Spring/Fall Transitions (BST/BFT) - nMDS crossover, visual examination
#       by C.A. Morgan
#   - Winter - 6 weeks prior to BST *for each year* (varies by year)
#   - Summer - 6 weeks post-BST

# Must run prior: 
#     - 'Lilly_etal_NCC_CopePhys_S3_SppProps.R'


library(tidyverse)
library(ggplot2)
library(reshape2)

dtinv <- 6 # Set the 'number of weeks' interval to calculate Winter and Summer 
            # dates


# ###  Step 1: INTERVAL DATES  ###
# First, input BST dates
bstdts <- as.Date(c("1996-07-01","1997-05-01","1998-07-16",
                    "1999-05-01","2000-04-01","2001-03-16",
                    "2002-04-16","2003-06-01","2004-05-16",
                    "2005-08-16","2006-05-16","2007-03-16",
                    "2008-03-01","2009-03-01","2010-06-16",
                    "2011-03-16","2012-05-01","2013-04-01",
                    "2014-04-01","2017-07-01","2018-05-16",
                    "2019-06-01","2020-04-01"),format="%Y-%m-%d")

# Second, calculate Winter and Summer date vectors by subtracting/adding 6 weeks
windts <- bstdts-(dtinv*7) # 6 weeks * 7 days/week 
sumdts <- bstdts+((dtinv+2)*7)



# ###  Step 2: INTERVAL COMMUNITIES  ###
# Get Year of each date
bstyrs <- format(bstdts,format="%Y")


# Get rows of each community
# WINTER
win_comm <- data.frame(matrix(nrow=length(windts),ncol=ncol(mstrdata)))
win_diffs <- data.frame(matrix(nrow=length(windts),ncol = 1))

for(p in 1:length(windts)){
  did = which.min(abs(dtmstr-windts[p]))
  dval = dtmstr[did]-windts[p]
  
  # If the closest matching date is empty (no samples), move to next closest
  #     date, but *farther* away from BST (i.e., more in winter)
  if(sum(is.na(mstrdata[did,]))==ncol(mstrdata)){
    did = did-1
    # If that one is also empty, try sample on other side of original winter
    #     date (closer to BST)
    if(sum(is.na(mstrdata[did,]))==ncol(mstrdata)){
      did = did+2
    }
  }
  
  win_comm[p,] = mstrdata[did,]
}

colnames(win_comm) <- colnames(mstrdata)
rownames(win_comm) <- bstyrs
# rownames(win_diffs) <- winyrs


# BST
bst_comm <- data.frame(matrix(nrow=length(bstdts),ncol=ncol(mstrdata)))
# bst_diffs <- data.frame(matrix(nrow=length(bstdts),ncol = 1))

for(p in 1:length(bstdts)){
  did = which.min(abs(dtmstr-bstdts[p]))
  
  # This one is pretty straightforward: If the BST sample is empty... well,
  #     so be it (but it shouldn't be... because the BST is sort of necessarily
  #     defined by a present sample)
  bst_comm[p,] = mstrdata[did,]
}
colnames(bst_comm) <- colnames(mstrdata)
rownames(bst_comm) <- bstyrs



# SUMMER
sum_comm <- data.frame(matrix(nrow=length(sumdts),ncol=ncol(mstrdata)))
# sum_diffs <- data.frame(matrix(nrow=length(sumdts),ncol=1))

for(p in 1:length(sumdts)){
  did = which.min(abs(dtmstr-sumdts[p]))
  dval = dtmstr[did]-sumdts[p]
  
  # If the closest matching date is empty (no samples), move to next closest
  #     date, but *farther* away from BST (i.e., farther ahead in *summer*)
  if(sum(is.na(mstrdata[did,]))==ncol(mstrdata)){
    did = did+1
    # If that one is also empty, try sample on other side of original summer
    #     date (closer to BST)
    if(sum(is.na(mstrdata[did,]))==ncol(mstrdata)){
      did = did-2
    }
  }
  
  sum_comm[p,] = mstrdata[did,]
}

colnames(sum_comm) <- colnames(mstrdata)
rownames(sum_comm) <- bstyrs



######################################
# Second, calculate PSIs between years

# First, get *subset* of years used (i.e., no 2015, 2016)
yrssub <- as.numeric(bstyrs)

# Winter
winpcts = (win_comm/rowSums(win_comm))*100
win_psi = data.frame(matrix(nrow=length(yrssub),ncol=length(yrssub)))

for(yr in 1:length(yrssub)){
  yidx = which(yrssub %in% yrssub[yr])
  yrpsis = vector() 
  for(yc in 1:length(yrssub)){
    yrdfs = winpcts[yidx,]-winpcts[yc,]
    ypsi = 100-0.5*sum(abs(yrdfs))
    yrpsis = c(yrpsis,ypsi)
  }
  win_psi[yr,] = yrpsis
}
rownames(win_psi) <- yrssub
colnames(win_psi) <- yrssub
# # Calculate 'avg PSI' for each year (aka the average of its PSIs with all other years)
# win_avg_psi <- colMeans(win_psi)
# winpsi_tib <- as_tibble(win_psi)


# BST
bstpcts = (bst_comm/rowSums(bst_comm))*100
bst_psi = data.frame(matrix(nrow=length(yrssub),ncol=length(yrssub)))

for(yr in 1:length(yrssub)){
  yidx = which(yrssub %in% yrssub[yr])
  yrpsis = vector() 
  for(yc in 1:length(yrssub)){
    yrdfs = bstpcts[yidx,]-bstpcts[yc,]
    ypsi = 100-0.5*sum(abs(yrdfs))
    yrpsis = c(yrpsis,ypsi)
  }
  bst_psi[yr,] = yrpsis
}
rownames(bst_psi) <- yrssub
colnames(bst_psi) <- yrssub



# Summer
sumpcts = (sum_comm/rowSums(sum_comm))*100
sum_psi = data.frame(matrix(nrow=length(yrssub),ncol=length(yrssub)))

for(yr in 1:length(yrssub)){
  yidx = which(yrssub %in% yrssub[yr])
  yrpsis = vector() 
  for(yc in 1:length(yrssub)){
    yrdfs = sumpcts[yidx,]-sumpcts[yc,]
    ypsi = 100-0.5*sum(abs(yrdfs))
    yrpsis = c(yrpsis,ypsi)
  }
  sum_psi[yr,] = yrpsis
}
rownames(sum_psi) <- yrssub
colnames(sum_psi) <- yrssub




###############################
# ### STEP 3: Plot PSIs by year
## Color Scheme -> four categories
symc = c(15,16,17,18,19,3,4,8)
colel = "orangered"
colwm = "orange2"
colcd = "skyblue"
colla = "royalblue"
colnu = "grey50"

# Lay out symbols and colors
psisyms = c(symc[1],symc[1],symc[1],symc[1],symc[2],symc[1],symc[2],
            symc[2],symc[2],symc[3],symc[3],symc[3],symc[3],symc[4],
            symc[2],symc[4],symc[5],symc[4],symc[5],symc[4],symc[3],
            symc[6],symc[6],symc[7],symc[7])
psicols = c(colwm,colnu,colel,colla,colla,colcd,colcd,
            colwm,colnu,colwm,colcd,colnu,colla,colcd,
            colel,colla,colcd,colnu,colnu,colwm,colel,
            colnu,colcd,colnu,colcd)

# Remove the '2015' and '2016' colors & symbols
psisyms2 = psisyms[-(21:20)]
psicols2 = psicols[-(21:20)]



# ### PLOT 1 - WINTER PSI
# First, melt DF down to single columns (from square matrix)
win_psi_mlt <- as_tibble(melt(as.matrix(win_psi), id = colnamesas.matrix(win_psi)))
colnames(win_psi_mlt) <- c("Year", "Comp_yr", "PSI")
win_psi_srt <- win_psi_mlt |>
  arrange(Year) |>
  mutate(PSI = replace(PSI, PSI == 100, 0))


plt05_1 <- ggplot(data = win_psi_srt, aes(x = Year, y = PSI)) +  # For some reason, have to plot the GAM first...
  geom_point(aes(x = Year, y = PSI, color = factor(Comp_yr), shape = factor(Comp_yr))) + 
  
  ylim(c(0,100)) +
  scale_color_manual(name = 'Year',
                     labels = yrssub,
                     values = psicols2) + 
  scale_shape_manual(name = 'Year',
                     labels = yrssub,
                     values = psisyms2) +
  scale_x_continuous(breaks = seq(1996,2020,1),
                     labels = c("1996","","","",
                                "2000","","","",
                                "2004","","","",
                                "2008","","","",
                                "2012","","","",
                                "2016","","","","2020")) +

  theme(axis.text = element_text(colour = "black", size = 12),
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "bottom", 
        axis.title = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "grey50", fill = NA, size = 0.4),
        legend.key=element_blank()) +
  guides(color = guide_legend(ncol = 7),
         shape = guide_legend(ncol = 7))

# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Plots_v4/P5_1_Win_PSI.png", plot = plt05_1, width = 2000, height = 1600, units = 'px')




# ### PLOT 2 - BST PSI
# First, melt DF down to single columns (from square matrix)
bst_psi_mlt <- as_tibble(melt(as.matrix(bst_psi), id = colnames(as.matrix(bst_psi))))
colnames(bst_psi_mlt) <- c("Year", "Comp_yr", "PSI")
bst_psi_srt <- bst_psi_mlt |>
  arrange(Year) |>
  mutate(PSI = replace(PSI, PSI == 100, 0))


plt05_2 <- ggplot(data = bst_psi_srt, aes(x = Year, y = PSI)) +  # For some reason, have to plot the GAM first...
  geom_point(aes(x = Year, y = PSI, color = factor(Comp_yr), shape = factor(Comp_yr))) + 
  
  ylim(c(0,100)) +
  scale_color_manual(name = 'Year',
                     labels = yrssub,
                     values = psicols2) + 
  scale_shape_manual(name = 'Year',
                     labels = yrssub,
                     values = psisyms2) +
  scale_x_continuous(breaks = seq(1996,2020,1),
                     labels = c("1996","","","",
                                "2000","","","",
                                "2004","","","",
                                "2008","","","",
                                "2012","","","",
                                "2016","","","","2020")) +
  
  theme(axis.text = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "bottom", 
        axis.title = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "grey50", fill = NA, size = 0.4),
        legend.key=element_blank()) +
  guides(color = guide_legend(ncol = 7),
         shape = guide_legend(ncol = 7))

# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Plots_v4/P5_2_BST_PSI.png", plot = plt05_2, width = 2000, height = 1600, units = 'px')




# ### PLOT 3 - SUMMER PSI
# First, melt DF down to single columns (from square matrix)
sum_psi_mlt <- as_tibble(melt(as.matrix(sum_psi), id = colnames(as.matrix(sum_psi))))
colnames(sum_psi_mlt) <- c("Year", "Comp_yr", "PSI")
sum_psi_srt <- sum_psi_mlt |>
  arrange(Year) |>
  mutate(PSI = replace(PSI, PSI == 100, 0))


plt05_3 <- ggplot(data = sum_psi_srt, aes(x = Year, y = PSI)) +  # For some reason, have to plot the GAM first...
  geom_point(aes(x = Year, y = PSI, color = factor(Comp_yr), shape = factor(Comp_yr))) + 
  
  ylim(c(0,100)) +
  scale_color_manual(name = 'Year',
                     labels = yrssub,
                     values = psicols2) + 
  scale_shape_manual(name = 'Year',
                     labels = yrssub,
                     values = psisyms2) +
  scale_x_continuous(breaks = seq(1996,2020,1),
                     labels = c("1996","","","",
                                "2000","","","",
                                "2004","","","",
                                "2008","","","",
                                "2012","","","",
                                "2016","","","","2020")) +
  
  theme(axis.text = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "bottom", 
        axis.title = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "grey50", fill = NA, size = 0.4),
        legend.key=element_blank()) +
  guides(color = guide_legend(ncol = 7),
         shape = guide_legend(ncol = 7))

# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Plots_v4/P5_3_Sum_8wks_PSI.png", plot = plt05_3, width = 2000, height = 1600, units = 'px')


