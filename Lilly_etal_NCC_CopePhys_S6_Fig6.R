########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 6b - Plots for Fig. 6 (ALF (raw, cumu) vs. BST date & mag)
##  Laura E. Lilly
##  Updated: 1 Aug 2023
########################################
# Compare:
#   - BST date vs. alongshore flow (instantaneous, cumulative)
#   - Summer PSI (BST mag) vs. alongshore flow (instantaneous, cumulative)

# Seasonal dates were same as for Cluster Analysis:
#   - Biological Spring/Fall Transitions (BST/BFT) - nMDS crossover, visual examination
#       by C.A. Morgan
#   - Winter - 6 weeks prior to BST *for each year* (varies by year)
#   - Summer - 8 weeks after BST

# ## NOTE: 
# Must run 'Lilly_etal_NCC_CopePhys_S5a_PSISeason.R' *prior* to this script



library(tidyverse)

######## Input copepod species file
### File 1: 'Raw' ALF transition dates - eyeballed by LEL
alffl_raw = read.csv(paste0('Phys_Inds/NH10_ALF_TrnsDts_Raw.csv'))

# ### File 2: 'Cumu' ALF transition dates - calced by BTC
alffl_cumu = read.csv(paste0('Phys_Inds/NH10_ALF_TrnsDts_Cumu.csv'))



# # ### Step 1: Convert Dates -> Yearday
# modys = c(31,28,31,30,31,30,31,31,30,31,30,31) # Number of days in each month -> to multiply by
# 
# # BST dates --> convert to yearday
# bstyrdy = vector()
# for(w in 1:length(sprdts)){
#   mno = month(sprdts[w])-1 # Subtract 1 because you only want number of *whole* months prior
#   if(mno > 0){
#     mdsum = sum(modys[1:mno])
#   } else if (mno == 0){
#     mdsum = 0
#   }
#   dsum = mdsum+day(sprdts[w])
#   bstyrdy = c(bstyrdy,dsum)
# }
# yrlbls = year(sprdts)


######## ALF dates
### Part 1: 'Raw' file - need to calculate 'yrdy'
alfyrdy_raw = vector()
for(a in 1:length(alffl_raw$Year)){
  
  mno = alffl_raw$Mo_spr[a]-1 # Subtract 1 because you only want number of *whole* months prior
  if(is.na(mno)==TRUE){
    mdsum = NA
  } else if(mno > 0){
    mdsum = sum(modys[1:mno])
  } else if(mno == 0){
    mdsum = 0
  }
  dsum = mdsum+alffl_raw$Dy_spr[a]
  alfyrdy_raw = c(alfyrdy_raw,dsum)
  
}
alfylbls_raw = alffl_raw$Year


### Part 2: 'Cumu' file - just grab 'yrdy' variable
# Just get rows for 1996-2020
alfylbls_cumu = alffl_cumu$Year[match(1996,alffl_cumu$Year):nrow(alffl_cumu)]
alfyrdy_cumu = alffl_cumu$Spr_YrDy[match(1996,alffl_cumu$Year):nrow(alffl_cumu)]


# Create table of BST & ALF yrdys 'raw' -> to match up by year
alfdf_raw = data.frame(alfylbls_raw,alfyrdy_raw)
colnames(alfdf_raw) = c("Year","ALF_yrdy_raw")
bstdf = data.frame(yrlbls,bstyrdy)
colnames(bstdf) = c("Year","BST_yrdy_raw")
# Merge all tables together by 'Year': ALF, BST_date, Summer_PSI
alf_bst_df_raw = merge(x=alfdf_raw,y=bstdf,by="Year",all.x=TRUE)

# Create table of BST & ALF yrdys 'cumu' -> to match up by year
alfdf_cumu = data.frame(alfylbls_cumu,alfyrdy_cumu)
colnames(alfdf_cumu) = c("Year","ALF_yrdy_cumu")
bstdf = data.frame(yrlbls,bstyrdy)
colnames(bstdf) = c("Year","BST_yrdy")
# Merge all tables together by 'Year': ALF, BST_date, Summer_PSI
alf_bst_df_cumu = merge(x=alfdf_cumu,y=bstdf,by="Year",all.x=TRUE)



# ### Step 2: Calculate avg PSI for each year within a season
# ALT #1 (*NOT* an 'OPTION): Just use 'normal' BST that already contains 2015, 2016 -> but
#     assign it the same vector name so I can switch out to Option 1 if needed
bstyrdy_plt = bstyrdy
# # ALT #2 (*NOT* an 'OPTION): Add 'placeholder' values for 2015, 2016 --> if I'm using a
# #     BST dates vector *without* 2015, 2016
# bstyrdy_plt = append(bstyrdy,c(50,50),after=19) # Add values of '50' so they fall on x-axis


# # WINTER (each year compared to all other years)
# win_col_psi = colMeans(win_psi)
# # win_col_noAnoms = win_col_psi[-(21:20)]

# SUMMER
sum_col_psi = colMeans(sum_psi)
# sum_col_noAnoms = sum_col_psi[-(21:20)]



# ### Step 3: CORRELATIONS

# FIG 6A (RAW): ALF_south_trans year-day vs. BST yearday
# NOTE: *CHANGE* text labels 22-25 -> 20-23 to remove spacing for 2015, 2016
dev.new()
par(pty='s')

fig6a <- plot(alf_bst_df_raw[,2],alf_bst_df_raw[,3],ann=FALSE,
     xlim=c(60,160),ylim=c(50,250),type='p',pch=16,cex.axis=1.7)
# Add year-labels - color-coded
text(alf_bst_df_raw[c(2,9,12,18,19,22,24),2],
     alf_bst_df_raw[c(2,9,12,18,19,22,24),3]+7,
     alf_bst_df_raw[c(2,9,12,18,19,22,24),1],col='grey50') # Neutral
text(alf_bst_df_raw[c(3,15),2],
     alf_bst_df_raw[c(3,15),3]+7,
     alf_bst_df_raw[c(3,15),1],col='orangered') #El Niño
text(alf_bst_df_raw[c(1,8,10),2],
     alf_bst_df_raw[c(1,8,10),3]+7,
     alf_bst_df_raw[c(1,8,10),1],col='orange') # Warm years
text(alf_bst_df_raw[c(6,7,11,14,17,23,25),2],
     alf_bst_df_raw[c(6,7,11,14,17,23,25),3]+7,
     alf_bst_df_raw[c(6,7,11,14,17,23,25),1],col='skyblue2') # cool years
text(alf_bst_df_raw[c(4,5,13,16),2],
     alf_bst_df_raw[c(4,5,13,16),3]+7,
     alf_bst_df_raw[c(4,5,13,16),1],col='royalblue3') # La Niña

# Test for a linear model
scttr6a = lm(alf_bst_df_raw[,3] ~ alf_bst_df_raw[,2])
abline(scttr6a)
summary(scttr6a)




# FIG 6B (CUMU): ALF_south_trans year-day vs. BST yearday
# NOTE: *CHANGE* text labels 22-25 -> 20-23 to remove spacing for 2015, 2016
dev.new()
par(pty='s')

fig6b <- plot(alf_bst_df_cumu[,2],alf_bst_df_cumu[,3],ann=FALSE,
     xlim=c(60,160),ylim=c(50,250),type='p',pch=16,cex.axis=1.7)
# Add year-labels - color-coded
text(alf_bst_df_cumu[c(2,9,12,18,19,22,24),2],
     alf_bst_df_cumu[c(2,9,12,18,19,22,24),3]+7,
     alf_bst_df_cumu[c(2,9,12,18,19,22,24),1],col='grey50') # Neutral
text(alf_bst_df_cumu[c(3,15),2],
     alf_bst_df_cumu[c(3,15),3]+7,
     alf_bst_df_cumu[c(3,15),1],col='orangered') #El Niño
text(alf_bst_df_cumu[c(1,8,10),2],
     alf_bst_df_cumu[c(1,8,10),3]+7,
     alf_bst_df_cumu[c(1,8,10),1],col='orange') # Warm years
text(alf_bst_df_cumu[c(6,7,11,14,17,23,25),2],
     alf_bst_df_cumu[c(6,7,11,14,17,23,25),3]+7,
     alf_bst_df_cumu[c(6,7,11,14,17,23,25),1],col='skyblue2') # cool years
text(alf_bst_df_cumu[c(4,5,13,16),2],
     alf_bst_df_cumu[c(4,5,13,16),3]+7,
     alf_bst_df_cumu[c(4,5,13,16),1],col='royalblue3') # La Nina

# Test for a linear model
scttr6b = lm(alf_bst_df_cumu[,3] ~ alf_bst_df_cumu[,2])
abline(scttr6b)
summary(scttr6b)



#########
# FIG. 6C (RAW): ALF_yearday vs. Summer PSI
dev.new()
par(pty='s')

fig6c <- plot(alf_bst_df_raw[,2],sum_col_psi,ann=FALSE,
     xlim=c(60,160),ylim = c(50,75),type='p',pch=16,cex.axis=1.7)
text(alf_bst_df_raw[c(2,9,12,18,19,22,24),2],
     sum_col_psi[c(2,9,12,18,19,22,24)]-0.7,
     alf_bst_df_raw[c(2,9,12,18,19,22,24),1],col='grey50') # Neutral
text(alf_bst_df_raw[c(3,15,21),2],
     sum_col_psi[c(3,15,21)]-0.7,
     alf_bst_df_raw[c(3,15,21),1],col='orangered') #El Niño
text(alf_bst_df_raw[c(1,8,10,20),2],
     sum_col_psi[c(1,8,10,20)]-0.7,
     alf_bst_df_raw[c(1,8,10,20),1],col='orange') # Warm years
text(alf_bst_df_raw[c(6,7,14,17),2],
     sum_col_psi[c(6,7,14,17)]-0.7,
     alf_bst_df_raw[c(6,7,14,17),1],col='skyblue2') # cool years
text(alf_bst_df_raw[c(11,23,25),2]-5,
     sum_col_psi[c(11,23,25)],
     alf_bst_df_raw[c(11,23,25),1],col='skyblue2') # cool years
text(alf_bst_df_raw[c(4,13,16),2],
     sum_col_psi[c(4,13,16)]-0.7,
     alf_bst_df_raw[c(4,13,16),1],col='royalblue3') # La Niña
text(alf_bst_df_raw[c(5),2]+5,
     sum_col_psi[c(5)],
     alf_bst_df_raw[c(5),1],col='royalblue3') # La Niña

# Test for a linear model
scttr6c = lm(sum_col_psi ~ alf_bst_df_raw[,2])
abline(scttr6c)
summary(scttr6c)




# FIG. 6D (CUMU): ALF_yearday vs. Summer PSI
dev.new()
par(pty='s')

fig6d <- plot(alf_bst_df_cumu[,2],sum_col_psi,ann=FALSE,
     xlim=c(60,160),ylim = c(50,75),type='p',pch=16,cex.axis=1.7)
text(alf_bst_df_cumu[c(2,9,12,18,19,22,24),2],
     sum_col_psi[c(2,9,12,18,19,22,24)]-0.7,
     alf_bst_df_cumu[c(2,9,12,18,19,22,24),1],col='grey50') # Neutral
text(alf_bst_df_cumu[c(3,15,21),2],
     sum_col_psi[c(3,15,21)]-0.7,
     alf_bst_df_cumu[c(3,15,21),1],col='orangered') #El Niño
text(alf_bst_df_cumu[c(1,8,10,20),2],
     sum_col_psi[c(1,8,10,20)]-0.7,
     alf_bst_df_cumu[c(1,8,10,20),1],col='orange') # Warm years
text(alf_bst_df_cumu[c(6,7,14,17),2],
     sum_col_psi[c(6,7,14,17)]-0.7,
     alf_bst_df_cumu[c(6,7,14,17),1],col='skyblue2') # cool years
text(alf_bst_df_cumu[c(11,23,25),2]-5,
     sum_col_psi[c(11,23,25)],
     alf_bst_df_cumu[c(11,23,25),1],col='skyblue2') # cool years
text(alf_bst_df_cumu[c(4,13,16),2],
     sum_col_psi[c(4,13,16)]-0.7,
     alf_bst_df_cumu[c(4,13,16),1],col='royalblue3') # La Niña
text(alf_bst_df_cumu[c(5),2]+5,
     sum_col_psi[c(5)],
     alf_bst_df_cumu[c(5),1],col='royalblue3') # La Niña

# Test for a linear model
scttr6d = lm(sum_col_psi ~ alf_bst_df_cumu[,2])
abline(scttr6d)
summary(scttr6d)



# FIG. 6D, v2 (CUMU): ALF_yearday vs. Summer PSI - NO 1998, 2015
alf_bst_df_cumu2 <- alf_bst_df_cumu[-c(20,3),]
sum_col_psi2 <- sum_col_psi[-c(20,3)]

dev.new()
par(pty='s')
fig6d2 <- plot(alf_bst_df_cumu2,sum_col_psi2,ann=FALSE,
              xlim=c(60,160),ylim = c(50,75),type='p',pch=16,cex.axis=1.7)
text(alf_bst_df_cumu[c(2,9,12,18,19,22,24),2],
     sum_col_psi[c(2,9,12,18,19,22,24)]-0.7,
     alf_bst_df_cumu[c(2,9,12,18,19,22,24),1],col='grey50') # Neutral
text(alf_bst_df_cumu[c(3,15,21),2],
     sum_col_psi[c(3,15,21)]-0.7,
     alf_bst_df_cumu[c(3,15,21),1],col='orangered') #El Niño
text(alf_bst_df_cumu[c(1,8,10,20),2],
     sum_col_psi[c(1,8,10,20)]-0.7,
     alf_bst_df_cumu[c(1,8,10,20),1],col='orange') # Warm years
text(alf_bst_df_cumu[c(6,7,14,17),2],
     sum_col_psi[c(6,7,14,17)]-0.7,
     alf_bst_df_cumu[c(6,7,14,17),1],col='skyblue2') # cool years
text(alf_bst_df_cumu[c(11,23,25),2]-5,
     sum_col_psi[c(11,23,25)],
     alf_bst_df_cumu[c(11,23,25),1],col='skyblue2') # cool years
text(alf_bst_df_cumu[c(4,13,16),2],
     sum_col_psi[c(4,13,16)]-0.7,
     alf_bst_df_cumu[c(4,13,16),1],col='royalblue3') # La Niña
text(alf_bst_df_cumu[c(5),2]+5,
     sum_col_psi[c(5)],
     alf_bst_df_cumu[c(5),1],col='royalblue3') # La Niña

# Test for a linear model
scttr6d2 = lm(sum_col_psi2 ~ alf_bst_df_cumu2)
abline(scttr6d2)
summary(scttr6d2)
