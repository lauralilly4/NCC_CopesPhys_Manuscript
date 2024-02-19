########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 6a - Plots for Fig. 5 (Winter PSI vs. BST date & mag)
##  Laura E. Lilly
##  Updated: 1 Aug 2023
########################################
# Compare:
#   - BST date vs. Winter PSI
#   - Summer PSI (BST mag) vs. Winter PSI

# Seasonal dates were same as for Cluster Analysis:
#   - Biological Spring/Fall Transitions (BST/BFT) - nMDS crossover, visual examination
#       by C.A. Morgan
#   - Winter - 6 weeks prior to BST *for each year* (varies by year)
#   - Summer - 8 weeks after BST

# ## NOTE #1: 
# Must run 'Lilly_etal_NCC_CopePhys_S5a_PSISeason.R' *prior* to this script


library(tidyverse)


# ### Step 1: Convert Dates -> Yearday
modys = c(31,28,31,30,31,30,31,31,30,31,30,31) # Number of days in each month -> to multiply by

# BST dates --> convert to yearday
bstyrdy = vector()
for(w in 1:length(sprdts)){
  mno = month(sprdts[w])-1 # Subtract 1 because you only want number of *whole* months prior
  if(mno > 0){
    mdsum = sum(modys[1:mno])
  } else if (mno == 0){
    mdsum = 0
  }
  dsum = mdsum+day(sprdts[w])
  bstyrdy = c(bstyrdy,dsum)
}
yrlbls = year(sprdts)


# ### Step 2: Calculate avg PSI for each year within a season
# ALT #1: Just use 'normal' BST that already contains 2015, 2016 -> but
#     assign it the same vector name so I can switch out to Option 1 if needed
bstyrdy_plt = bstyrdy
# # ALT #2: Add 'placeholder' values for 2015, 2016 --> if I'm using a
# #     BST dates vector *without* 2015, 2016
# bstyrdy_plt = append(bstyrdy,c(50,50),after=19) # Add values of '50' so they fall on x-axis


# WINTER (each year compared to all other years)
win_col_psi = colMeans(win_psi)
# win_col_noAnoms = win_col_psi[-(21:20)]

# SUMMER
sum_col_psi = colMeans(sum_psi)
# sum_col_noAnoms = sum_col_psi[-(21:20)]



# ### Step 3: CORRELATIONS

# FIG. 5A:  Winter PSI avg. vs. BST yearday
# # ## OPTION 1: NO 2015, 2016 (b/c no BSTs)
# dev.new()
# par(pty='s')
# 
# plot(pre_avg_noAnoms,bstyrdy,xlim=rev(range(40,65)),ann=FALSE,
#      ylim=c(50,250),type='p',pch=16,cex.axis=1.7)
# # Add year-labels - color-coded
# text(pre_avg_noAnoms[c(2,9,12,18,19,20,22)],
#      bstyrdy[c(2,9,12,18,19,20,22)]+7,
#      yrlbls[c(2,9,12,18,19,20,22)],col='grey50') # Neutral
# text(pre_avg_noAnoms[c(3,15)],
#      bstyrdy[c(3,15)]+7,
#      yrlbls[c(3,15)],col='orangered') #El Nino
# text(pre_avg_noAnoms[c(1,8,10)],
#      bstyrdy[c(1,8,10)]+7,
#      yrlbls[c(1,8,10)],col='orange') # Warm years
# text(pre_avg_noAnoms[c(6,7,11,14,17,21,23)],
#      bstyrdy[c(6,7,11,14,17,21,23)]+7,
#      yrlbls[c(6,7,11,14,17,21,23)],col='skyblue2') # cool years
# text(pre_avg_noAnoms[c(4,5,13,16)],
#      bstyrdy[c(4,5,13,16)]+7,
#      yrlbls[c(4,5,13,16)],col='royalblue3') # La Niña
# 
# # Test for a linear model
# scttr = lm(bstyrdy ~ pre_avg_noAnoms)
# abline(scttr)
# text(43,215,'R-sq. = 0.29',cex=1.2)
# text(43,205,'p-val < 0.01',cex=1.2)


# ## OPTION 2: WITH 2015, 2016 (on x-axis)
dev.new()
par(pty='s')

plot(win_col_psi,bstyrdy_plt,xlim=c(42,68),ylim=c(50,250),type='p',pch=16)
# Add year-labels - color-coded
text(win_col_psi[c(2,9,12,18,19,22,24)],
     bstyrdy_plt[c(2,9,12,18,19,22,24)]+7,
     yrlbls[c(2,9,12,18,19,22,24)],col='grey50') # Neutral
text(win_col_psi[c(3,15)],
     bstyrdy_plt[c(3,15)]+7,
     yrlbls[c(3,15)],col='orangered') #El Nino
text(win_col_psi[c(21)]+1,
     bstyrdy_plt[c(21)]-7,
     yrlbls[c(21)],col='orangered') # El Nino, pt. 2 -> Label 2016 BELOW
text(win_col_psi[c(1,8,10)],
     bstyrdy_plt[c(1,8,10)]+7,
     yrlbls[c(1,8,10)],col='orange') # Warm years
text(win_col_psi[c(20)]-1,
     bstyrdy_plt[c(20)]-7,
     yrlbls[c(20)],col='orange') # Warm years, pt. 2 -> Label 2015 BELOW
text(win_col_psi[c(6,7,11,17,23,25)],
     bstyrdy_plt[c(6,7,11,17,23,25)]+7,
     yrlbls[c(6,7,11,17,23,25)],col='skyblue2') # cool years
text(win_col_psi[c(14)]-1.5,
     bstyrdy_plt[c(14)],
     yrlbls[c(14)],col='skyblue2') # cool years, pt. 2
text(win_col_psi[c(4,13,16)],
     bstyrdy_plt[c(4,13,16)]+7,
     yrlbls[c(4,13,16)],col='royalblue3') # La Niña
text(win_col_psi[c(5)]-1.5,
     bstyrdy_plt[c(5)],
     yrlbls[c(5)],col='royalblue3') # La Nina, pt. 2

# Test for a linear model
scttr5a = lm(bstyrdy_plt ~ win_col_psi)
abline(scttr5a)
# summary(scttr5a)





#########
# FIG. 5B:
dev.new()
par(pty='s')

plot(win_col_psi,sum_col_psi,xlim=c(42,68),ylim = c(35,75),type='p',pch=16)
# Add year-labels - color-coded
text(win_col_psi[c(2,9,12,18,19,22,24)]+1.5,
     sum_col_psi[c(2,9,12,18,19,22,24)],
     yrlbls[c(2,9,12,18,19,22,24)],col='grey50') # Neutral
text(win_col_psi[c(3,21)]+1.5,
     sum_col_psi[c(3,21)],
     yrlbls[c(3,21)],col='orangered') #El Nino
text(win_col_psi[c(15)]-1.5,
     sum_col_psi[c(15)],
     yrlbls[c(15)],col='orangered') # El Nino, pt. 2 -> Label 2016 BELOW
text(win_col_psi[c(1,8,10)]+1.5,
     sum_col_psi[c(1,8,10)],
     yrlbls[c(1,8,10)],col='orange') # Warm years
text(win_col_psi[c(20)]+1.5,
     sum_col_psi[c(20)],
     yrlbls[c(20)],col='orange') # Warm years, pt. 2 -> Label 2015 BELOW
text(win_col_psi[c(6,7,11,14,23)]+0.5,
     sum_col_psi[c(6,7,11,14,23)]-1,
     yrlbls[c(6,7,11,14,23)],col='skyblue2') # cool years
text(win_col_psi[c(17)]-1,
     sum_col_psi[c(17)]+1,
     yrlbls[c(17)],col='skyblue2') # cool years, pt. 2 -> 2012
text(win_col_psi[c(25)]+1,
     sum_col_psi[c(25)]+1,
     yrlbls[c(25)],col='skyblue2') # cool years, pt. 3 -> 2020
text(win_col_psi[c(5,13,16)],
     sum_col_psi[c(5,13,16)]-1,
     yrlbls[c(5,13,16)],col='royalblue3') # La Niña
text(win_col_psi[c(4)]-1.5,
     sum_col_psi[c(4)],
     yrlbls[c(4)],col='royalblue3') # La Niña

# Test for a linear model
scttr5b = lm(sum_col_psi ~ win_col_psi)
abline(scttr5b)


