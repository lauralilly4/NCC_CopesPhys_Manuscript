########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 5a: Calculate Percent Similarity Index (PSI) between years - by season
##  Laura E. Lilly
##  Updated: 14 Jun 2023
########################################
# Percent Similarity Index (PSI) of interannual community variability
#   within a season
# Seasonal dates were same as for Cluster Analysis:
#   - Biological Spring/Fall Transitions (BST/BFT) - nMDS crossover, visual examination
#       by C.A. Morgan
#   - Winter - 6 weeks prior to BST *for each year* (varies by year)
#   - Summer - 8 weeks after BST

# NOTE: Must run 'Lilly_etal_NCC_CopePhys_S4_ClusterSeason.R' *prior* 
#     to this script



`%!in%` <- Negate(`%in%`)


# ### Step 1: For each season's community (e.g., 'win_comm' from 'S4' code), 
#         Calculate PSIs between years

# Winter
winpcts = (win_comm/rowSums(win_comm))*100
win_psi = data.frame(matrix(nrow=length(dtyrs),ncol=length(dtyrs)))

for(yr in 1:length(dtyrs)){
  yidx = which(dtyrs %in% dtyrs[yr])
  yrpsis = vector() 
  for(yc in 1:length(dtyrs)){
    yrdfs = winpcts[yidx,]-winpcts[yc,]
    ypsi = 100-0.5*sum(abs(yrdfs))
    yrpsis = c(yrpsis,ypsi)
  }
  win_psi[yr,] = yrpsis
}
rownames(win_psi) = dtyrs
colnames(win_psi) = dtyrs

# Calculate 'avg PSI' for each year (aka the average of its PSIs with all other years)
win_avg_psi = colMeans(win_psi)


# Spring
sprpcts = (spr_comm/rowSums(spr_comm))*100
spr_psi = data.frame(matrix(nrow=length(dtyrs),ncol=length(dtyrs)))

for(yr in 1:length(dtyrs)){
  yidx = which(dtyrs %in% dtyrs[yr])
  yrpsis = vector() 
  for(yc in 1:length(dtyrs)){
    yrdfs = sprpcts[yidx,]-sprpcts[yc,]
    ypsi = 100-0.5*sum(abs(yrdfs))
    yrpsis = c(yrpsis,ypsi)
  }
  spr_psi[yr,] = yrpsis
}
rownames(spr_psi) = dtyrs
colnames(spr_psi) = dtyrs

# Calculate 'avg PSI' for each year (aka the average of its PSIs with all other years)
spr_avg_psi = colMeans(spr_psi)



# Summer
sumpcts = (sum_comm/rowSums(sum_comm))*100
sum_psi = data.frame(matrix(nrow=length(dtyrs),ncol=length(dtyrs)))

for(yr in 1:length(dtyrs)){
  yidx = which(dtyrs %in% dtyrs[yr])
  yrpsis = vector() 
  for(yc in 1:length(dtyrs)){
    yrdfs = sumpcts[yidx,]-sumpcts[yc,]
    ypsi = 100-0.5*sum(abs(yrdfs))
    yrpsis = c(yrpsis,ypsi)
  }
  sum_psi[yr,] = yrpsis
}
rownames(sum_psi) = dtyrs
colnames(sum_psi) = dtyrs

# Calculate 'avg PSI' for each year (aka the average of its PSIs with all other years)
sum_avg_psi = colMeans(sum_psi)


# Fall
falpcts = (fal_comm/rowSums(fal_comm))*100
fal_psi = data.frame(matrix(nrow=length(dtyrs),ncol=length(dtyrs)))

for(yr in 1:length(dtyrs)){
  yidx = which(dtyrs %in% dtyrs[yr])
  yrpsis = vector() 
  for(yc in 1:length(dtyrs)){
    yrdfs = falpcts[yidx,]-falpcts[yc,]
    ypsi = 100-0.5*sum(abs(yrdfs))
    yrpsis = c(yrpsis,ypsi)
  }
  fal_psi[yr,] = yrpsis
}
rownames(fal_psi) = dtyrs
colnames(fal_psi) = dtyrs

# Calculate 'avg PSI' for each year (aka the average of its PSIs with all other years)
fal_avg_psi = colMeans(fal_psi)



# ### STEP 3: Plot PSIs by year
## Color Scheme -> four categories
symc = c(15,16,17,18,19,3,4,8)
colel = "orangered"
colwm = "orange2"
colcd = "skyblue"
colla = "royalblue"
colnu = "grey50"

# Years: 96,97,98,99,00,01,02,
#        03,04,05,06,07,08,09,
#        10,11,12,13,14,15,16,
#        17,18,19,20
# ### V2
psisyms = c(symc[1],symc[1],symc[1],symc[1],symc[2],symc[1],symc[2],
            symc[2],symc[2],symc[3],symc[3],symc[3],symc[3],symc[4],
            symc[2],symc[4],symc[5],symc[4],symc[5],symc[4],symc[3],
            symc[6],symc[6],symc[7],symc[7])
psicols = c(colwm,colnu,colel,colla,colla,colcd,colcd,
            colwm,colnu,colwm,colcd,colnu,colla,colcd,
            colel,colla,colcd,colnu,colnu,colwm,colel,
            colnu,colcd,colnu,colcd)

# # Remove the '2015' and '2016' colors & symbols
# psisyms2 = psisyms[-(21:20)]
# psicols2 = psicols[-(21:20)]


# Plot WINTER
dev.new()
par(mar=c(6.5,5,3,4))

plot(as.numeric(dtyrs[1]),win_psi[1,2],typ='p',xaxt='n',xlab='',ylab='',bty='n',pch=psisyms[2],col=psicols[2],xlim=c(1996,2020),ylim=c(0,100),cex.axis=1.7)
for(p in 1:nrow(win_psi)){
  for(c in 1:ncol(win_psi)){
    if(is.na(win_psi[p,c]) == TRUE){
      next
    } else if(win_psi[p,c] == 100){
      points(as.numeric(dtyrs[p]),-3,pch=psisyms[c],col=psicols[c],cex=1.2)
    } else {
      points(as.numeric(dtyrs[p]),win_psi[p,c],pch=psisyms[c],col=psicols[c],cex=1)
    }
  }
}
axis(side=1,at=as.numeric(dtyrs),line=1,labels=dtyrs,cex.axis=1.7,las=2)


# Plot SPRING
dev.new()
par(mar=c(6.5,5,3,4))

plot(as.numeric(dtyrs[1]),spr_psi[1,2],typ='p',xaxt='n',xlab='',ylab='',bty='n',pch=psisyms[2],col=psicols[2],xlim=c(1996,2020),ylim=c(0,100),cex.axis=1.7)
for(p in 1:nrow(spr_psi)){
  for(c in 1:ncol(spr_psi)){
    if(is.na(spr_psi[p,c]) == TRUE){
      next
    } else if(spr_psi[p,c] == 100){
      points(as.numeric(dtyrs[p]),-3,pch=psisyms[c],col=psicols[c],cex=1.2)
    } else {
      points(as.numeric(dtyrs[p]),spr_psi[p,c],pch=psisyms[c],col=psicols[c],cex=1)
    }
  }
}
axis(side=1,at=as.numeric(dtyrs),line=1,labels=dtyrs,cex.axis=1.7,las=2)


# Plot SUMMER
dev.new()
par(mar=c(6.5,5,3,4))

plot(as.numeric(dtyrs[1]),sum_psi[1,2],typ='p',xaxt='n',xlab='',ylab='',bty='n',pch=psisyms[2],col=psicols[2],xlim=c(1996,2020),ylim=c(0,100),cex.axis=1.7)
for(p in 1:nrow(sum_psi)){
  for(c in 1:ncol(sum_psi)){
    if(is.na(sum_psi[p,c]) == TRUE){
      next
    } else if(sum_psi[p,c] == 100){
      points(as.numeric(dtyrs[p]),-3,pch=psisyms[c],col=psicols[c],cex=1.2)
    } else {
      points(as.numeric(dtyrs[p]),sum_psi[p,c],pch=psisyms[c],col=psicols[c],cex=1)
    }
  }
}
axis(side=1,at=as.numeric(dtyrs),line=1,labels=dtyrs,cex.axis=1.7,las=2)


# Plot FALL
dev.new()
par(mar=c(6.5,5,3,4))

plot(as.numeric(dtyrs[1]),fal_psi[1,2],typ='p',xaxt='n',xlab='',ylab='',bty='n',pch=psisyms[2],col=psicols[2],xlim=c(1996,2020),ylim=c(0,100),cex.axis=1.7)
for(p in 1:nrow(fal_psi)){
  for(c in 1:ncol(fal_psi)){
    if(is.na(fal_psi[p,c]) == TRUE){
      next
    } else if(fal_psi[p,c] == 100){
      points(as.numeric(dtyrs[p]),-3,pch=psisyms[c],col=psicols[c],cex=1.2)
    } else {
      points(as.numeric(dtyrs[p]),fal_psi[p,c],pch=psisyms[c],col=psicols[c],cex=1)
    }
  }
}
axis(side=1,at=as.numeric(dtyrs),line=1,labels=dtyrs,cex.axis=1.7,las=2)


