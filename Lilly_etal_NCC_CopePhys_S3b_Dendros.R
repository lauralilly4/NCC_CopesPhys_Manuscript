########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 3b: Plot Clusters (dendrograms) of community composition 
#       similarity by season -> OLD figure
##  Laura E. Lilly
##  Updated: 10 Mar 2024
########################################
# Cluster plots (dendrograms) of interannual community variability
#   within a season
# Seasonal dates were determined by:
#   - Biological Spring/Fall Transitions - nMDS crossover, visual examination
#       by C.A. Morgan

# OLD setup for dates: variable increments for Winter vs. Summer, 
#   with multi-step process to find a date with values (i.e., if no
#   values are -6/+8 weeks, go forward/back 2 weeks)
#   - Winter - 6 weeks prior to BST (varies by year)
#   - Summer - 8 weeks after BST

# NOTE: Must run 'Lilly_etal_NCC_CopePhys_S3_SppProps.R' *prior* to this script


library(dendextend)
library(vegan)
`%!in%` <- Negate(`%in%`)


# ### Step 0: Create vectors of dates: 
# Spring (BST) - based on C.A. Morgan: nMDS, cluster, visual
# NOTE: 2015, 2016 dates are *artificial*, based on average BST dates of
#   all other years, because BST did not occur in those years
sprdts <- as.Date(c("1996-07-01","1997-05-01","1998-07-16",
                    "1999-05-01","2000-04-01","2001-03-16",
                    "2002-04-16","2003-06-01","2004-05-16",
                    "2005-08-16","2006-05-16","2007-03-16",
                    "2008-03-01","2009-03-01","2010-06-16",
                    "2011-03-16","2012-05-01","2013-04-01",
                    "2014-04-01","2015-05-01","2016-05-01",
                    "2017-07-01","2018-05-16","2019-06-01",
                    "2020-04-01"),format="%Y-%m-%d")

# Fall (BFT)
faldts <-  as.Date(c("1996-10-01","1997-08-16","1998-08-01",
                     "1999-11-01","2000-10-16","2001-11-01",
                     "2002-11-01","2003-10-01","2004-10-01",
                     "2005-09-16","2006-10-16","2007-10-01",
                     "2008-10-16","2009-12-01","2010-11-16",
                     "2011-09-16","2012-10-16","2013-09-16",
                     "2014-09-16","2015-10-01","2016-10-01",
                     "2017-10-01","2018-11-01","2019-09-01",
                     "2020-09-16"),format="%Y-%m-%d")
                    

# Get Year of each date
dtyrs <- format(sprdts,format="%Y")


# ### Step 1: Get rows of each community
# WINTER
win_comm <- data.frame(matrix(nrow=length(sprdts),ncol=ncol(mstrdata)))
win_dts <- as.Date(NA,length(sprdts))
for(p in 1:length(sprdts)){
  prid = which(dtmstr %in% sprdts[p])-3
  if(sum(is.na(mstrdata[prid,]))==ncol(mstrdata)){
    prid = which(dtmstr %in% sprdts[p])-4 #Caveat for years that don't have a sample 3 periods before
    if(sum(is.na(mstrdata[prid,]))==ncol(mstrdata)){
      prid = which(dtmstr %in% sprdts[p])-2 # Caveat for years of no -3,-4 samples (only 2008)
    }
  }
  win_comm[p,] = mstrdata[prid,]
  win_dts[p] = dtmstr[prid]
} 
colnames(win_comm) <- colnames(mstrdata)
rownames(win_comm) <- dtyrs

win_comm[is.na(win_comm)] <- 0
winclust <- win_comm[rowSums(win_comm[])>0,] # Keep only rows with *some* non-zeroes
# winclust2 <- winclust[-(21:20),] # REMOVE 2015, 2016 for WINTER


# Spring
spr_comm <- data.frame(matrix(nrow=length(sprdts),ncol=ncol(mstrdata)))
spr_dts <- as.Date(NA,length(sprdts))
for(s in 1:length(sprdts)){
  dtid = which(dtmstr %in% sprdts[s])
  spr_comm[s,] = mstrdata[dtid,]
  spr_dts[s] = dtmstr[dtid]
} 
colnames(spr_comm) <- colnames(mstrdata)
rownames(spr_comm) <- dtyrs

spr_comm[is.na(spr_comm)] <- 0
sprclust <- spr_comm[rowSums(spr_comm[])>0,] # Keep only rows with *some* non-zeroes


# SUMMER
sum_comm <- data.frame(matrix(nrow=length(sprdts),ncol=ncol(mstrdata)))
sum_dts <- as.Date(NA,length(sprdts))
for(p in 1:length(sprdts)){
  prid = which(dtmstr %in% sprdts[p])+4
  if(sum(is.na(mstrdata[prid,]))==ncol(mstrdata)){
    prid = which(dtmstr %in% sprdts[p])+5 # Caveat for years that don't have a sample 5 periods after
    if(sum(is.na(mstrdata[prid,]))==ncol(mstrdata)){
      prid = which(dtmstr %in% sprdts[p])+3 # Caveat for years of no +4,+5 samples (only 2008)
    }
  }
  sum_comm[p,] = mstrdata[prid,]
  sum_dts[p] = dtmstr[prid]
} 
colnames(sum_comm) <- colnames(mstrdata)
rownames(sum_comm) <- dtyrs

sum_comm[is.na(sum_comm)] <- 0
sumclust <- sum_comm[rowSums(sum_comm[])>0,] # Keep only rows with *some* non-zeroes


# FALL
fal_comm <- data.frame(matrix(nrow=length(faldts),ncol=ncol(mstrdata)))
fal_dts <- as.Date(NA,length(faldts))
for(f in 1:length(faldts)){
  dtid = which(dtmstr %in% faldts[f])
  fal_comm[f,] = mstrdata[dtid,]
  fal_dts[f] = dtmstr[dtid]
} 
colnames(fal_comm) <- colnames(mstrdata)
rownames(fal_comm) <- dtyrs

fal_comm[is.na(fal_comm)] <- 0
falclust <- fal_comm[rowSums(fal_comm[])>0,] # Keep only rows with *some* non-zeroes




### STEP 3: Calculate cluster distances
# NOTE: input matrix: rows = samples, columns = variables -> must TRANSFORM
#       if you want to group/compare *species* across biweeklys

# # METHOD 2: Bray-Curtis similarities -> Do NOT scale!
windist <- vegdist(winclust, method = 'bray')
# windist2 <- vegdist(winclust2, method = 'bray')
sprdist <- vegdist(sprclust, method = 'bray')
sumdist <- vegdist(sumclust, method = 'bray')
faldist <- vegdist(falclust, method = 'bray')


# Calculate cluster analysis from distance matrix ('distmat')
winclustavg <- hclust(windist, method = 'average')
# winclustavg2 <- hclust(windist2, method = 'average')
sprclustavg <- hclust(sprdist, method = 'average')
sumclustavg <- hclust(sumdist, method = 'average')
falclustavg <- hclust(faldist, method = 'average')



### STEP 4: Plot Dendrogram
# WINTER

## Option 1: w/ 2015, 2016
windendavg = as.dendrogram(winclustavg)
windendavg = set_labels(windendavg,as.character(rownames(winclust)[winclustavg$order]))
labels_cex(windendavg) = 1.4
windendavg = color_labels(windendavg,col="grey50",labels=labels(windendavg)[c(2,6,8,10,17,19,20)])
windendavg = color_labels(windendavg,col="orangered",labels=labels(windendavg)[c(11,14,15)])
windendavg = color_labels(windendavg,col="orange",labels=labels(windendavg)[c(7,9,12,13)])
windendavg = color_labels(windendavg,col="skyblue2",labels=labels(windendavg)[c(1,3,4,5,18,21,24)])
windendavg = color_labels(windendavg,col="royalblue2",labels=labels(windendavg)[c(16,22,23,25)])
dev.new(width=8,height=20)
par(mar=c(12,5,3,4))
plot(windendavg,cex.axis=1.3,ylim=c(0,0.8))

# # ## Option 2: NO 2015, 2016
# windendavg <- as.dendrogram(winclustavg2)
# windendavg <- set_labels(windendavg,as.character(rownames(winclust2)[winclustavg2$order]))
# labels_cex(windendavg) = 1.4
# windendavg <- color_labels(windendavg,col="grey50",labels=labels(windendavg)[c(2,10,11,12,19,20,23)])
# windendavg <- color_labels(windendavg,col="orangered",labels=labels(windendavg)[c(6,7)])
# windendavg <- color_labels(windendavg,col="orange",labels=labels(windendavg)[c(5,21,22)])
# windendavg <- color_labels(windendavg,col="skyblue2",labels=labels(windendavg)[c(1,3,4,8,16,17,18)])
# windendavg <- color_labels(windendavg,col="royalblue2",labels=labels(windendavg)[c(9,13,14,15)])
# dev.new(width=8,height=20)
# par(mar=c(12,5,3,4))
# plot(windendavg,cex.axis=1.3,ylim=c(0,0.8))


# ## SPRING (BST)
sprdendavg <- as.dendrogram(sprclustavg)
sprdendavg <- set_labels(sprdendavg,as.character(rownames(sprclust)[sprclustavg$order]))
labels_cex(sprdendavg) = 1.4
sprdendavg <- color_labels(sprdendavg,col="grey50",labels=labels(sprdendavg)[c(6,9,13,17,21,22,25)])
sprdendavg <- color_labels(sprdendavg,col="orangered",labels=labels(sprdendavg)[c(1,4,20)])
sprdendavg <- color_labels(sprdendavg,col="orange",labels=labels(sprdendavg)[c(2,3,10,15)])
sprdendavg <- color_labels(sprdendavg,col="skyblue2",labels=labels(sprdendavg)[c(5,7,14,16,18,19,23)])
sprdendavg <- color_labels(sprdendavg,col="royalblue2",labels=labels(sprdendavg)[c(8,11,12,24)])
dev.new(width=8,height=20)
par(mar=c(12,5,3,4))
plot(sprdendavg,cex.axis=1.3,ylim=c(0,0.8))


# ## SUMMER
sumdendavg <- as.dendrogram(sumclustavg)
sumdendavg <- set_labels(sumdendavg,as.character(rownames(sumclust)[sumclustavg$order]))
labels_cex(sumdendavg) = 1.4
sumdendavg <- color_labels(sumdendavg,col="grey50",labels=labels(sumdendavg)[c(2,3,4,18,21,22,23)])
sumdendavg <- color_labels(sumdendavg,col="orangered",labels=labels(sumdendavg)[c(1,5,25)])
sumdendavg <- color_labels(sumdendavg,col="orange",labels=labels(sumdendavg)[c(8,14,15,24)])
sumdendavg <- color_labels(sumdendavg,col="skyblue2",labels=labels(sumdendavg)[c(6,7,12,13,16,17,20)])
sumdendavg <- color_labels(sumdendavg,col="royalblue2",labels=labels(sumdendavg)[c(9,10,11,19)])
dev.new(width=8,height=20)
par(mar=c(12,5,3,4))
plot(sumdendavg,cex.axis=1.3,ylim=c(0,0.8))


# ## FALL (BFT)
faldendavg <- as.dendrogram(falclustavg)
faldendavg <- set_labels(faldendavg,as.character(rownames(falclust)[falclustavg$order]))
labels_cex(faldendavg) = 1.4
faldendavg <- color_labels(faldendavg,col="grey50",labels=labels(faldendavg)[c(1,2,5,14,15,18,21,23)])
faldendavg <- color_labels(faldendavg,col="orangered",labels=labels(faldendavg)[c(4,10,22)])
faldendavg <- color_labels(faldendavg,col="orange",labels=labels(faldendavg)[c(3,11,17,24)])
faldendavg <- color_labels(faldendavg,col="skyblue2",labels=labels(faldendavg)[c(6,8,9,12,13,19,20)])
faldendavg <- color_labels(faldendavg,col="royalblue2",labels=labels(faldendavg)[c(1,7,16,25)])
dev.new(width=8,height=20)
par(mar=c(12,5,3,4))
plot(faldendavg,cex.axis=1.3,ylim=c(0,0.8))

# # ## FALL (BFT), v2: NO 2015, 2016
# faldendavg <- as.dendrogram(falclustavg)
# faldendavg <- set_labels(faldendavg,as.character(rownames(falclust)[falclustavg$order]))
# labels_cex(faldendavg) = 1.4
# faldendavg <- color_labels(faldendavg,col="grey50",labels=labels(faldendavg)[c(2,10,11,14,17,19,23)])
# faldendavg <- color_labels(faldendavg,col="orangered",labels=labels(faldendavg)[c(1,18)])
# faldendavg <- color_labels(faldendavg,col="orange",labels=labels(faldendavg)[c(7,13,20)])
# faldendavg <- color_labels(faldendavg,col="skyblue2",labels=labels(faldendavg)[c(3,5,6,8,9,15,16)])
# faldendavg <- color_labels(faldendavg,col="royalblue2",labels=labels(faldendavg)[c(4,12,21,22)])
# dev.new(width=8,height=20)
# par(mar=c(12,5,3,4))
# plot(faldendavg,cex.axis=1.3,ylim=c(0,0.8))

