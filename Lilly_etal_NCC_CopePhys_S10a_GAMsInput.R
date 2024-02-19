########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 5: Variable Coefficient GAMs, pt. 1
##  Laura E. Lilly
##  Updated: 18 May 2023
########################################
# Step 5.1: Input physical datasets and format to a standard date-frame


library(lubridate)
library(dplyr)
library(mgcv)
library(TimeWarp)

# Load datafiles - Copepod nMDS score, Physical variables
scrfl <- read.csv('NH05_Cope_biom_MDSscore_v4_CAM_RawDts.csv')
wndfl <- read.csv('Phys_Inds/NH10_ALF_TrnsDts_Cumu.csv')
alffl2 <- read.csv('Phys_Inds/NH10_ALF_FlowMagHourly.csv')
romsfl1 <- read.csv('Phys_Inds/ROMS_MGJ_wcra_vars_44.65N_124.2W_1980-2010_monthly.csv')
romsfl2 <- read.csv('Phys_Inds/ROMS_MGJ_wcra_vars_44.65N_124.2W_2011-2021_monthly.csv')
sshfl <- read.csv('Phys_Inds/NH05_SSH_from_AVISO.csv')
bvfl <- read.csv('Phys_Inds/NH05_BV_from_NHL.csv')
cutifl <- read.csv('Phys_Inds/CUTI_monthly.csv')
beutifl <- read.csv('Phys_Inds/BEUTI_monthly.csv')


scridx <- readline("Which nMDS dim? 1, 2    ") # Select nMDS dimension


####################################
# ### Step 0a: Assign variable names & reconfigure (as needed)
# Copepod nMDS: Mon, Dy
nmdsyrall <- scrfl$Year
nmdsyrunq <- unique(nmdsyrall)
nmdsmoall <- scrfl$Mon
nmdsdyall <- scrfl$Day
nmdsdtarr <- scrfl[,1:3] # dataframe w/ extra col: 'Day' <- 1 or 16
nmdsdts <- as.Date(with(nmdsdtarr,paste(nmdsdtarr$Year,nmdsdtarr$Mon,nmdsdtarr$Day,sep="-")),"%Y-%m-%d")
nmdsvar <- scrfl[,as.numeric(scridx)+3] # nMDS scores

# ALF_SprTr
alfyrall <- wndfl[,1]
alfsprall <- wndfl[,2]  # Spring Transition
alffalall <- wndfl[,3]  # Fall Transition

# ## ALF_mag -> v2 (Hourly)
magyrall <- alffl2[,1]
magmoall <- alffl2[,2]
magdyall <- alffl2[,3]
maghrall <- alffl2[,4]
magalngall <- alffl2[,5]
magacrsall <- alffl2[,6]
magdtarr <- data.frame(cbind(magyrall,magmoall,magdyall))
colnames(magdtarr) <- c("Year","Mon","Day")
magdts <- as.Date(with(magdtarr,paste(magdtarr$Year,magdtarr$Mon,magdtarr$Day,sep="-")),"%Y-%m-%d")

# S2: Calculate *daily-avg* values from all hourly values
magallarr <- data.frame(magdts,magalngall,magacrsall)

magdlyalng <- magallarr %>%
  mutate(magdts = as.Date(magdts,format = "%Y-%m-%d")) %>% 
  group_by(magdts) %>%
  summarise(AvgDyAlng = mean(magalngall))

magdlyacrs <- magallarr %>%
  mutate(magdts = as.Date(magdts,format = "%Y-%m-%d")) %>% 
  group_by(magdts) %>%
  summarise(AvgDyAcrs = mean(magacrsall))


# Upwelling indices (from 2 ROMS files)
uiyrall <- c(romsfl1$year,romsfl2$year)
uimoall <- c(romsfl1$month,romsfl2$month)
uisstall <- c(romsfl1$SST..C.,romsfl2$SST..C.)
uiildall <- c(romsfl1$isothermal.layer.depth..m.,romsfl2$isothermal.layer.depth..m.)
uidtarr1 <- cbind(romsfl1[,1:2],rep(1,nrow(romsfl1)))
uidtarr2 <- cbind(romsfl2[,1:2],rep(1,nrow(romsfl2))) 
colnames(uidtarr1) <- c("Year","Mon","Day")
colnames(uidtarr2) <- c("Year","Mon","Day")
uidtarr <- rbind(uidtarr1,uidtarr2)
uidts <- as.Date(with(uidtarr,paste(uidtarr$Year,uidtarr$Mon,uidtarr$Day,sep="-")),"%Y-%m-%d")

# SSH
sshyrall <- sshfl$Year
sshmoall <- sshfl$Mon
sshdyall <- sshfl$Day
sshvarall <- sshfl$sshtimser

# BV
bvyrall <- bvfl$Year
bvmoall <- bvfl$Mon
bvdyall <- bvfl$Day
bvvarall <- bvfl$bvcut

# CUTI
cutiyrall <- cutifl$year
cutimoall <- cutifl$month
cuti44all <- cutifl$X44N
cuti45all <- cutifl$X45N
cutivarall <- rowMeans(cbind(cuti44all,cuti45all))

# BEUTI
beutiyrall <- beutifl$year
beutimoall <- beutifl$month
beuti44all <- beutifl$X44N
beuti45all <- beutifl$X45N
beutivarall <- rowMeans(cbind(beuti44all,beuti45all))


####################################
# ### Step 0b: Get subset of each Phys 
#     Index that matches years of Cope 
#     nMDS timeseries
# ALF_SprTr
alfyr <- alfyrall[which(alfyrall %in% nmdsyrall)]
alfspr <- alfsprall[which(alfyrall %in% nmdsyrall)]

# ALF_mag
magyr <- magyrall[which(magyrall %in% nmdsyrall)]
magmo <- magmoall[which(magyrall %in% nmdsyrall)]
magalng <- magalngall[which(magyrall %in% nmdsyrall)]

# Upwelling vars
uiyr <- uiyrall[which(uiyrall %in% nmdsyrall)]
uimo <- uimoall[which(uiyrall %in% nmdsyrall)]
uisst <- uisstall[which(uiyrall %in% nmdsyrall)]
uiild <- uiildall[which(uiyrall %in% nmdsyrall)]

# SSH
sshyr <- sshyrall[which(sshyrall %in% nmdsyrall)]
sshmo <- sshmoall[which(sshyrall %in% nmdsyrall)]
sshdy <- sshdyall[which(sshyrall %in% nmdsyrall)]
sshvar <- sshvarall[which(sshyrall %in% nmdsyrall)]
sshdtsmat <- data.frame(sshyr,sshmo,sshdy)
sshdts <- as.Date(with(sshdtsmat,paste(sshdtsmat$sshyr,sshdtsmat$sshmo,sshdtsmat$sshdy,sep="-")),"%Y-%m-%d")
sshidsall <- data.frame(sshyr,sshmo,sshdy,sshvar)

# BV
bvyr <- bvyrall[which(bvyrall %in% nmdsyrall)][-c(186,290)]
bvmo <- bvmoall[which(bvyrall %in% nmdsyrall)][-c(186,290)]
bvdy <- bvdyall[which(bvyrall %in% nmdsyrall)][-c(186,290)]
bvvar <- bvvarall[which(bvyrall %in% nmdsyrall)][-c(186,290)]
bvdtsmat <- data.frame(bvyr,bvmo,bvdy)
bvdts <- as.Date(with(bvdtsmat,paste(bvdtsmat$bvyr,bvdtsmat$bvmo,bvdtsmat$bvdy,sep="-")),"%Y-%m-%d")
bvidsall <- data.frame(bvyr,bvmo,bvdy,bvvar)

# CUTI
cutiyr <- cutiyrall[which(cutiyrall %in% nmdsyrall)]
cutimo <- cutimoall[which(cutiyrall %in% nmdsyrall)]
cutivar <- cutivarall[which(cutiyrall %in% nmdsyrall)]

# BEUTI
beutiyr <- beutiyrall[which(beutiyrall %in% nmdsyrall)]
beutimo <- beutimoall[which(beutiyrall %in% nmdsyrall)]
beutivar <- beutivarall[which(beutiyrall %in% nmdsyrall)]


####################################
# ### Step 1: Convert 'Dates' of nMDS 
#     scores (Copes) -> 'Yearday' values
modys <- c(31,28,31,30,31,30,31,31,30,31,30,31) # Number of days in 
                                              # each month -> to multiply by

nmdsyrdy <- vector()
for(m in 1:length(nmdsmoall)){
  mno = nmdsmoall[m]-1 # Subtract 1 - only want number of *whole* months prior
  if(mno > 0){
    mdsum = sum(modys[1:mno])
  } else if (mno == 0){
    mdsum = 0
  }
  dsum = mdsum+nmdsdyall[m]
  nmdsyrdy = c(nmdsyrdy,dsum)
}


####################################
# ### Step 2: Repeat monthly dates of 
#     sPhys Inds -> to daily resolution
#     - ALF_SprTr -> repeat to *yearly* res
#     - Upwelling Inds, CUTI, BEUTI
#     - NOT SSH and BV -> already have daily-res values

# ALF_SprTr -> *Yearly*
# Set up repeat sequence of each (yearly) Spring & Fall Trans Date (ALF/Wind)
alfarr <- data.frame(alfyr,alfspr) # Combine 'Years' & 'SprTr_dts' into one DF

alfsprrep <- data.frame(matrix(nrow=length(alfyr)*365,ncol=2))

for(a in 1:length(alfyr)){
  arep = matrix(rep(alfarr[a,],each=365),nrow=365)
  aidx = (365*(a-1))+1 # Create aidx for matrix indexing
  
  for(af in 1:365){
    afid = aidx+af-1
    alfsprrep[afid,] = arep[af,]
  }
}


# Upwelling Inds
uiyrrep <- vector()
uimorep <- vector()
uidyrep <- vector()
uisstrep <- vector()
uiildrep <- vector()

for(um in 1:length(uimo)){
  repno = modys[uimo[um]]
  
  uiyrrep = c(uiyrrep,rep(uiyr[um],repno)) 
  uimorep = c(uimorep,rep(uimo[um],repno)) 
  uidyrep = c(uidyrep,seq(1,repno,1))
  uisstrep = c(uisstrep,rep(uisst[um],repno)) 
  uiildrep = c(uiildrep,rep(uiild[um],repno)) 
}

uidtsmat <- data.frame(uiyrrep,uimorep,uidyrep)
uidts <- as.Date(with(uidtsmat,paste(uidtsmat$uiyrrep,uidtsmat$uimorep,uidtsmat$uidyrep,sep="-")),"%Y-%m-%d")

uiindsall <- data.frame(cbind(uiyrrep,uimorep,uidyrep,uisstrep,uiildrep))


# CUTI
cutiyrrep <- vector()
cutimorep <- vector()
cutidyrep <- vector()
cutivarrep <- vector()

for(c in 1:length(cutimo)){
  repno = modys[cutimo[c]]
  cutiyrrep = c(cutiyrrep,rep(cutiyr[c],repno)) 
  cutimorep = c(cutimorep,rep(cutimo[c],repno)) 
  cutidyrep = c(cutidyrep,seq(1,modys[cutimo[c]],1))
  cutivarrep = c(cutivarrep,rep(cutivar[c],repno)) 
}

cutidtsmat <- data.frame(cutiyrrep,cutimorep,cutidyrep)
cutidts <- as.Date(with(cutidtsmat,paste(cutidtsmat$cutiyrrep,cutidtsmat$cutimorep,cutidtsmat$cutidyrep,sep="-")),"%Y-%m-%d")
cutiidsall <- data.frame(cbind(cutiyrrep,cutimorep,cutidyrep,cutivarrep))


# BEUTI
beutiyrrep <- vector()
beutimorep <- vector()
beutidyrep <- vector()
beutivarrep <- vector()

for(b in 1:length(beutimo)){
  repno = modys[beutimo[b]]
  beutiyrrep = c(beutiyrrep,rep(beutiyr[b],repno)) 
  beutimorep = c(beutimorep,rep(beutimo[b],repno)) 
  beutidyrep = c(beutidyrep,seq(1,modys[beutimo[b]],1))
  beutivarrep = c(beutivarrep,rep(beutivar[b],repno)) 
}
beutidtsmat <- data.frame(beutiyrrep,beutimorep,beutidyrep)
beutidts <- as.Date(with(beutidtsmat,paste(beutidtsmat$beutiyrrep,beutidtsmat$beutimorep,beutidtsmat$beutidyrep,sep="-")),"%Y-%m-%d")
beutiidsall <- data.frame(cbind(beutiyrrep,beutimorep,beutidyrep,beutivarrep))


########################################
# ### Step 3: Match Phys Inds to nMDS dates
# ALF_SprTr
alf_sprtr <- vector()

for(al in 1:length(nmdsdts)){
  alid = which(alfsprrep[,1] %in% nmdsyrall[al])
  alf_sprtr = rbind(alf_sprtr,alfsprrep[alid[1],])
}

# UIs
uimtch <- vector() # 1-Yr, 2-Mo, 3-Dy, 4-SST, 5-ILD
for(u in 1:length(nmdsdts)){
  uid = dateMatch(nmdsdts[u],uidts,how=("nearest"))
  uimtch = rbind(uimtch,uiindsall[uid[1],])
}

# SSH
sshmtch <- vector()
for(s in 1:length(nmdsdts)){
  sid = dateMatch(nmdsdts[s],sshdts, how=("nearest"))
  sshmtch = rbind(sshmtch,sshidsall[sid[1],])
}

# BV
bvmtch <- vector()
for(b in 1:length(nmdsdts)){
  bid = dateMatch(nmdsdts[b],bvdts, how=("nearest"))
  bvmtch = rbind(bvmtch,bvidsall[bid[1],])
}

# CUTI
cutimtch <- vector()
for(cu in 1:length(nmdsdts)){
  cid = dateMatch(nmdsdts[cu],cutidts, how=("nearest"))
  cutimtch = rbind(cutimtch,cutiidsall[cid[1],])
}

# BEUTI
beutimtch <- vector()
for(bu in 1:length(nmdsdts)){
  bid = dateMatch(nmdsdts[bu],beutidts, how=("nearest"))
  beutimtch = rbind(beutimtch,beutiidsall[bid[1],])
}


########################################
# ### Step 4: *Cumul* ALF_mag 30 days 
#     prior to Cope nMDS
nocumul <- 14  # Number of cumul. days to sum
magdlyall <- data.frame(magdlyalng,magdlyacrs[,2])

magcumuall <- data.frame(matrix(nrow=length(nmdsdts),ncol=2))
for(n in 1:length(nmdsdts)){
  mdtid = which(magdlyall[,1] %in% nmdsdts[n])
  if(length(mdtid) == 0 || mdtid<nocumul){
    next
  } else{
    alngcum = sum(magdlyall[(mdtid-(nocumul-1)):mdtid,2],na.rm=TRUE)
    acrscum = sum(magdlyall[(mdtid-(nocumul-1)):mdtid,3],na.rm=TRUE)
    magcumuall[n,1] = alngcum
    magcumuall[n,2] = acrscum 
  }
}
magcumudts <- data.frame(nmdsdts,magcumuall)
colnames(magcumudts) <- c("Date","alng_mag","acrs_mag")


