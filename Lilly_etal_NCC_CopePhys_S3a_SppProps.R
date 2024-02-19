########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 3: Calculate proportions of copepod spp. for each biweekly period
##  Laura E. Lilly
##  Updated: 18 May 2023
########################################
# From raw values for each copepod spp. --> calculate proportions of spp
#   for visual examination and subsequent analyses
# NOTE: Proportions are based on *log-transformed* density values
#   because otherwise 2-3 spp. completely dominate community
# Do NOT need to run other scripts prior to this one


library(lubridate)


# ### Input copepod species file
copeflin <- read.csv('NH05_CopeDens_log_subSpp_1996_2020.csv')

copefl <- cbind(copeflin[,1:4],(10^(copeflin[,5:ncol(copeflin)]-1))-0.1)


######### Part 1: Set up 24-week timeseries (biweekly) #########
### Step 1: Check each month-day and re-categorize 
# into new vector as 1 or 16
dy2 <- vector()
for(d in 1:length(copefl$Day)){
  if (copefl[d,2] %in% seq(1,15)){
    d2 = 1
  } else if (copefl[d,2] %in% seq(16,31)){
    d2 = 16
  }
  dy2 = rbind(dy2,d2)    
}

# Add 'Dy2' column next to 'Day'
copetbl <- cbind(copefl[,3],copefl[,1:2],dy2,copefl[,5:ncol(copefl)]) # dataframe w/ extra col: 'Day' = 1 or 16
colnames(copetbl[,1:2]) <- c("Year","Mon")
copedtsnew <- as.Date(with(copetbl,paste(copetbl[,1],copetbl$Mon,copetbl$dy2,sep="-")),"%Y-%m-%d")
copes_nodts <- data.frame(copetbl[5:ncol(copetbl)]) # dataframe w/ NO dates


### Step 1.2: Average across all repeat samples for each date
copedtunq <- unique(copedtsnew)

copeunq <- data.frame()
for(dl in 1:length(copedtunq)){
  uids = which(copedtsnew %in% copedtunq[dl])
  dtavg = colMeans(copes_nodts[uids,],na.rm=TRUE)
  copeunq = rbind(copeunq,dtavg)
}

copeallunq <- data.frame(copedtunq,copeunq)
names(copeallunq) <- c("Date",colnames(copes_nodts))



### Step 1.3: Create timeseries of biweekly values for *every* month
yrsunq <- unique(copefl$Year)
dymstrall <- rep(c(1,16),(copetbl[nrow(copetbl),4]-copetbl[1,4]+1)*12) # Sequence of Days: 1, 16, 1, 16, etc.
momstr <- rep(seq(1,12,1),(yrsunq[length(yrsunq)]-yrsunq[1]+1))
momstrall <- rep(momstr,each=2) # Sequence of Months, repeated 2x each: 5,5,6,6,etc.
yrmstrall <- rep(seq(yrsunq[1],yrsunq[length(yrsunq)],1),each=24)


# Combine columns into dataframe
dtcmb <- data.frame(yrmstrall,momstrall,dymstrall)
# Convert columns -> "master" list of biweekly dates
dtmstr <- as.Date(with(dtcmb,paste(yrmstrall,momstrall,dymstrall,sep="-")),"%Y-%m-%d")


# Create empty data array to match DateMaster
mstrdata <- data.frame(matrix(nrow=length(dtmstr),ncol=ncol(copetbl)-4))

for(dn in 1:nrow(copeallunq)){
  midx = which(dtmstr %in% as.Date(copeallunq[dn,1],format="%y-%m-%d"))
  mstrdata[midx,] = copeallunq[dn,2:ncol(copeallunq)]
}
names(mstrdata) = names(copeallunq[2:ncol(copeallunq)])



# ### Step 2: Convert all raw values -> proportions
sppsums <- rowSums(mstrdata,na.rm=TRUE)
propsmat <- mstrdata/sppsums
propsmat[is.na(propsmat)] <- 0


### Step 3: Get yearlong cutout of props 
yrin <- readline(prompt = "Year to plot? (1996-2020)   ")
chnkidx <- which(format(dtmstr,format="%Y") %in% as.numeric(yrin))
pltdts <- dtmstr[chnkidx]
chnkmat <- propsmat[chnkidx,]

# ### Change all NA values -> zero
chnkzros <- chnkmat
chnkzros[is.na(chnkzros)] <- 0
chnktrns <- t(chnkzros)


## Step 4 - Plot proportions
# ### Define color set and plot proportions
colprps <- c("cyan","green2","yellow3","magenta","grey","navyblue","wheat","lavender","white",
            "tomato1","blue2","green4","lemonchiffon","red","peachpuff","saddlebrown","grey70",
            "blue4","olivedrab1","orange3","darkgreen","red1","springgreen2","bisque2","thistle3",
            "salmon2","royalblue3","plum2","mistyrose","orange1","saddlebrown",
            "skyblue1","red4","black","yellow1","plum4")


dev.new(width=8,height=20)
par(mar=c(6.5,5,3,4))
prpplt <- barplot(chnktrns,col=colprps,xaxt='n')
axis(side=1,at=prpplt,labels=format(pltdts,format="%m-%d"),cex.axis=1.7,las=2)

