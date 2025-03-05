#########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 2 - Species-level climatologies -> to use in Fig. 2
##  Laura E. Lilly
##  Updated: 3 Mar 2025
########################################
# Calculate species-level means and CIs to use for plotting Fig. 2


library(dplyr)
library(tidyverse)
library(lubridate)
# library(prettyGraphics)

# ### Input copepod species file
copefl = read.csv(paste0('Biol_files_2025/NH05_CopeDens_log_subSpp_1996_2020_from_CAM.csv'))
# meanfl = read.csv(paste0('NH05_CopeSpp_Climatology_Means_v2_PetersonGroups.csv'))
# cifl = read.csv(paste0('NH05_CopeSpp_Climatology_CIs_v2_PetersonGroups.csv'))

######### BOXPLOTS #########

########## OTHER SPECIES ##########
# # OPTION 1: Select species - single
# copespp = readline("Which species? [PSEUDO,CALMAR,ACALON,CENABD,METR,CALPAC,ACATON,OITSIM,NEOPLU]  ")
# sppcut = as.data.frame(cbind(copefl$Mon,copefl$Day,copefl$Year,copefl[,copespp]))
# colnames(sppcut) = c(colnames(copefl[1:3]),copespp)
########## OTHER SPECIES ##########


########## START OF CLAUSO ########## 
# # OPTION 2: Select all Clauso spp. -> to combine
# # Clausocalanus spp., Cl. arcuicornis, Cl. pergens, Cl. parapergens,
# #   Cl. paululus, Cl. lividus
clausospp_nms <- c("CLASO","CLAPER","CLAARC", "CLAPAR", "CLAPAU",
                   "CLALIV")
clauso_cut <- copefl |>
  select(Mon,Day,Year,all_of(clausospp_nms)) |>
  mutate(Date = as.Date(paste(as.numeric(Year),as.numeric(Mon),as.numeric(Day), sep = "-"), format = "%Y-%m-%d")) |>
  mutate(CLASO = (10^(CLASO-1))-0.1,
         CLAPER = (10^(CLAPER-1))-0.1,
         CLAARC = (10^(CLAARC-1))-0.1,
         CLAPAR = (10^(CLAPAR-1))-0.1,
         CLAPAU = (10^(CLAPAU-1))-0.1,
         CLALIV = (10^(CLALIV-1))-0.1) |>
  select(Date,clausospp_nms)


# Group by date and summarise across species
Clauso_sums <- clauso_cut |>
  pivot_longer(!Date, names_to = "Clauso_spp", values_to = "Dens") |>
  group_by(Date) |>
  summarise(Clauso_tot = sum(Dens))

sppcut <- Clauso_sums |>
  mutate(Year = year(Date),
         Mon = month(Date),
         Day = day(Date)) |>
  select(Mon,Day,Year,Clauso_tot)

########## END OF CLAUSO ########## 



# ### Step 1: Average all values within the same year-month (e.g., 1996-4)
yrsunq = unique(sppcut$Year)
mosunq = unique(sppcut$Mon)

yrmoavgs = data.frame(matrix(ncol=ncol(sppcut),nrow=length(yrsunq)*length(mosunq)))
di = 0

for(sy in 1:length(yrsunq)){
  for(sm in 1:length(mosunq)){
    di = di+1
    rowcut = sppcut[which(sppcut$Year %in% yrsunq[sy] & sppcut$Mon %in% mosunq[sm]),]
    if(nrow(rowcut) == 0){
      rowavg = c(mosunq[sm],1,yrsunq[sy],NA)
    } else if(nrow(rowcut) == 1) {
      rowavg = c(rowcut[1:3],rowcut[4])
    } else {
      rowavg = c(rowcut[1,1:3],mean(rowcut[,4],na.rm=TRUE))
    }
    yrmoavgs[di,] = rowavg
  }
}

colnames(yrmoavgs) = colnames(sppcut)
##### FOR CLAUSO SUM ONLY: log-transform means to then use for all calcs
yrmoavgs <- yrmoavgs |>
  mutate(avg_log10 = log10(Clauso_tot+0.1)+1) |>
  select(Mon,Day,Year,avg_log10)


# Then calculate yearly avgs (across all months) -> for climatologies of *means*
moavgs = data.frame(matrix(ncol=ncol(yrmoavgs)-2,nrow=length(mosunq)))
cimarg = data.frame(matrix(ncol=ncol(yrmoavgs)-2,nrow=length(mosunq)))
mosort = sort(mosunq)

for(mm in 1:length(mosort)){
  mavg = mean(yrmoavgs[(yrmoavgs[,1] %in% mosort[mm]),4],na.rm=TRUE)
  if(length(mavg) == 0){
    mavg = NA
    cimargin = NA
  } else{
  cimargin = qt(0.975,df=length(yrmoavgs[(yrmoavgs[,1] %in% mosort[mm]),4])-1)*(sd(yrmoavgs[(yrmoavgs[,1] %in% mosort[mm]),4],na.rm=TRUE))/sqrt(length(yrmoavgs[(yrmoavgs[,1] %in% mosort[mm]),4]))
  }
  moavgs[mm,] = c(mosort[mm],mavg)
  cimarg[mm,] = c(mosort[mm],cimargin)
}



# ### Calculate MIN, MAX, and DIFF between monthly means for each spp.
sppmin = min(meanfl[,which(colnames(meanfl) == paste0(copespp,"_mean"))])
sppmax = max(meanfl[,which(colnames(meanfl) == paste0(copespp,"_mean"))])
sppdif = sppmax/sppmin




# 
# ######################################
# # # ### Step 2: Plot boxplots by month
# 
# 
# ### CORRELATIONS
# crspp1 = meanfl[,8]
# crspp2 = meanfl[,9]
# cor.test(crspp1,crspp2,method="spearman")
# 
# ### Calculate CI values: hi and lo
# cilo = meanfl-cifl
# cihi = meanfl+cifl
# 
# corrcols = c("purple1","plum3","skyblue3","navyblue","royalblue","orangered1","gold1","orange2","yellowgreen","green2","forestgreen")
# 
# 
# 
# # ## COOL spp only
# colcl = c("purple2","royalblue","blue3","blue2")
# # colcl = rep("grey25",4) # GREY-only option
# symcl = c(15,16,17,5)
# dev.new(width=8,height=20)
# par(mar=c(6.5,5,3,4))
# 
# # First year -> to "establish" plot
# cis1 = data.frame(meanfl[,2],cilo[,2],cihi[,2])
# colnames(cis1) = c("fit","lowerCI","upperCI")
# plot(meanfl$Month,meanfl$PSEUDO_mean,type='l',col=colcl[1],ylim=c(0,6),cex=1,xlab='',ylab='',xaxt='n',lwd=1,cex.axis=1.8,cex.main=1.5,frame.plot=FALSE)
# axis(1,labels=FALSE,tick=FALSE)
# points(meanfl$Month,meanfl$PSEUDO_mean,pch=symcl[1],col=colcl[1],cex=1.5)
# add_error_envelope(meanfl$Month,ci = cis1,type = "poly",
#                    add_ci = list(col = scales::alpha(colcl[1], 0.1), border = FALSE),
#                    add_fit = list(col = colcl[1], lwd = 2, lty = 1))
# 
# # Subsequent years
# for(pu in 3:5){
#   cidf = data.frame(meanfl[,pu],cilo[,pu],cihi[,pu])
#   colnames(cidf) = c("fit","lowerCI","upperCI")
#   if(pu <= 4){
#     add_error_envelope(meanfl$Month,ci = cidf,type = "poly",
#                        add_ci = list(col = scales::alpha(colcl[pu-1], 0.1), border = FALSE),
#                        add_fit = list(col = colcl[pu-1], lwd = 2, lty = 1))
#   } else {
#     add_error_envelope(meanfl$Month,ci = cidf,type = "poly",
#                        add_ci = list(col = scales::alpha(colcl[pu-1], 0.1), border = FALSE),
#                        add_fit = list(col = colcl[pu-1], lwd = 1.8, lty = 5))
#   }
#   points(meanfl$Month,meanfl[,pu],pch=symcl[pu-1],col=colcl[pu-1],cex=1.5)
# }
# 
# axis(side=1,at=seq(1,12,1),labels=seq(1,12,1),cex.axis=1.7,las=1)
# nmscl = c("Pseudocalanus spp","C. marshallae","A. longiremis","Cp. abdominalis")
# 
# 
# 
# # ## WARM spp only
# colwm = rep(c("red3","orangered","orange3","goldenrod"),3)
# symwm = c(16,17,18,0,1,2,4,5,6,8)
# dev.new(width=8,height=20)
# par(mar=c(6.5,5,3,4))
# 
# # First year -> to "establish" plot
# cis1 = data.frame(meanfl[,7],cilo[,7],cihi[,7])
# colnames(cis1) = c("fit","lowerCI","upperCI")
# plot(meanfl$Month,meanfl$ACATON_mean,type='l',col=colwm[1],ylim=c(0,6),cex=1,xlab='',ylab='',xaxt='n',lwd=1,cex.axis=1.8,cex.main=1.5,frame.plot=FALSE)
# axis(1,labels=FALSE,tick=FALSE)
# points(meanfl$Month,meanfl$ACATON_mean,pch=symwm[1],col=colwm[1],cex=1.5)
# add_error_envelope(meanfl$Month,ci = cis1,type = "poly",
#                    add_ci = list(col = scales::alpha(colwm[1], 0.1), border = FALSE),
#                    add_fit = list(col = colwm[1], lwd = 1.5, lty = 1))
# 
# # Subsequent years
# for(pu in c(8:15,17)){
#   cidf = data.frame(meanfl[,pu],cilo[,pu],cihi[,pu])
#   colnames(cidf) = c("fit","lowerCI","upperCI")
#   if(pu <= 15){
#     add_error_envelope(meanfl$Month,ci = cidf,type = "poly",
#                      add_ci = list(col = scales::alpha(colwm[pu-6], 0.08), border = FALSE),
#                      add_fit = list(col = colwm[pu-6], lwd = 1.5, lty = 1))
#     points(meanfl$Month,meanfl[,pu],pch=symwm[pu-6],col=colwm[pu-6],cex=1.5)
#   } else {
#     add_error_envelope(meanfl$Month,ci = cidf,type = "poly",
#                      add_ci = list(col = scales::alpha(colwm[pu-7], 0.08), border = FALSE),
#                      add_fit = list(col = colwm[pu-7], lwd = 1.5, lty = 5))
#     points(meanfl$Month,meanfl[,pu],pch=symwm[pu-7],col=colwm[pu-7],cex=1.5)
#   }
# 
# }
# 
# axis(side=1,at=seq(1,12,1),labels=seq(1,12,1),cex.axis=1.7,las=1)
# nmswm = c("A. tonsa","C. pacificus","Co. styliremis","Co. tenuis","Clasocalanus spp.",
#           "Corycaeus anglicus","Ct. vanus","Mc. tenuicornis","Paracalanus spp",
#           "Cl. arcuicornis")
# 
