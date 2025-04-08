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
  mutate(CLASO_raw = (10^(CLASO-1))-0.1,
         CLAPER_raw = (10^(CLAPER-1))-0.1,
         CLAARC_raw = (10^(CLAARC-1))-0.1,
         CLAPAR_raw = (10^(CLAPAR-1))-0.1,
         CLAPAU_raw = (10^(CLAPAU-1))-0.1,
         CLALIV_raw = (10^(CLALIV-1))-0.1) |>
  select(Date,CLASO_raw,CLAPER_raw,CLAARC_raw,CLAPAR_raw,CLAPAU_raw,CLALIV_raw)


# Group by date and summarise across species
Clauso_sums <- clauso_cut |>
  pivot_longer(!Date, names_to = "Clauso_spp", values_to = "Dens") |>
  group_by(Date) |>
  summarise(Clauso_tot = sum(Dens, na.rm = TRUE))


# new_files <- gsub('_samples\\.txt', '', files)

sppcut_raw <- Clauso_sums |>
  mutate(Year = year(Date),
         Mon = month(Date),
         Day = day(Date)) |>
  select(Mon,Day,Year,Clauso_tot)

##### Log-transform means to then use for all calcs
sppcut <- sppcut_raw |>
  mutate(avg_log10 = log10(Clauso_tot+0.1)+1,
         avg_log10_num = as.numeric(avg_log10)) |>
  select(Mon,Day,Year,avg_log10_num)

########## END OF CLAUSO ########## 



# ### Step 1: Average all values within the same year-month (e.g., 1996-4)
yrsunq = unique(sppcut$Year)
mosunq = sort(unique(sppcut$Mon))

yrmoavgs = data.frame(matrix(ncol=ncol(sppcut),nrow=length(yrsunq)*length(mosunq)))
di = 0

for(sy in 1:length(yrsunq)){
  for(sm in 1:length(mosunq)){
    di = di+1
    rowcut = sppcut[which(sppcut$Year %in% yrsunq[sy] & sppcut$Mon %in% mosunq[sm]),]
    if(nrow(rowcut) == 0){
      mnvl = NA
      rowavg = c(mosunq[sm],1,yrsunq[sy],mnvl)
    } else if(nrow(rowcut) == 1) {
      mnvl = rowcut[4]
      rowavg = c(rowcut[1:3],mnvl)
    } else {
      mnvl = mean(rowcut$avg_log10_num, na.rm = TRUE)
      rowavg = c(rowcut[1,1:3],mnvl)
    }
    yrmoavgs[di,] = rowavg
  }
}

colnames(yrmoavgs) = colnames(sppcut)



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


# ### Copy values over to 'Climatology_Means' and 'Climatology_CIs'


# ### Calculate MIN, MAX, and DIFF between monthly means for each spp.
sppmin = min(meanfl[,which(colnames(meanfl) == paste0(copespp,"_mean"))])
sppmax = max(meanfl[,which(colnames(meanfl) == paste0(copespp,"_mean"))])
sppdif = sppmax/sppmin

