########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 6: Variable Coefficient GAMs, pt. 1
##  Laura E. Lilly
##  Updated: 29 Jul 2024
########################################
# Load all physical datasets -> Then check correlation matrix to determine
#     acceptable combos for GAMS
library(lubridate)
library(dplyr)
library(mgcv)
library(plyr)
library(corrplot)


# ### Load physical files: 
#   Along & Across flows, Upwelling Indices (for SST, ILD), SSH,
#   Buoyancy Variability, CUTI, and BEUTI
flowfl <- read.csv('Phys_files_2024/NH10_AlongAcrossFlows_daily.csv')
uifl1 <- read.csv('Phys_files_2024/ROMS_MGJ_wcra_vars_44.65N_124.2W_1980-2010_monthly.csv')
uifl2 <- read.csv('Phys_files_2024/ROMS_MGJ_wcra_vars_44.65N_124.2W_2011-2021_monthly.csv')
sshfl <- read.csv('Phys_files_2024/NH05_SSH_from_AVISO.csv')
bvfl <- read.csv('Phys_files_2024/NH05_BV_from_NHL.csv')
cutifl <- read.csv('Phys_files_2024/CUTI_monthly.csv')
beutifl <- read.csv('Phys_files_2024/BEUTI_monthly.csv')


# scridx <- readline("Which nMDS dim? 1, 2    ") # Select nMDS dimension

################# Variables Configuration #################
# # First, create a daily timeseries -> to use for monthly-resolution variables
# daily_df <- tibble(Date = seq(as.Date(paste(1996,01,01,sep = "-")), 
#                               as.Date(paste(2020,12,31,sep = "-")), 
#                               "days")) |> 
#   mutate(month = month(Date),
#          .before = Date)

# Create Days-of-Month vector -> to use for monthly repeats
modys <- c(31,28,31,30,31,30,31,31,30,31,30,31)


# # 0) Reconfigure nMDS file to have *dates* (not date pieces)
# nmds_df <- scrfl |>
#   mutate(Date = as.Date(paste(Year,Mon,Day,sep = "-"),format = "%Y-%m-%d")) |>
#   select(Date,contains(scridx))

# 1a,b) Load Across & Along flows into DF
flow_df <- flowfl |>
  mutate(Date = as.POSIXct(Date,format="%d-%B-%Y %H:%M:%S"))

# 1c,d) Load SST and ILD from ROMS Upwelling Indices
upwl_df <- data.frame(rbind(uifl1,uifl2)) |>
  mutate(Date = as.Date(paste(year,month,1,sep = "-"),format = "%Y-%m-%d"),
         SST = `SST..C.`,
         ILD = `isothermal.layer.depth..m.`,
         month = month(Date)) |>
  select(Date,month,SST,ILD)

# Same for-loop as for CUTI and BEUTI
upwlyrrep <- vector()
upwlmorep <- vector()
upwldyrep <- vector()
upwlsstrep <- vector()
upwlildrep <- vector()

for(c in 1:length(upwl_df$month)){
  repno = modys[upwl_df$month[c]]
  
  upwlyrrep = c(upwlyrrep,rep(year(upwl_df$Date)[c],repno))
  upwlmorep = c(upwlmorep,rep(upwl_df$month[c],repno))
  upwldyrep = c(upwldyrep,seq(1,modys[upwl_df$month[c]],1))
  upwlsstrep = c(upwlsstrep,rep(upwl_df$SST[c],repno))
  upwlildrep = c(upwlildrep,rep(upwl_df$ILD[c],repno))
}
# Daily-resolution DF of CUTI
upwl_daily <- data.frame(Date = as.Date(paste(upwlyrrep,upwlmorep,upwldyrep,sep = "-")),
                         SST = upwlsstrep,
                         ILD = upwlildrep)


# 1e) Load SSH file
ssh_df <- sshfl |>
  mutate(Date = as.Date(paste(Year,Mon,Day,sep = "-"),format = "%Y-%m-%d"),
         SSH = sshtimser) |>
  select(Date,SSH)

# 1f) Load Buoyancy Freq. file
bv_df <- bvfl |>
  mutate(Date = as.Date(paste(Year,Mon,Day,sep = "-"),format = "%Y-%m-%d"),
         BV = bvcut) |>
  select(Date,BV)

# 1g) Load CUTI file and select latitude (44N and 45N -> then average)
cuti_df <- cutifl |>
  mutate(Date = as.Date(paste(year,month,1,sep = "-"),format = "%Y-%m-%d"),
         N44_5 = rowMeans(cbind(X44N,X45N), na.rm = TRUE),
         month = month(Date)) |>
  select(Date,month,N44_5) 

# Run for-loop to repeat each monthly value for the number of days in that month.
#     (It's a beast...but it works!)
cutiyrrep <- vector()
cutimorep <- vector()
cutidyrep <- vector()
cutivarrep <- vector()

for(c in 1:length(cuti_df$month)){
  repno = modys[cuti_df$month[c]]
  
  cutiyrrep = c(cutiyrrep,rep(year(cuti_df$Date)[c],repno))
  cutimorep = c(cutimorep,rep(cuti_df$month[c],repno))
  cutidyrep = c(cutidyrep,seq(1,modys[cuti_df$month[c]],1))
  cutivarrep = c(cutivarrep,rep(cuti_df$N44_5[c],repno))
}
# Daily-resolution DF of CUTI
cuti_daily <- data.frame(Date = as.Date(paste(cutiyrrep,cutimorep,cutidyrep,sep = "-")),
                         cuti44_5 = cutivarrep)



# 1h) Same thing for BEUTI
beut_df <- beutifl |>
  mutate(Date = as.Date(paste(year,month,1,sep = "-"),format = "%Y-%m-%d"),
         N44_5 = rowMeans(cbind(X44N,X45N), na.rm = TRUE),
         month = month(Date)) |>
  select(Date,month,N44_5) 

# For-loop for daily reps
beutyrrep <- vector()
beutmorep <- vector()
beutdyrep <- vector()
beutvarrep <- vector()

for(c in 1:length(beut_df$month)){
  repno = modys[beut_df$month[c]]
  
  beutyrrep = c(beutyrrep,rep(year(beut_df$Date)[c],repno))
  beutmorep = c(beutmorep,rep(beut_df$month[c],repno))
  beutdyrep = c(beutdyrep,seq(1,modys[beut_df$month[c]],1))
  beutvarrep = c(beutvarrep,rep(beut_df$N44_5[c],repno))
}
# Daily-resolution DF of CUTI
beut_daily <- data.frame(Date = as.Date(paste(beutyrrep,beutmorep,beutdyrep,sep = "-")),
                         beut44_5 = beutvarrep)


########### Combine variables and check correlations ############
phys_df_1 <- join_all(list(upwl_daily,ssh_df,bv_df,cuti_daily,beut_daily), 
                    by = 'Date', type = 'left')
phys_df <- left_join(phys_df_1, flow_df, by = 'Date') |>
  filter(Date >= as.Date("1996-01-01"),
         Date <= as.Date("2020-12-31"))
colnames(phys_df)[6:9] <- c("CUTI","BEUTI","AlongFlow","AcrossFlow")


varscorr = cor(phys_df[,2:ncol(phys_df)], use = 'complete.obs')
corrplot(varscorr,method='number',type='upper',is.corr=FALSE,order='original', col.lim = c(-1,1))
pairs(varscorr)


phys_real <- phys_df[complete.cases(phys_df[,2:ncol(phys_df)]),]

