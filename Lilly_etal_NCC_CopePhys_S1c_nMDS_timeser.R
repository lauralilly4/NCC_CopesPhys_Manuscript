########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 1c: Fig. -> Plot long-term timeseries of nMDS1 & 2 scores
##  Laura E. Lilly
##  Updated: 30 Mar 2025
########################################
# From saved file ('Lilly_etal_NCC_CopePhys_S1_nMDS)
# NOTE: Calculated seasonal cycle based on *monthly* avgs
#   because, while data were daily, they were unevenly sampled
#   throughout months across years

library(tidyverse)
library(mgcv)
library(tidygam)
library(imputeTS)


colsel <- readline("Which nMDS or spp?  ")

# ### Open datafile
### Option 1 - Cope spp.
scrfl <- read.csv('Biol_files_2025/NH05_CopeDens_log_subSpp_1996_2020_from_CAM.csv')
# ### Step 0: Convert file days -> R Dates & select species column
scrs_tib <- as_tibble(scrfl) |>
  mutate(Samp_Date = make_date(Year,Mon,Day))


# # ### Option 2 - nMDS
# scrfl <- read.csv('Biol_files_2025/NH05_Cope_nMDSscore_log10dens_Samples_k2_Jan2025.csv')
# colnames(scrfl) <- c("Samp_Date","NMDS1","NMDS2")
# ### Step 0: Convert file days -> R Dates & select species column
# scrs_tib <- as_tibble(scrfl) |>
#   mutate(Samp_Date = as.Date(Samp_Date, format = "%m/%d/%Y"),
#          NMDS1 = -(NMDS1),
#          NMDS2 = -(NMDS2))


# Select only column of intersts (plus Date)
scrs_cut <- scrs_tib[,c("Samp_Date",colsel)]
colnames(scrs_cut) <- c("Date","value")


# ### Step 1: Create daily-resolution time structure
stdt <- as.Date(paste(min(year(scrs_cut$Date)),01,01,sep="-"),format="%Y-%m-%d")
endt <- as.Date(paste(max(year(scrs_cut$Date)),12,31,sep="-"),format="%Y-%m-%d")
mstrdt <- seq(stdt,endt,by="day")

# Remove Leap Days -> to make every timeseries the same length
lpdt <- format(as.Date(paste(2012,02,29,sep='-')),"%m-%d")
mstrmd <- format(mstrdt,"%m-%d")
lpidx <- which(mstrmd %in% lpdt)
mstrdt <- mstrdt[-lpidx]


# ### Step 2: Slot actual copepod data points into larger timeseries
# ***If a sample exists for a Leap Day, pause to deal with it
mstrdayall <- as_tibble(data.frame(mstrdt,matrix(nrow=length(mstrdt),ncol=ncol(scrs_cut)-1)))
sppidx <- vector(length=nrow(scrs_cut))

for(dn in 1:nrow(scrs_cut)){
  mstrvec = mstrdayall[,1, drop = TRUE]
  scrval = scrs_cut[dn,1, drop = TRUE]
  
  midx = which(mstrvec == scrval)
  if(length(midx) == 0){
    print(scrs_cut[dn,1])
    continue = FALSE
  } else {
    continue = TRUE
  }
  mstrdayall[midx,] = scrs_cut[dn,]
  sppidx[dn] = midx
}
names(mstrdayall) =  names(scrs_cut)
mstrdayall$DOY <- yday(mstrdayall$Date)


# # Calculate monthly avg across all dates -> then remove seasonal cycle
# mstrdayall$month <- month(mstrdayall$Date)
# mstrdayall$year <- year(mstrdayall$Date)

# # Calculate "seasonal cycle" -> month avgs across years
# mon_avgs <- aggregate(value ~ month, 
#                         data = mstrdayall, 
#                         mean)




# # ### OPTION 1 - Trend Removal: Subtract month avg from each daily value for that month
# # ### ALT: Calculate a GAM and subtract that (basically the GAMs in Fig. 1)
# gm1 <- gam(value ~ s(DOY, bs = 'cc'), data = mstrdayall)
# model_gm1 <- predict_gam(gm1) |>
#   mutate(Date = as.Date(DOY, origin = "1996-01-01"),
#          value = fit)

# mstr_detrend <- data.frame(matrix(ncol = 1, nrow = nrow(mstrdayall)))
# 
# for(m in 1:nrow(mstrdayall)){
#   # mon = mstrdayall$month[m]
#   # ### 'Monthly avg' method
#   # monsub = mon_avgs[mon_avgs$month == mon,2]
#   # ### GAM method
#   monsub = model_gm1$value[which.min(abs(model_gm1$DOY-yday(mstrdayall$Date[m])))]
#   mondtr = mstrdayall$value[m]-monsub
#   
#   mstr_detrend[m,1] = mondtr
# }
# # Add Dates column to detrended vector
# mstr_detrend$Date = mstrdayall$Date
# colnames(mstr_detrend) <- c("value","Date")



# #### OPTION 2: Do actual decomposition (using R functions)
# Step 1: Combine any daily samples from each month to get a monthly average value
#   for each month -> to have a more complete timeseriess
mstr_grp <- mstrdayall |>
  mutate(Year = year(Date),
         Month = month(Date),
         Day = day(Date)) |>
  filter(Date >= as.Date("1996-05-01")) |>
  group_by(Year,Month) |>
  summarise(yrmo_avg = mean(value, na.rm = TRUE)) |>
  mutate(Date_mo = make_date(Year,Month,1))


ts_intrp <- na_interpolation(mstr_grp)
ts_fill <- ts(ts_intrp$yrmo_avg, start = c(1996,5), 
              end = c(2020,12),
              frequency = 12)
ts_decomp <- decompose(ts_fill)
lm_trend <- lm(ts_decomp$trend ~ mstr_grp$Date_mo)
summary(lm_trend)

# Combine four trend components into "long" DF, and add trendline for plotting
pts_df <- data.frame(mstr_grp$Date_mo,ts_decomp$x,ts_decomp$seasonal,
                     ts_decomp$trend,ts_decomp$random)
colnames(pts_df) <- c("Date","Data","Seasonal Cycle","Trend", "Residual")

long_df <- pts_df |>
  pivot_longer(cols = !(Date),
               values_to = "Value",
               names_to = "Type") |>
  arrange(Type)

trend_seq <- c(rep(NA,3*nrow(pts_df)),rep(NA,6),lm_trend$fitted.values, 
               rep(NA,6))

trend_df <- data.frame(long_df,trend_seq) |>
  mutate(LM_true = case_when(!is.na(trend_seq) ~ trend_seq,
                             TRUE ~ NA_real_))






#########################################################
# ### Step 4a: Plot entire 25-yr timeseries -> for trends

# ### PLOT 4.0: Plot raw timeseries of data -> to check against decomp
pltx0a <- ggplot(data = mstrdayall, aes(x = Date, y = value)) +
  geom_point(aes(x = Date, y = value))


# ### PLOT 4.1 - v2: Decomposed timeseries
pltx1a <- ggplot(trend_df, aes(x = Date, y = Value, group = Type)) +
  geom_line(aes(x = Date, y = Value, group = Type)) + 
  geom_smooth(aes(x = Date, y = trend_seq)) +

  annotate("text", label = "Adj. R^2 = 0.01", x = trend_df$Date[2000],
           y = -0.75, size = 4, colour = "black") +
  annotate("text", label = "p < 0.05", x = trend_df$Date[2000],
           y = -0.9, size = 4, colour = "black") +

  facet_wrap(~Type, scales = "free") + 

  theme(axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 12, colour = "black",
                                   angle = 90, vjust = 0.5, hjust = 1),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
  )

# ggsave(paste0("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/PX1_",colsel,"_LoadsTimeser_v2_decomp.png"), plot = pltx1a, width = 2400, height = 1600, units = 'px')








# # ### OLD VERSION
# # 4.1: NMDS1
# lm_1 <- lm(value ~ Date, data = mstr_detrend)
# summary(lm_1)
# 
# 
# pltx1a <- ggplot(mstr_detrend, aes(x = Date, y = value)) +
#   geom_point(aes(x = Date, y = value), size = 0.7, color = 'grey40') +
#   geom_smooth(method = "lm", se = FALSE) +
#   
#   annotate("text", label = "Adj. R^2 = 0.05", x = mstr_detrend$Date[2000], 
#            y = -0.75, size = 4, colour = "black") +
#   annotate("text", label = "p < 0.05", x = mstr_detrend$Date[2000], 
#            y = -0.9, size = 4, colour = "black") +
#   
#   scale_x_continuous(breaks = seq.Date(as.Date(paste(min(year(mstr_detrend$Date)),01,01,sep = "-")),
#                                        as.Date(paste(max(year(mstr_detrend$Date)),01,01,sep = "-")),365),
#                      labels = c("1996","","","","2000","","","",
#                                 "2004","","","","2008","","","",
#                                 "2012","","","","2016","","","","2020"),
#                      limits = c(as.Date(paste(min(year(mstr_detrend$Date)),01,01,sep = "-")),
#                                 as.Date(paste(max(year(mstr_detrend$Date)),12,31,sep = "-"))),
#                      name = "") +
#   scale_y_continuous(breaks = c(-1,-0.5,0,0.5,1),
#                      labels = c("-1","","0","","1"),
#                      limits = c(-1,1),
#                      name = "") + 
#   
#   theme(axis.title.y = element_text(size = 14, colour = "black"),
#         axis.text.y = element_text(size = 14, colour = "black"),
#         axis.title.x = element_text(size = 14, colour = "black"),
#         axis.text.x = element_text(size = 14, colour = "black", 
#                                    angle = 90, vjust = 0.5, hjust = 1),
#         panel.background = element_blank(),
#         axis.line.x = element_line(),
#         axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
#   )
# 
# # ggsave(paste0("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/PX1_",colsel,"_LoadsTimeser_detrend_neg.png"), plot = pltx1a, width = 2000, height = 1000, units = 'px')
# 
# 
