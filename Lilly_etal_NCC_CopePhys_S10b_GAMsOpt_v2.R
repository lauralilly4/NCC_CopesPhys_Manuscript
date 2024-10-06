########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 5b: Variable Coefficient GAMs - Evaluate for optimal GAM
##  Laura E. Lilly
##  Updated: 6 Oct 2024
########################################
# Step 5.2: After determining which physical variables can 
#     be run in the same GAM (no variables that correlate with each 
#     other), test all possible combinations of variables for each VC GAM
#     (uncomment individual variables to test different combinations).
# Variable Coefficient GAMs: These equations include a 'by' term, which
#   tests whether a habitat variable has a fluctuating influence on species
#   abundance throughout the year (e.g., 'Temperature *by* date-time)
# NOTE: RUN/EVALUATE EACH GAM ONE-BY-ONE (i.e., do not auto-run the
#   entire code)
# Good GAM resource: https://r.qcbs.ca/workshop08/workshop08-en/workshop08-en.html#51

# RUN AFTER: 'Lilly_etal_NCC_CopePhys_S10a_GAMsIn_v2.R'

library(mgcViz)



############## File load & combine ##############
# ### CHOOSE & comment out

# # # ### Copepod nMDS scores
scrfl <- read.csv('NH05_Cope_biom_MDSscore_v4_CAM_RawDts.csv')
scridx <- readline("Which nMDS dim? 1, 2    ") # Select nMDS dimension
# Reconfigure nMDS file to have *dates* (not date pieces)
cope_df <- scrfl |>
  mutate(Date = as.Date(paste(Year,Mon,Day,sep = "-"),format = "%Y-%m-%d")) |>
  select(Date,contains(scridx))

# # ### OR
# # ### Copepod species
# copefl = read.csv('NH05_CopeDens_log_subSpp_1996_2020.csv')
# sppnm = readline("Which spp short name?   ")
# cope_df <- copefl |>
#   mutate(Date = as.Date(paste(Year,Mon,Day,sep = "-"),format = "%Y-%m-%d")) |>
#   select(Date, sppnm)

# Then combine Copepod DF with Phys Vars DF
gam_df <- left_join(phys_df,cope_df,by = 'Date') |>
  mutate(DOY = yday(Date),
         Year = year(Date))
colnames(gam_df)[10] <- "copespp"



################  RUN GAMs ################
# ### GAM 1: y ~ SST + SSH + BV + CUTI + Across_flow + 'by' terms
gam11 = gam(copespp ~ 
              s(DOY,bs='cc') + 
              s(SST, k = 4) + 
              s(SSH, k = 4) + 
              # s(BV, k = 4) + 
              # s(CUTI, k = 4) +
              s(AcrossFlow, k = 4) + 
              # s(DOY, by = SST, bs='cc') + 
              # s(DOY, by = SSH, bs='cc') + 
              # s(DOY, by = BV, bs='cc') + 
              # s(DOY, by = CUTI, bs='cc') + 
              s(DOY, by = AcrossFlow, bs='cc') + 
              factor(Year), 
            data = gam_df,
            # gamma=1.2
            )
summary(gam11)
AIC(gam11)
plot(gam11)


# ### GAM 2: y ~ SST + SSH + BV + BEUTI + Across_flow + 'by' terms
gam21 = gam(copespp ~ 
              s(DOY,bs='cc') + 
              s(SST, k = 4) + 
              s(SSH, k = 4) + 
              # s(BV, k = 4) + 
              # s(BEUTI, k = 4) +
              s(AcrossFlow, k = 4) + 
              # s(DOY, by = SST, bs='cc') + 
              # s(DOY, by = SSH, bs='cc') + 
              # s(DOY, by = BV, bs='cc') + 
              # s(DOY, by = BEUTI, bs='cc') + 
              s(DOY, by = AcrossFlow, bs='cc') + 
              factor(Year), 
            data = gam_df,
            # gamma=1.2
)
summary(gam21)
AIC(gam21)
plot(gam21)

gam21v <- getViz(gam21)
plt21_1 <- plot(sm(gam21v,1)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt21_2 <- plot(sm(gam21v,2)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt21_3 <- plot(sm(gam21v,3)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())



# ### GAM 3: y ~ SST + SSH + BV + AlongFlow + Across_flow + 'by' terms
gam31 = gam(copespp ~ 
              s(DOY,bs='cc') + 
              # s(SST, k = 4) + 
              s(SSH, k = 4) + 
              # s(BV, k = 4) + 
              s(AlongFlow, k = 4) +
              s(AcrossFlow, k = 4) + 
              # s(DOY, by = SST, bs='cc') + 
              # s(DOY, by = SSH, bs='cc') + 
              # s(DOY, by = BV, bs='cc') + 
              # s(DOY, by = AlongFlow, bs='cc') + 
              s(DOY, by = AcrossFlow, bs='cc') + 
              factor(Year), 
            data = gam_df,
            # gamma=1.2
)
summary(gam31)
AIC(gam31)
plot(gam31)

gam31v <- getViz(gam31)
plt31_1 <- plot(sm(gam31v,1)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt31_2 <- plot(sm(gam31v,2)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt31_3 <- plot(sm(gam31v,3)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())



# ### GAM 4: y ~ ILD + SSH + BV + AlongFlow + Across_flow + 'by' terms
gam41 = gam(copespp ~ 
              s(DOY,bs='cc') + 
              s(ILD, k = 4) + 
              s(SSH, k = 4) + 
              # s(BV, k = 4) + 
              s(AlongFlow, k = 4) +
              s(AcrossFlow, k = 4) + 
              # s(DOY, by = ILD, bs='cc') + 
              # s(DOY, by = SSH, bs='cc') + 
              # s(DOY, by = BV, bs='cc') + 
              # s(DOY, by = AlongFlow, bs='cc') + 
              s(DOY, by = AcrossFlow, bs='cc') + 
              factor(Year), 
            data = gam_df,
            # gamma=1.2
)
summary(gam41)
AIC(gam41)
plot(gam41)

gam41v <- getViz(gam41)
plt41_1 <- plot(sm(gam41v,1)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt41_2 <- plot(sm(gam41v,2)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt41_3 <- plot(sm(gam41v,3)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())








############################################################
################### OLD - DON'T USE ########################
############################################################

# ### GAM 1.3: y ~ ILD + SSH + BV + Along_flow + Across_flow
gam13 = gam(copespp ~ 
              s(DOY,bs='cc') + 
              s(ILD, k = 4) + 
              s(SSH, k = 4) + 
              # s(BV, k = 4) + 
              # s(AlongFlow, k = 4) + 
              s(AcrossFlow, k = 4) + 
              factor(Year), 
            data = gam_df,
            # gamma=1.2
)
summary(gam13)
AIC(gam13)
plot(gam13)
gam13v <- getViz(gam13)

plt13_1 <- plot(sm(gam13v,1)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt13_2 <- plot(sm(gam13v,2)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt13_3 <- plot(sm(gam13v,3)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt13_4 <- plot(sm(gam13v,4)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())



############# STEP 2: 'BY' GAMs
# Only evaluate 'by' for five terms: SST, SSH, CUTI, Along_flow, Across_flow

# ### GAM 2.1.1: y ~ SST + SSH + BV + CUTI + Along_flow + 
#                   Across_flow + (by = SST)
gam211 = gam(copespp ~ 
              s(DOY, bs = 'cc') + 
              # s(SST, k = 4) + 
              # s(SSH, k = 4) + 
              # s(BV, k = 4) + 
              # s(CUTI, k = 4) +
              # s(AlongFlow, k = 4) + 
              # s(AcrossFlow, k = 4) + 
              s(DOY, by = SST, bs='cc') + 
              factor(Year), 
            data = gam_df,
            # gamma=1.2
)
summary(gam211)
AIC(gam211)
gam211v <- getViz(gam211)

plt211_1 <- plot(sm(gam211v,1)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt211_2 <- plot(sm(gam211v,2)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt211_3 <- plot(sm(gam211v,3)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt211_4 <- plot(sm(gam211v,4)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())



# ### GAM 2.1.2: y ~ SST + SSH + BV + CUTI + Along_flow + 
#                   Across_flow + (by = SSH)
gam212 = gam(copespp ~ 
               s(DOY, bs = 'cc') + 
               # s(SST, k = 4) + 
               # s(SSH, k = 4) + 
               # s(BV, k = 4) + 
               # s(CUTI, k = 4) +
               # s(AlongFlow, k = 4) + 
               s(AcrossFlow, k = 4) + 
               s(DOY, by = SSH,bs='cc') + 
               factor(Year), 
             data = gam_df,
             # gamma=1.2
)
summary(gam212)
AIC(gam212)
# plot(gam212)
gam212v <- getViz(gam212)

plt212_1 <- plot(sm(gam212v,1)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt212_2 <- plot(sm(gam212v,2)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt212_3 <- plot(sm(gam212v,3)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())



# ### GAM 2.1.3: y ~ SST + SSH + BV + CUTI + Along_flow + 
#                   Across_flow + (by = CUTI)
gam213 = gam(copespp ~ 
               s(DOY, bs = 'cc') + 
               # s(SST, k = 4) + 
               s(SSH, k = 4) + 
               s(BV, k = 4) + 
               # s(CUTI, k = 4) +
               # s(AlongFlow, k = 4) + 
               s(AcrossFlow, k = 4) + 
               s(DOY, by = CUTI,bs='cc') + 
               factor(Year), 
             data = gam_df,
             # gamma=1.2
)
summary(gam213)
AIC(gam213)
# plot(gam213)
gam213v <- getViz(gam213)

plt213_1 <- plot(sm(gam213v,1)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt213_2 <- plot(sm(gam213v,2)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt213_3 <- plot(sm(gam213v,3)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt213_4 <- plot(sm(gam213v,4)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt213_5 <- plot(sm(gam213v,5)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())



# ### GAM 2.1.4: y ~ SST + SSH + BV + CUTI + Along_flow + 
#                   Across_flow + (by = Along_flow)
gam214 = gam(copespp ~ 
               s(DOY, bs = 'cc') + 
               s(SST, k = 4) + 
               # s(SSH, k = 4) + 
               # s(BV, k = 4) + 
               # s(CUTI, k = 4) +
               # s(AlongFlow, k = 4) + 
               # s(AcrossFlow, k = 4) + 
               s(DOY, by = AlongFlow,bs='cc') + 
               factor(Year), 
             data = gam_df,
             # gamma=1.2
)
summary(gam214)
AIC(gam214)
# plot(gam214)
gam214v <- getViz(gam214)

plt214_1 <- plot(sm(gam214v,1)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt214_2 <- plot(sm(gam214v,2)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt214_3 <- plot(sm(gam214v,3)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())
plt214_4 <- plot(sm(gam214v,4)) + 
  theme(axis.text = element_text(size = 32),
        axis.title = element_blank())




# ### GAM 2.1.5: y ~ SST + SSH + BV + CUTI + Along_flow + 
#                   Across_flow + (by = Across_flow)
gam215 = gam(copespp ~ 
               s(DOY, bs = 'cc') + 
               # s(SST, k = 4) + 
               # s(SSH, k = 4) + 
               # s(BV, k = 4) + 
               # s(CUTI, k = 4) +
               # s(AlongFlow, k = 4) + 
               # s(AcrossFlow, k = 4) + 
               s(DOY, by = AcrossFlow,bs='cc') + 
               factor(Year), 
             data = gam_df,
             # gamma=1.2
)
summary(gam215)
AIC(gam215)
plot(gam215)

