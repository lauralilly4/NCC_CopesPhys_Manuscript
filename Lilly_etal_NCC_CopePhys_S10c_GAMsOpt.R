########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 5: Variable Coefficient GAMs, pt. 3
##              - Evaluate for optimal GAM
##  Laura E. Lilly
##  Updated: 18 May 2023
########################################
# Step 5.3: Step 5.3: After determining which physical variables can 
#     be run in the same GAM (NO variables that correlate with each 
#     other), test all possible combinations of variables for each VC GAM
#     (uncomment individual variables to test different combinations).

# Variable Coefficient GAMs: These equations include a 'by' term, which
#   tests whether a habitat variable has a fluctuating influence on species
#   abundance throughout the year (e.g., 'Temperature *by* date-time)

# NOTE: RUN/EVALUATE EACH GAM ONE-BY-ONE (i.e., do not auto-run the
#   entire code)



####################################
# ### Step 3a: Evaluate all potential GAM combinations 
#               - *REGULAR GAMs*
gam1 = gam(spprls ~ 
             s(spyrdyrls,bs='cc') 
           # + s(magalng_rla,k=4)
           # + s(magacrs_rla,k=4)
           # + s(uisst_rla,k=4) 
           + s(ssh_rla,k=4) 
           # + s(bv_rla,k=4) 
           # + s(uiild_rla,k=4) 
           # + s(cuti_rla,k=4) 
           # + s(beuti_rla,k=4)
           + factor(alf_trspr_rla)
           # + factor(yr_rla)
           ,gamma=1.2
)
summary(gam1)
AIC(gam1)
plot(gam1)



####################################
# ### Step 3b: Evaluate all potential GAM combinations 
#               - *VC GAMs*

# ### 2.1: VC variable = Alongshore magnitude
vgam21 = gam(spprls ~
               s(spyrdyrls,bs='cc')
             + s(magalng_rla,k=4)
             # + s(magacrs_rla,k=4)
             # + s(uisst_rla,k=4)
             # + s(ssh_rla,k=4)
             # + s(uibv_rla,k=4)
             # + s(uiild_rla,k=4)
             # + s(cuti_rla,k=4)
             # + s(beuti_rla,k=4)
             + s(spyrdyrls,by=magalng_rla,bs='cc')
             # + factor(alf_trspr_rla)
             + factor(yr_rla)
             ,gamma=1.2
)
summary(vgam21)
AIC(vgam21)
plot(vgam21)


# ### 2.2: VC variable = Across-shore magnitude
vgam22 = gam(spprls ~
               s(spyrdyrls,bs='cc')
             + s(magalng_rla,k=4)
             # + s(magacrs_rla,k=4)
             # + s(uisst_rla,k=4)
             # + s(ssh_rla,k=4)
             # + s(bv_rla,k=4)
             # + s(uiild_rla,k=4)
             # + s(cuti_rla,k=4)
             # + s(beuti_rla,k=4)
             + s(spyrdyrls,by=magacrs_rla,bs='cc')
             # + factor(alf_trspr_rla)
             + factor(yr_rla)
             ,gamma=1.2
)
summary(vgam22)
AIC(vgam22)
plot(vgam22)


# ### 2.3: VC variable = SST
vgam23 = gam(spprls ~ 
               #  s(spyrdyrls) 
               # + s(magalng_rla,k=4)
               # + s(magacrs_rla,k=4)
               + s(uisst_rla,k=4) 
             # + s(ssh_rla,k=4) 
             # + s(bv_rla,k=4) 
             # + s(uiild_rla,k=4) 
             # + s(cuti_rla,k=4) 
             # + s(beuti_rla,k=4) 
             + s(spyrdyrls,by=uisst_rla,bs='cc')
             + factor(alf_trspr_rla)
             + factor(yr_rla)
             ,gamma=1.2
)
summary(vgam23)
AIC(vgam23)
plot(vgam23)


# ### 2.4: VC variable = SSH
vgam24 = gam(spprls ~ 
               s(spyrdyrls,bs='cc') 
             # + s(magalng_rla,k=4)
             + s(magacrs_rla,k=4)
             # + s(uisst_rla,k=4) 
             # + s(ssh_rla,k=4) 
             # + s(bv_rla,k=4) 
             # + s(uiild_rla,k=4) 
             # + s(cuti_rla,k=4) 
             # + s(beuti_rla,k=4)
             + s(spyrdyrls,by=ssh_rla,bs='cc')
             # + factor(alf_trspr_rla)
             + factor(yr_rla)
             ,gamma=1.2
)
summary(vgam24)
AIC(vgam24)
plot(vgam24)


# ### 2.5: VC variable = CUTI
vgam25 = gam(spprls ~ 
               s(spyrdyrls)
             # + s(magalng_rla,k=4)
             # + s(magacrs_rla,k=4)
             # + s(uisst_rla,k=4) 
             # + s(ssh_rla,k=4) 
             # + s(bv_rla,k=4) 
             # + s(uiild_rla,k=4) 
             # + s(cuti_rla,k=4) 
             # + s(beuti_rla,k=4)
             + s(spyrdyrls,by=cuti_rla)
             # + factor(alf_trspr_rla)
             + factor(yr_rla)
             ,gamma=1.2
)
summary(vgam25)
AIC(vgam25)
plot(vgam25)