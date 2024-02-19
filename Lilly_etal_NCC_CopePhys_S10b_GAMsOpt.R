########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 5: Variable Coefficient GAMs, pt. 2
##            - Align phys. vars. & test for correlation
##  Laura E. Lilly
##  Updated: 18 May 2023
########################################
# Step 5.2: 1) Align all physical variables to the 'complete cases' (no NA)
#     of one variable. 2) Test for correlation between physical variables
# NOTE: IF correlation exists between two physical variables, do NOT 
#     use both of those variables in the same GAM


library(corrplot)

# ### Input copepod species file
copefl = read.csv('NH05_CopeDens_log_subSpp_1996_2020.csv')


####################################
# ### Step 0: Assign dates vector
spmo = copefl[,1]
spdy = copefl[,2]

modys = c(31,28,31,30,31,30,31,31,30,31,30,31)

spyrdy = vector()
for(m in 1:length(spmo)){
  mno = spmo[m]-1 # Subtract 1 because you only want number of *whole* months prior
  if(mno > 0){
    mdsum = sum(modys[1:mno])
  } else if (mno == 0){
    mdsum = 0
  }
  dsum = mdsum+spdy[m]
  spyrdy = c(spyrdy,dsum)
}


####################################
# ### Step 1: Align all variables by 'complete cases'
### Aligned by: 'cumumag'
magalng_rla = magcumudts[complete.cases(magcumudts) & magcumudts[,2] != 0,2]
magacrs_rla = magcumudts[complete.cases(magcumudts) & magcumudts[,2] != 0,3]
alf_trspr_rla = alf_sprtr[complete.cases(magcumudts) & magcumudts[,2] != 0,2]
nmdsrla = nmdsvar[complete.cases(magcumudts) & magcumudts[,2] != 0]
uicompa = uimtch[complete.cases(magcumudts) & magcumudts[,2] != 0,]
cutiarr_rla = cutimtch[complete.cases(magcumudts) & magcumudts[,2] != 0,]
beutiarr_rla = beutimtch[complete.cases(magcumudts) & magcumudts[,2] != 0,]
ssharr_rla = sshmtch[complete.cases(magcumudts) & magcumudts[,2] != 0,]
bvarr_rla = bvmtch[complete.cases(magcumudts) & magcumudts[,2] != 0,]

# Give names to individual UI variables
uisst_rla = uicompa$uisstrep
uiild_rla = uicompa$uiildrep
ssh_rla = ssharr_rla$sshvar
bv_rla = bvarr_rla$bvvar
cuti_rla = cutiarr_rla$cutivarrep
beuti_rla = beutiarr_rla$beutivarrep
yr_rla = uicompa$uiyrrep


# Assign year-days vector to nMDS scores (dim=1) vector
copespp = readline("Which species? [PSEUDO,OITSIM,CALMAR,ACALON,CENABD,CTNCAL,CLAARC,PARA,CALPAC,METR,ACATON]  ")
sppdens = copefl[,copespp]

spp_yrdy = cbind(copefl[,3],spyrdy,sppdens)
colnames(spp_yrdy) = c("Year","Yearday","Species")

spplt = spp_yrdy[,3] # Just assign a variable name to this vector, because I will use it a lot!
spyrdyplt = spp_yrdy[,2]
spyrsall = spp_yrdy[,1]
spyrsunq = unique(spyrsall)

sppdtarr = copefl[,1:3] # dataframe w/ extra col: 'Day' = 1 or 16
sppdts = as.Date(with(sppdtarr,paste(sppdtarr$Year,sppdtarr$Mon,sppdtarr$Day,sep="-")),"%Y-%m-%d")


# Get only 'yrdy' values corresponding to *non-NA* 'wndsprmtch'
spprls = spplt[complete.cases(magcumudts) & magcumudts[,2] != 0] 
spyrdyrls = spyrdyplt[complete.cases(magcumudts) & magcumudts[,2] != 0] 



####################################
# ### Step 2: Check covariability of all physical indices

# Option 1: Correlation table
varsmat = cbind(alf_trspr_rla,uisst_rla,ssh_rla,bv_rla,cuti_rla,beuti_rla)
varscorr = cor(varsmat,use='complete.obs')
corrplot(varscorr,method='number',type='upper',is.corr=FALSE,order='original')
pairs(varscorr)

# Option 2: Manual correlations (one vs. one) between phys inds
pind1 = beuti_rla # Index 1: CHANGE for each test
pind2 = uiild_rla # INDEX 2: CHANGE for each test
cor.test(pind1,pind2,method="pearson")
plot(pind1,pind2)



