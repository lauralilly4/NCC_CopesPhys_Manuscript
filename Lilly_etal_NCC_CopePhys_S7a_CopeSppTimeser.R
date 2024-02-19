###########    NOAA/OSU post-doc    ###########
###   Step 2.8e: Copepod spp. timeseries  ###
###############################################
# PURPOSE: 
#     - Plot long-term timeseries for copepod species and community
#     - Calculate long-term trends (or not)



library(prettyGraphics)
# Change WD to "R codes" folder
setwd('~/../../Documents/OSU_NOAA_Codes/R_codes')


# ### Input copepod species file
copefl = read.csv(paste0('../Datasets/Zoops_NHL_csv_v2_fromCAM/NH05_CopeDens_log_subSpp_1996_2020.csv'))


# Select species
copespp = readline("Which species? [PSEUDO,CALMAR,ACALON,CENABD,OITSIM,CTNCAL,PARA,CLAARC,CALPAC,METR,ACATON]  ")
sppcut = copefl[,copespp]
# colnames(sppcut) = c(colnames(copefl[1:3]),copespp)

sppdts_df = data.frame(copefl$Mon,copefl$Day,copefl$Year)
sppdts = as.Date(with(sppdts_df,paste(copefl$Year,copefl$Mon,copefl$Day,sep="-")),format="%Y-%m-%d")

sppcut[sppcut == 0] = NA


############  PLOTS  ############
dev.new()
par(mar=c(6.5,5,3,4))

plot(sppdts,sppcut, type='l',lwd=2,col='skyblue3',main = copespp,xlim=c(sppdts[1],sppdts[length(sppdts)]),ylim=c(0,6),cex.axis=1.7)

trln = lm(sppcut ~ sppdts)
summary(trln)
abline(trln)
# text(as.Date("2015-01-01"),5.5,"p-val = 0.72")
