#########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 2 - Species-level climatologies (Fig. 2 -> now Fig. 3
##  Laura E. Lilly
##  Updated: 15 May 2024
########################################
#   - Box plots (1996-2020)
#   - Spp-level climatologies


library(tidyverse)
library(reshape2)


# ### Input copepod species file
# copefl = read.csv(paste0('NH05_CopeDens_log_subSpp_1996_2020.csv'))
meanfl = read.csv(paste0('NH05_CopeSpp_Climatology_Means_v2_PetersonGroups.csv'))
cifl = read.csv(paste0('NH05_CopeSpp_Climatology_CIs_v2_PetersonGroups.csv'))



######### CALCULATE BOXPLOTS #########
# Select species
copespp = readline("Which species? [PSEUDO,CALMAR,ACALON,CENABD,METR,CALPAC,ACATON,OITSIM,NEOPLU]  ")
sppcut = as.data.frame(cbind(copefl$Mon,copefl$Day,copefl$Year,copefl[,copespp]))
colnames(sppcut) = c(colnames(copefl[1:3]),copespp)


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


# Then calculate yearly avgs (across all months) -> for climatologies of *means*
moavgs = data.frame(matrix(ncol=ncol(yrmoavgs)-2,nrow=length(mosunq)))
cimarg = data.frame(matrix(ncol=ncol(yrmoavgs)-2,nrow=length(mosunq)))
mosort = sort(mosunq)

for(mm in 1:length(mosort)){
  mavg = mean(yrmoavgs[(yrmoavgs[,1] %in% mosort[mm]),4],na.rm=TRUE)
  cimargin = qt(0.975,df=length(yrmoavgs[(yrmoavgs[,1] %in% mosort[mm]),4])-1)*(sd(yrmoavgs[(yrmoavgs[,1] %in% mosort[mm]),4],na.rm=TRUE))/sqrt(length(yrmoavgs[(yrmoavgs[,1] %in% mosort[mm]),4]))
  
  moavgs[mm,] = c(mosort[mm],mavg)
  cimarg[mm,] = c(mosort[mm],cimargin)
}


# ### Calculate MIN, MAX, and DIFF between monthly means for each spp.
sppmin = min(meanfl[,which(colnames(meanfl) == paste0(copespp,"_mean"))])
sppmax = max(meanfl[,which(colnames(meanfl) == paste0(copespp,"_mean"))])
sppdif = sppmax/sppmin


############# CALCS for CLIMATOLOGIES ############# 
# First, get Cool and Warm DFs -> to melt
# Mean
coolmns_df <- meanfl |>
  select(Month,PSEUDO_mean,CALMAR_mean,ACALON_mean,CENABD_mean) |>
  melt(id = "Month")

# CIs -> get two columns (min and max)
cismin <- cbind(meanfl[,1],meanfl[,2:ncol(meanfl)]-cifl[,2:ncol(cifl)])
colnames(cismin)[1] <- "Month"
cismax <- cbind(meanfl[,1],meanfl[,2:ncol(meanfl)]+cifl[,2:ncol(cifl)])
colnames(cismax)[1] <- "Month"

coolcis_min <- cismin |>
  select(Month,PSEUDO_mean,CALMAR_mean,ACALON_mean,CENABD_mean) |>
  melt(id = "Month")
  # rename(CI_min = value)
coolcis_max <- cismax |>
  select(Month,PSEUDO_mean,CALMAR_mean,ACALON_mean,CENABD_mean) |>
  melt(id = "Month") |>
  rename(CI_max = value)

coolcis_df <- cbind(coolcis_min,coolcis_max$CI_max)
colnames(coolcis_df)[4] <- "CI_max"


warmmns_df <- meanfl |>
 select(Month,ACATON_mean,CALPAC_mean,CALSTY_mean,CALTENU_mean,CLASO_mean,
        CORANG_mean,CTNCAL_mean,MESOCALTEN_mean,PARA_mean,CLAARC_mean) |>
 melt(id = "Month")

warmcis_min <- cismin |>
  select(Month,ACATON_mean,CALPAC_mean,CALSTY_mean,CALTENU_mean,CLASO_mean,
         CORANG_mean,CTNCAL_mean,MESOCALTEN_mean,PARA_mean,CLAARC_mean) |>
  melt(id = "Month")
warmcis_max <- cismax |>
  select(Month,ACATON_mean,CALPAC_mean,CALSTY_mean,CALTENU_mean,CLASO_mean,
         CORANG_mean,CTNCAL_mean,MESOCALTEN_mean,PARA_mean,CLAARC_mean) |>
  melt(id = "Month") |>
  rename(CI_max = value)

warmcis_df <- cbind(warmcis_min,warmcis_max$CI_max)
colnames(warmcis_df)[4] <- "CI_max"






############# PLOT CLIMATOLOGIES ###########
colcl = c("purple2","royalblue","blue3","blue2")
symcl = c(0,1,2,17)

# PLOT - Cool spp.
plt03_1 <- ggplot(data = coolmns_df, aes(x = Month, y = value, group = variable)) +

  geom_ribbon(data = coolcis_df, aes(x = Month, ymin = value, ymax = CI_max, group = variable, fill = variable), alpha = 0.1) +
  geom_line(aes(group = variable, color = variable)) + 
  geom_point(aes(group = variable, color = variable, shape = variable)) + 
  
  ylab("Log10(density)") + 
  scale_x_continuous(breaks = seq(month.name),
                     labels = month.name) + 
  scale_color_manual(name = "Species",
                     values = colcl,
                     labels = c("Pseudocalanus spp.", "C. marshallae", "A. longiremis","C. abdominalis")) +
  scale_fill_manual(name = "Species",
                   values = colcl,
                   labels = c("Pseudocalanus spp.", "C. marshallae", "A. longiremis","C. abdominalis")) + 
  scale_shape_manual(name = "Species",
                     values = symcl,
                     labels = c("Pseudocalanus spp.", "C. marshallae", "A. longiremis","C. abdominalis")) +
  
    theme(axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black", 
                                   angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top",
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
  )
# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/P3a_CoolSppClims.png", plot = plt03_1, width = 2000, height = 1600, units = 'px')




# PLOT - Warm spp.
colwm = rep(c("red3","orangered","orange3","goldenrod"),3)
symwm = c(16,17,18,0,1,2,4,5,6,8)

plt03_2 <- ggplot(data = warmmns_df, aes(x = Month, y = value, group = variable)) +
  
  geom_ribbon(data = warmcis_df, aes(x = Month, ymin = value, ymax = CI_max, group = variable, fill = variable), alpha = 0.1) +
  geom_line(aes(group = variable, color = variable)) + 
  geom_point(aes(group = variable, color = variable, shape = variable)) + 
  
  ylab("Log10(density)") + 
  scale_x_continuous(breaks = seq(month.name),
                     labels = month.name) + 
  scale_color_manual(name = "Species",
                     values = colwm,
                     labels = c("A. tonsa","C. pacificus","Co. styliremis","Co. tenuis",
                                "Clasocalanus spp.","Cr. anglicus","Ct. vanus","Mc. tenuicornis",
                                "Paracalanus spp.","Cl. arcuicornis")) +
  scale_fill_manual(name = "Species",
                    values = colwm,
                    labels = c("A. tonsa","C. pacificus","Co. styliremis","Co. tenuis",
                               "Clasocalanus spp.","Cr. anglicus","Ct. vanus","Mc. tenuicornis",
                               "Paracalanus spp.","Cl. arcuicornis")) +
  scale_shape_manual(name = "Species",
                     values = symwm,
                     labels = c("A. tonsa","C. pacificus","Co. styliremis","Co. tenuis",
                                "Clasocalanus spp.","Cr. anglicus","Ct. vanus","Mc. tenuicornis",
                                "Paracalanus spp.","Cl. arcuicornis")) +
  
  theme(axis.title.y = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black", 
                                   angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "top",
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
  )
# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/P3b_WarmSppClims.png", plot = plt03_2, width = 2000, height = 1600, units = 'px')
