########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 1a: Run nMDS of copepod spp.
##  Laura E. Lilly
##  Updated: 10 Mar 2024
########################################
# Data are from Newport Hydro Line (NHL):
#   - *Year-day dates*
#   - 36 major spp -> already log(0.1+dens)+1 transformed by C.A. Morgan
#   - includes data from JSOES cruises


# Load environment
library(tidyverse)
library(ggplot2)
library(vegan)  # For nMDS


# Load datafile
nhfl <- read.csv('Biol_files_2025/NH05_CopeDens_log_subSpp_1996_2020_from_CAM.csv')


# ### Step 1: Combine date columns and create new dataframe
dtfl <- data.frame(nhfl[,3],nhfl[,1],nhfl[,2])
colnames(dtfl) <- c(colnames(nhfl[3]),colnames(nhfl[1]),colnames(nhfl[2]))
dtvec <- as.Date(with(dtfl,paste(dtfl$Year,dtfl$Mon,dtfl$Day,sep="-")),"%Y-%m-%d")
sppsub <- data.frame(dtvec,nhfl[,5:ncol(nhfl)]) 


# ### Step 2: Run nMDS
# Get only rows (dates) with *some* non-zero spp
sppnonzro <- sppsub[rowSums(sppsub[,2:ncol(sppsub)])>0,2:ncol(sppsub)]

# Get corresponding 'Dates' vector (to reapply after nMDS)
dtsnonzro <- sppsub[rowSums(sppsub[,2:ncol(sppsub)])>0,1]


##### Average all dulicate values for each date 
dtsnonzrounq <- unique(dtsnonzro)

sppdtavg <- data.frame(matrix(nrow=length(dtsnonzrounq),ncol=ncol(sppnonzro)))
for(d in 1:length(dtsnonzrounq)){
  dtid <- which(dtsnonzro %in% dtsnonzrounq[d])
  sppavg <- colMeans(sppnonzro[dtid,])
  sppdtavg[d,] <- sppavg
}
rownames(sppdtavg) <- dtsnonzrounq
colnames(sppdtavg) <- colnames(sppnonzro)


# # Run nMDS calculations
# First, set up spp. data
sppnmslist <- c("Eucalanus","A.hudsoni","Pseudocalanus",
                        "O.similis","Ct.vanus","Cp.abdominalis",
                        "C.pacificus","C.marshallae","A.tonsa","A.longiremis",
                        "Paracalanus","O.spinirostris","Metridia","Scolecithricella",
                        "Oncaea","M.pusillus","Sc.minor","T.discaudatus",
                        "Ep.longipedata","Cr.anglicus","Clausocalanus",
                        "Cl.pergens","Co.tenuis","Co.styliremis",
                        "A.danae","Cl.arcuicornis","Cl.parapergens",
                        "M.tenuicornis","N.plumchrus","Microsetella",
                        "Lucicutia","Cl.paululus","Cd.bipinnata",
                        "Cl.lividus","Co.pavo","Co.pavon")

# Assign "Warm", "Cool", or "neutral" color to each spp name
sppnmscols <- c("grey50","grey50","purple3",
                "grey50","orange2","purple3",
                "orange2","purple3","orange2","purple3",
                "orange2","grey50","grey50","grey50",
                "grey50","grey50","grey50","grey50",
                "grey50","orange2","orange2",
                "grey50","orange2","orange2",
                "grey50","orange2","grey50",
                "orange2","grey50","grey50",
                "grey50","grey50","grey50",
                "grey50","grey50","grey50")

colnames(sppdtavg) <- sppnmslist

# Run nMDS with two different dimensions applied: k = 2, k = 3
spp_nmds2 <- metaMDS(sppdtavg,trymax=100,k=2) 
# spp_nmds3 <- metaMDS(sppdtavg,trymax=100,k=3) 



# ### Convert to DFs -> for ggplot (can't believe I'm saying this)

# k = 2
# First, nMDS scores (d = 1, 2) by *Species* -> for nMDS plots
# -> Add the color-codings for species names (for plotting)
sppnmds2df <- as.data.frame(scores(spp_nmds2$species)) |>
  mutate(Spp_col = sppnmscols)
# Second, get nMDS scores by *Sample Date* (to save for Fig. 1, Fig. 3 on
#   timeseries of community shifts, etc.)
sampnmds2df <- data.frame(scores(spp_nmds2$sites))


# # k = 3
# sppnmds3df <- as.data.frame(scores(spp_nmds3$species)) |>
#   mutate(Spp_col = sppnmscols)
# 
# sampnmds3df <- as.data.frame(scores(spp_nmds3$sites))



###################################################
# ### Step 4: Plot nMDS ordinations and save scores


# ### Dim:  k = 2
# Plot 4a: Stress plot
# png("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Plots_v4/P0_nMDS_StressPlot_k2.png", width = 1000, height = 800, units = 'px')
dev.new(width=8,height=20)
stressplot(spp_nmds2)
dev.off()


# Plot 4b: nMDS scores by species
plt1a <- ggplot(sppnmds2df, aes(x = MDS1, y = MDS2)) + 
  geom_point(size = 2, aes(color = Spp_col)) + 
  geom_text(aes(label = rownames(sppnmds2df), color = Spp_col), size = 3, nudge_y = 0.05, show.legend = FALSE) +
  
  labs(x = "NMDS1", colour = "Spp. group", y = "NMDS2", shape = "Type") + 
  
  scale_colour_manual(values = c("grey30","orange2","purple3"),
                      labels = c("Neutral","Warm","Cool")) +
  
  theme(axis.text = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "right", 
        axis.title = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        legend.key=element_blank())

# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/P1_nMDS_Ord_k2.png", plot = plt1a, width = 2000, height = 1600, units = 'px')


# # Save nMDS scores -> by Species
# write.csv(sppnmds2df,'NH05_CopeDens_log10_nMDSscr_Spp_k2.csv')
# # Save nMDS scores -> by Sample Dates
# write.csv(sampnmds2df,'NH05_CopeDens_log10_nMDSscr_Samp_k2.csv')





# ### Dim:  k = 3
# Plot 5a: Stress plot
# png("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/P0_nMDS_StressPlot_k3.png", width = 1000, height = 800, units = 'px')
# dev.new(width=8,height=20)
stressplot(spp_nmds3)
dev.off()


# Plot 5b: nMDS scores by species
plt1b <- ggplot(sppnmds3df, aes(x = MDS1, y = MDS2)) + 
  geom_point(size = 2, aes(color = Spp_col)) + 
  geom_text(aes(label = rownames(sppnmds3df), color = Spp_col), size = 3, nudge_y = 0.05) +
  
  labs(x = "NMDS1", colour = "Spp. group", y = "NMDS2", shape = "Type") + 
  
  scale_colour_manual(values = c("grey30","orange2","purple3"),
                      labels = c("Neutral","Cool","Warm")) +
  
  theme(axis.text = element_text(colour = "black", size = 12),
        legend.text = element_text(size = 12, colour ="black"),
        legend.position = "right", 
        axis.title = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.key=element_blank())

# ggsave("../../../OSU_NOAA_postdoc/Project1_SeasonalUpwelling/Figures/Plots_v4/P1_nMDS_Ord_k3.png", plot = plt1b, width = 2000, height = 1600, units = 'px')


# # Save nMDS scores -> by Species
# write.csv(sppnmds3df,'NH05_CopeDens_log10_nMDSscr_Spp_k3.csv')
# # Save nMDS scores -> by Sample Dates
# write.csv(sampnmds3df,'NH05_CopeDens_log10_nMDSscr_Samp_k3.csv')





