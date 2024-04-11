########################################
##  NOAA/OSU post-doc (NCC Copepods)
##  Step 8a (formerly 6b), v2 - Get dates of Inst and Cumu ALF transitions
##  Plots for Fig. 8 (ALF (raw, cumu) vs. BST date & mag)
##  Laura E. Lilly
##  Updated: 8 Apr 2024
########################################


library(tidyverse)
library(lubridate)
library(reshape2)



########## Data In ########## 
flwsfl <- read.csv('NH10_Flows_Inst_Cumu.csv')

# Dataset cleanup and date conversion
flws_data <- flwsfl |>
    mutate_all(~ifelse(is.nan(.), NA, .)) |>
    mutate(Date = as.Date(Date, format = "%d-%b-%Y"),
      Year = year(Date),
         .after = Date)


########## Reshape DFs and get info #########
# ### DF #1: Instantaneous flow
# First, reconfigure dates
inst_flw <- flws_data |>
  mutate(dt_mody = format(Date, format = "%m-%d"),
         # Will have to remove Leap Days if I want to add a repeat Date seq for plotting
         # plt_dt = rep(seq.Date(as.Date("1999-01-01"), 
         #                  as.Date("1999-12-31"),
         #                  by = "day")
         ) |>
  select(dt_mody,Year,Inst_flow)

inst_flw_yrly <- inst_flw |>
  pivot_wider(names_from = Year, values_from = Inst_flow)

# Second, get date of first negative value for each year
inst_neg1 <- data.frame(matrix(nrow = 1, ncol = ncol(inst_flw_yrly)-1))

for(i in 2:ncol(inst_flw_yrly)){
  idx1 = which(inst_flw_yrly[,i] < 0)[1]
  dt1 = inst_flw_yrly[idx1,1]
  
  inst_neg1[1,i-1] <- dt1
}
colnames(inst_neg1) <- unique(inst_flw$Year)





#############################
# ### DF #2: Cumulative flow
cumu_flw <- flws_data |>
  mutate(dt_mody = format(Date, format = "%m-%d")) |>
  select(dt_mody,Year,Cumu_flow)







########## Plots ##########

# ### Plot 1: Instantaneous flow - year-chunks
plt1 <- ggplot(data = inst_flw, aes(x = dt_mody, y = Inst_flow)) + 
  
  geom_line(aes(color = factor(Year), group = Year)) + 
  geom_segment(aes(x = dt_mody[1], xend = dt_mody[length(dt_mody)],
                   y = 0, yend = 0), color = "black", lwd = 1) + 
  
  ylab("Flow (m/s)") + 
  xlab("Date") 
  # scale_x_date(date_labels = "%b", date_breaks = "1 month")  




# ### Plot 2: Instantaneous flow - year-chunks
plt2 <- ggplot(data = cumu_flw, aes(x = dt_mody, y = Cumu_flow)) + 
  
  geom_line(aes(color = factor(Year), group = Year)) + 
  geom_segment(aes(x = dt_mody[1], xend = dt_mody[length(dt_mody)],
                   y = 0, yend = 0), color = "black", lwd = 1) + 
  
  ylab("Flow (m/s)") + 
  xlab("Date")



