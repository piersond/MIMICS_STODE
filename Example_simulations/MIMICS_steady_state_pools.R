library(dplyr)
library(ggplot2)
library(DT)
library(purrr)
library(Metrics)

## Set working drive
setwd("C:/GitHub/MIMICS_STODE")


#----------------------------------------
# Load necessary model components
#----------------------------------------

# Load reverse Michaelis-Menton MIMICS function
source("functions/RXEQ.R")

# Load function to calculate model input variables
source("functions/calc_Tpars.R")

# Load function to calc MIMICS steady state pools
source("functions/MIMICS_calc_steady_state_pools.R")

# Set MIMICS parameters
#source("parameters/MIMICS_parameters_sandbox_20231129.R")  #--> Default "Sandbox" - Wieder et al. 2015
source("parameters/MIMICS_parameters_ReynoldsCreek_20240130.R")  #--> For Reynolds Creek. See Pierson et al. 2022


#-------------------------------------------------------
# Load forcing data: 
#  - requires columns: SITE, ANPP, CLAY, TSOI and LIG_N or LIG + CN
#  - optional columns: SOC, GWC, MAT
#-------------------------------------------------------

# Forcing data for forcing_data sites
#forcing_data <- read.csv("example_simulations/data/LTER_SITE_1.csv", as.is=T)
forcing_data <- read.csv("example_simulations/data/ReynoldsCreek_SOC_HiRes.csv", as.is=T, check.names = T)
colnames(forcing_data)[1] <- "Set"

#forcing_data$CLAY <- forcing_data$estCLAY
#forcing_data$ANPP <- forcing_data$ANPP


#---------------------------------------------------------------------------
# Example: Find steady state pools for a single site
#---------------------------------------------------------------------------

# Get MIMICS steady state (SS) pools for a single site in the forcing_data forcing dataset
MIMout_single_raw <- MIMICS_SS(forcing_data[1,])
MIMICS_ss_tbl <- MIMICS_SS_format(MIMout_single_raw)


#---------------------------------------------------------------------------
# Example: Find steady state pools for a forcing dataset
#---------------------------------------------------------------------------

MIMruns <- forcing_data %>% split(1:nrow(forcing_data)) %>% map(~ MIMICS_SS(df=.))
MIMICS_ss_dataset <- lapply(MIMruns, MIMICS_SS_format) %>% bind_rows()


#---------------------------------
# Plot field vs MIMICS SOC
#---------------------------------
plot_data <- MIMICS_ss_dataset

# Calc SOC vs. MIMSOC r2 and RMSE
r2_test <- cor.test(plot_data$SOC, plot_data$MIMSOC)
r_val <- round(as.numeric(unlist(r2_test ['estimate'])),2)
lb2 <- paste("R^2 == ", r_val)
rmse <- round(rmse(plot_data$SOC, plot_data$MIMSOC),2)

# Plot SOC vs. MIMSOC
ggplot(plot_data, aes(x=MIMSOC, y=SOC, color=TSOI)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  geom_point(size=4, alpha=0.8) +
  geom_text(aes(label=paste0(SITE)),hjust=-0.2, vjust=0.2) +
  annotate("text", label = lb2, x = 2, y = 8.5, size = 4, colour = "black", parse=T) +
  annotate("text", label = paste0("RMSE = ", rmse), x = 2, y = 7.4, size = 4, colour = "black") +
  ylim(0,10) + xlim(0,10) +
  theme_minimal()

