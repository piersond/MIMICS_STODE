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

# Load function to simulate litter bag decomposition
source("functions/MIMICS_sim_litterbag.R")

# Set MIMICS parameters
source("parameters/MIMICS_parameters_sandbox_20231129.R")


#--------------------------------------
# Datasets for example simulations
#--------------------------------------

# Forcing data for LTER sites
LTER <- read.csv("example_simulations/data/LTER_SITE_1.csv", as.is=T)

# Import and calc LiDET litter bag data 
LiDET_BAGS <- data.frame(TYPE = c('TRAEf', 'PIREf','THPLf','ACSAf','QUPRf','DRGLf'),
                         BAG_LIG = c(16.2, 19.2, 26.7, 15.9, 23.5, 10.9),
                         BAG_N = c(0.38, 0.59, 0.62, 0.81, 1.03, 1.97), 
                         BAG_CN = c(133.3,92.7, 83.1, 61.8, 50.5, 24.2))

LiDET_BAGS$CALC_N <- (1 / LiDET_BAGS$BAG_CN) / 2.5 * 100   
LiDET_BAGS$CALC_MET <- 0.85 - 0.013 * LiDET_BAGS$BAG_LIG/LiDET_BAGS$CALC_N

BAG_init_size <- 100
BAGS <- LiDET_BAGS %>% select(TYPE, CALC_MET)
BAGS$BAG_LITm <- BAG_init_size * BAGS$CALC_MET  
BAGS$BAG_LITs <- BAG_init_size * (1-BAGS$CALC_MET)  


#---------------------------------------------
# Example: Run single litter bag decomp simulation
#---------------------------------------------

# Use an LTER site
#forcing_input <- LTER[6,]

# Build your own example
forcing_input <- data.frame(Site = "TEST SITE",
                            ANPP = 744,
                            TSOI = 15,
                            MAT = 25,
                            CLAY = 15,
                            LIG = 21,
                            N = 1.02,
                            CN = 49,
                            lig_N = 20.588,
                            grav.moisture = 35)

# Select a litter type from the LiDET example litter bags
BAG_input <- BAGS[6,]

# Run litterbag decomp simulation
BAG_out <- MIMICS_LITBAG(forcing_input, BAG_input, dailyInput=NA, nspin_yrs=10, nspin_days=150, litadd_day=125, verbose=T)

# Calc MIMICS C pools
BAG_out$MIMSOC <- rowSums(BAG_out[,6:12])
BAG_out$MIMMIC <- rowSums(BAG_out[,8:9])  
BAG_out$MIMLIT <- rowSums(BAG_out[,6:7])  
BAG_out$MICtoSOC <- BAG_out$MIMMIC/BAG_out$MIMSOC

# Plot litterbag decomp over time
ggplot(BAG_out) + 
  geom_line(aes(y=LITBAGm, x=DAY, colour="LITm"), size=1, linetype="dashed") +
  geom_line(aes(y=LITBAGs, x=DAY, colour="LITs"), size=1, linetype="dashed") +
  geom_line(aes(y=LITBAGs+LITBAGm, x=DAY, colour="Total"), size=1) +
  scale_colour_manual(values = c("Dark Green","Orange","Black")) +
  ylab("Litter Bag C") +
  xlab("Day") +
  labs(title = "MIMICS Litter Bag Simulation",
       subtitle = paste0("Site: ", forcing_input$Site, "\nLitter type: ", BAG_input[1]),
       colour = "LIT Pool") +
  theme_bw()


#-------------------------------------------------------------
# Example: Run multiple litter bags thru decomp simulation
#-------------------------------------------------------------

BAGS_out <- BAGS %>% split(1:nrow(BAGS)) %>% map(~ MIMICS_LITBAG(litBAG=.,
                                                                 forcing_df=LTER[6,],
                                                                 nspin_yrs=20,
                                                                 nspin_days=0,
                                                                 litadd_day=100,
                                                                 verbose=T)) %>% bind_rows()
ggplot() +
  geom_line(data=BAGS_out, aes(y=LITBAGs+LITBAGm, x=DAY, color=Litter_Type), linewidth=1, alpha=0.5) +
  ylab("Litter Bag C") +
  xlab("Day") +
  labs(color="LiDET\nLitter Type") +
  ggtitle("MIMICS Litter Bag Simulation") +
  theme_bw()

