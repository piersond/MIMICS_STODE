library(purrr)
library(furrr)
library(dplyr)
library(tibble)
library(future)

#library(ggplot2)
#library(Metrics)

## Set working drive
setwd("C:/GitHub/MIMICS_STODE")

#----------------------------------------
# Load model components
#----------------------------------------

# Load reverse Michaelis-Menton MIMICS function
source("functions/RXEQ.R")

# Load function to calculate model input variables
source("functions/calc_Tpars.R")

# Load function to simulate litter bag decomposition
source("functions/MIMICS_sim_litterbag.R")

# Set MIMICS parameters
source("parameters/MIMICS_parameters_sandbox_20231129.R")

#------------------------------------------------------------

# Load a litter bag for this example
source("C:/github/MIMICS_STODE/Utilities/load_LiDET_litterbags.R")
litBAG_in <- BAGS[3,]

# Load MIMICS_SS_MC.R output 
ss_mc <- readRDS("MC_SS_output.rds")

# i = 1
# ss_mc[[i]][1]  # C pools
# ss_mc[[i]][2]  # Tpars
# ss_mc[[i]][3]  # forcings
# ss_mc[[i]][4]  # CO2 r, k
# ss_mc[[i]][5]  # Parameter multipliers

# FOR TESTING USE ONLY
#----------------------------------------------
source("Utilities/temp_curve_builder.R")
# Create the temperature vector
daily_TSOI <- generate_temp_curve(winter_low = 5, summer_high = 25, mean_temp = 15)

# Now build out the daily inputs, keeping all same except TSOI
dailyInputs <- data.frame(ANPP = rep(750, 365), 
                          CLAY = rep(30, 365), 
                          TSOI = daily_TSOI,
                          MAT = rep(15, 365), 
                          LIG_N = rep(20, 365),
                          CN = rep(30, 365),  
                          LIG = rep(15, 365),
                          GWC = rep(40, 365), 
                          theta_frzn = rep(0, 365))

rMIMICS_LITBAG <- function(MIMICS_SS_output, litterbag, daily_climate = NA){
  
  # Randomize parameters here
  Vslope_x <- MIMICS_SS_output[['Vslope_x']]
  Vslope <<- rep(0.063, 6) * Vslope_x
  
  output <- MIMICS_LITBAG(forcing_df=as.data.frame(MIMICS_SS_output[[3]]), 
                          SS_output=MIMICS_SS_output, 
                          litBAG=litterbag, 
                          dailyInput=daily_climate, 
                          loop_dailyInput=TRUE, 
                          nspin_yrs=3, 
                          nspin_days=100, 
                          litadd_day=100, 
                          verbose=F) %>% 
    tail(1) #<-- ONLY SAVING LAST ROW OF LITBAG SIM. REPLACE AS NEEDED !!!

  # Save parameter multipliers to output
  output <- as.data.frame(output)
  output$Vslope_x <- Vslope_x
  
  return(as.data.frame(output))
}

# Single use example
LITbag_out1 <- rMIMICS_LITBAG(ss_mc[[1]], litBAG_in, dailyInputs)


# Set up parallelism
no_cores <- availableCores() - 1 
plan(multisession, gc = TRUE, workers = no_cores) 
start_time <- Sys.time()

# Run it
litbag_results <- ss_mc %>% 
  future_map(~rMIMICS_LITBAG(.,litterbag = litBAG_in), 
             .progress = TRUE,
             .options = furrr_options(globals = TRUE)) %>% 
                      bind_rows()

# Release cores
plan(sequential)
print(paste0("Wall time: ", as.character(Sys.time() - start_time)))

### ~270 sec/1400 litterbag simulation on home CPU => 0.193 s/run ...3x speed on Derecho => 0.065 s/run => 55000 runs in ~1 hour 

# QC plot
library(ggplot2)
ggplot(litbag_results, aes(x=Vslope_x, y=LITBAGm+LITBAGs, color=SITE)) + geom_point() +theme_bw()
