library(furrr)
library(dplyr)

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

# Load function to apply random parameters and run MIMICS steady state simulation
source("functions/MC_parameterization/MIMICS_repeat.R")

# Set MIMICS parameters
source("parameters/MIMICS_parameters_sandbox_20231129.R")

# Set default parameters for MC to remember
source("functions/MC_parameterization/set_parameter_defaults.R")

#-------------------------------------------------------
# Load forcing data
forcing_data <- read.csv("example_simulations/data/LTER_SITE_1.csv", as.is=T)

# Prep forcing data
colnames(forcing_data) <- toupper(colnames(forcing_data))
forcing_data$FID <- seq(1,nrow(forcing_data),1)

####################################
# Use the brute force MIMICS ftn
####################################

# Set desired number of random parameter runs
MIM_runs <- 10000

### Create random parameter dataframe
#------------------------------------------------------------------------------
#  !!! NOTE, when adding or removing parameters, must setup coressponding 
#      parameter values here, in MIMrepeat() and in "set_parameter_defaults.R

rand_params <- data.frame(Vslope_x = runif(MIM_runs,1, 1),   #0.5, 2
                          Vint_x = runif(MIM_runs, 1, 1)#,  #0.8, 1.3
                          # Kslope_x = runif(MIM_runs, 0.2, 1.2),#,  #0.5, 2  
                          # Kint_x = runif(MIM_runs, 0.2, 1.2)#, #0.5,2 
                          # Vslope_MOD_x1 = runif(MIM_runs, 0.7, 1.3),
                          # Vslope_MOD_x2 = runif(MIM_runs, 0.7, 1.3),
                          # Vint_MOD_x1 = runif(MIM_runs, 0.7, 1.3),
                          # Vint_MOD_x2 = runif(MIM_runs, 0.7, 1.3)
                          # Tau_x = runif(MIM_runs, 0.3, 3),  
                          # CUE_x = runif(MIM_runs, 0.5, 1.5),  
                          # desorb_x = runif(MIM_runs, 0.001, 0.3),  
                          # fPHYS_x = runif(MIM_runs, 0.4, 4),
                          # fW_p1_x = runif(MIM_runs, 0.5, 1.5),
                          # fW_p2_x = runif(MIM_runs, 0.5, 1.5)
  
)

# Add ID# for each parameter set 
rand_params$run_num <- seq(1,MIM_runs,1)


#------------------------------------
# Begin parallel processing
#------------------------------------

# Set number of cores to use
no_cores <- availableCores() - 1
plan(multisession, gc = TRUE, workers = no_cores) # TRY USING 'multicore' on Derecho

# Run MIMICS!

print(paste0("Starting ", MIM_runs, " runs"))
print(paste0("Start time: ", Sys.time()))

start_time <- Sys.time()
MC_MIMICS <- rand_params %>% split(1:nrow(rand_params)) %>% future_map(~MIMrepeat(forcing_df = forcing_data ,rparams = .), .progress=TRUE) %>% bind_rows()

#check_df <- MC_MIMICS %>% filter(DAY == 200) %>% filter(SITE == "BART") %>% filter(run_num == 1)

wall_time <- Sys.time() - start_time
print(paste0("Wall time: ", as.character(wall_time)))


# Release CPU cores
plan(sequential)
nbrOfWorkers()

# Clean up memory
gc()

## Join random parameters to MIMICS output table
MC_MIMICS <- MC_MIMICS %>% left_join(rand_params)

##########################################
# Save MC output data
##########################################
#saveRDS(MC_MIMICS, paste0("MC_output/MSBio_MC_Tsens_", as.character(MIM_runs), "_", format(Sys.time(), "%Y%m%d_%H%M%S"),  ".rds"))

