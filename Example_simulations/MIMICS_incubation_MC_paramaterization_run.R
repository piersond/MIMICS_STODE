### MIMICS MC for litterbag simulations

########################################
# Load R packages
########################################
library(dplyr)
library(rootSolve)
library(purrr)
library(furrr)

setwd("C:/github/MIMICS_STODE")

########################################
# Load MIMICS data and ftns
########################################
source("Parameters/MIMICS_parameters_MSBio_Incubation_v2.R") #MSBio litter incubation parameters
source("functions/RXEQ.R")
source("functions/calc_Tpars.R")
source("functions/MIMICS_sim_incubation.R")
source("functions/MC_parameterization/MIMICS_incubation_repeat.R")
source("functions/MC_parameterization/set_parameter_defaults.R")


########################################
# Load forcing data
########################################
#load site data
data <- read.csv("example_simulations/Data/MSBio_MIM_forcings.csv")

####################################
# Use the brute force MIMICS ftn
####################################

# Set desired number of random parameter runs
MIM_runs <- 100

### Create random parameter dataframe
  #!!! NOTE - Must change variable in Incubation_MIMICS_repeat as well
rand_params <- data.frame(
  Vslope_x = runif(MIM_runs, 0.5, 2),
  Vint_x = runif(MIM_runs, 0.8, 1.3),
  Kslope_x = runif(MIM_runs, 0.5, 2),
  Kint_x = runif(MIM_runs, 0.5, 2)#,
  # vMOD_x = runif(MIM_runs, 0.5, 2),
  # kMOD_x = runif(MIM_runs, 0.5, 2),
  # Tau_x = runif(MIM_runs, 0.3, 3),
  # Tau_r = runif(MIM_runs, 0.3, 3),
  # Tau_k = runif(MIM_runs, 0.3, 3),
  # CUE_x = runif(MIM_runs, 0.5, 1.4), 
  # CUE_r = runif(MIM_runs, 0.5, 1.4),
  # CUE_k = runif(MIM_runs, 0.5, 1.4)
)

rand_params$run_num <- seq(1,MIM_runs,1)

# Set number of cores to use
no_cores <- availableCores() - 1
plan(multicore, gc = TRUE, workers = no_cores) #vs. multisession

# Run MIMICS!

print(paste0("Starting ", MIM_runs, " parameterization runs"))
print(paste0("Start time: ", Sys.time()))

start_time <- Sys.time()


# Take each row of the random parameters (a single run number) and run the MIMrepeat function for each row
MC_MIMICS <- rand_params %>% split(1:nrow(rand_params)) %>% future_map(~INC_MIMrepeat(forcing_df = data, rparams = .), .progress=TRUE) %>% 
  bind_rows() 


wall_time <- Sys.time() - start_time
print(paste0("Wall time: ", as.character(wall_time)))


# Release CPU cores
plan(sequential)
nbrOfWorkers()

# Clean up memory
gc()


#####################################
# Join MC output with forcing data
MC_output_tbl <- left_join(data, MC_MIMICS %>% select(-SITE), by="ID") %>% left_join(., rand_params, by="run_num")

##########################################
# Save MC output data - for use with computing clusters
##########################################
saveRDS(MC_output_tbl, paste0("MIM_INC-MC_", as.character(MIM_runs), "_", format(Sys.time(), "%Y%m%d_%H%M%S_"),  ".rds"))


#####################################
# Filter and plot best fit for MSBio data
# smry_df <- MC_output_tbl %>% mutate(diff = abs(CO2C_prop - CO2_prop_totC)) %>%
#                              group_by(run_num) %>% summarize(cost = sum(diff),
#                                                              Vslope_x = unique(Vslope_x)) 
# library(ggplot2)
# ggplot(smry_df, aes(x=Vslope_x, y=cost, color=cost)) + 
#   geom_point(size=3, alpha=0.7) +
#   scale_color_gradientn(colors = c("green", "red")) +
#   theme_minimal()
                              

#--- NOTES ---------------------
# 10 runs = 12.79 sec on 19 high-end cpu cores
# 100 runs = 2.09 minutes on 19 high-end cpu cores
