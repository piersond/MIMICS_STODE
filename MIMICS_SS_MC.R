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

# Load function to calc MIMICS steady state pools
source("functions/MIMICS_calc_steady_state_pools.R")

# Set MIMICS parameters
source("parameters/MIMICS_parameters_sandbox_20231129.R")

#------------------------------------------------------------

# Example data
forcing_data <- read.csv("example_simulations/data/LTER_SITE_1.csv", as.is=T)
MIM_runs = 100

# Function template for randomizing parameters AND
#  returning complete list output from MIMICS_SS()
randomize_MIMICS_SS <- function(df){
  
  # Randomize parameters here
  Vslope_x <- runif(1, 0.7, 1.3)
  Vslope <<- rep(0.063, 6) * Vslope_x
  
  output <- df %>% 
    rowwise() %>%
    mutate(MIMout = list(append(MIMICS_SS(cur_data()), 
                                c(Vslope_x=Vslope_x)))) # Store parameter multipliers here
  return(output)
}

# Parallel repeat run of dataset through randomized MIMICS steady-state
no_cores <- availableCores() - 1
plan(multisession, gc = TRUE, workers = no_cores) # Switch to 'multicore' on Derecho

print(paste0("Starting ", MIM_runs, " runs"))
print(paste0("Start time: ", Sys.time()))
start_time <- Sys.time()

MIM_repeat <- future_map(1:MIM_runs, ~ {
  forcing_data %>% randomize_MIMICS_SS()}) %>% 
  bind_rows() %>%   
  select(MIMout)

wall_time <- Sys.time() - start_time
print(paste0("Wall time: ", as.character(wall_time)))

# Release cores
plan(sequential)

# Turn tibble into large list of lists
MIM_repeat <- MIM_repeat[[1]]

saveRDS(MIM_repeat, "MC_SS_output.rds")

