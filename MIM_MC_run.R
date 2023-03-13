### MIMICS MC

setwd("C:/github/MIMICS_STODE")

########################################
# Load MIMICS data and ftns from Brute Forcing script
########################################
source("MIM_ftns/MIMICS_repeat_base.R")


########################################
# Load forcing data
########################################
data <- read.csv("Data/LTER_SITE_1.csv", as.is=T)


####################################
# Use the brute force MIMICS ftn
####################################

# Set desired number of random parameter runs
MIM_runs <- 1000

### Create random parameter dataframe
## Parameter range informed by range observed over 10+ MCMC analysis results
rand_params <- data.frame(Vslope_x = runif(MIM_runs, 0.4, 4),  
                          Vint_x = runif(MIM_runs, 0.3, 3),  
                          Kslope_x = runif(MIM_runs, 0.4, 4),  
                          Kint_x = runif(MIM_runs, 0.3, 3),  
                          Tau_x = runif(MIM_runs, 0.3, 3),  
                          CUE_x = runif(MIM_runs, 0.5, 1.5),  
                          desorb_x = runif(MIM_runs, 0.001, 0.3),  
                          fPHYS_x = runif(MIM_runs, 0.4, 4)  
                          )

rand_params$run_num <- seq(1,MIM_runs,1)

#### SET UP PARALLEL PROCESSING

# Set number of cores to use
no_cores <- availableCores() - 1
plan(multicore, gc = FALSE, workers = no_cores)

# Run MIMICS!

print(paste0("Starting ", MIM_runs, " runs"))
print(paste0("Start time: ", Sys.time()))

start_time <- Sys.time()
MC_MIMICS <- rand_params %>% split(1:nrow(rand_params)) %>% future_map(~ MIMrepeat(forcing_df = data, rparams = ., output_type = "summary"), .progress=TRUE) %>% bind_rows()
print(paste0("Task time: ", Sys.time() - start_time))

# Release CPU cores
plan(sequential)
nbrOfWorkers()

# Clean up memory
gc()


## Join parameters to MIMICS output table
MC_MIMICS <- MC_MIMICS %>% left_join(rand_params)


##########################################
# Save MC output data
##########################################
#saveRDS(MC_MIMICS, paste0("MC/Output/", "MC_MIMICS_data-r", as.character(MIM_runs), "_", format(Sys.time(), "%Y%m%d_%H%M%S_"),  ".rds"))




