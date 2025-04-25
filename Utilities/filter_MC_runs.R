library(dplyr)
library(Metrics)

getwd()

# Load the MC output .rds file
mc_rds_path <- "C:/github/MIMICS_MSBio/MSBio_redux/Derecho_MC/MSBio_MC_Tsens_12800_20240408_125630.rds"
MC <- readRDS(mc_rds_path)

colnames(MC)

# Calculate cost
    # e.g., for MSBio, calculate difference between the CO2C_prop and CO2_prop_totC
MC$cost <- abs(MC$CO2C_PROP - MC$CO2_prop_totC) 

# Summarize by each MC run and calculate mean cost and RMSE
smry <- MC %>% #filter(moisture.trt != 20) %>%  #<---!!! NOTE THE FILTER HERE
            group_by(run_num) %>% summarize(n = n(),
                                               cost = mean(cost),
                                                RMSE = rmse(CO2_prop_totC, CO2C_PROP))

#!!!
# ADTL STEP NEEDED FOR MIMICS_SS
#--> Add filter for plausible ranges for pools etc.
#!!!

# Find the best n runs using either cost or RMSE
best_n = 1
#best_run <- as.numeric(smry %>% slice_min(order_by = RMSE, n = best_n) %>% select(run_num))
best_run <- as.numeric(smry %>% slice_min(order_by = cost, n = best_n) %>% select(run_num))

# Select best fit parameters
best_fit <- MC %>% filter(run_num == best_run)

# Isolate parameter multipliers for the best run
colnames(best_fit)
best_run_params <- best_fit[33:38] %>% distinct()
