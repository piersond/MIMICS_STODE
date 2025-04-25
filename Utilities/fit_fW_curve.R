library(dplyr)
library(ggplot2)

setwd("C:/github/MIMICS_STODE")

# Pull in raw data from incubation experiment
lab_data <- read.csv("Example_simulations/Data/MSBio_moisture_control_data.csv")

plot_fW_curve <- function(fW_p1, fW_p2, o_threshold, o_spread, d_threshold, d_strength){

  # Build fW curve from 0-100 grav moisture
  df <- data.frame(theta = 0, fW = 0, fW_mod = 0, fW_mod2 = 0)
  for(i in 1:100){
    theta_liq  <- i/100
    theta_frzn = 0
    air_filled_porosity = max(0.0, 1.0-theta_liq-theta_frzn)
  
    f <- function(x, p1, p2) {x^p1 * (1-x)^p2}
    fW_p3 <- optimize(f, interval=c(0.01,1), p1=fW_p1, p2=fW_p2, maximum = T)$objective
  
    #fW = (theta_liq^3 * air_filled_porosity^2.5)/0.022600567942709
    fW = (theta_liq^fW_p1 * air_filled_porosity^fW_p2)/fW_p3
 
    # optimum breadth modifier
    fW_mod <- ifelse(fW > o_threshold, o_threshold + (fW-o_threshold)/o_spread, fW)
    
    # dry threshold adjustment
    fW_mod <- ifelse(theta_liq < d_threshold, fW_mod - (d_threshold-theta_liq) * d_strength, fW_mod)

    # zero adjustment
    fW_mod <- ifelse(fW_mod < 0, 0, fW_mod)
        
    df[i+1,1] <- i
    df[i+1,2] <- fW
    df[i+1,3] <- fW_mod
  }
  
  # Plot fW curve vs. lab data
  ggplot(lab_data %>% group_by(gravm) %>% summarise(cumlCO2 = mean(cumlCO2)),
                                                    aes(y=cumlCO2, x=gravm)) + geom_point(aes(group=gravm), size=3) +
    theme_bw() +
    geom_line(data=df, aes(x=theta, y=fW*max(lab_data$cumlCO2))) +
    geom_line(data=df, aes(x=theta, y=fW_mod*max(lab_data$cumlCO2)), linetype="dashed") +
    #geom_line(aes(group=site, color=site)) +
    ylim(0, max(lab_data$cumlCO2))
}

#fW_p1, fW_p2, o_threshold, o_spread, d_threshold, d_strength
plot_fW_curve(0.8, 2, 0.85, 1.5, 0.15, 4)


# create a function to fit the fW curve to the lab data

try_fW_curve <- function(n){

  #!!! HARD CODED DATAFRAME !!!
  target_data <- lab_data %>% 
                 filter(gravm < 12) %>%
                 group_by(gravm) %>% 
                 summarise(mean_cumlCO2 = mean(cumlCO2))

  # Set random curve parameters
  fW_p1 = runif(1, 0.5, 3)
  fW_p2 = runif(1, 1, 5)
  f <- function(x, p1, p2) {(x^p1 * (1-x)^p2)}
  fW_p3 <- optimize(f, interval=c(0.01,1), p1=fW_p1, p2=fW_p2, maximum = T)$objective
  
  # Get fW curve points at moisture intervals
  df <- NULL
  for(i in 1:nrow(target_data)){
    theta_liq  <- unique(target_data$gravm)[i]/100
    theta_frzn = 0
    air_filled_porosity = max(0.0, 1.0-theta_liq-theta_frzn)
    fW = (theta_liq^fW_p1 * air_filled_porosity^fW_p2)/fW_p3 
    
    # append to list
    df$fW[i] = fW
  }
  
  threshold = runif(1, 0.8, 0.94)
  spread = runif(1, 1.2, 2)
  df$fW_mod <- ifelse(df$fW > threshold, threshold + (df$fW-threshold)/spread , df$fW)
  
  trial_data = target_data %>% mutate(fW = df$fW * max(target_data$mean_cumlCO2)) %>%
               mutate(fW_mod = df$fW_mod * max(target_data$mean_cumlCO2))
  
  # get root mean square error
  rmse = sqrt(mean((trial_data$mean_cumlCO2 - trial_data$fW)^2))
  rmse_mod = sqrt(mean((trial_data$mean_cumlCO2 - trial_data$fW_mod)^2))
  
  return(data.frame(n = n,
                    p1 = fW_p1,
                    p2 = fW_p2,
                    p3 = fW_p3,
                    threshold = threshold,
                    spread = spread,
                    RMSE = rmse,
                    RMSE_mod = rmse_mod))  
}

#######################################################
### Parallelize fitting the fW curve to the lab data ###

# In parallel, try fitting the fW curve to the lab data 1000 times
library(purrr)
library(furrr)
furrr_options(seed=TRUE)

# Build a list with 1000 elements
run_n = 4000
run_list <- lapply(seq(1, run_n), function(x) x)

print(paste0("Start time: ", Sys.time()))
start_time <- Sys.time()

fW_trials <- future_map(run_list, try_fW_curve, .progress = TRUE) %>% bind_rows() 

wall_time <- Sys.time() - start_time
print(paste0("Wall time: ", as.character(wall_time)))

# Release CPU cores
plan(sequential)
nbrOfWorkers()

# Clean up memory
gc()


#######################################################
### PLOT BEST FIT ###

# Find the best fit
best_fit <- fW_trials %>% filter(RMSE == min(RMSE))

# Set best fit curve parameters
fW_p1 = best_fit$p1
fW_p2 = best_fit$p2
f <- function(x, p1, p2) {x^p1 * (1-x)^p2}
fW_p3 <- optimize(f, interval=c(0.01,1), p1=fW_p1, p2=fW_p2, maximum = T)$objective

# Get fW curve points at moisture intervals
fW_curve = NULL
for(i in 1:100){
  theta_liq  <- i/100
  theta_frzn = 0
  air_filled_porosity = max(0.0, 1.0-theta_liq-theta_frzn)
  fW = (theta_liq^fW_p1 * air_filled_porosity^fW_p2)/fW_p3

  fW_curve <- rbind(fW_curve, data.frame(theta = theta_liq, fW = fW))
}

# Plot best fit fW curve vs. lab data
ggplot(lab_data, aes(y=cumlCO2, x=gravm)) + geom_point(aes(group=moisture.trt, color=site), size=3) +
  theme_bw() +
  geom_line(data=fW_curve, aes(x=theta*100, y=fW*max(lab_data$cumlCO2))) +
  #geom_line(aes(group=site, color=site)) +
  ggtitle(paste0("Best fit fW curve: \np1 = ", fW_p1, "\np2 = ", fW_p2, "\np3 = ", fW_p3, "\nRMSE = ", best_fit$RMSE))

# Save best fit parameters
#write.csv(best_fit, "best_fit_fW_curve_params.csv", row.names=FALSE)

# Save best fit curve
#write.csv(df, "best_fit_fW_curve.csv", row.names=FALSE)

# Save best fit curve plot
#ggsave("best_fit_fW_curve.png")


### Plot best fit for RMSE_mod
best_fit_mod <- fW_trials %>% filter(RMSE_mod == min(RMSE_mod))

# Set best fit curve parameters
fW_p1 = best_fit_mod$p1
fW_p2 = best_fit_mod$p2
f <- function(x, p1, p2) {x^p1 * (1-x)^p2}
fW_p3 <- optimize(f, interval=c(0.01,1), p1=fW_p1, p2=fW_p2, maximum = T)$objective

# Get fW curve points at moisture intervals
fW_mod_curve = NULL
for(i in 1:100){
  theta_liq  <- i/100
  theta_frzn = 0
  air_filled_porosity = max(0.0, 1.0-theta_liq-theta_frzn)
  fW = (theta_liq^fW_p1 * air_filled_porosity^fW_p2)/fW_p3

  fW_mod_curve <- rbind(fW_mod_curve, data.frame(theta = theta_liq, fW_mod = fW))
}

# Apply mod
threshold = best_fit_mod$threshold
spread = best_fit_mod$spread
fW_curve$fW_mod <- ifelse(fW_curve$fW > threshold, threshold + (fW_curve$fW-threshold)/spread , fW_curve$fW)

# Plot best fit fW curve vs. lab data
ggplot(lab_data, aes(y=cumlCO2, x=gravm)) + geom_point(aes(group=moisture.trt, color=site), size=3) +
  theme_bw() +
  geom_line(data=fW_curve, aes(x=theta*100, y=fW*max(lab_data$cumlCO2))) +
  geom_line(data=fW_mod_curve, aes(x=theta*100, y=fW_mod*max(lab_data$cumlCO2)), linetype="dashed") +
  #geom_line(aes(group=site, color=site)) +
  ggtitle(paste0("Best fit fW curve: \np1 = ", fW_p1, "\np2 = ", fW_p2, "\np3 = ", fW_p3, "\nRMSE = ", best_fit_mod$RMSE_mod))


