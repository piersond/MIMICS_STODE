#rm(list = ls())

## Set working drive
#setwd("/Users/wwieder/Will/git_repos_local/MIMICS_STODE")

#Libraries
library(rootSolve)
library(boot)
library(dplyr)
library(purrr)
library(ggplot2)
library(Metrics)
library(DT)

# Bring in RXEQ function
source("RXEQ/RXEQ_ftn.R")

# calculate Tpars via R script
source("calc_Tpars.R")

# Set MIMICS parameters via R script
source("Parameters/MIMICS_parameters_sandbox_20231129.R")

###############################################
# MIMICS single point function
#> Input dataframe ("df") must contain columns: SITE, ANPP, fCLAY, TSOI, LIG_N or LIG + CN // Optional: MAT, GWC  
#>> With output required for running forward
###############################################
MIMICS_SS <- function(df){
  
  #DEBUG
  #df = LTER[1,]
  
  # Convert all column names to upper case
  colnames(df) <- toupper(colnames(df))
  
  ### Setup a var to collect run notes
  note <- ""
  
  ### Bring in forcing ANPP value, convert gDW to gC
  ANPP <- df$ANPP/2
  
  ### Bring in CLAY value, convert from percent to decimal
  fCLAY <- df$CLAY/100
  
  ### Bring in TSOI value
  TSOI <- df$TSOI
  
  ### Bring in lig:N forcing data
  LIG_N <- df$LIG_N
  LIG <- df$LIG
  CN <- df$CN

  # Bring in soil moisture information
  theta_liq  <- df$GWC/100  # GWC = Gravimetric water content
  theta_frzn <- 0           # Not used here. TODO Needs validation. 

  ### Bring in mean annual temperature data
  MAT <- df$MAT  
  
  # Use TSOI if 'MAT' column is not in forcing data
  if(is.null(MAT)){
    MAT <- TSOI
  }

  ############################################################
  # MIMICS MODEL CODE STARTS HERE
  ############################################################
  # function calculates fMET with LIG_N if provided in input data.
  Tpars <- calc_Tpars_Conly(ANPP=ANPP, fCLAY=fCLAY, TSOI=TSOI, MAT=MAT,     
                            CN=CN, LIG=LIG, LIG_N=LIG_N,
                            theta_liq=theta_liq, theta_frzn=theta_frzn)  #<---- IMPORTANT USER OPTIONS HERE
                            #fWmethod=df$fWmethod, historic=df$historic, #<--- Moved to parameter file
                            #fixed_fMET=df$fixed_fMET, tauMethod=df$tauMethod)
    
  # Create arrays to hold output
  lit     <- Tpars$I
  mic     <- Tpars$I
  som     <- rep(NA, 3) 
  som[1]  <- Tpars$I[1]
  som[2]  <- Tpars$I[2]
  som[3]  <- Tpars$I[1] 
  CO2     <- rep(0, 2) 

  Ty    <- c(LIT_1 = lit[1], LIT_2 = lit[2], 
              MIC_1 = mic[1], MIC_2 = mic[2], 
              SOM_1 = som[1], SOM_2 = som[2], SOM_3 = som[3])
  
  ## Set global parameters to ensure passing of variables to stode function
  .GlobalEnv$VMAX <- Tpars$VMAX
  .GlobalEnv$KM <- Tpars$KM
  .GlobalEnv$fPHYS <- Tpars$fPHYS
  .GlobalEnv$fCHEM <- Tpars$fCHEM
  .GlobalEnv$fAVAI <- Tpars$fAVAI
  .GlobalEnv$I <- Tpars$I
  .GlobalEnv$tau <- Tpars$tau
  .GlobalEnv$LITmin <- Tpars$LITmin
  .GlobalEnv$SOMmin <- Tpars$SOMmin
  .GlobalEnv$MICtrn <- Tpars$MICtrn
  .GlobalEnv$desorb <- Tpars$desorb
  .GlobalEnv$DEsorb <- Tpars$DEsorb
  .GlobalEnv$OXIDAT <- Tpars$OXIDAT
  .GlobalEnv$beta <- Tpars$beta
  
  # ------------RUN THE MODEL-------------
  test  <- stode(y = Ty, time = 1e7, fun = RXEQ, parms = Tpars, positive = TRUE)
  
  ### Calc and get MIMICS output 
  MIMLIT    <- (test[[1]][[1]]+test[[1]][[2]]) * depth *1e4 / 1e6 #convert kgC/m2 from mgC/cm3 (0-30 cm) 
  MIMMIC    <- (test[[1]][[3]]+test[[1]][[4]]) * depth *1e4 / 1e6
  MIM_CO    <-  test[[1]][[3]]/test[[1]][[4]]
  MIMSOC    <- sum(test[[1]]) * depth * 1e4 / 1e6   
  
  table <- as.numeric(test[[1]])
  
  MIMout <- list()
  MIMout[[1]] <- as.numeric(test[[1]]) * depth *1e4 / 1e6   # MIMICS steady state C pools
  MIMout[[2]] <- Tpars                                      # MIMICS parameters used for/from simulation
  MIMout[[3]] <- df                                         # Forcing dataframe
  MIMout[[4]] <- as.numeric(test[[2]]) * depth *1e4 / 1e6   # CO2-r & CO2-k pools 
  
  return(MIMout)
}


#################################################
# Helper ftn to summarize MIMICS_SS() output
#################################################

MIMICS_SS_format <- function(MIMICS_SS_output) {
  
  MIMout_C_pools <- as.data.frame(t(as.data.frame(MIMICS_SS_output[[1]]))) 
  colnames(MIMout_C_pools) <- c("LITm", "LITs", "MICr", "MICk", "SOMp", "SOMc", "SOMa")
 
  MIMout_C_pools$MIMSOC <- rowSums(MIMout_C_pools[,1:7])
  MIMout_C_pools$MIMLIT <- MIMout_C_pools$LITm + MIMout_C_pools$LITs    
  MIMout_C_pools$MIMMIC <- MIMout_C_pools$MICr + MIMout_C_pools$MICk
  MIMout_C_pools$MIC_ratio <- MIMout_C_pools$MICr/MIMout_C_pools$MICk
  
  MIMout_CO2_pools <- as.data.frame(t(as.data.frame(MIMICS_SS_output[[4]])))
  colnames(MIMout_CO2_pools) <- c("CO2r", "CO2k")
  
  SOMp_desorb_yr <- as.numeric(MIMICS_SS_output[[2]]['desorb']) * 24 * 365
  MIMout_C_pools$SOMp_Turnover_yrs <- MIMout_C_pools$SOMp/SOMp_desorb_yr
  
  MIMout_single_tbl <- cbind(as.data.frame(MIMICS_SS_output[[3]]), MIMout_C_pools %>%
                               select(MIMSOC, MIMLIT, MIMMIC, MIC_ratio, LITm, LITs, MICr, 
                                      MICk, SOMa, SOMc, SOMp, SOMp_Turnover_yrs),
                             MIMout_CO2_pools)
  return(MIMout_single_tbl)
}


#------------------------------------------------------------------------------
#--- EXAMPLE USE ---
#------------------------------------------------------------------------------

#>> Set working drive to the script directory (repo main)

#-------------------------------------------------------
# Load LTER forcing dataset for example simulations
#-------------------------------------------------------

# Forcing data for LTER sites
LTER <- read.csv("Data/LTER_SITE_1.csv", as.is=T)


#---------------------------------------------------------------------------
# Example: Find steady state pools for a single site
#---------------------------------------------------------------------------

# Get MIMICS steady state (SS) pools for a single site in the LTER forcing dataset
MIMout_single_raw <- MIMICS_SS(LTER[1,])
MIMICS_ss_tbl <- MIMICS_SS_format(MIMout_single_raw)

# View table
#datatable(MIMICS_ss_tbl %>% mutate_if(is.numeric, round, digits=3))

#---------------------------------------------------------------------------
# Example: Find steady state pools for a forcing dataset
#---------------------------------------------------------------------------

MIMruns <- LTER %>% split(1:nrow(LTER)) %>% map(~ MIMICS_SS(df=.))
MIMruns_format <- lapply(MIMruns, MIMICS_SS_format) %>% bind_rows()
MIMICS_ss_dataset <- LTER[,1:2] %>% cbind(MIMruns_format %>% select(-SITE, -SOC))  # Bind data info columns to MIMICS output

# View table
datatable(MIMICS_ss_dataset %>% mutate_if(is.numeric, round, digits=3))

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
  geom_text(aes(label=paste0(Site)),hjust=-0.2, vjust=0.2) +
  annotate("text", label = lb2, x = 2, y = 8.5, size = 4, colour = "black", parse=T) +
  annotate("text", label = paste0("RMSE = ", rmse), x = 2, y = 7.4, size = 4, colour = "black") +
  ylim(0,10) + xlim(0,10) +
  theme_minimal()
 
# Check against sandbox
# # Run the sandbox MIMICS script, then run...
# # plot_data$Wieder_2015_MIMSOC <- MIMSOC
# # ggplot(plot_data, aes(x=MIMSOC, y=Wieder_2015_MIMSOC)) + geom_point(size=3) +
# #   geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
# #   xlab("Pierson_STODE\nMIMSOC") +
# #   ylab("Wieder et al. 2015 - Sandbox\nMIMSOC") +
# #   theme_minimal()

# quick way to set title for plot below 
if (.GlobalEnv$beta==1) {
  gg_title = 'tauMethod = NPP' 
} else {
  gg_title = 'tauMethod = beta' 
}

ggplot(plot_data, aes(x=ANPP, y=MIMMIC, color=TSOI)) +
  geom_point(size=4, alpha=0.8) +
  ggtitle(gg_title) +
  theme_minimal()
