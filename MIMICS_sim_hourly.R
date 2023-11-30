## Set working drive
setwd("C:/github/MIMICS_STODE")

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

# Bring in MIMICS steady state ftn
source("MIMICS_steady_state_pools.R")

# Set MIMICS parameters via R script
source("Parameters/MIMICS_parameters_sandbox_20231129.R")

#######################################################
# Ftn to run MIMICS from steady state at hourly time-step
#######################################################
MIMICS_hourly <- function(MIMss, nday){
  
  #DEBUG
  #MIMss = MIMout_single  
  #nday = 365 * 30
  
  #Grab forcings
  Site   <- as.character(MIMss[[3]][1])
  # ANPP   <- MIMss[[3]][2]
  # tsoi   <- MIMss[[3]][3]
  # clay   <- MIMss[[3]][4]/100
  # lig_N  <- MIMss[[3]][5]
  
  # Set run time vars
  day    <- seq(1,nday,1)
  year   <- day/365
  doy    <- 1
  
  #Init pool arrays
  LIT    <- array(NA, dim = c(2,nday))
  MIC    <- array(NA, dim = c(2,nday))
  SOM    <- array(NA, dim = c(3,nday))
  CO2    <- array(NA, dim = c(2,nday))
  
  #Init pools
  LIT_1    <- (MIMss[[1]][1] / depth) / (1e4 / 1e6)    
  LIT_2    <- (MIMss[[1]][2] / depth) / (1e4 / 1e6)    
  MIC_1    <- (MIMss[[1]][3] / depth) / (1e4 / 1e6)    
  MIC_2    <- (MIMss[[1]][4] / depth) / (1e4 / 1e6)    
  SOM_1    <- (MIMss[[1]][5] / depth) / (1e4 / 1e6)    
  SOM_2    <- (MIMss[[1]][6] / depth) / (1e4 / 1e6)    
  SOM_3    <- (MIMss[[1]][7] / depth) / (1e4 / 1e6)    
  
  # Init rates
  I        <- MIMss[[2]][1:2]
  VMAX     <- MIMss[[2]][3:8]
  KM       <- MIMss[[2]][9:14]
  CUE      <- MIMss[[2]][15:18]
  fPHYS    <-  MIMss[[2]][19:20]
  fCHEM    <-  MIMss[[2]][21:22]
  fAVAI    <-  MIMss[[2]][23:24]
  fI       <-  MIMss[[2]][25:26]
  tau      <- MIMss[[2]][27:28]
  desorb   <- MIMss[[2]][37]
  #DEsorb   <- MIMss[[2]][38]
  #OXIDAT   <- MIMss[[2]][39]
  KO       <- MIMss[[2]][40:41]
  
  sim_year = 0
  
  # Run MIMICS hourly for ndays 
  for (d in 1:nday)  { 
    for (h in 1:24)   {
      
      #Fluxes at each time step
      LITmin  <- rep(NA, dim=4)
      MICtrn  <- rep(NA, dim=6)
      SOMmin  <- rep(NA, dim=2)
      DEsorb  <- rep(NA, dim=1)
      OXIDAT  <- rep(NA, dim=1)
      
      #----- Reverse MM version ----------
      
      #Flows to and from MIC_1
      LITmin[1] = MIC_1 * VMAX[1] * LIT_1 / (KM[1] + MIC_1)   #MIC_1 decomp of MET lit
      LITmin[2] = MIC_1 * VMAX[2] * LIT_2 / (KM[2] + MIC_1)   #MIC_1 decomp of STRUC lit
      SOMmin[1] = MIC_1 * VMAX[3] * SOM_3 / (KM[3] + MIC_1)   #Decomp of SOMa by MIC_1
      
      MICtrn[1] = MIC_1 * tau[1]  * fPHYS[1]                  #MIC_1 turnover to SOMp
      MICtrn[2] = MIC_1 * tau[1]  * fCHEM[1]                  #MIC_1 turnover to SOMc  
      MICtrn[3] = MIC_1 * tau[1]  * fAVAI[1]                  #MIC_1 turnover to SOMa 
      
      #Flows to and from MIC_2
      LITmin[3] = MIC_2 * VMAX[4] * LIT_1 / (KM[4] + MIC_2)   #decomp of MET litter
      LITmin[4] = MIC_2 * VMAX[5] * LIT_2 / (KM[5] + MIC_2)   #decomp of SRUCTURAL litter
      SOMmin[2] = MIC_2 * VMAX[6] * SOM_3 / (KM[6] + MIC_2)   #decomp of PHYSICAL SOM by MIC_1
      
      MICtrn[4] = MIC_2 * tau[2]  * fPHYS[2]                  #MIC_2 turnover to SOMp 
      MICtrn[5] = MIC_2 * tau[2]  * fCHEM[2]                  #MIC_2 turnover to SOMc 
      MICtrn[6] = MIC_2 * tau[2]  * fAVAI[2]                  #MIC_2 turnover to SOMa  
      
      
      DEsorb    = SOM_1 * desorb  #* (MIC_1 + MIC_2)  #desorbtion of PHYS to AVAIL (function of fCLAY)
      OXIDAT    = ((MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + MIC_2)) +
                     (MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + MIC_1)))  #oxidation of C to A
      
      
      LIT_1 = LIT_1 + I[1]*(1-FI[1]) - LITmin[1] - LITmin[3] 
      MIC_1 = MIC_1 + CUE[1]*(LITmin[1]+ SOMmin[1]) + CUE[2]*(LITmin[2]) - sum(MICtrn[1:3])
      SOM_1 = SOM_1 + I[1]*FI[1] + MICtrn[1] + MICtrn[4]- DEsorb 
      CO2_1  = (1-CUE[1])*(LITmin[1]+ SOMmin[1]) + (1-CUE[2])*(LITmin[2])
      
      LIT_2 = LIT_2 + I[2] * (1-FI[2]) - LITmin[2] - LITmin[4]
      MIC_2 = MIC_2 + CUE[3]*(LITmin[3]+ SOMmin[2]) + CUE[4]*(LITmin[4]) - sum(MICtrn[4:6])  
      SOM_2 = SOM_2 + I[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT
      CO2_2  = (1-CUE[3])*(LITmin[3]+ SOMmin[2]) + (1-CUE[4])*(LITmin[4])
      
      SOM_3 = SOM_3 + MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2]	
      
      #write out daily results
      if (h == 24) {
        LIT[1,d] <- LIT_1
        LIT[2,d] <- LIT_2
        MIC[1,d] <- MIC_1
        MIC[2,d] <- MIC_2
        SOM[1,d] <- SOM_1
        SOM[2,d] <- SOM_2
        SOM[3,d] <- SOM_3
        CO2[1,d] <- CO2_1
        CO2[2,d] <- CO2_2
        
        #advancy day of year counter
        if (doy == 365) {
          doy <- 1
          sim_year = sim_year + 1
          print(paste0("Finished MIMICS simulation year ", sim_year))
        } else {
          doy <- doy + 1
        } 
      }	
    }		#close hour loop
  }		#close daily loop
  
  #return(list(LIT, MIC, SOM))
  SS_MIMout <- rbind(as.data.frame(LIT), 
                     as.data.frame(MIC), 
                     as.data.frame(SOM),
                     as.data.frame(CO2))
  
  SS_MIMout <- SS_MIMout * depth *1e4 / 1e6 
  SS_MIMout <- as.data.frame(t(SS_MIMout))
  colnames(SS_MIMout) <- c("LITm", "LITs", "MICr", "MICk", "SOMp", "SOMc", "SOMa", "CO2r", "CO2k")
  SS_MIMout$MIMSOC <- rowSums(SS_MIMout[,1:7])
  SS_MIMout <- cbind(data.frame(SITE = as.character(Site),
                                DAY=seq(1:nrow(SS_MIMout))),
                     SS_MIMout)
  
  return(SS_MIMout)
}



#------------------------------------------------------------------------------
#--- EXAMPLE USE ---
#------------------------------------------------------------------------------

#-------------------------------------------------------
# Load LTER forcing dataset for example simulations
#-------------------------------------------------------

# Forcing data for LTER sites
LTER <- read.csv("Data/LTER_SITE_1.csv", as.is=T)

#---------------------------------------------------------------------------
# Example: Run MIMICS simulation forward from steady-state for 30 years
#---------------------------------------------------------------------------

# Get MIMICS steady state (SS) pools for a single site in the LTER forcing dataset
MIMout_single <- MIMICS_SS(LTER[6,])
MIMout_tbl <- MIMICS_SS_format(MIMout_single)

# R simulation forward for specified # of days
sim_days <- 30*365
SS_MIMrun <- MIMICS_hourly(MIMout_single, sim_days)

# Plot MIMICS total SOM to verify steady state (expect flat blue line)
ggplot(SS_MIMrun) +
  geom_line(aes(y=MIMSOC, x=DAY), color="blue", linewidth=1) +
  ylim(min(SS_MIMrun$MIMSOC)-1,max(SS_MIMrun$MIMSOC)+1) +
  ggtitle("Check MIMICS steady-state C pools over time") +
  theme_bw()
