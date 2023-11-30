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

###############################
# Litter bag simulation ftn
###############################

MIMICS_LITBAG <- function(forcing_df, litBAG, nspin_yrs=10, nspin_days=200, litadd_day=143, verbose=T){
  
  #DEBUG
  #forcing_df <- LTER[6,]
  #litBAG <- BAGS[1,]
  #nspin_yrs <- 10
  #nspin_days <- 200
  #litadd_day <- 143
  
  if(verbose){
    print("-------------------------------------------------------")
    print(paste0("Starting ", forcing_df$Site, " - ", litBAG$TYPE))
    print("-------------------------------------------------------")
  }
  
  # Get MIMICS steady_state output
  MIMss <- MIMICS_SS(forcing_df)
  
  nday   <- 365 * nspin_yrs + nspin_days
  day    <- seq(1,nday,1)
  year   <- (day-litadd_day)/365
  doy    <- 1
  
  #Init arrays to store daily output data
  LIT    <- array(NA, dim = c(2,nday))
  LITBAG <- array(NA, dim = c(2,nday))
  MIC    <- array(NA, dim = c(2,nday))
  SOM    <- array(NA, dim = c(3,nday))
  
  #Init MIMICS pools
  LIT_1    <- (MIMss[[1]][1] / depth) / (1e4 / 1e6)  
  LIT_2    <- (MIMss[[1]][2] / depth) / (1e4 / 1e6)  
  LITbag_1 <- 0
  LITbag_2 <- 0
  MIC_1    <- (MIMss[[1]][3] / depth) / (1e4 / 1e6)  
  MIC_2    <- (MIMss[[1]][4] / depth) / (1e4 / 1e6)  
  MICbag_1 <- (MIMss[[1]][3] / depth) / (1e4 / 1e6)  
  MICbag_2 <- (MIMss[[1]][4] / depth) / (1e4 / 1e6)  
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
  DEsorb   <- MIMss[[2]][38]
  OXIDAT   <- MIMss[[2]][39]
  KO       <- MIMss[[2]][40:41]
  
  sim_year = 0
  
  for (d in 1:nday)  { 
    for (h in 1:24)   {
      #Fluxes at each time step
      LITmin  <- rep(NA, dim=4)
      LITbag  <- rep(NA, dim=4)
      MICtrn  <- rep(NA, dim=6)
      SOMmin  <- rep(NA, dim=2)
      DEsorb  <- rep(NA, dim=1)
      OXIDAT  <- rep(NA, dim=1)
      
      # --- Reverse MM version ---
      
      #Flows to and from MIC_1
      LITmin[1] = MIC_1 * VMAX[1] * LIT_1 / (KM[1] + MIC_1)   #MIC_1 decomp of MET lit
      LITmin[2] = MIC_1 * VMAX[2] * LIT_2 / (KM[2] + MIC_1)   #MIC_1 decomp of STRUC lit
      LITbag[1] = MIC_1 * VMAX[1] * LITbag_1 / (KM[1] + MICbag_1)   #MIC_1 mineralization of METABOLIC litter
      LITbag[2] = MIC_1 * VMAX[2] * LITbag_2 / (KM[2] + MICbag_2)   #MIC_1 mineralization of STRUC litter
      SOMmin[1] = MIC_1 * VMAX[3] * SOM_3 / (KM[3] + MIC_1)   #Decomp of SOMa by MIC_1
      
      MICtrn[1] = MIC_1 * tau[1]  * fPHYS[1]                  #MIC_1 turnover to SOMp
      MICtrn[2] = MIC_1 * tau[1]  * fCHEM[1]                  #MIC_1 turnover to SOMc  
      MICtrn[3] = MIC_1 * tau[1]  * fAVAI[1]                  #MIC_1 turnover to SOMa 
      
      #Flows to and from MIC_2
      LITmin[3] = MIC_2 * VMAX[4] * LIT_1 / (KM[4] + MIC_2)   #decomp of MET litter
      LITmin[4] = MIC_2 * VMAX[5] * LIT_2 / (KM[5] + MIC_2)   #decomp of SRUCTURAL litter
      LITbag[3] = MIC_2 * VMAX[4] * LITbag_1 / (KM[4] + MICbag_1)   #mineralization of MET litter
      LITbag[4] = MIC_2 * VMAX[5] * LITbag_2 / (KM[5] + MICbag_2)   #mineralization of SRUCTURAL litter
      SOMmin[2] = MIC_2 * VMAX[6] * SOM_3 / (KM[6] + MIC_2)   #decomp of PHYSICAL SOM by MIC_1
      
      MICtrn[4] = MIC_2 * tau[2]  * fPHYS[2]                  #MIC_2 turnover to SOMp 
      MICtrn[5] = MIC_2 * tau[2]  * fCHEM[2]                  #MIC_2 turnover to SOMc 
      MICtrn[6] = MIC_2 * tau[2]  * fAVAI[2]                  #MIC_2 turnover to SOMa  
      
      
      DEsorb    = SOM_1 * desorb  #* (MIC_1 + MIC_2)  #desorbtion of PHYS to AVAIL (function of fCLAY)
      OXIDAT    = ((MIC_2 * VMAX[5] * SOM_2 / (KO[2]*KM[5] + MIC_2)) +
                     (MIC_1 * VMAX[2] * SOM_2 / (KO[1]*KM[2] + MIC_1)))  #oxidation of C to A
      
      LIT_1 = LIT_1 + I[1]*(1-FI[1]) - LITmin[1] - LITmin[3]
      LITbag_1 <- LITbag_1 - LITbag[1] - LITbag[3] #+ I[1]*(1-FI[1])
      MIC_1 = MIC_1 + CUE[1]*(LITmin[1]+ SOMmin[1]) + CUE[2]*(LITmin[2]) - sum(MICtrn[1:3])
      #MICbag_1
      SOM_1 = SOM_1 + I[1]*FI[1] + MICtrn[1] + MICtrn[4]- DEsorb 
      
      LIT_2 = LIT_2 + I[2] * (1-FI[2]) - LITmin[2] - LITmin[4]
      LITbag_2 <- LITbag_2 - LITbag[2] - LITbag[4] #+ I[2] * (1-FI[2])
      MIC_2 = MIC_2 + CUE[3]*(LITmin[3]+ SOMmin[2]) + CUE[4]*(LITmin[4]) - sum(MICtrn[4:6])  
      #MICbag_2
      SOM_2 = SOM_2 + I[2]*FI[2] + MICtrn[2] + MICtrn[5] - OXIDAT
      
      SOM_3 = SOM_3 + MICtrn[3] + MICtrn[6] + DEsorb + OXIDAT - SOMmin[1] - SOMmin[2]	#add litter bag on Oct 1
      
      if (d == litadd_day)   {
        if (h == 24)  {
          LITbag_1 <- LITbag_1 + as.numeric(litBAG[3])
          LITbag_2 <- LITbag_2 + as.numeric(litBAG[4])
        }
      }
      
      #write out daily results
      if (h == 24) {
        LIT[1,d] <- LIT_1
        LIT[2,d] <- LIT_2
        LITBAG[1,d]  <- LITbag_1
        LITBAG[2,d]  <- LITbag_2
        MIC[1,d] <- MIC_1
        MIC[2,d] <- MIC_2
        SOM[1,d] <- SOM_1
        SOM[2,d] <- SOM_2
        SOM[3,d] <- SOM_3
        
        
        #View(as.data.frame(LITBAG))
        
        
        #advance day of year counter
        if (doy == 365) {
          doy <- 1
          sim_year = sim_year + 1
          if(verbose){
            print(paste0("Finished MIMICS simulation year ", sim_year))
          }
        } else {
          doy <- doy + 1
        }                         
      }	   						
    }		#close hour loop
  }		#close daily loop
  #return(list(LIT, LITBAG, MIC, SOM))
  
  LITBAG_out <- rbind(as.data.frame(LITBAG), 
                      as.data.frame(LIT), 
                      as.data.frame(MIC),
                      as.data.frame(SOM))
  
  LITBAG_out <- LITBAG_out * depth *1e4 / 1e6 
  LITBAG_out <- as.data.frame(t(LITBAG_out))
  colnames(LITBAG_out) <- c("LITBAGm", "LITBAGs", "LITm", "LITs", "MICr", "MICk", "SOMp", "SOMc", "SOMa")
  LITBAG_out <- cbind(data.frame(SITE = forcing_df$Site,
                                 Litter_Type = as.character(litBAG[1]),
                                 DAY=seq(1:nrow(LITBAG_out))), LITBAG_out)
  
  return(LITBAG_out)
}


#------------------------------------------------------------------------------
#--- EXAMPLE USE ---
#------------------------------------------------------------------------------

#--------------------------------------
# Datasets for example simulations
#--------------------------------------

# Forcing data for LTER sites
LTER <- read.csv("Data/LTER_SITE_1.csv", as.is=T)

# LiDET litter bags
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
BAG_out <- MIMICS_LITBAG(forcing_input, BAG_input, nspin_yrs=10, nspin_days=150, litadd_day=125, verbose=T)

# Calc MIMICS C pools
BAG_out$MIMSOC <- rowSums(BAG_out[,6:12])
BAG_out$MIMMIC <- rowSums(BAG_out[,8:9])  
BAG_out$MIMLIT <- rowSums(BAG_out[,6:7])  
BAG_out$MICtoSOC <- BAG_out$MIMMIC/BAG_out$MIMSOC

# Plot litterbag decomp over time
ggplot(BAG_out) + 
  geom_line(aes(y=LITBAGm, x=DAY, colour="LITm"), linewidth=1, linetype="dashed") +
  geom_line(aes(y=LITBAGs, x=DAY, colour="LITs"), linewidth=1, linetype="dashed") +
  geom_line(aes(y=LITBAGs+LITBAGm, x=DAY, colour="Total"), linewidth=1) +
  scale_colour_manual(values = c("Dark Green","Orange","Black")) +
  ylab("Litter Bag C") +
  xlab("Day") +
  labs(title = "MIMICS Litter Bag Simulation",
       subtitle = paste0("Site: ", forcing_input$Site, "\nLitter type: ", BAG_input[1]),
       colour = "LIT Pool") +
  theme_bw()


#---------------------------------------------
# Example: Run multiple litter bags thru decomp simulation
#---------------------------------------------

# BAGS_out <- BAGS %>% split(1:nrow(BAGS)) %>% map(~ MIMICS_LITBAG(litBAG=., 
#                                                                  forcing_df=LTER[6,],
#                                                                  nspin_yrs=20, 
#                                                                  nspin_days=0, 
#                                                                  litadd_day=100, 
#                                                                  verbose=T)) %>% bind_rows()
# ggplot() + 
#   geom_line(data=BAGS_out, aes(y=LITBAGs+LITBAGm, x=DAY, color=Litter_Type), linewidth=1, alpha=0.5) +
#   ylab("Litter Bag C") +
#   xlab("Day") +
#   labs(color="LiDET\nLitter Type") +
#   ggtitle("MIMICS Litter Bag Simulation") +
#   theme_bw()
