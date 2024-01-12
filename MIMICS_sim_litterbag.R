rm(list = ls()) #Clears the working environment

## Set working drive
#setwd("/Users/wwieder/Will/git_repos_local/MIMICS_STODE")
#setwd("C:/github/MIMICS_STODE")

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

MIMICS_LITBAG <- function(forcing_df, litBAG, dailyInput=NA, nspin_yrs=10, nspin_days=200, litadd_day=143, verbose=T){
  
  #DEBUG
  #forcing_df <- LTER[6,]
  #litBAG <- BAGS[1,]
  #nspin_yrs <- 10
  #nspin_days <- 200
  #litadd_day <- 143
  
  if(verbose){
    print("-------------------------------------------------------")
    print(paste0("Starting ", forcing_df$Site, " - ", litBAG[1]))
    print("-------------------------------------------------------")
  }
  
  # Get MIMICS steady_state output
  MIMss <- MIMICS_SS(forcing_df)

  # Create dataframe to store MIMICS pools over timesteps
  MIMfwd = MIMss[[1]]
  MIMfwd = (MIMfwd / depth) / (1e4 / 1e6)  #convert from gC m-2 to mgC cm-3
  
  # Get Tpars from ss simulation  
  Tpars = MIMss[[2]]

  #Init arrays to store daily output data
  nday   <- 365 * nspin_yrs + nspin_days
  day    <- seq(1,nday,1)
  year   <- (day-litadd_day)/365
  doy    <- 1
    
  LIT    <- array(NA, dim = c(2,nday))
  LITBAG <- array(NA, dim = c(2,nday))
  MIC    <- array(NA, dim = c(2,nday))
  SOM    <- array(NA, dim = c(3,nday))
  LITbag_1 = 0. # initial litter bag 
  LITbag_2 = 0

  sim_year = 0
  i = 1

  for (d in 1:nday)  { 
    # For recalc of Tpars from daily forcing data
    if (!is.na(dailyInput)) {
      # Set daily Tpars from daily state variables in "dailyInput" dataframe (added in function arguments)
      Tpars_mod = calc_Tpars(fWmethod = 0, historic=FALSE, LiDET_fMET=FALSE, #<-- ENSURE these settings match the ss run
                             ANPP = dailyInput$ANPP[d], 
                             fCLAY = dailyInput$fCLAY[d], 
                             TSOI = dailyInput$TSOI[d], 
                             MAT = dailyInput$MAT[d], 
                             LIG_N = dailyInput$LIG_N[d],
                             CN = dailyInput$CN[d],  # Only needed if LIG_N not supplied 
                             LIG = dailyInput$LIG[d],  # Only needed if LIG_N not supplied 
                             theta_liq = dailyInput$GWC[d], 
                             theta_frzn = 0) # Change to column name if frozen water content is available 
    } else {
      # Use ss Tpars (i.e., same forcing variables for each sim day)
      Tpars_mod = Tpars
    }  
   
    # Run simulation at hourly timestep
    for (h in 1:24)   {
      Ty <- c( LIT_1 = MIMfwd[1], LIT_2 = MIMfwd[2], 
                 MIC_1 = MIMfwd[3], MIC_2 = MIMfwd[4], 
                 SOM_1 = MIMfwd[5], SOM_2 = MIMfwd[6], 
                 SOM_3 = MIMfwd[7])
      
      # Run MIMICS simulation step
      step = RXEQ(t=NA, y=Ty, pars=Tpars_mod)
      
      # Litter bag decomp simulation step
      LITbag  <- rep(NA,4)
      LITbag[1] <- Ty[[3]] * Tpars_mod$VMAX[1] * LITbag_1 / (Tpars_mod$KM[1] + Ty[[3]])   #MIC_1 mineralization of METABOLIC litter
      LITbag[2] <- Ty[[3]] * Tpars_mod$VMAX[2] * LITbag_2 / (Tpars_mod$KM[2] + Ty[[3]])   #MIC_1 mineralization of STRUC litter
      LITbag[3] <- Ty[[4]] * Tpars_mod$VMAX[4] * LITbag_1 / (Tpars_mod$KM[4] + Ty[[4]])   #mineralization of MET litter
      LITbag[4] <- Ty[[4]] * Tpars_mod$VMAX[5] * LITbag_2 / (Tpars_mod$KM[5] + Ty[[4]])   #mineralization of SRUCTURAL litter

      # Update MIMICS pools
      #---------------------------------------------
      MIMfwd = MIMfwd + unlist(step[[1]]) # MIMICS pools
      
      LITbag_1 <- LITbag_1 - LITbag[1] - LITbag[3] # Litter bag pools
      LITbag_2 <- LITbag_2 - LITbag[2] - LITbag[4] 
      

      # add litter on correct day
      if (d == litadd_day)   {
        if (h == 24)  {
          LITbag_1 <- LITbag_1 + as.numeric(litBAG[3]) # Metabolic pool add
          LITbag_2 <- LITbag_2 + as.numeric(litBAG[4]) # Structural pool add
        }
      }
          
      #write out daily results
      if (h == 24) {
        LIT[1,d] <- MIMfwd[1]
        LIT[2,d] <- MIMfwd[2]
        LITBAG[1,d]  <- LITbag_1
        LITBAG[2,d]  <- LITbag_2
        MIC[1,d] <- MIMfwd[3]
        MIC[2,d] <- MIMfwd[4]
        SOM[1,d] <- MIMfwd[5]
        SOM[2,d] <- MIMfwd[6]
        SOM[3,d] <- MIMfwd[7]
          
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
      } # close h=24 loop	   						
    }		#close hour loop
  }		#close daily loop

  # Compile and format output
  LITBAG_out <- rbind(as.data.frame(LITBAG), 
                      as.data.frame(LIT), 
                      as.data.frame(MIC),
                      as.data.frame(SOM))
  
  LITBAG_out <- LITBAG_out * depth * 1e4 / 1e6 # Soil depth and unit conversion 
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

# Import and calc LiDET litter bag data 
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
BAG_out <- MIMICS_LITBAG(forcing_input, BAG_input, dailyInput=NA, nspin_yrs=10, nspin_days=150, litadd_day=125, verbose=T)

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

BAGS_out <- BAGS %>% split(1:nrow(BAGS)) %>% map(~ MIMICS_LITBAG(litBAG=.,
                                                                 forcing_df=LTER[6,],
                                                                 nspin_yrs=20,
                                                                 nspin_days=0,
                                                                 litadd_day=100,
                                                                 verbose=T)) %>% bind_rows()
ggplot() +
  geom_line(data=BAGS_out, aes(y=LITBAGs+LITBAGm, x=DAY, color=Litter_Type), linewidth=1, alpha=0.5) +
  ylab("Litter Bag C") +
  xlab("Day") +
  labs(color="LiDET\nLitter Type") +
  ggtitle("MIMICS Litter Bag Simulation") +
  theme_bw()
