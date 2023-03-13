## Set working drive
setwd("C:/github/MIMICS_STODE")

#Libraries
library(rootSolve)
library(boot)
library(dplyr)
library(purrr)
library(furrr)
library(ggplot2)
library(Metrics)

#bring in RXEQ function
source("MIM_ftns/RXEQ_ftn.R")
source("MIM_ftns/MIMICS_set_parameters.R")

###########################################
# MIMICS single point function
###########################################
MIMICS1 <- function(df){
  
  #DEBUG
  #df = data[1,]
  
  ### Setup a var to collect run notes
  note <- ""
  
  ###Bring in forcing ANPP value
  ANPP <- df$ANPP/2
  
  ### Bring in CLAY value, convert from percent to decimal
  fCLAY <- df$CLAY/100
  
  ### Bring in TSOI value
  TSOI <- df$TSOI
  
  ### Bring in lig:N forcing data
  lig_N <- df$lig_N
  
  ###########################################################
  ### Set fMET equation
  ###########################################################
  ## Option A: Defualt fMET equation using lig:N values
  #fMET <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * lig_N) 
  
  ## Option B: LTER "SHORTCUT" fMET value (average from LiDET)
  fMET <- 0.3846423
  
  ############################################################
  # MIMICS MODEL CODE STARTS HERE
  ############################################################
  
  # Calc litter input rate
  EST_LIT <- (ANPP / (365*24)) * 1e3 / 1e4
  #print(EST_LIT)# gC/m2/h (from gC/m2/y) then mgC/cm2/h(from gC/m2/h) 
  
  # ------------ caclulate parameters ---------------
  Vmax     <- exp(TSOI * Vslope + Vint) * aV 
  Km       <- exp(TSOI * Kslope + Kint) * aK
  
  #ANPP strongly correlated with MAP
  Tau_MOD1 <- sqrt(ANPP/Tau_MOD[1])         
  Tau_MOD2 <- Tau_MOD[4]                        
  Tau_MOD1[Tau_MOD1 < Tau_MOD[2]] <- Tau_MOD[2]
  Tau_MOD1[Tau_MOD1 > Tau_MOD[3]] <- Tau_MOD[3] 
  
  tau <- c(tau_r[1]*exp(tau_r[2]*fMET), 
           tau_K[1]*exp(tau_K[2]*fMET))   
  tau <- tau * Tau_MOD1 * Tau_MOD2 * Tau_MULT 
  
  fPHYS    <- c(fPHYS_r[1] * exp(fPHYS_r[2]*fCLAY), 
                fPHYS_K[1] * exp(fPHYS_K[2]*fCLAY)) 	            
  fCHEM    <- c(fCHEM_r[1] * exp(fCHEM_r[2]*fMET) * fCHEM_r[3], 
                fCHEM_K[1] * exp(fCHEM_K[2]*fMET) * fCHEM_K[3]) 	
  fAVAI    <- 1 - (fPHYS + fCHEM)
  
  desorb   <- fSOM_p[1] * exp(fSOM_p[2]*(fCLAY))                  
  
  desorb <- desorb * desorb_MULT
  fPHYS <- fPHYS * fPHYS_MULT
  
  pSCALAR  <- PHYS_scalar[1] * exp(PHYS_scalar[2]*(sqrt(fCLAY)))  #Scalar for texture effects on SOMp
  
  v_MOD    <- vMOD  
  k_MOD    <- kMOD 
  k_MOD[3] <- k_MOD[3] * pSCALAR    
  k_MOD[6] <- k_MOD[6] * pSCALAR    
  
  VMAX     <- Vmax * v_MOD 
  KM       <- Km / k_MOD
  
  #----------initialize pools---------------
  I       <- array(NA, dim=2)             
  I[1]    <- (EST_LIT / depth) * fMET     
  I[2]    <- (EST_LIT / depth) * (1-fMET)
  lit     <- I   
  mic     <- I  
  som     <- rep(NA, 3) 
  som[1]  <- I[1]
  som[2]  <- I[2]
  som[3]  <- I[1] 
  LITmin  <- rep(NA, dim=4)
  MICtrn  <- c(NA,NA,NA,NA,NA,NA)
  SOMmin  <- rep(NA, dim=2)
  DEsorb  <- rep(NA, dim=1)
  OXIDAT  <- rep(NA, dim=1)
  
  #Calculate RXEQ pools  
  Tpars <- c( I = I, VMAX = VMAX, KM = KM, CUE = CUE, 
              fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
              tau = tau, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
              desorb = desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO)
  
  Ty    <- c( LIT_1 = lit[1], LIT_2 = lit[2], 
              MIC_1 = mic[1], MIC_2 = mic[2], 
              SOM_1 = som[1], SOM_2 = som[2], SOM_3 = som[3])
  
  ## Set global parameters to allow pass to stode function
  .GlobalEnv$VMAX <- VMAX
  .GlobalEnv$KM <- KM
  .GlobalEnv$fPHYS <- fPHYS
  .GlobalEnv$fCHEM <- fCHEM
  .GlobalEnv$fAVAI <- fAVAI
  .GlobalEnv$I <- I
  .GlobalEnv$tau <- tau
  .GlobalEnv$LITmin <- LITmin
  .GlobalEnv$SOMmin <- SOMmin
  .GlobalEnv$MICtrn <- MICtrn
  .GlobalEnv$desorb <- desorb
  .GlobalEnv$DEsorb <- DEsorb
  .GlobalEnv$OXIDAT <- OXIDAT
  
  # ------------RUN THE MODEL-------------
  test  <- stode(y = Ty, time = 1e6, fun = RXEQ, parms = Tpars, positive = TRUE)
  
  ### Calc and get MIMICS output 
  MIMLIT    <- (test[[1]][[1]]+test[[1]][[2]])  * depth *1e4 / 1e6 #convert kgC/m2 from mgC/cm3 (0-30 cm) 
  MIMMIC    <- (test[[1]][[3]]+test[[1]][[4]])  * depth *1e4 / 1e6
  MIM_CO    <-  test[[1]][[3]]/test[[1]][[4]]
  MIMSOC    <- sum(test[[1]])  * depth *1e4 / 1e6   
  
  table <- as.numeric(test[[1]])
  
  MIMout <- data.frame(Site = df$Site,
                       fCLAY = fCLAY,
                       TSOI = TSOI,
                       ANPP = ANPP,
                       LIGN = lig_N,
                       EST_LIT = EST_LIT,
                       MIMSOC = MIMSOC,
                       MIMMIC = MIMMIC,
                       MIMLIT = MIMLIT,
                       MIM_CO = MIM_CO,
                       desorb = as.numeric(desorb),
                       SOMpTOv = 1/(as.numeric(desorb)*24*365), #convert from per hr to per yr
                       LITm = table[1] * depth *1e4 / 1e6, #convert kgC/m2 from mgC/cm3 (0-30 cm) 
                       LITs = table[2] * depth *1e4 / 1e6,
                       MICr = table[3] * depth *1e4 / 1e6,
                       MICK = table[4] * depth *1e4 / 1e6,
                       SOMp = table[5] * depth *1e4 / 1e6,
                       SOMc = table[6] * depth *1e4 / 1e6,
                       SOMa = table[7] * depth *1e4 / 1e6,
                       #CO2r = table[8],
                       #CO2r = table[9],
                       JITn = "",
                       DEBUG = note)
  return(MIMout)
}

#####################
# Example use of 
#####################

# ###############################################
# #> Single point run
# ###############################################
# data <- data.frame(Site = 1,
#                    ANPP = 141,
#                    TSOI = -7.0,
#                    CLAY = 5,
#                    lig_N = 0)
# 
# MIMout_single <- MIMICS1(data[1,])
# 
# 
###############################################
#>  Dataset run from .csv
###############################################
# data <- read.csv("Data/LTER_SITE_1.csv", as.is=T)
# 
# MIMrun <- data %>% split(1:nrow(data)) %>% map(~ MIMICS1(df=.)) %>% bind_rows()
# MIMrun <- data[,1:2] %>% cbind(MIMrun %>% select(-Site))  # Bind data info columns to MIMICS output
# 
# 
# #################################################
# #> Plot SOC vs MIMSOC
# #################################################
# 
# plot_data <- MIMrun
# 
# #calc SOMp turnover time
# plot_data$desorb_yr <- plot_data$desorb*24*365
# plot_data$SOMpTO <- plot_data$SOMp/plot_data$desorb_yr
# 
# # Calc SOC vs. MIMSOC r2 and RMSE
# r2_test <- cor.test(MIMrun$SOC, MIMrun$MIMSOC)
# r_val <- round(as.numeric(unlist(r2_test ['estimate'])),2)
# lb2 <- paste("R^2 == ", r_val)
# 
# rmse <- round(rmse(plot_data$SOC, plot_data$MIMSOC),2)
# 
# # Plot SOC vs. MIMSOC
# ggplot(plot_data, aes(x=MIMSOC, y=SOC, color=TSOI)) +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
#   geom_point(size=4, alpha=0.8) +
#   geom_text(aes(label=paste0(Site)),hjust=-0.2, vjust=0.2) +
#   annotate("text", label = lb2, x = 2, y = 8.5, size = 4, colour = "black", parse=T) +
#   annotate("text", label = paste0("RMSE = ", rmse), x = 2, y = 7.4, size = 4, colour = "black") +
#   ylim(0,10) + xlim(0,10) +
#   theme_minimal()

