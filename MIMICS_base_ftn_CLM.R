## Set working drive
setwd("C:/github/MIMICS_STODE")

#Libraries
library(rootSolve)
library(boot)
library(dplyr)
library(purrr)
library(ggplot2)
library(Metrics)

#bring in RXEQ function
source("RXEQ_ftn.R")

# Here are the MIMICS parameter's we're using in CLM.  
# The modifiers (aV and vmod) are multiplied here in a single product (mimics_vmod)

#// mimics_cn_k = 10 ;
#// mimics_cn_mod_num = 0.4 ;
#// mimics_cn_r = 6 ;
#// mimics_densdep = 1 ;
#?? mimics_desorp = 1.05e-06, -2 ;
#?? mimics_desorpQ10 = 1 ;
#++ mimics_fmet = 0.75, 0.85, 0.013, 40 ;
#== mimics_fphys_k = 0.02, 0.8 ;
#== mimics_fphys_r = 0.03, 1.3 ;
#== mimics_kint = 3.19, 3.19, 3.19, 3.19, 3.19, 3.19, 0, 0 ;
#++ mimics_kmod = 0.001953125, 0.0078125, 0.00390625, 0.0078125, 0.00390625, 0.002604167, 0, 0 ;
#++ mimics_ko_k = 4 ;
#++ mimics_ko_r = 4 ;
#++ mimics_kslope = 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0, 0 ;
#++ mimics_mge = 0.5, 0.25, 0.5, 0.7, 0.35, 0.7, 0, 0 ;
#// mimics_nue_into_mic = 0.85 ;
#++ mimics_p_scalar = 0.8, -3 ;
#?? mimics_t_soi_ref = 25 ;
#== mimics_tau_k = 0.00024, 0.1 ;
#?? mimics_tau_mod_factor = 0.01 ;
#?? mimics_tau_mod_max = 1.3 ;
#?? mimics_tau_mod_min = 0.6 ;
#== mimics_tau_r = 0.00052, 0.3 ;
#== mimics_vint = 5.47, 5.47, 5.47, 5.47, 5.47, 5.47, 0, 0 ;
#++ mimics_vmod = 1.25e-07, 2.5e-08, 1.25e-07, 3.75e-08, 3.75e-08, 2.5e-08, 0, 0 ;
#== mimics_vslope = 0.063, 0.063, 0.063, 0.063, 0.063, 0.063, 0, 0 ;

########################################
# Set MIMICS default parameters
########################################
Vslope  <- rep(0.063, 6)
Vint    <- rep(5.47, 6)
#aV      <- rep(0.000008, 6)  
aV <- 1

#Kslope  <- rep(c(0.025, 0.035, 0.025),2)
Kslope  <- rep(c(0.02, 6))

Kint    <- rep(3.19, 6)
aK      <- rep(10, 6)
#vMOD    <- c(10, 2, 10, 3, 3, 2)
vMOD <- c(1.25e-07, 2.5e-08, 1.25e-07, 3.75e-08, 3.75e-08, 2.5e-08, 0, 0)

#kMOD    <- c(8, 2, 4, 2, 4, 6)
kMOD = c(0.001953125, 0.0078125, 0.00390625, 0.0078125, 0.00390625, 0.002604167, 0, 0)

#KO      <- c(6, 6)
KO      <- c(4, 4)

#CUE     <- c(0.55, 0.25, 0.75, 0.35)
CUE <- c(0.5, 0.25, 0.5, 0.7, 0.35, 0.7, 0, 0)

tau_r   <- c(0.00052, 0.3)
tau_K   <- c(0.00024, 0.1)
#Tau_MOD <- c(100, 0.8, 1.2, 2)
Tau_MOD <- c(0.1, 1.3, 0.6, 1)

Tau_MULT <- 1
fPHYS_r <- c(0.3, 1.3)
fPHYS_K <- c(0.2, 0.8)
fCHEM_r <- c(0.1, -3, 1)
fCHEM_K <- c(0.3, -3, 1)
fSOM_p  <- c(0.000015, -1.5)

#PHYS_scalar <- c(2, -2, NA, NA, NA, NA)
PHYS_scalar <- c(0.8, -3)

FI      <- c(0.05, 0.05)

#fmet_p <- c(1, 0.85, 0.013)
fmet_p = c(0.75, 0.85, 0.013, 40)

depth <- 20 # set soil depth
h2y        <- 24*365
MICROtoECO <- depth * 1e4 * 1e-3  # mgC/cm3 to g/m2

#Set default multipliers
Tau_MULT = 1
desorb_MULT = 1
fPHYS_MULT = 1


########################################
# Apply parameter multipliers
########################################
# Vslope = Vslope * 1.693578
# Vint = Vint * 0.633318
# Kslope = Kslope * 1.782366
# Kint = Kint * 0.3609913
# CUE = CUE * 1
# Tau_MULT = 1
# desorb_MULT = 2.3635554
# fPHYS_MULT = 2.0716163

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
                       DEBUG = note)
  return(MIMout)
}

#####################
# Example use of 
#####################

###############################################
#> Single point run
###############################################
data <- data.frame(Site = 1,
                   ANPP = 141,
                   TSOI = -7.0,
                   CLAY = 5,
                   lig_N = 0)

MIMout_single <- MIMICS1(data[1,])


###############################################
#>  Dataset run from .csv
###############################################
data <- read.csv("LTER_SITE_1.csv", as.is=T)

MIMrun <- data %>% split(1:nrow(data)) %>% map(~ MIMICS1(df=.)) %>% bind_rows()
MIMrun <- data[,1:2] %>% cbind(MIMrun %>% select(-Site))  # Bind data info columns to MIMICS output


#################################################
#> Plot SOC vs MIMSOC
#################################################

plot_data <- MIMrun

#calc SOMp turnover time
plot_data$desorb_yr <- plot_data$desorb*24*365
plot_data$SOMpTO <- plot_data$SOMp/plot_data$desorb_yr

# Calc SOC vs. MIMSOC r2 and RMSE
r2_test <- cor.test(MIMrun$SOC, MIMrun$MIMSOC)
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

