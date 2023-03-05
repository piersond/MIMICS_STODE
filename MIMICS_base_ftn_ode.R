## Set working drive
setwd("C:/github/MIMICS_STODE")

#Libraries
library(rootSolve)
library(deSolve)

library(boot)
library(dplyr)
library(purrr)
library(ggplot2)
library(Metrics)

#bring in RXEQ function
source("RXEQ_ftn.R")

########################################
# Set MIMICS default parameters
########################################
Vslope  <- rep(0.063, 6)
Vint    <- rep(5.47, 6)
aV      <- rep(0.000008, 6)  
Kslope  <- rep(c(0.025, 0.035, 0.025),2)
Kint    <- rep(3.19, 6)
aK      <- rep(10, 6)
vMOD    <- c(10, 2, 10, 3, 3, 2)
kMOD    <- c(8, 2, 4, 2, 4, 6)
KO      <- c(6, 6)
CUE     <- c(0.55, 0.25, 0.75, 0.35)
tau_r   <- c(0.00052, 0.3)
tau_K   <- c(0.00024, 0.1)
Tau_MOD <- c(100, 0.8, 1.2, 2)
Tau_MULT <- 1
fPHYS_r <- c(0.3, 1.3)
fPHYS_K <- c(0.2, 0.8)
fCHEM_r <- c(0.1, -3, 1)
fCHEM_K <- c(0.3, -3, 1)
fSOM_p  <- c(0.000015, -1.5)
PHYS_scalar <- c(2, -2, NA, NA, NA, NA)
FI      <- c(0.05, 0.05)
fmet_p <- c(1, 0.85, 0.013)
depth <- 30 # set soil depth
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
MIMICS1 <- function(df, solver="ODE", step=10000, end=1e7) {
  
  #DEBUG
  #data <- read.csv("LTER_SITE_1.csv", as.is=T)
  #df = data[5,]
  
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
  
  if(solver == "STODE") {
  
    # ------------RUN THE MODEL-------------
    test  <- stode(y = Ty, time = end, fun = RXEQ, parms = Tpars, positive = TRUE)
    
    # ### Calc and get MIMICS output 
    MIMLIT    <- (test[[1]][[1]]+test[[1]][[2]])  * depth *1e4 / 1e6 #convert kgC/m2 from mgC/cm3 (0-30 cm)
    MIMMIC    <- (test[[1]][[3]]+test[[1]][[4]])  * depth *1e4 / 1e6
    MIM_CO    <-  test[[1]][[3]]/test[[1]][[4]]
    MIMSOC    <- sum(test[[1]])  * depth *1e4 / 1e6
  
    table <- as.numeric(test[[1]])
    table2 <- as.numeric(test[[2]])
  
    MIMout <- data.frame(Site = df$Site,
                         fCLAY = fCLAY,
                         TSOI = TSOI,
                         ANPP = ANPP,
                         LIGN = lig_N,
                         EST_LIT = EST_LIT *1e4 / 1e6,
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
                         CO2r = table2[1] * depth *1e4 / 1e6,
                         CO2K = table2[2] * depth *1e4 / 1e6,
                         DEBUG = note)
    return(MIMout)
  }
  
  if(solver == "ODE") {
    
    time_steps = seq(1,end,step)
    n_steps = length(time_steps)
    
    t_sol_ss = ode(
      y = Ty,
      times = time_steps,
      func = RXEQ,
      parms = c(), 
      method = "ode45")
  
    sol = as.data.frame(t_sol_ss)
    colnames(sol) = c("time", "LITm", "LITs", "MICr", "MICK", "SOMp", "SOMc", "SOMa", "CO2r", "CO2K")
    
    forcing = data.frame(Site = rep(df$Site,  n_steps ),
                          fCLAY = rep(fCLAY,  n_steps ),
                          TSOI = rep(TSOI,  n_steps ),
                          ANPP = rep(ANPP,  n_steps ),
                          LIGN = rep(lig_N,  n_steps ),
                          EST_LIT = rep(EST_LIT *1e4 / 1e6,  n_steps ))
    
    # Convert units, * depth *1e4 / 1e6
    for(i in 2:ncol(sol)) {
      sol[,i] <- sol[,i] * depth *1e4 / 1e6
    }
    
    sol$MIMSOC = rowSums(sol[, c(2:8)])
    sol$MIMMIC = rowSums(sol[, c(4:5)])
    sol$MIMLIT = rowSums(sol[, c(2:3)])
    sol$MIM_CO = sol[,4]/sol[,5]
    sol$desorb = as.numeric(desorb)
    sol$SOMpTOv = 1/(as.numeric(desorb)*24*365)
    sol$DEBUG = note
    
    MIMout <- cbind(forcing, sol)
    return(MIMout)
  }
}


#######################################
# Example use
#######################################
data <- read.csv("LTER_SITE_1.csv", as.is=T)
df = data[5,]

# Single site run
MIMout_single <- MIMICS1(df, solver="STODE", step=10000, end=1e7)
MIMout_spinup <- MIMICS1(df, solver="ODE", step=10000, end=1e7)

# Dataset run
MIMout <- data %>% split(1:nrow(data)) %>% map(~ MIMICS1(df=., solver="ODE", step=10000, end=1e7)) %>% bind_rows()

#######################################
# Example plot
#######################################
p <- ggplot(MIMout, aes(x=time, y=MIMSOC, color=Site)) + geom_line(size=1) +
  xlab("Time step") +
  ggtitle("MIMICS ODE spin up for LTER sites") +
  theme_minimal()
p

# Save plot
#png(file = "MIMSOC-spinup.png", height = 500, width=1000)
#p
#dev.off()

