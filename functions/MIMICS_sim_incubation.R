# Required libraries
library(rootSolve)

###############################################
# MIMICS single point function
#> Input dataframe ("df") must contain columns: SITE, ANPP, fCLAY, TSOI, LIG_N or LIG + CN // Optional: MAT, GWC, W_SCALAR 
#>> With output required for running forward
###############################################
MIMICS_INCUBATION <- function(df, days=105, step="daily", output_type=2){
  
  #DEBUG
  # df = MSB_data[1,]
  # days = 105
  # step = "hourly"
  # output_type = 1
  
  # Convert all column names to upper case
  colnames(df) <- toupper(colnames(df))
  
  ### Setup a var to collect run notes
  note <- ""
  
  ### Bring in forcing ANPP value, convert gDW to gC
  #ANPP <- df$ANPP/2
  ANPP <- 0 #<--- 0 turns off additional C inputs
  
  ### Bring in CLAY value, convert from percent to decimal
  fCLAY <- df$CLAY/100
  
  ### Bring in TSOI value
  TSOI <- df$TINC  #<--- ***Use incubation temperature instead of TSOI***
  
  ### Bring in lig:N forcing data
  LIG_N <- df$LIG_N
  LIG <- df$LIG
  CN <- df$CN

  # Bring in soil moisture information
  theta_liq  <- df$VWC #GWC/100  #<--- ***Using VWC, not GWC. What's best?
  theta_frzn <- 0 

  ### Bring in mean annual temperature data
  MAT <- df$MAT  
  
  # Use TSOI if 'MAT' column is not in forcing data
  if(is.null(MAT)){
    MAT <- TSOI
    print("!!! Warning: Using TSOI for MAT")
  }
  
  #Bring in W_SCALAR if present
  W_SCALAR = df$W_SCALAR

  
  ############################################################
  # MIMICS SIMULATION STARTS HERE
  ############################################################
  
  # function calculates fMET with LIG_N if provided in input data.
  Tpars <- calc_Tpars_Conly(ANPP=ANPP, fCLAY=fCLAY, TSOI=TSOI, MAT=MAT,     
                            CN=CN, LIG=LIG, LIG_N=LIG_N,
                            theta_liq=theta_liq, theta_frzn=theta_frzn, W_SCALAR = W_SCALAR) 
  
  # Set fMET
  fMET  <- fmet_p[1] * (fmet_p[2] - fmet_p[3] * LIG_N)   
  
  # Open matrices to store model output
  LITmin  <- rep(NA, dim=4)
  MICtrn  <- rep(NA, dim=6)
  SOMmin  <- rep(NA, dim=2)
  DEsorb  <- rep(NA, dim=1)
  OXIDAT  <- rep(NA, dim=1)
  
  # Initialize model pools and fluxes
  LIT_INIT_TOT = 1000 # Init litter pool for incubation 
  I        <- rep(0,2)
  LIT_1    <- LIT_INIT_TOT * fMET      
  LIT_2    <- LIT_INIT_TOT * (1-fMET)
  MIC_1    <- 0.01
  MIC_2    <- 0.01
  SOM_1    <- 0    #<--- How much is coming back MIC from SOM? Maybe that's correct for an incubation?
  SOM_2    <- 0
  SOM_3    <- 0
  CO2_1    <- 0
  CO2_2    <- 0
  
  # Create vector of parameter values
  # tpars <- c(I = I, VMAX = VMAX, KM = KM, CUE = CUE, 
  #            fPHYS = fPHYS, fCHEM = fCHEM, fAVAI = fAVAI, FI = FI, 
  #            tau   = tau, LITmin = LITmin, SOMmin = SOMmin, MICtrn = MICtrn, 
  #            desorb= desorb, DEsorb = DEsorb, OXIDAT = OXIDAT, KO = KO)


  # Create dataframe to store simulation output
  if(step == "hourly"){
    MIMout <- data.frame(ID = rep(df$ID, 24*days),
                         SITE = rep(df$SITE, 24*days),
                         DAY = rep(NA, 24*days),
                         HOUR = rep(NA, 24*days),
                         LITm = rep(NA, 24*days),
                         LITs = rep(NA, 24*days),
                         MICr = rep(NA, 24*days),
                         MICK = rep(NA, 24*days),
                         SOMp = rep(NA, 24*days),
                         SOMc = rep(NA, 24*days),
                         SOMa = rep(NA, 24*days),
                         CO2_MICr = rep(NA, 24*days),
                         CO2_MICK = rep(NA, 24*days),
                         CO2_prop_totC= rep(NA, 24*days))
  } else if(step == "daily"){
      MIMout <- data.frame(ID = rep(df$ID, days),
                           SITE = rep(df$SITE, days),
                           DAY = rep(NA, days),
                           HOUR = rep(NA, days),
                           LITm = rep(NA, days),
                           LITs = rep(NA, days),
                           MICr = rep(NA, days),
                           MICK = rep(NA, days),
                           SOMp = rep(NA, days),
                           SOMc = rep(NA, days),
                           SOMa = rep(NA, days),
                           CO2_MICr = rep(NA, days),
                           CO2_MICK = rep(NA, days),
                           CO2_prop_totC= rep(NA, days))
  }
  
  ### BEGIN MODEL LOOP ###
  for (d in 1:days) {
    for (h in 1:24) {
      if(step == "hourly"){
        update <- RXEQ(y = c(LIT_1 = LIT_1, LIT_2 = LIT_2, 
                             MIC_1 = MIC_1, MIC_2 = MIC_2, 
                             SOM_1 = SOM_1, SOM_2 = SOM_2, 
                             SOM_3 = SOM_3),
                       pars = Tpars)
        
        # Update C pools
        LIT_1  <- LIT_1 + update[[1]][1]
        LIT_2  <- LIT_2 + update[[1]][2]
        MIC_1  <- MIC_1 + update[[1]][3]
        MIC_2  <- MIC_2 + update[[1]][4]
        SOM_1  <- SOM_1 + update[[1]][5]
        SOM_2  <- SOM_2 + update[[1]][6]
        SOM_3  <- SOM_3 + update[[1]][7]
        CO2_1  <- CO2_1 + update[[2]][1]
        CO2_2  <- CO2_2 + update[[2]][2]        
        
        # Set iteration number for indexing
        iter <- (d*24)-24+h
        
        # Store hourly output       
        MIMout$DAY[iter] <- d
        MIMout$HOUR[iter] <- h
        MIMout$LITm[iter] <- LIT_1
        MIMout$LITs[iter] <- LIT_2
        MIMout$MICr[iter] <- MIC_1
        MIMout$MICK[iter] <- MIC_2
        MIMout$SOMp[iter] <- SOM_1
        MIMout$SOMc[iter] <- SOM_2
        MIMout$SOMa[iter] <- SOM_3
        MIMout$CO2_MICr[iter] <- CO2_1
        MIMout$CO2_MICK[iter] <- CO2_2
        MIMout$CO2_prop_totC[iter] <- (MIMout$CO2_MICr[iter] + MIMout$CO2_MICK[iter]) / LIT_INIT_TOT
      
      } else if(step == "daily"){
          if(h == 24) {
            # Get model output from RXEQ ftn
            update <- RXEQ(y = c(LIT_1 = LIT_1, LIT_2 = LIT_2, 
                                     MIC_1 = MIC_1, MIC_2 = MIC_2, 
                                     SOM_1 = SOM_1, SOM_2 = SOM_2, 
                                     SOM_3 = SOM_3),
                               pars = Tpars)
            update[[1]] = update[[1]] * 24
            update[[2]] = update[[2]] * 24
        
            # Update C pools
            LIT_1  <- LIT_1 + update[[1]][1]
            LIT_2  <- LIT_2 + update[[1]][2]
            MIC_1  <- MIC_1 + update[[1]][3]
            MIC_2  <- MIC_2 + update[[1]][4]
            SOM_1  <- SOM_1 + update[[1]][5]
            SOM_2  <- SOM_2 + update[[1]][6]
            SOM_3  <- SOM_3 + update[[1]][7]
            CO2_1  <- CO2_1 + update[[2]][1]
            CO2_2  <- CO2_2 + update[[2]][2]
            
            # Store daily output
            MIMout$DAY[d] <- d
            MIMout$HOUR[d] <- h
            MIMout$LITm[d] <- LIT_1
            MIMout$LITs[d] <- LIT_2
            MIMout$MICr[d] <- MIC_1
            MIMout$MICK[d] <- MIC_2
            MIMout$SOMp[d] <- SOM_1
            MIMout$SOMc[d] <- SOM_2
            MIMout$SOMa[d] <- SOM_3
            MIMout$CO2_MICr[d] <- CO2_1
            MIMout$CO2_MICK[d] <- CO2_2
            MIMout$CO2_prop_totC[d] <- (MIMout$CO2_MICr[d] + MIMout$CO2_MICK[d]) / LIT_INIT_TOT
          }
      }
    }
  }
  #############################
  ### SET SIMULATION OUTPUT ###
  #############################
  
  if(output_type == 1) {
    return(MIMout)  
  } else {
    ftn_output <- MIMout %>% filter(DAY == days) 
    return(ftn_output)    
  }
}


