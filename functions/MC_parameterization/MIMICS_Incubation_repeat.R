###########################################
# MIMICS repeat run function for litterbag simulations 
###########################################

INC_MIMrepeat <- function(forcing_df, rparams) {
  
  # # Set global model parameters
  Vslope <<- Vslope_default * rparams$Vslope_x[1]
  Vint <<- Vint_default * rparams$Vint_x[1]
  Kslope <<- Kslope_default * rparams$Kslope_x[1]
  Kint <<- Kint_default * rparams$Kint_x[1]
  # Tau_MULT <<- Tau_MULT_default * rparams$Tau_x[1] #change where Tau_MULT comes in in MIMICS_INC_daily to make this just influence r or K
  # Tau_MULT.r <<- Tau_MULT_default * rparams$Tau_r[1]
  # Tau_MULT.k <<- Tau_MULT_default * rparams$Tau_k[1]
  # CUE <<- CUE_default * rparams$CUE_x[1]
  # CUE <<- c(CUE_default[1] * rparams$CUE_x[1], CUE_default[2] * rparams$CUE_x[1], CUE_default[3], CUE_default[4])  #add indexing to just get r or K-selected
  # CUE <<- c(CUE_default[1] * rparams$CUE_r[1], CUE_default[2] * rparams$CUE_r[1], CUE_default[3] * rparams$CUE_k[1], CUE_default[4] * rparams$CUE_k[1]) #seperate for r and K
  # vMOD <<- vMOD_default * rparams$vMOD_x[1]
  # kMOD <<- kMOD_default * rparams$kMOD_x[1]
  
  #full run of forcing data csv
  MIMrun <- forcing_df %>% split(1:nrow(forcing_df)) %>% map(MIMICS_INCUBATION) %>% bind_rows() 

  #add run number
  MIMrun$run_num <- rparams$run_num[1]
  
  # return output
  return(MIMrun)
}

