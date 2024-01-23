###########################################
# MIMICS repeat run function
###########################################

MIMrepeat <- function(forcing_df, rparams) {

  # Set global model parameters
  Vslope <<- Vslope_default * rparams$Vslope_x[1]
  Vint <<- Vint_default * rparams$Vint_x[1]
  #Kslope <<- Kslope_default * rparams$Kslope_x[1]
  #Kint <<- Kint_default * rparams$Kint_x[1]
    
  #Tau_MULT <<- Tau_MULT_default * rparams$Tau_x[1]
  #CUE <<- CUE_default * rparams$CUE_x[1]
  
  #full run of forcing data csv
  MIMrun <- forcing_df %>% split(1:nrow(forcing_df)) %>% map(MIMICS_SS_pools) %>% bind_rows() 

  #add run number
  MIMrun$run_num <- rparams$run_num[1]
  
  #-----------------------------------------
  ### ADD DATA SUMMARY HERE IF NEEDED
  
  
  #-----------------------------------------
  
  # return MIMrun
  return(MIMrun)
}
