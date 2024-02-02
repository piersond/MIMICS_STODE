########################################
# Set MIMICS parameters
########################################

# --- From pset .csv ----
#Vslope_x	    Vint_x	    Kslope_x	  Kint_x	    Tau_x	      CUE_x	      desorb_x	  fPHYS_x
#3.853107615	1.16299943	0.839868442	1.546675755	0.455164391	0.544205741	0.078479435	0.336044881


Vslope  <- rep(0.063, 6) * 3.853107615  #<----MULTIPLIER
Vint    <- rep(5.47, 6) * 1.16299943 #<----MULTIPLIER
aV      <- rep(0.000008, 6)  
Kslope  <- rep(c(0.025, 0.035, 0.025),2) * 0.839868442  #<----MULTIPLIER
Kint    <- rep(3.19, 6) * 1.546675755  #<----MULTIPLIER
aK      <- rep(10, 6)
vMOD    <- c(10, 2, 10, 3, 3, 2)
kMOD    <- c(8, 2, 4, 2, 4, 6)
KO      <- c(6, 6)
CUE     <- c(0.55, 0.25, 0.75, 0.35) * 0.544205741 #<----MULTIPLIER
tau_r   <- c(0.00052, 0.3)
tau_K   <- c(0.00024, 0.1)
Tau_MOD <- c(100, 0.8, 1.2, 2)
beta    <- c(1.5) # only used if tauMethod='beta'
fPHYS_r <- c(0.3, 1.3)
fPHYS_K <- c(0.2, 0.8)
fCHEM_r <- c(0.1, -3, 1)
fCHEM_K <- c(0.3, -3, 1)
fSOM_p  <- c(1.5e-5, -1.5)
PHYS_scalar <- c(2, -2, NA, NA, NA, NA)
FI      <- c(0.05, 0.05)
fmet_p <- c(1, 0.85, 0.013)

#Set default multipliers
tau_MULT <- 0.455164391 #<----MULTIPLIER
desorb_MULT <- 0.078479435  #<----MULTIPLIER
fPHYS_MULT <- 0.336044881  #<----MULTIPLIER

depth <- 30 # set soil depth
h2y        <- 24*365
MICROtoECO <- depth * 1e4 * 1e-3  # mgC/cm3 to g/m2

# fW coeeficient for Pierson-CORPSE moisture control
fW_p1 <- 1.212580 #* 0.6867031  # MSBio new
fW_p2 <- 2.748028 #* 0.6300376  # MSBio new

#Set default methods
fWmethod <- 0      #0-> fW=1, 1->CORPSE, 2->Calibrated, 3->water scalar from other model
historic <- FALSE   #modify Vmax based on historic MAT
fixed_fMET <- FALSE #calculate fMET based on litter chemistry
tauMethod <- 'NPP'  #'NPP' and 'beta' accepted




