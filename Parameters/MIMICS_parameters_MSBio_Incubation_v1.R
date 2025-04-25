########################################
# Set MIMICS parameters
########################################
Vslope  <- rep(0.063, 6) * 0.7585359 #<-- CALIBRATION FACTOR
Vint    <- rep(5.47, 6) * 1.0470140  #<-- CALIBRATION FACTOR
#aV      <- rep(0.000008, 6)
aV      <- rep(0.000000075, 6) #***UNIQUE TO INCUBATION. WHY USING THIS ADJ? COMES FROM OLD~SULMAN INC SCRIPT
#Kslope  <- rep(c(0.025, 0.035, 0.025),2
Kslope  <- rep(0.02, 6) * 0.6955175  #<-- UNIQUE + CALIBRATION FACTOR
Kint    <- rep(3.19, 6) * 0.7543530  #<-- CALIBRATION FACTOR
#aK      <- rep(10, 6)
aK      <- rep(0.15625, 6)  #***UNIQUE TO INCUBATION. WHY?
#vMOD    <- c(10, 2, 10, 3, 3, 2)
vMOD    <- c(2, 0.4, 2, 0.6, 0.6, 0.4) #***UNIQUE TO INCUBATION. WHY?
kMOD    <- c(8, 2, 4, 2, 4, 6)
KO      <- c(6, 6)
#CUE     <- c(0.55, 0.25, 0.75, 0.35)
CUE     <- c(0.5, 0.25, 0.7, 0.35) #***UNIQUE TO INCUBATION. WHY?
tau_r   <- c(0.00052, 0.3)
tau_K   <- c(0.00024, 0.1)
#Tau_MOD <- c(100, 0.8, 1.2, 2)
Tau_MOD <- c(100, 0.6, 1.3, 3.5)  #***UNIQUE TO INCUBATION. WHY? 
beta    <- c(1.1) # only used if tauMethod='beta'
fPHYS_r <- c(0, 0) # Turning off SOM_P partitioning
fPHYS_K <- c(0, 0) # Turning off SOM_P partitioning
fCHEM_r <- c(0.1, -3, 1)
fCHEM_K <- c(0.3, -3, 1)
fSOM_p  <- c(1.5e-5, -1.5)
PHYS_scalar <- c(2, -2, NA, NA, NA, NA)
FI      <- c(0.05, 0.05)
fmet_p <- c(1, 0.85, 0.013)

depth <- 20 # set soil depth
h2y        <- 24*365
MICROtoECO <- depth * 1e4 * 1e-3  # mgC/cm3 to g/m2

# fW coeeficient for Pierson-CORPSE moisture control
fW_p1 <- 0.5809101 #* 0.6867031  # MSBio new
fW_p2 <- 0.8590644 #* 0.6300376  # MSBio new

# Historic MAT based kinetic mod coefficients
#Vslope_MOD = c(1.069356e-02, 5.797809e-01)   #<-- OLD MAT DEPENDENCE PARAMETERS
#Vint_MOD = c( -3.169280e-04, 1.050568)       #<-- OLD MAT DEPENDENCE PARAMETERS
Vh_MOD = c(0.00104, 0.0228) #<-- NEW MAT DEPENDENCE PARAMETERS

#Set default methods
fWmethod=2      #0-> fW=1, 1->CORPSE, 2->Calibrated, 3->water scalar from other model
historic=TRUE  #modify Vmax based on historic MAT
fixed_fMET=FALSE #calculate fMET based on litter chemistry
tauMethod='beta'  #'NPP' and 'beta' accepted

#Set default multipliers
tau_MULT = 1
desorb_MULT = 1
fPHYS_MULT = 1


