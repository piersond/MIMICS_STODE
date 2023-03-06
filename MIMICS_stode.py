import numpy as np
import pandas as pd
import math
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt
import time
from RXEQ import RXEQ


# Set default parameter values
Vslope = np.tile(np.array([0.063]), 6).astype(float)
Vint = np.tile(np.array([5.47]), 6).astype(float)
aV = np.tile(np.array([0.000008]), 6).astype(float)
Kslope = np.tile(np.array([0.025, 0.035, 0.025]), 2).astype(float)
Kint = np.tile(np.array([3.19]), 6).astype(float)
aK      = np.repeat(np.array([10]), 6).astype(float)
vMOD    = np.array([10, 2, 10, 3, 3, 2]).astype(float)
kMOD    = np.array([8, 2, 4, 2, 4, 6]).astype(float)
KO      = np.array([6, 6]).astype(float)
CUE     = np.array([0.55, 0.25, 0.75, 0.35]).astype(float)
tau_r   = np.array([0.00052, 0.3]).astype(float)
tau_K   = np.array([0.00024, 0.1]).astype(float)
Tau_MOD = np.array([100, 0.8, 1.2, 2]).astype(float)
Tau_MULT = 1
fPHYS_r = np.array([0.3, 1.3]).astype(float)
fPHYS_K = np.array([0.2, 0.8]).astype(float)
fCHEM_r = np.array([0.1, -3, 1]).astype(float)
fCHEM_K = np.array([0.3, -3, 1]).astype(float)
fSOM_p  = np.array([0.000015, -1.5]).astype(float)
PHYS_scalar = np.array([2, -2, np.nan, np.nan, np.nan, np.nan]).astype(float)
FI = np.array([0.05, 0.05]).astype(float)
fmet_p = np.array([1, 0.85, 0.013]).astype(float)
depth = 30 
h2y = 24 * 365
MICROtoECO = depth * 1e4 * 1e-3  # mgC/cm3 to g/m2

#Set default multipliers
Tau_MULT = 1
desorb_MULT = 1
fPHYS_MULT = 1

# Load forcing data
df = pd.read_csv('Data/LTER_SITE_1.csv', delimiter=',')
print(df["Site"][0])

# Begin MIMICS model
ANPP = df["ANPP"][0]/2
lig_N = df["lig_N"][0] 
fCLAY = df["CLAY2"][0]/100 
TSOI = df["MAT"][0]
 
fMET = 0.3846423 #mean(lit_fMET)

# Calc litter input rate
EST_LIT = (ANPP / (365*24)) * 1e3 / 1e4

### Calculate parameters ###
Vmax = np.exp(TSOI * Vslope + Vint) * aV
Km = np.exp(TSOI * Kslope + Kint) * aK

# Calculate Tau
Tau_MOD1 = math.sqrt(ANPP/Tau_MOD[0])         
Tau_MOD2 = Tau_MOD[3]                        

if(Tau_MOD1 < Tau_MOD[1]):
    Tau_MOD1 = Tau_MOD[1]

if(Tau_MOD1 > Tau_MOD[2]):
    Tau_MOD1 = Tau_MOD[2]
 
tau = np.array([tau_r[0] * np.exp(tau_r[1]*fMET), tau_K[0] * np.exp(tau_K[1]*fMET)])   
tau = tau * Tau_MOD1 * Tau_MOD2 * Tau_MULT   

fPHYS = np.array([fPHYS_r[0] * np.exp(fPHYS_r[1]*fCLAY), fPHYS_K[0] * np.exp(fPHYS_K[1]*fCLAY)]) 
fCHEM = np.array([fCHEM_r[0] * np.exp(fCHEM_r[1]*fMET) * fCHEM_r[2], fCHEM_K[0] * np.exp(fCHEM_K[1]*fMET) * fCHEM_K[2]]) 	
fAVAI = 1 - (fPHYS + fCHEM)
desorb = fSOM_p[0] * np.exp(fSOM_p[1]*(fCLAY))                  

desorb = desorb * desorb_MULT
fPHYS = fPHYS * fPHYS_MULT

pSCALAR = PHYS_scalar[0] * np.exp(PHYS_scalar[1]*(math.sqrt(fCLAY)))
v_MOD = vMOD  
k_MOD = kMOD 
k_MOD[2] = k_MOD[2] * pSCALAR    
k_MOD[5] = k_MOD[5] * pSCALAR    

VMAX = Vmax * v_MOD 
KM = Km / k_MOD

print(VMAX)

### Initialize pools ###          
I = np.array([(EST_LIT / depth) * fMET, (EST_LIT / depth) * (1-fMET)]).astype(float)
lit = I   
mic = I  
som = np.array([I[0], I[1], I[0]])
LITmin = np.tile(np.array([np.nan]), 4).astype(float)
MICtrn = np.tile(np.array([np.nan]), 6).astype(float)
SOMmin = np.tile(np.array([np.nan]), 2).astype(float)
DEsorb = np.tile(np.array([np.nan]), 1).astype(float)
OXIDAT = np.tile(np.array([np.nan]), 1).astype(float)


### MIMICS SOLVE OVER TIMESERIES
y0 = [lit[0], lit[1], mic[0], mic[1], som[0], som[1], som[2]]
t = np.linspace(0, 5000000, 5000)
stode_in = (I, VMAX, KM, CUE, fPHYS, fCHEM, fAVAI, FI, tau, LITmin, SOMmin, MICtrn, desorb, DEsorb, OXIDAT, KO)
start_time = time.time()   
sol = odeint(RXEQ, y0, t, args=stode_in, printmessg=1, hmax=5000)
print("--- %s seconds ---" % (time.time() - start_time))

#Print C pools for last point in the timeseries
print(sol[len(sol)-1,:]* depth *1e4 / 1e6)
MIMout = sol[len(sol)-1,:]

MIMSOC = sum(MIMout) * depth * 1e4 / 1e6
MIMLIT = (MIMout[0] + MIMout[1]) * depth * 1e4 / 1e6
MIMMIC = (MIMout[2] + MIMout[3]) * depth * 1e4 / 1e6
SOMp = MIMout[4] * depth * 1e4 / 1e6
SOMc = MIMout[5] * depth * 1e4 / 1e6
SOMa = MIMout[6] * depth * 1e4 / 1e6

print("MIMSOC: " + str(MIMSOC))
print("MIMLIT: " + str(MIMLIT))
print("MIMMIC: " + str(MIMMIC))
print("SOMp: " + str(SOMp))
print("SOMc: " + str(SOMc))
print("SOMa: " + str(SOMa))


### MIMICS SOLVE FOR SINGLE TIMEPOINT (>3x faster)
t_eval = [1e7]
t_span = [0, 1e7]
# start_time = time.time()  
sol2 = solve_ivp(lambda t,y: RXEQ(y, 1, I, VMAX, KM, CUE, fPHYS, fCHEM, fAVAI, FI, tau, LITmin, SOMmin, MICtrn, desorb, DEsorb, OXIDAT, KO), 
                 t_span, y0, method='RK45', t_eval=t_eval)
# print("--- %s seconds ---" % (time.time() - start_time))

MIMout = sol2.y.transpose()[0]

MIMSOC = sum(MIMout) * depth * 1e4 / 1e6
MIMLIT = (MIMout[0] + MIMout[1]) * depth * 1e4 / 1e6
MIMMIC = (MIMout[2] + MIMout[3]) * depth * 1e4 / 1e6
SOMp = MIMout[4] * depth * 1e4 / 1e6
SOMc = MIMout[5] * depth * 1e4 / 1e6
SOMa = MIMout[6] * depth * 1e4 / 1e6

print("MIMSOC: " + str(MIMSOC))
print("MIMLIT: " + str(MIMLIT))
print("MIMMIC: " + str(MIMMIC))
print("SOMp: " + str(SOMp))
print("SOMc: " + str(SOMc))
print("SOMa: " + str(SOMa))


### Plot C pools timeseries
# plt.plot(t, sol[:, 0] * depth *1e4 / 1e6, 'b', label='LITm')
# plt.plot(t, sol[:, 1] * depth *1e4 / 1e6, 'g', label='LITs')
# plt.plot(t, sol[:, 2] * depth *1e4 / 1e6, 'r', label='MICr')
# plt.plot(t, sol[:, 3] * depth *1e4 / 1e6, 'c', label='MICK')
# plt.plot(t, sol[:, 4] * depth *1e4 / 1e6, 'm', label='SOMp')
# plt.plot(t, sol[:, 5] * depth *1e4 / 1e6, 'y', label='SOMc')
# plt.plot(t, sol[:, 6] * depth *1e4 / 1e6, 'k', label='SOMa')
# plt.legend(loc='best')
# plt.xlabel('t')
# plt.ylim([0, 3])
# plt.grid()
# plt.show()