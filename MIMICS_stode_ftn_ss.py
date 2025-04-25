import numpy as np
import pandas as pd
import math
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt
import dask.dataframe as dd
from dask.multiprocessing import get
import time
from multiprocesspandas import applyparallel
from RXEQ import RXEQ
import matplotlib


# Set default parameter values
Vslope = np.tile(np.array([0.063]), 6).astype(float)  #MCMC
Vint = np.tile(np.array([5.47]), 6).astype(float)     #MCMC
aV = np.tile(np.array([0.000008]), 6).astype(float)
Kslope = np.tile(np.array([0.025, 0.035, 0.025]), 2).astype(float)  #MCMC
Kint = np.tile(np.array([3.19]), 6).astype(float)   #MCMC
aK      = np.repeat(np.array([10]), 6).astype(float)
vMOD    = np.array([10, 2, 10, 3, 3, 2]).astype(float)
kMOD    = np.array([8, 2, 4, 2, 4, 6]).astype(float)
KO      = np.array([6, 6]).astype(float)
CUE     = np.array([0.55, 0.25, 0.75, 0.35]).astype(float) #MCMC
tau_r   = np.array([0.00052, 0.3]).astype(float)
tau_K   = np.array([0.00024, 0.1]).astype(float)
Tau_MOD = np.array([100, 0.8, 1.2, 2]).astype(float)
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
Tau_MULT = 1  #MCMC
desorb_MULT = 1  #MCMC
fPHYS_MULT = 1  #MCMC

#-----------------------------------------
# SET PARAMETER FIT MULTIPLIERS
Vslope = Vslope * 3.853107615
Vint = Vint * 1.16299943
Kslope = Kslope * 0.839868442
Kint = Kint * 1.546675755 
CUE = CUE * 0.544205741
Tau_MULT = 0.455164391 
desorb_MULT = 0.078479435  
fPHYS_MULT = 0.336044881  
#------------------------------------------


# Define MIMICS model functions
### MIMICS steady state solver
def MIMICS_SS(df, multi=True): 
    
    if multi:
        ANPP = (df["pGPP"] + 400)/2
        lig_N = df["lig_N"]
        fCLAY = df["CLAY"]/100 
        TSOI = np.tile(np.array([df["TSOI"]]), 6).astype(float)
    else:
        ANPP = (df.pGPP.values[0] + 400)/2
        lig_N = df.lig_N.values[0] 
        fCLAY = df.CLAY.values[0]/100 
        TSOI = np.tile(np.array([df.TSOI.values[0]]), 6).astype(float)
        
    fMET = 0.3846423
    
    #lig = df["LIG"]/100
    #Nnew = 1/df["CN"]/2.5                  	                       
    #fMET = fmet_p[0] * (fmet_p[1] - fmet_p[2] * lig / Nnew) 

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
    
    tau = np.array([tau_r[0] * np.exp(tau_r[1]*fMET), 
                    tau_K[0] * np.exp(tau_K[1]*fMET)])   
    tau = tau * Tau_MOD1 * Tau_MOD2 * Tau_MULT   

    fPHYS = np.array([fPHYS_r[0] * np.exp(fPHYS_r[1]*fCLAY), 
                      fPHYS_K[0] * np.exp(fPHYS_K[1]*fCLAY)]) 
    fCHEM = np.array([fCHEM_r[0] * np.exp(fCHEM_r[1]*fMET) * fCHEM_r[2], 
                      fCHEM_K[0] * np.exp(fCHEM_K[1]*fMET) * fCHEM_K[2]]) 	
    
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

    if len(VMAX) == 1:
        VMAX = VMAX[0]
        KM = KM[0]
        fAVAI = fAVAI[0]

    ### Initialize pools ###          
    I = np.array([(EST_LIT / depth) * fMET, 
                  (EST_LIT / depth) * (1-fMET)]).astype(float)
    lit = I   
    mic = I  
    som = np.array([I[0], I[1], I[0]])
    LITmin = np.tile(np.array([np.nan]), 4).astype(float)
    MICtrn = np.tile(np.array([np.nan]), 6).astype(float)
    SOMmin = np.tile(np.array([np.nan]), 2).astype(float)
    DEsorb = np.tile(np.array([np.nan]), 1).astype(float)
    OXIDAT = np.tile(np.array([np.nan]), 1).astype(float)

    y0_arr = np.array([lit[0], lit[1], mic[0], mic[1], som[0], som[1], som[2]]).ravel()
    y0 = y0_arr.tolist() 
    
    t_eval = [1e6]
    t_span = [0, 1e6]

    sol = solve_ivp(lambda t,y: RXEQ(y, 1, I, VMAX, KM, CUE, fPHYS, fCHEM, fAVAI, FI, tau, LITmin, SOMmin, MICtrn, desorb, DEsorb, OXIDAT, KO), 
                    t_span, y0, method='RK45', t_eval=t_eval)
    
    #return(sol.y.transpose()[0] * depth *1e4 / 1e6)    
    result = sol.y.transpose()[0] * depth * 1e4 / 1e6
    return pd.Series(result, index=['LITm', 'LITs', 'MICr', 'MICk', 'SOMp', 'SOMc', 'SOMa'])

### MIMICS steady state pool calculation and format function
def MIMICS_SS_format(mimout_c_pools):
    mimout_c_pools["MIMSOC"] = mimout_c_pools[["LITm", "LITs", "MICr", "MICk", "SOMp", "SOMc", "SOMa"]].sum(axis=1)
    mimout_c_pools["MIMLIT"] = mimout_c_pools["LITm"] + mimout_c_pools["LITs"]
    mimout_c_pools["MIMMIC"] = mimout_c_pools["MICr"] + mimout_c_pools["MICk"]
    mimout_c_pools["MIC_ratio"] = mimout_c_pools["MICr"] / mimout_c_pools["MICk"]
    return mimout_c_pools


#############################
### EXAMPLE SIMULATIONS
#############################

df = pd.read_csv('Data/MC_cross-val_data.csv', delimiter=',')

# Add unique ID
df['UID'] = df.index + 1

# Filter set == "Validation"
df = df[df["Set"] == "VAL"]

### Method 1: Regular
# start_time = time.time()
# MIMout_set = df.apply(lambda row: MIMICS_SS(row, multi=True), axis=1) #lambda df: df.assign(result=df.col_1 * df.col_2))
# print("--- %s seconds ---" % (time.time() - start_time)) 
# print(MIMout_set)

### Method 2: With multiprocessing (!!! ERROR IN ARRAY SHAPE NEEDS TO BE FIXED to use multi=True)
if __name__ == '__main__':
    tstart = time.time()   
    MIMout_raw = df.groupby(["UID"]).apply_parallel(lambda row: MIMICS_SS(row, multi=False), 
                                                       num_processes=7) #<--- Adjust based on CPU power
    print(time.time() - tstart)
    print(MIMout_raw)
    
    # Save output
    MIMout_raw.to_csv('MIMICS_SS_raw_output.csv', index=False)
    
    # Format and save MIMICS C pools
    MIMout_ss_pools = MIMICS_SS_format(MIMout_raw)
    
    MIMout_ss_pools.to_csv('MIMICS_SS_pools_output.csv', index=False)


    #########################
    # SOC QA plot using matplotlib
    #########################
    
    # Get r^2 and RMSE values for the plot
    # Calculate r^2
    r2 = np.corrcoef(df["SOC"], MIMout_ss_pools["MIMSOC"])[0, 1] ** 2
    # Calculate RMSE
    rmse = np.sqrt(np.mean((df["SOC"] - MIMout_ss_pools["MIMSOC"]) ** 2))
    
    # Plot MIMICS SOC vs. measured SOC
    plt.figure(figsize=(8, 8))  
    plt.scatter(df["SOC"], MIMout_ss_pools["MIMSOC"], color='black', alpha=0.5)
    plt.xlabel('Measured SOC (g/m2)')
    plt.ylabel('MIMICS SOC (g/m2)')
    plt.title('MIMICS SOC vs. Measured SOC')
    
    # Add r^2 and RMSE to the plot
    plt.text(0.1, 0.9, f'R² = {r2:.2f}', transform=plt.gca().transAxes)
    plt.text(0.1, 0.85, f'RMSE = {rmse:.2f}', transform=plt.gca().transAxes)
    
    # Add 1:1 line
    plt.plot([0, 100], [0, 100], color='green', linestyle='--')
    
    # Set limits and grid
    plt.xlim(0, 15)
    plt.ylim(0, 15)
    plt.grid()
    plt.show()
    
    
    
    