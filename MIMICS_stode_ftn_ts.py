import numpy as np
import pandas as pd
import math
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import dask.dataframe as dd
from dask.multiprocessing import get
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

def MIMICS_TS(df, arr=False): 
     
     # Begin MIMICS model
   
    if not arr:
        Site = df["Site"]
        ANPP = df["ANPP"]/2
        lig_N = df["lig_N"]
        fCLAY = df["CLAY2"]/100 
        TSOI = np.tile(np.array([df["MAT"]]), 6).astype(float)
    else:
        Site = df.Site.values[0]
        ANPP = df.ANPP.values[0]/2
        lig_N = df.lig_N.values[0] 
        fCLAY = df.CLAY2.values[0]/100 
        TSOI = np.tile(np.array([df.MAT.values[0]]), 6).astype(float)
        
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

    if len(VMAX) == 1:
        VMAX = VMAX[0]
        KM = KM[0]
        fAVAI = fAVAI[0]

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
    t = np.linspace(0, 10000000, 10000)
    stode_in = (I, VMAX, KM, CUE, fPHYS, fCHEM, fAVAI, FI, tau, LITmin, SOMmin, MICtrn, desorb, DEsorb, OXIDAT, KO)
    sol = odeint(RXEQ, y0, t, args=stode_in, hmax=5000)
    
    return(sol)    


####################
### EXAMPLES 
####################

df = pd.read_csv('Data/LTER_SITE_1.csv', delimiter=',')

# ### Method 1: Apply function each row in dataframe (84.5 sec for RCrk_all)
# start_time = time.time()
# MIMout_set = df.apply(lambda row: MIMICS_TS(row, arr=False), axis=1) 
# print("--- %s seconds ---" % (time.time() - start_time)) 
# #print(MIMout_set)

# MIMout1 = MIMout_set[0]
# MIMpools = MIMout1[len(MIMout1)-1,:] 

# MIMSOC = sum(MIMpools) * depth * 1e4 / 1e6
# MIMLIT = (MIMpools[0] + MIMpools[1]) * depth * 1e4 / 1e6
# MIMMIC = (MIMpools[2] + MIMpools[3]) * depth * 1e4 / 1e6
# SOMp = MIMpools[4] * depth * 1e4 / 1e6
# SOMc = MIMpools[5] * depth * 1e4 / 1e6
# SOMa = MIMpools[6] * depth * 1e4 / 1e6
# print("MIMSOC: " + str(MIMSOC))
# print("MIMLIT: " + str(MIMLIT))
# print("MIMMIC: " + str(MIMMIC))
# print("SOMp: " + str(SOMp))
# print("SOMc: " + str(SOMc))
# print("SOMa: " + str(SOMa))

# # Create time series plot of MIMICS C pools
# t = np.linspace(0, 10000000, 10000)
# plt.plot(t, MIMout1[:, 0], 'b', label='LITm')
# plt.plot(t, MIMout1[:, 1], 'g', label='LITs')
# plt.plot(t, MIMout1[:, 2], 'r', label='MICr')
# plt.plot(t, MIMout1[:, 3], 'c', label='MICK')
# plt.plot(t, MIMout1[:, 4], 'm', label='SOMp')
# plt.plot(t, MIMout1[:, 5], 'y', label='SOMc')
# plt.plot(t, MIMout1[:, 6], 'k', label='SOMa')
# plt.legend(loc='best')
# plt.xlabel('t')
# plt.ylim([0, 5])
# plt.grid()
# plt.show()

### Method 2: Use Dask to multiprocess (50.2 sec using 10 partitions for RCrk_all)
# from multiprocessing import Process, freeze_support
# if __name__ == '__main__':
#     freeze_support()  # needed for Windows

#     ddf = dd.from_pandas(df, npartitions=10) #<--- Adjust based on CPU power
#     start_time = time.time()                              
#     MIMout_set = ddf.map_partitions(lambda df: df.apply((lambda row: MIMICS_TS(row, arr=False)), 
#                                                         axis=1)).compute(scheduler='processes') 
#     print("--- %s seconds ---" % (time.time() - start_time)) 
#     print(MIMout_set)


# ### Method 3: Use multiprocess (26 sec using 7 processes for RCrk_all)
# import multiprocessing as mp
# from multiprocessing import  Pool
# from functools import partial
# import numpy as np

# def parallelize(data, func, num_of_processes=8):
#     data_split = np.array_split(data, num_of_processes)
#     pool = Pool(num_of_processes)
#     data = pd.concat(pool.map(func, data_split))
#     pool.close()
#     pool.join()
#     return data

# def run_on_subset(func, data_subset):
#     return data_subset.apply(func, axis=1)

# def parallelize_on_rows(data, func, num_of_processes=7):
#     return parallelize(data, partial(run_on_subset, func), num_of_processes)

# if __name__ == '__main__':
#     tstart = time.time()   
#     MIMout = parallelize_on_rows(df, MIMICS_TS)
#     print(time.time() - tstart)
#     ### Returns the data, but still need to figure out how to parse it

 
### Method 4: Use multiprocesspandas (26 sec using 7 processes for RCrk_all)
        ### ~Same time savings as multiprocess in functions, but implimentation is far less complex
from multiprocesspandas import applyparallel
if __name__ == '__main__':
    #tstart = time.time()
    MIMout_set = df.groupby(["Site"]).apply_parallel(lambda row: MIMICS_TS(row, arr=True), 
                                                       num_processes=7) #<--- Adjust based on CPU power
    #print(time.time() - tstart)
    
    Site_names = MIMout_set.index.values
    
    MIMout1 = MIMout_set.to_numpy()[9][0] # Grabs C pool timeseries' for first site (row) in dataset 
    #print(MIMout1)

    # Create time series plot of MIMICS C pools
    t = np.linspace(0, 10000000, 10000)

    plt.plot(t, MIMout1[:, 0], 'b', label='LITm')
    plt.plot(t, MIMout1[:, 1], 'g', label='LITs')
    plt.plot(t, MIMout1[:, 2], 'r', label='MICr')
    plt.plot(t, MIMout1[:, 3], 'c', label='MICK')
    plt.plot(t, MIMout1[:, 4], 'm', label='SOMp')
    plt.plot(t, MIMout1[:, 5], 'y', label='SOMc')
    plt.plot(t, MIMout1[:, 6], 'k', label='SOMa')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.ylim([0, 5])
    plt.grid()
    plt.show()

    MIMpools = MIMout1[len(MIMout1)-1,:] 
    #print(MIMpools)
    
    Site = Site_names[0]
    MIMSOC = sum(MIMpools) * depth * 1e4 / 1e6
    MIMLIT = (MIMpools[0] + MIMpools[1]) * depth * 1e4 / 1e6
    MIMMIC = (MIMpools[2] + MIMpools[3]) * depth * 1e4 / 1e6
    SOMp = MIMpools[4] * depth * 1e4 / 1e6
    SOMc = MIMpools[5] * depth * 1e4 / 1e6
    SOMa = MIMpools[6] * depth * 1e4 / 1e6

    print("Site: " + str(Site))
    print("MIMSOC: " + str(MIMSOC))
    print("MIMLIT: " + str(MIMLIT))
    print("MIMMIC: " + str(MIMMIC))
    print("SOMp: " + str(SOMp))
    print("SOMc: " + str(SOMc))
    print("SOMa: " + str(SOMa))
    
    table = pd.DataFrame(columns=['Site','MIMSOC','MIMLIT','MIMMIC', 'SOMp', 'SOMc', 'SOMa'])
    for i, j in enumerate(Site_names):
        
        MIMrun = MIMout_set.to_numpy()[i][0] 
        MIMpools = MIMrun[len(MIMrun)-1,:] 
        
        table = table.append({
            "Site": j,
            "MIMSOC": sum(MIMpools) * depth * 1e4 / 1e6,
            "MIMLIT": (MIMpools[0] + MIMpools[1]) * depth * 1e4 / 1e6,
            "MIMMIC": (MIMpools[2] + MIMpools[3]) * depth * 1e4 / 1e6,
            "SOMp": MIMpools[4] * depth * 1e4 / 1e6,
            "SOMc": MIMpools[5] * depth * 1e4 / 1e6,
            "SOMa": MIMpools[6] * depth * 1e4 / 1e6
            }, ignore_index=True)
        
    print(table)
    
    
