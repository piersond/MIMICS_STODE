def RXEQ(y, t, I, VMAX, KM, CUE, fPHYS, fCHEM, fAVAI, FI, tau, LITmin, SOMmin, MICtrn, desorb, DEsorb, OXIDAT, KO):
  
    # Prevent negative pool values
    LIT_1 = y[0] if y[0] > 0 else 0.000001
    LIT_2 = y[1] if y[1] > 0 else 0.000001
    MIC_1 = y[2] if y[2] > 0 else 0.000001
    MIC_2 = y[3] if y[3] > 0 else 0.000001
    SOM_1 = y[4] if y[4] > 0 else 0.000001
    SOM_2 = y[5] if y[5] > 0 else 0.000001
    SOM_3 = y[6] if y[6] > 0 else 0.000001
    
    # Flows to and from MIC_1
    LITmin[0] = MIC_1 * VMAX[0] * LIT_1 / (KM[0] + MIC_1)   #MIC_1 decomp of MET lit
    LITmin[1] = MIC_1 * VMAX[1] * LIT_2 / (KM[1] + MIC_1)   #MIC_1 decomp of STRUC lit
    MICtrn[0] = MIC_1 * tau[0]  * fPHYS[0]                  #MIC_1 turnover to PHYSICAL SOM 
    MICtrn[1] = MIC_1 * tau[0]  * fCHEM[0]                  #MIC_1 turnover to CHEMICAL SOM  
    MICtrn[2] = MIC_1 * tau[0]  * fAVAI[0]                  #MIC_1 turnover to AVAILABLE SOM  
    SOMmin[0] = MIC_1 * VMAX[2] * SOM_3 / (KM[2] + MIC_1)   #decomp of SOMa by MIC_1

    # Flows to and from MIC_2
    LITmin[2] = MIC_2 * VMAX[3] * LIT_1 / (KM[3] + MIC_2)   #decomp of MET litter
    LITmin[3] = MIC_2 * VMAX[4] * LIT_2 / (KM[4] + MIC_2)   #decomp of SRUCTURAL litter
    MICtrn[3] = MIC_2 * tau[1]  * fPHYS[1]                  #MIC_2 turnover to PHYSICAL  SOM 
    MICtrn[4] = MIC_2 * tau[1]  * fCHEM[1]                  #MIC_2 turnover to CHEMICAL  SOM  
    MICtrn[5] = MIC_2 * tau[1]  * fAVAI[1]                  #MIC_2 turnover to AVAILABLE SOM  
    SOMmin[1] = MIC_2 * VMAX[5] * SOM_3 / (KM[5] + MIC_2)   #decomp of SOMa by MIC_2

    DEsorb    = SOM_1 * desorb  #* (MIC_1 + MIC_2)  	#desorbtion of PHYS to AVAIL (function of fCLAY)
    OXIDAT    = ((MIC_1 * VMAX[1] * SOM_2 / (KO[0]*KM[1] + MIC_1)) +
                (MIC_2 * VMAX[4] * SOM_2 / (KO[1]*KM[4] + MIC_2)) )  #oxidation of C to A

    dLIT_1 = I[0]*(1-FI[0]) - LITmin[0] - LITmin[2]

    dMIC_1 = CUE[0]*(LITmin[0]+ SOMmin[0]) + CUE[1]*(LITmin[1]) - (MICtrn[[0]] + MICtrn[[1]] + MICtrn[[2]])#sum(MICtrn[1:3])
    dSOM_1 = I[0]*FI[0] + MICtrn[0] + MICtrn[3]- DEsorb

    dLIT_2 = I[1] * (1-FI[1]) - LITmin[1] - LITmin[3]
    dMIC_2 = CUE[2]*(LITmin[2]+ SOMmin[1]) + CUE[3]*(LITmin[3]) - (MICtrn[[3]] + MICtrn[[4]] + MICtrn[[5]])#sum(MICtrn[4:6])
    dSOM_2 = I[1]*FI[1] + MICtrn[1] + MICtrn[4] - OXIDAT

    dSOM_3 = MICtrn[2] + MICtrn[5] + DEsorb + OXIDAT - SOMmin[0] - SOMmin[1]

    return(dLIT_1, dLIT_2, dMIC_1, dMIC_2, dSOM_1, dSOM_2, dSOM_3)

