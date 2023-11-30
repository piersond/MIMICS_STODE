# Plot Vmax curve
#----------------------------------------------------------------------
#----------------------------------------------------------------------

# Simulate data
#----------------------------------------------------------------------
#aV      <- rep(0.000000075, 6)
#Vslope_MOD = c(0.0171 , 0.548)
#Vint_MOD = c(-0.00429, 1.165)

aV = aV      <- rep(0.000008, 6) 
Vslope_MOD = c(1.069356e-02, 5.797809e-01)   # MSBio new
Vint_MOD = c( -3.169280e-04, 1.050568)        # MSBio new

df <- NULL
groups = c("DEFAULT", "MAT = 10", "MAT = 15", "MAT = 20", "MAT = 25")
MATs = c(0, 10, 15, 20, 25)

for(j in 1:5) {
  for(i in 5:35) {
    
    TSOI = i
    MAT = MATs[j]
    
    Vslope  <- rep(0.063, 6)
    Vint    <- rep(5.47, 6)
    
    if(MAT != 0) {
      #Vslope =Vslope * (Vslope_MOD[1]*MAT + Vslope_MOD[2])
      #Vint = Vint * (Vint_MOD[1]*MAT + Vint_MOD[2])
      
      #DEBUG
      Vslope =Vslope + (MAT*0.00104) #  0.0024 vs 0.003
      Vint = Vint - (MAT*0.0228) #(Vint_MOD[1]*MAT + Vint_MOD[2])
    }
    
    
    Vmax     <- exp(TSOI * Vslope + Vint) * aV 
    
    output <- data.frame(group = groups[j],
                         iter = i, 
                         TSOI = TSOI, 
                         MAT = MAT, 
                         Vslope = Vslope, 
                         Vint = Vint, 
                         Vmax = Vmax)
    
    df <- rbind(df, output)
  }
}

# Plot
#---------------------------------------------------------------------

library(dplyr)
library(ggplot2)

cols <- c("DEFAULT" = "black", 
          "MAT = 10" = "#FA6F01", 
          "MAT = 15" = "#F55301", 
          "MAT = 20" = "#F03801", 
          "MAT = 25" = "#EB1C01")

ggplot(df %>% filter(group != "DEFAULT"), aes(x=TSOI, y=Vmax, color=group)) + geom_line() +
  theme_bw() +
  scale_colour_manual(values = cols) +
  geom_line(data = df %>% filter(group == "DEFAULT"), aes(x=TSOI, y=Vmax), size=1) +
  labs(color="") +
  ylab("Vmax [1]") +
  xlim(10,30) +
  #ylim(0,0.02) +
  geom_vline(xintercept=c(15,25), alpha = 0.5, linetype="dashed", color="dark blue")


# Get MAT proportional change to Vmax 
#-------------------------------------------

cold <- df %>% filter(TSOI == 15) %>% unique()
warm <- df %>% filter(TSOI == 25) %>% unique()

cold[2,7]/cold[5,7] # Looking for 1.1158
warm[2,7]/warm[5,7] # Looking for 0.95
