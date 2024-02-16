library(dplyr)
library(ggplot2)
library(DT)
library(purrr)
library(Metrics)

## Set working drive
setwd("C:/GitHub/MIMICS_STODE")

#----------------------------------------
# Load necessary model components
#----------------------------------------

# Load reverse Michaelis-Menton MIMICS function
source("functions/RXEQ.R")

# Load function to calculate model input variables
source("functions/calc_Tpars.R")

# Set MIMICS parameters
source("parameters/MIMICS_parameters_MSBio_Incubation.R")  #--> Default "Sandbox" - Wieder et al. 2015

# Bring in MIMICS incubation simulation function
source("functions/MIMICS_sim_incubation.R")

#-----------------------------------------

# Bring in MSBio forcing data
MSB_data <- read.csv("Example_simulations/Data/MSBio_MIM_forcings.csv")
forcing_df <- MSB_data

#---------------------------------------------------------
# Run single incubation
#---------------------------------------------------------
#MIMinc_out <- MIMICS_INCUBATION(forcing_df[1,], days=105, step="hourly", output_type = 1)


#---------------------------------------------------------
#>> Run all rows in forcing dataset
#---------------------------------------------------------
forcing_df <- MSB_data
MIMrun <- forcing_df %>% split(1:nrow(forcing_df)) %>% 
            map(~MIMICS_INCUBATION(df=., days=105, step="daily", output_type = 2)) %>% 
            bind_rows() 

MC_MIMICS <- MIMrun %>% left_join(forcing_df %>% select(-SITE), by="ID")


##########################################################
# Plot lab vs MIMICS incubation cumulative respiration
##########################################################

plot_df <- MC_MIMICS #%>% filter(moisture.trt != 20)

MC_CO2Cp_cor <- round(cor(plot_df$CO2_prop_totC, plot_df$CO2C_prop), 4)
MC_CO2Cp_rmse <- round(rmse(plot_df$CO2_prop_totC, plot_df$CO2C_prop), 4)
MC_CO2Cp_fit <- lm(plot_df$CO2_prop_totC ~ plot_df$CO2C_prop)

# Extract coefficient values
intercept <- coef(MC_CO2Cp_fit)[1]
slope <- coef(MC_CO2Cp_fit)[2]

# Print the equation
line_eqn = paste0("y = ", round(slope, 6), "x + ", round(intercept, 4))

# Plot the relationship
plt <- ggplot(plot_df,
              aes(y=CO2_prop_totC, x=CO2C_prop, color=factor(moisture.trt))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth=1, alpha=0.6) +
  geom_point(size=3) +
  #geom_smooth(method = "lm", se = FALSE, alpha=0.8) +
  geom_line(stat="smooth",method = "lm", formula = y ~ x,
            linewidth = 1,
            #linetype ="dashed",
            alpha = 0.3) +
  ylab("MIMICS (CO2C / total C)") +
  xlab("Laboratory (CO2C / total C)") +
  labs(#title = "Calibrated kinetic terms & fW curve\nCumulative respiration after 105 day incubation
  #      ",#\nUsing site MAT and moisture controls on decomposition rate", #Default MIMICS parameters",
        subtitle = paste0("r^2 = ", MC_CO2Cp_cor, "   ", "RMSE = ", MC_CO2Cp_rmse),
        caption = "Dashed line = 1:1",
        color = "Moisture") +
  #ylim(0, 0.13) + xlim(0, 0.13) +
  theme_minimal()
plt
