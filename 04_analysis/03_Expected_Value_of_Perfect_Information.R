###########################################################
### Expected Value of Perfect Information Analysis      ###
###########################################################

### Packages
library(voi)

## Run the model
source("04_analysis/01_model_run.R")

## Expected Value of Perfect Information - single WTP
# Specify willingness to pay
wtp_fix = 20000
# Extract net benefit for this particular WTP
nb <- m_net_benefit[ , , wtp_seq == wtp_fix]
# Calculate EVPI from net benefit
evpi(nb)

## Baseline Cost-Effectiveness Formatting
# The output from the cost-effectiveness model should be formatted for 
# the voi package.
# Use BCEA package to create a cost-effectiveness object.
chemotherapy_output <- list(e = m_costs_effects[, "Effects", ],
                            c = m_costs_effects[, "Costs", ],
                            k = seq(0, 50000, length.out = 501))

## Expected Value of Perfect Information
# Calculate
EVPI <- evpi(chemotherapy_output)
# WTP = 20000
EVPI$evpi[EVPI$k == wtp_fix]

# Plot
pdf("06_figs/EVPI.pdf")
plot(EVPI,
     xlab = "Willingness-to-Pay",
     ylab = "EVPI",
     main = "Expected Value of Perfect Information",
     type = "l")
dev.off()

ELC <- array(NA, dim = c(length(chemotherapy_output$k), 2))
for(i in 1:length(chemotherapy_output$k)){
  NB <- chemotherapy_output$e * chemotherapy_output$k[i] - chemotherapy_output$c
  loss <- NB - apply(NB, 1, max)
  
ELC[i, ] <- - apply(loss, 2, mean)
}

full_ELC <- c(ELC[,1], ELC[, 2])
int_type <- as.factor(c(rep("Standard of Care", 501), 
                        rep("Novel Intervention", 501)))
int_type <- relevel(int_type, ref = "Standard of Care")
full_wtp <- c(chemotherapy_output$k, chemotherapy_output$k)
ELC_plot <- data.frame(ELC = full_ELC, Intervention = int_type, 
                       k = full_wtp)

library(ggplot2)
ggplot(data = ELC_plot, aes(x = k, y = ELC, color = Intervention, linetype = Intervention)) + 
  geom_line(size = 1.1) +
  theme_bw() +
  xlab("Willingness to Pay") + 
  ylab("Expected Loss")
