################################################################################
#### Model Inputs from Assumptions for the Chemotherapy Model
################################################################################

## Time Horizon
time_horizon <- 50

## Parameters for the PSA distribution of the risk reduction of side effects
rr_side_effects_mu <- 0.67
rr_side_effects_sd <- 0.2

## Parameters for the PSA distribution of the recovery time for patients who are
## not hospitalised
t_recovery_home_mu <- 5.3
t_recovery_home_sd <- sqrt(1.7)

## Parameters for the PSA distribution of the recovery time for patients who are
## hospitalised
t_recovery_hosp_mu <- 7.6
t_recovery_hosp_sd <- sqrt(3.5)

## Parameters for the PSA distribution of the costs of treating patients at home
c_home_care_mu <- 830
c_home_care_sd <- sqrt(150)

## Parameters for the PSA distribution of the costs of treating patients in
## hospital
c_hospital_mu <- 2400
c_hospital_sd <- sqrt(1880)

## Parameters for the PSA distribution of the one-off cost of death
c_death_mu <- 1710
c_death_sd <- sqrt(760)


## Parameters for the PSA distribution of the utility for recovered patients
u_recovery_mu <- 0.98
u_recovery_sd <- sqrt(0.001)

## Parameters for the PSA distribution of the utility of patients who are treated
## at home. 
u_home_care_mu <- 0.7
u_home_care_sd <- sqrt(0.02)

## Parameters for the PSA distribution of the utility of treating patients in 
## hospital
u_hospital_mu <- 0.3
u_hospital_sd <- sqrt(0.03)

## Drug costs
c_treatment_1 <- 120
c_treatment_2 <- 10300
