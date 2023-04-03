######################
### Plots for Book ###
######################
### Packages
library(voi)
library(ggplot2)
library(dplyr)
library(R2jags)

## Run the model
source("04_analysis/02_baseline_model_output.R")
source("06_figs/01_plotting_functions.R")

## Baseline Cost-Effectiveness Formatting
# The output from the cost-effectiveness model should be formatted for 
# the voi package.
# Use BCEA package to create a cost-effectiveness object.
chemotherapy_output <- list(e = m_costs_effects[, "Effects", ],
                            c = m_costs_effects[, "Costs", ],
                            k = seq(0, 50000, length.out = 501))

# Using a bespoke analysis function - trial only updates odds ratio.
# Data generation function
OR_datagen_fn <- function(inputs, n = 500){
  p_side_effects_t1 <- inputs[, "p_side_effects_t1"]
  logor_side_effects <- inputs[, "logor_side_effects"]
  
  # Odds for side effects for treatment 1
  odds_side_effects_t1 <- p_side_effects_t1 / (1 - p_side_effects_t1)
  # Odds for side effects on treatment 2
  odds_side_effects_t2 <- odds_side_effects_t1 * exp(logor_side_effects)
  # Probability of side effects under treatment 2
  p_side_effects_t2    <- odds_side_effects_t2 / (1 + odds_side_effects_t2)
  
  # Data generation
  X1 <- rbinom(length(p_side_effects_t1), n, p_side_effects_t1)
  X2 <- rbinom(length(p_side_effects_t2), n, p_side_effects_t2)
  data_save <- data.frame(X1 = X1, X2 = X2)
  return(data_save)
}

# Analysis function based on JAGS
OR_analysis_fn <- function(data, args, pars){
  X1 <- data$X1
  X2 <- data$X2
  
  data_jags <- list(X1 = X1,
                    X2 = X2,
                    n = args$n,
                    n_side_effects = args$n_side_effects,
                    n_patients = args$n_patients,
                    logor_side_effects_mu = args$logor_side_effects_mu,
                    logor_side_effects_sd = args$logor_side_effects_sd)
  
  LogOR_trial <- function(){
    # Probability of side effects under treatment 1
    p_side_effects_t1 ~ dbeta(1 + n_side_effects, 
                              1 + n_patients - n_side_effects)
    
    # Log odds of side effects on treatment 2
    logor_side_effects ~ dnorm(logor_side_effects_mu, logor_side_effects_sd)
    # Odds of side effects on treatment 1
    odds_side_effects_t1 <- p_side_effects_t1 / (1 - p_side_effects_t1)
    # Odds for side effects on treatment 2
    odds_side_effects_t2 <- odds_side_effects_t1 * exp(logor_side_effects)
    
    # Probability of side effects under treatment 2
    p_side_effects_t2    <- odds_side_effects_t2 / (1 + odds_side_effects_t2)
    
    X1 ~ dbin(p_side_effects_t1, n)
    X2 ~ dbin(p_side_effects_t2, n)
  }
  
  filein <- file.path(tempdir(),fileext="datmodel.txt")
  R2OpenBUGS::write.model(LogOR_trial,filein)
  
  # Perform the MCMC simulation with OpenBUGS.
  # Close OpenBUGS once it has finished (if debug is set to TRUE)
  bugs.data <- jags(
    data =  data_jags,
    parameters.to.save = pars,
    model.file = filein, 
    n.chains = 1, 
    n.iter = args$n.iter, 
    n.thin = 1, 
    n.burnin = 250, progress.bar = "none") 
  
  
  return(data.frame(logor_side_effects = bugs.data$BUGSoutput$sims.matrix[, pars[1]]))
}

# EVSI calculation using the momemt matching method.
evsi_OR <- evsi(outputs = chemotherapy_output,
                inputs = m_params,
                pars = c("logor_side_effects"),
                pars_datagen = c("p_side_effects_t1", "logor_side_effects"),
                n = seq(50, 1500, by = 50),
                method = "mm",
                datagen_fn = OR_datagen_fn,
                model_fn = calculate_costs_effects,
                analysis_args = list(n_side_effects = n_side_effects,
                                     n_patients = n_patients,
                                     n = 500,
                                     logor_side_effects_mu = logor_side_effects_mu,
                                     logor_side_effects_sd = logor_side_effects_sd,
                                     n.iter = 5250),
                analysis_fn = OR_analysis_fn, 
                par_fn = generate_psa_parameters)

### Format for plotting
plotting_1 <- evsi.plot.adapt(chemotherapy_output, m_params, c("logor_side_effects"), evsi_OR, "gam")

pdf("06_figs/EVSIwtpN.pdf")
evsi.wtp.plot(plotting_1)
dev.off()

pdf("06_figs/EVSIwtp.pdf")
evsi.wtp.plot(plotting_1, N = 250)
dev.off()

pdf("06_figs/EVSIbyN.pdf")
evsi.ss.plot(plotting_1)
dev.off()

pdf("06_figs/ENBS.pdf")
evsi.enbs.plot(plotting_1, c(5e6, 1e7), c(28000,42000), 
               k = 20000, Pop = 46000, Time = 10)
dev.off()
optim.ss(plotting_1, mean(c(5e6, 1e7)), mean(c(28000,42000)), 
         k = 20000, Pop = 46000, Time = 10)

pdf("06_figs/coss.pdf")
coss(plotting_1, c(5e6, 1e7), c(28000,42000), Pop = 46000, Time = 5)
dev.off()

pdf("06_figs/ENBS-pop.pdf")
evsi.prob.plot(plotting_1, setup = c(5e6, 1e7), pp = c(28000,42000), k = 20000, 
               N = 398, Pop = c(0,60000), Time = c(0,10))
dev.off()

## Side Effects
# Using a bespoke analysis function - trial only updates odds ratio.
# Data generation function

OR_datagen_fn <- function(inputs, n = 500){
  p_side_effects_t1 <- inputs[, "p_side_effects_t1"]
  p_side_effects_t2 <- inputs[, "p_side_effects_t2"]
  X1 <- rbinom(length(p_side_effects_t1), n, p_side_effects_t1)
  X2 <- rbinom(length(p_side_effects_t2), n, p_side_effects_t2)
  # Create odds ratio as summary statistic
  OR <- (n - X2) / X2 / ((n - X1) / X1)
  data_save <- data.frame(OR = OR)
  return(data_save)
}

# EVSI calculation using GAM regression.
evsi_OR <- evsi(outputs = chemotherapy_output,
                inputs = m_params,
                pars = c("p_side_effects_t1", "p_side_effects_t2"),
                n = seq(50, 1500, by = 50),
                method = "gam",
                datagen_fn = OR_datagen_fn,
                par_fn = generate_psa_parameters)


plotting_2 <- evsi.plot.adapt(chemotherapy_output, m_params, c("logor_side_effects"), evsi_OR, "gam")

pdf("06_figs/EVSI_WTP_SE.pdf")
evsi.wtp.plot(plotting_2)
dev.off()

pop.adjust <- 46000 * (1 / (1 + 0.035)^3)
pdf("06_figs/ENBS_SE.pdf")
evsi.enbs.plot(plotting, c(1260000, 1400000), 2 * c(1560.55, 1600), 
               k = 20000, Pop = pop.adjust, Time = 7)
dev.off()

optim.ss(plotting, c(1260000, 1400000), 2 * c(1560.55, 1600), 
         k = 20000, Pop = pop.adjust, Time = 7)

## Full Study
## Moment Matching
# Data generation function
full_datagen_fn <- function(inputs, n = 500){
  p_side_effects_t1 <- inputs[, "p_side_effects_t1"]
  logor_side_effects <- inputs[, "logor_side_effects"]
  # Odds for side effects for treatment 1
  odds_side_effects_t1 <- p_side_effects_t1 / (1 - p_side_effects_t1)
  # Odds for side effects on treatment 2
  odds_side_effects_t2 <- odds_side_effects_t1 * exp(logor_side_effects)
  # Probability of side effects under treatment 2
  p_side_effects_t2 <- odds_side_effects_t2 / (1 + odds_side_effects_t2)
  p_hospitalised_total <- inputs[, "p_hospitalised_total"]
  p_died <- inputs[, "p_died"]
  lambda_home <- inputs[, "lambda_home"]
  lambda_hosp <- inputs[, "lambda_hosp"]
  rate_recover_hosp <- -log(1 -lambda_hosp )
  rate_recover_home <- -log(1 - lambda_home)
  
  X1 <- X2 <- X_hosp <- X_dead <- N_recover_home <-  N_recover_hospital <- vector("numeric", length = dim(inputs)[1])
  T_home <- T_hosp <- matrix(NA, nrow = dim(inputs)[1], ncol = 2 * n)
  for(i in 1:dim(inputs)[1]){
    # Simulate the number of patients with side effects
    X1[i] <- rbinom(1, n, p_side_effects_t1[i])
    X2[i] <- rbinom(1, n, p_side_effects_t2[i])
    
    # Simulate the number of patients hospitalised 
    X_hosp[i] <- rbinom(1, X1[i] + X2[i], p_hospitalised_total[i])
    # Simulate the number of patients die
    X_dead[i] <- rbinom(1, X_hosp[i], p_died[i])
    
    ## Simulate recovery times for patients
    N_recover_home[i] <- X1[i] + X2[i] - X_hosp[i]
    if(N_recover_home[i] > 0){
      T_home[i, 1:N_recover_home[i]] <- rexp(N_recover_home[i], rate_recover_home[i])
    }
    N_recover_hospital[i] <- X_hosp[i] - X_dead[i]
    if(N_recover_hospital[i] > 0){
      T_hosp[i, 1:N_recover_hospital[i]] <- rexp(N_recover_hospital[i], rate_recover_hosp[i])
    }
    
  }
  
  data_save_dat <- data.frame(cbind(X1 = X1, X2 = X2, 
                                    X_hosp = X_hosp, X_dead = X_dead,
                                    N_recover_home = N_recover_home,
                                    N_recover_hospital = N_recover_hospital,
                                    T_home = T_home, T_hosp = T_hosp))
  return(data_save_dat)
}

## Regression Based
full_datagen_fn_RB <- function(inputs, n = 500){
  # Load the data
  p_side_effects_t1 <- inputs[, "p_side_effects_t1"]
  logor_side_effects <- inputs[, "logor_side_effects"]
  # Odds for side effects for treatment 1
  odds_side_effects_t1 <- p_side_effects_t1 / (1 - p_side_effects_t1)
  # Odds for side effects on treatment 2
  odds_side_effects_t2 <- odds_side_effects_t1 * exp(logor_side_effects)
  # Probability of side effects under treatment 2
  p_side_effects_t2 <- odds_side_effects_t2 / (1 + odds_side_effects_t2)
  p_hospitalised_total <- inputs[, "p_hospitalised_total"]
  p_died <- inputs[, "p_died"]
  lambda_home <- inputs[, "lambda_home"]
  lambda_hosp <- inputs[, "lambda_hosp"]
  rate_recover_hosp <- -log(1 -lambda_hosp )
  rate_recover_home <- -log(1 - lambda_home)
  
  X1 <- X2 <- X_hosp <- X_dead <- N_recover_home <- 
    N_recover_hospital <- vector("numeric", length = dim(inputs)[1])
  T_home <- T_hosp <- matrix(NA, nrow = dim(inputs)[1], ncol = 2 * n)
  for(i in 1:dim(inputs)[1]){
    # Simulate the number of patients with side effects
    X1[i] <- rbinom(1, n, p_side_effects_t1[i])
    X2[i] <- rbinom(1, n, p_side_effects_t2[i])
    
    # Simulate the number of patients hospitalised 
    X_hosp[i] <- rbinom(1, X1[i] + X2[i], p_hospitalised_total[i])
    # Simulate the number of patients die
    X_dead[i] <- rbinom(1, X_hosp[i], p_died[i])
    
    ## Simulate recovery times for patients
    N_recover_home[i] <- X1[i] + X2[i] - X_hosp[i]
    if(N_recover_home[i] > 0){
      T_home[i, 1:N_recover_home[i]] <- rexp(N_recover_home[i], rate_recover_home[i])
    }
    N_recover_hospital[i] <- X_hosp[i] - X_dead[i]
    if(N_recover_hospital[i] > 0){
      T_hosp[i, 1:N_recover_hospital[i]] <- rexp(N_recover_hospital[i], rate_recover_hosp[i])
    }
  }
  
  OR <- (n - X2) / X2 / ((n - X1) / X1)
  p_hosp <- X_hosp / (X1 + X2)
  p_dead <- X_dead / X_hosp
  T_home_sum <- rowSums(T_home, na.rm = TRUE)
  T_hosp_sum <- rowSums(T_hosp, na.rm = TRUE)
  
  data_save_dat <- data.frame(cbind(OR = OR, 
                                    p_hosp = p_hosp, p_dead = p_dead,
                                    T_home_sum = T_home_sum, 
                                    T_hosp_sum = T_hosp_sum))
  return(data_save_dat)
}

# EVSI calculation using MARS regression - large number of parameters.
evsi_OR_allout <- evsi(outputs = chemotherapy_output,
                       inputs = m_params,
                       pars = c("p_side_effects_t1", "logor_side_effects",
                                "p_hospitalised_total", "p_died",
                                "lambda_home", "lambda_hosp"),
                       n = seq(50, 1500, by = 50),
                       method = "earth",
                       datagen_fn = full_datagen_fn_RB,
                       par_fn = generate_psa_parameters)

plotting_3 <- evsi.plot.adapt(chemotherapy_output, m_params, c("logor_side_effects",
                                                             "p_hospitalised_total", "p_died",
                                                             "lambda_home", "lambda_hosp"), 
                            evsi_OR_allout, "earth")

pop.adjust <- 46000 * (1 / (1 + 0.035)^3)
pdf("06_figs/ENBS_SEFU.pdf")
evsi.enbs.plot(plotting_3, c(1260000, 1400000), 2 * c(1560.55, 1600), 
               k = 20000, Pop = pop.adjust, Time = 7)
dev.off()

optim.ss(plotting_3, c(1260000, 1400000), 2 * c(1560.55, 1600), 
         k = 20000, Pop = pop.adjust, Time = 7)

pdf("06_figs/COSS_SEFU.pdf")
coss(plotting_3, c(1260000, 1400000), 2 * c(1560.55, 1600), 
     Pop = pop.adjust, Time = 7)
dev.off()

pdf("06_figs/prob_SEFU.pdf")
evsi.prob.plot(plotting_3, setup = c(1260000, 1400000), pp = 2 * c(1560.55, 1600), k = 20000, 
               N = 1080, Pop = c(0,60000), Time = c(0,15))
dev.off()

# Analysis function based on JAGS
full_analysis_fn <- function(data, args, pars){
  ## Format Data - Adjust for 0 recovery times
  T_home <- NA
  if(data$N_recover_home > 0){
    T_home <- as.numeric(as.matrix(data[, (1:data$N_recover_home) + 6]))
  } 
  
  T_hosp <- NA
  if(data$N_recover_hospital > 0){
    T_hosp <- as.vector(as.matrix(data[, 
                                       (6 + 2 * args$n) + (1:data$N_recover_hospital)]))
  } 
  
  # Create the data list for JAGS
  data_jags <- list(X1 = data$X1,
                    X2 = data$X2,
                    X_hosp = data$X_hosp,
                    X_dead = data$X_dead,
                    T_home = T_home,
                    T_hosp = T_hosp,
                    N_recover_home = ifelse(data$N_recover_home > 0,
                                            data$N_recover_home, 
                                            1),
                    N_recover_hosp = ifelse(data$N_recover_hospital > 0,
                                            data$N_recover_hospital,
                                            1),
                    n = args$n,
                    n_side_effects = args$n_side_effects,
                    n_patients = args$n_patients,
                    logor_side_effects_mu = args$logor_side_effects_mu,
                    logor_side_effects_sd = args$logor_side_effects_sd,
                    p_recovery_home_alpha = betaPar(args$p_recovery_home_mu, 
                                                    args$p_recovery_home_sd)$alpha,
                    p_recovery_home_beta = betaPar(args$p_recovery_home_mu, 
                                                   args$p_recovery_home_sd)$beta,
                    p_recovery_hosp_alpha = betaPar(args$p_recovery_hosp_mu, 
                                                    args$p_recovery_hosp_sd)$alpha,
                    p_recovery_hosp_beta = betaPar(args$p_recovery_hosp_mu,
                                                   args$p_recovery_hosp_sd)$beta,
                    n_died = args$n_died,
                    n_hospitalised = args$n_hospitalised)
  
  LogOR_addoutcomes_trial <- function(){
    
    ## Models for the data
    X1 ~ dbin(p_side_effects_t1, n)
    X2 ~ dbin(p_side_effects_t2, n)
    
    X_hosp ~ dbinom(p_hospitalised_total, X1 + X2)
    X_dead ~ dbin(p_died, X_hosp)
    
    rate_recover_home <- -log(1 - lambda_home)
    rate_recover_hosp <- -log(1 - lambda_hosp)
    
    for(i in 1:N_recover_home){
      T_home[i] ~ dexp(rate_recover_home)
    }
    
    for(i in 1:N_recover_hosp){
      T_hosp[i] ~ dexp(rate_recover_hosp)
    }
    
    # Probability of side effects under treatment 1
    p_side_effects_t1 ~ dbeta(1 + n_side_effects, 
                              1 + n_patients - n_side_effects)
    
    # Log odds of side effects on treatment 2
    logor_side_effects ~ dnorm(logor_side_effects_mu, logor_side_effects_sd)
    # Odds of side effects on treatment 1
    odds_side_effects_t1 <- p_side_effects_t1 / (1 - p_side_effects_t1)
    # Odds for side effects on treatment 2
    odds_side_effects_t2 <- odds_side_effects_t1 * exp(logor_side_effects)
    
    # Probability of side effects under treatment 2
    p_side_effects_t2    <- odds_side_effects_t2 / (1 + odds_side_effects_t2)
    
    ## Variables to define transition probabilities
    # Probability that a patient is hospitalised over the time horizon
    p_hospitalised_total ~ dbeta(1 + n_hospitalised, 
                                 1 + n_side_effects - n_hospitalised)
    # Probability that a patient dies over the time horizon given they were 
    # hospitalised
    p_died ~ dbeta(1 + n_died, 1 + n_hospitalised - n_died)
    # Lambda_home: Conditional probability that a patient recovers considering 
    # that they are not hospitalised
    lambda_home ~ dbeta(p_recovery_home_alpha, p_recovery_home_beta)
    # Lambda_hosp: Conditional probability that a patient recovers considering 
    # that they do not die
    lambda_hosp ~ dbeta(p_recovery_hosp_alpha, p_recovery_hosp_beta)
  }
  
  filein <- file.path(tempdir(),fileext="datmodel.txt")
  R2OpenBUGS::write.model(LogOR_addoutcomes_trial,filein)
  
  # Perform the MCMC simulation with OpenBUGS.
  # Close OpenBUGS once it has finished (if debug is set to TRUE)
  bugs.data <- jags(
    data =  data_jags,
    parameters.to.save = pars,
    model.file = filein, 
    n.chains = 1, 
    n.iter = args$n.iter, 
    n.thin = 1, 
    n.burnin = 250, progress.bar = "none") 
  
  return(data.frame(logor_side_effects = bugs.data$BUGSoutput$sims.matrix[, "logor_side_effects"],
                    p_hospitalised_total= bugs.data$BUGSoutput$sims.matrix[, "p_hospitalised_total"],
                    p_died = bugs.data$BUGSoutput$sims.matrix[, "p_died"],
                    lambda_home = bugs.data$BUGSoutput$sims.matrix[, "lambda_home"],
                    lambda_hosp = bugs.data$BUGSoutput$sims.matrix[, "lambda_hosp"]))
}

# EVSI calculation using the momemt matching method.
evsi_OR_allout_MM <- evsi(outputs = chemotherapy_output,
                          inputs = m_params,
                          pars = c("logor_side_effects",
                                   "p_hospitalised_total", "p_died",
                                   "lambda_home", "lambda_hosp"),
                          pars_datagen = c("p_side_effects_t1", 
                                           "logor_side_effects",
                                           "p_hospitalised_total", "p_died",
                                           "lambda_home", "lambda_hosp"),
                          n = seq(50, 1500, by = 50),
                          method = "mm",
                          datagen_fn = full_datagen_fn,
                          model_fn = calculate_costs_effects,
                          analysis_args = list(n_side_effects = n_side_effects,
                                               n_patients = n_patients,
                                               n = 50,
                                               logor_side_effects_mu = logor_side_effects_mu,
                                               logor_side_effects_sd = logor_side_effects_sd,
                                               betaPar = betaPar,
                                               p_recovery_home_mu = p_recovery_home_mu,
                                               p_recovery_home_sd = p_recovery_home_sd,
                                               p_recovery_hosp_mu = p_recovery_hosp_mu,
                                               p_recovery_hosp_sd = p_recovery_hosp_sd,
                                               n.iter = 5250,
                                               n_died = n_died,
                                               n_hospitalised = n_hospitalised,
                                               p_side_effects_t1 = m_params$p_side_effects_t1),
                          analysis_fn = full_analysis_fn, 
                          par_fn = generate_psa_parameters,
                          npreg_method = "earth")

plotting_4 <- evsi.plot.adapt(chemotherapy_output, m_params, c("logor_side_effects",
                                                             "p_hospitalised_total", "p_died",
                                                             "lambda_home", "lambda_hosp"), 
                            evsi_OR_allout_MM, "earth")

pop.adjust <- 46000 * (1 / (1 + 0.035)^3)

pdf("06_figs/ENBS_SEFU_MM.pdf")
evsi.enbs.plot(plotting_4, c(1260000, 1400000), 2 * c(1560.55, 1600), 
               k = 20000, Pop = pop.adjust, Time = 7)
dev.off()

optim.ss(plotting_4, c(1260000, 1400000), 2 * c(1560.55, 1600), 
         k = 20000, Pop = pop.adjust, Time = 7)

pdf("06_figs/COSS_SEFU_MM.pdf")
coss(plotting_4, c(1260000, 1400000), 2 * c(1560.55, 1600), 
     Pop = pop.adjust, Time = 7)
dev.off()

pdf("06_figs/prob_SEFU_MM.pdf")
evsi.prob.plot(plotting_4, setup = c(1260000, 1400000), pp = 2 * c(1560.55, 1600), k = 20000, 
               N = 1020, Pop = c(0,60000), Time = c(0,15))
dev.off()


### Utilities ###
# Data generation function
utility_datagen_fn <- function(inputs, n = 20, 
                               sd_recovery_fun = function(){return(runif(1, 0.000001, 0.15))},
                               sd_home_care_fun = function(){return(runif(1, 0.00001, 0.6))},
                               sd_hospital_fun = function(){return(runif(1, 0.00001, 0.4))}
){
  # Load the data
  X_home_care <- X_hospital <- X_recovery <- matrix(NA, nrow = dim(inputs)[1], ncol = n[1])
  X_recovery_mean1 <- X_recovery_mean2 <-  X_home_care_mean1 <- X_home_care_mean2 <-
    X_hospital_mean1 <- X_hospital_mean2 <- vector("numeric", dim(inputs)[1])
  for(i in 1:dim(inputs)[1]){
    set.seed(123 + i)
    m_recovery <- inputs[i, "u_recovery"]
    m_home_care <- inputs[i, "u_home_care"]
    m_hospital <- inputs[i, "u_hospital"]
    sd_recovery <- sd_recovery_fun()
    sd_home_care <- sd_home_care_fun()
    sd_hospital <- sd_hospital_fun()
    
    # Simulate the costs 
    X_recovery[i, ] <- truncnorm::rtruncnorm(n[1], mean = m_recovery, sd = sd_recovery,
                                             a = -Inf, b = 1)
    X_home_care[i, ] <- truncnorm::rtruncnorm(n[1], mean = m_home_care, sd = sd_home_care,
                                              a = -Inf, b = 1)
    X_hospital[i, ] <- truncnorm::rtruncnorm(n[1], mean = m_hospital, sd = sd_hospital,
                                             a = -Inf, b = 1)  }
  
  data_save_dat <- data.frame(cbind(X_recovery = X_recovery,
                                    X_home_care = X_home_care,
                                    X_hospital = X_hospital))
  return(data_save_dat)
}


# Analysis function based on JAGS
utility_analysis_fn <- function(data, args, pars){
  # Create the data list for JAGS
  data_jags <- list(X_recovery = as.vector(as.matrix(data[, (1:args$n)])),
                    X_home_care = as.vector(as.matrix(data[, args$n + (1:args$n)])),
                    X_hospital = as.vector(as.matrix(data[, 2*args$n + (1:args$n)])),
                    n = args$n,
                    alpha_recovery = args$betaPar(
                      args$u_recovery_mu,
                      args$u_recovery_sd
                    )$alpha,
                    beta_recovery = args$betaPar(
                      args$u_recovery_mu,
                      args$u_recovery_sd
                    )$beta,
                    alpha_home_care = args$betaPar(
                      args$u_home_care_mu,
                      args$u_home_care_sd
                    )$alpha,
                    beta_home_care = args$betaPar(
                      args$u_home_care_mu,
                      args$u_home_care_sd
                    )$beta,
                    alpha_hospital = args$betaPar(
                      args$u_hospital_mu,
                      args$u_hospital_sd
                    )$alpha,
                    beta_hospital = args$betaPar(
                      args$u_hospital_mu,
                      args$u_hospital_sd
                    )$beta)
  
  trial <- function(){
    for(i in 1:n){
      X_recovery[i] ~ dnorm(u_recovery, tau_recovery);T(, 1)
      X_home_care[i] ~ dnorm(u_home_care, tau_home_care);T(, 1)
      X_hospital[i] ~ dnorm(u_hospital, tau_hospital);T(, 1)
    }
    
    u_recovery ~ dbeta(alpha_recovery, beta_recovery)
    u_home_care ~ dbeta(alpha_home_care, beta_home_care)
    u_hospital ~ dbeta(alpha_hospital, beta_hospital)
    
    sd_recovery ~ dunif(0.000001, 0.15)
    sd_home_care ~ dunif(0.00001, 0.6)
    sd_hospital ~ dunif(0.00001, 0.4)
    
    tau_recovery <- 1 / sd_recovery ^ 2
    tau_home_care <- 1 / sd_home_care ^ 2
    tau_hospital <- 1 / sd_hospital ^ 2
  }
  
  filein <- file.path(tempdir(),fileext="datmodel.txt")
  R2OpenBUGS::write.model(trial,filein)
  
  # Perform the MCMC simulation with OpenBUGS.
  # Close OpenBUGS once it has finished (if debug is set to TRUE)
  bugs.data <- jags(
    data =  data_jags,
    parameters.to.save = pars,
    model.file = filein, 
    n.chains = 1, 
    n.iter = args$n.iter, 
    n.thin = 1, 
    n.burnin = 250, progress.bar = "none") 
  
  return(data.frame(u_recovery = bugs.data$BUGSoutput$sims.matrix[, "u_recovery"],
                    u_home_care = bugs.data$BUGSoutput$sims.matrix[, "u_home_care"],
                    u_hospital = bugs.data$BUGSoutput$sims.matrix[, "u_hospital"]))
}

# EVSI calculation using the momemt matching method.
evsi_utility <- evsi(outputs = chemotherapy_output,
                     inputs = m_params,
                     pars = c("u_recovery", "u_home_care", "u_hospital"),
                     n = seq(20, 300, by = 10),
                     method = "mm",
                     datagen_fn = utility_datagen_fn,
                     model_fn = calculate_costs_effects,
                     analysis_args = list(n = 20,
                                          betaPar = betaPar,
                                          u_recovery_mu = u_recovery_mu,
                                          u_recovery_sd = u_recovery_sd,
                                          u_home_care_mu = u_home_care_mu,
                                          u_home_care_sd = u_home_care_sd,
                                          u_hospital_mu = u_hospital_mu,
                                          u_hospital_sd = u_hospital_sd,
                                          n.iter = 5000),
                     analysis_fn = utility_analysis_fn, 
                     par_fn = generate_psa_parameters,
                     Q = 50)

plotting_6 <- evsi.plot.adapt(chemotherapy_output, m_params, c("u_recovery", "u_home_care", "u_hospital"), 
                            evsi_utility, "gam")

pdf("06_figs/EVSI_WTP_U_MM.pdf")
evsi.wtp.plot(plotting_6)
dev.off()
pop.adjust <- 46000 * (1 / (1 + 0.035)^2)
pdf("06_figs/ENBS_U_MM.pdf")
evsi.enbs.plot(plotting_6, c(90000, 95000), 3 * c(370-25, 370), 
               k = 20000, Pop = pop.adjust, Time = 8)
dev.off()


optim.ss(plotting_6, c(90000, 95000), 3 * c(370-25, 370), 
         k = 20000, Pop = pop.adjust, Time = 8)

pdf("06_figs/COSS_U_MM.pdf")
coss(plotting_6, c(90000, 95000), 3 * c(370-25, 370), 
     Pop = pop.adjust, Time = 8)
dev.off()


pdf("06_figs/prob_U_MM.pdf")
evsi.prob.plot(plotting_6, setup = c(90000, 95000), pp = 3 * c(370-25, 370), k = 20000, 
               N = 300, Pop = c(0,60000), Time = c(0,15))
dev.off()



#### STUDY 3: Long-term Survival ####
longterm_datagen_fn <- function(inputs, n = 46000){
  rate_longterm <- inputs[, "rate_longterm"]
  sum_of_surv <- rgamma(dim(inputs)[1], shape = 2 * n, scale = n / rate_longterm)
  return(data.frame(surv_sum = sum_of_surv))
}

# EVSI calculation using GAM regression.
evsi_longterm <- evsi(outputs = chemotherapy_output,
                      inputs = m_params,
                      pars = c("rate_longterm"),
                      n = 46000,
                      method = "gam",
                      datagen_fn = longterm_datagen_fn,
                      par_fn = generate_psa_parameters)

plotting_7 <- evsi.plot.adapt(chemotherapy_output, m_params, c("rate_longterm"), 
                              evsi_longterm, "gam")
pdf("06_figs/EVSI_LT.pdf")
evsi.wtp.plot(plotting_7)
dev.off()

evsi_longterm %>% filter(k %in% c(20000, 25000, 30000))


pop.adjust <- 46000 * (1 / (1 + 0.035)^0.5)
ENBS.fun(evsi_longterm %>% filter(k %in% c(20000, 25000, 30000)),
         cost=c(60000, 60000), Pop = pop.adjust, Time = 9.5, Dis = 0.035
)

optim.ss(plotting, c(60000, 60000), c(0, 0), 
         k = 30000, Pop = pop.adjust, Time = 9.5)

pdf("06_figs/EVSI_LT_Prob.pdf")
evsi.prob.plot(plotting, c(60000, 60000), c(0, 0),
               k = 20000, Pop = c(0, pop.adjust * 2), Time = c(0, 15))
dev.off()