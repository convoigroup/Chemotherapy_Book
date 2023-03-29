#######################################################################
### Expected Value of Sample Information Analysis - Moment Matching ###
#######################################################################
### Packages
library(voi)
library(ggplot2)
library(dplyr)
library(R2jags)

## Run the model
source("04_analysis/02_baseline_model_output.R")

## Baseline Cost-Effectiveness Formatting
# The output from the cost-effectiveness model should be formatted for 
# the voi package.
# Use BCEA package to create a cost-effectiveness object.
chemotherapy_output <- list(e = m_costs_effects[, "Effects", ],
                            c = m_costs_effects[, "Costs", ],
                            k = seq(0, 50000, length.out = 501))

## EVSI Calculations
#### STUDY 1: Randomised Trial for Log-Odds ratio ####
## Using the default trials in the voi package

# A randomised trial for binary outcomes requires a beta prior distribution
# Beta prior for standard care is set using the number of events
beta_params_t1 <- c(1 + n_side_effects, 
                    1 + n_patients - n_side_effects)
# Beta prior for the novel intervention is approximated from the mean and 
# standard deviation of the PA distribution for the probability of side effects.
beta_params_t2 <- betaPar(mean(m_params$p_side_effects_t2),
        sd(m_params$p_side_effects_t2))

# EVSI calculation with moment matching method
evsi_default <- evsi(outputs = chemotherapy_output,
                     inputs = m_params,
                     study = "trial_binary",
                     pars = c("p_side_effects_t1", "p_side_effects_t2"),
                     n = seq(50, 1500, by = 50),
                     method = "mm",
                     model_fn = calculate_costs_effects,
                     analysis_args = list(a1 = beta_params_t1[1],
                                          b1 = beta_params_t1[2],
                                          a2 = beta_params_t2$alpha,
                                          b2 = beta_params_t2$beta),
                     par_fn = generate_psa_parameters)

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

#### STUDY 2: Randomised Trial with Multiple Outcomes ####
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

  return(data.frame(p_side_effects_t1 = args$p_side_effects_t1,
                    p_side_effects_t2 = bugs.data$BUGSoutput$sims.matrix[, "p_side_effects_t2"],
                    p_hospitalised_total= bugs.data$BUGSoutput$sims.matrix[, "p_hospitalised_total"],
                    p_died = bugs.data$BUGSoutput$sims.matrix[, "p_died"],
                    lambda_home = bugs.data$BUGSoutput$sims.matrix[, "lambda_home"],
                    lambda_hosp = bugs.data$BUGSoutput$sims.matrix[, "lambda_hosp"]))
}

# EVSI calculation using the momemt matching method.
evsi_OR_allout_MM <- evsi(outputs = chemotherapy_output,
                inputs = m_params,
                pars = c("p_side_effects_t1", "logor_side_effects",
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

#### STUDY 3: Long-term Survival ####
longterm_datagen_fn <- function(inputs, n = 46000){
  rate_longterm <- inputs[, "rate_longterm"]
  sum_of_surv <- rgamma(dim(inputs)[1], shape = 2 * n, scale = n / rate_longterm)
  return(data.frame(surv_sum = sum_of_surv))
}
data <- longterm_datagen_fn(m_params)
pars <- "rate_longterm"

longterm_analysis_fn <- function(data, args, pars){
  # Load key function
  gammaPar <- args$gammaPar
  
  # Data list for JAGS
  data_jags <- list(
    surv_sum = data$surv_sum[1],
    alpha_rate = gammaPar(args$rate_longterm_mu, 
                          args$rate_longterm_sd)$alpha,
    beta_rate = gammaPar(args$rate_longterm_mu, 
                         args$rate_longterm_sd)$beta,
    n = args$n
  )
  
  longterm_jags <- function(){
    ## Models for the data
    surv_sum ~ dgamma(2 * n, rate_longterm / n)
    rate_longterm ~ dgamma(alpha_rate, beta_rate)
  }
  
  filein <- file.path(tempdir(),fileext="datmodel.txt")
  R2OpenBUGS::write.model(longterm_jags,filein)
  
  # Perform the MCMC simulation with OpenBUGS.
  # Close OpenBUGS once it has finished (if debug is set to TRUE)
  bugs.data <- jags(
    data =  data_jags,
    parameters.to.save = pars,
    model.file = filein, 
    n.chains = 3, 
    n.iter = args$n.iter, 
    n.thin = 1, 
    n.burnin = 250, progress.bar = "none") 
  
  return(data.frame(rate_longterm = bugs.data$BUGSoutput$sims.matrix[, "rate_longterm"]))
}


# EVSI calculation using the momemt matching method.
evsi_longterm <- evsi(outputs = chemotherapy_output,
                       inputs = m_params,
                       pars = c("rate_longterm"),
                       n = 46000,
                       method = "mm",
                       datagen_fn = longterm_datagen_fn,
                       model_fn = calculate_costs_effects,
                       analysis_args = list(n = 40000,
                                            gammaPar = gammaPar,
                                            rate_longterm_mu = rate_longterm_mu,
                                            rate_longterm_sd = rate_longterm_sd,
                                            n.iter = 2000),
                       analysis_fn = longterm_analysis_fn, 
                       par_fn = generate_psa_parameters)

plotting <- evsi.plot.adapt(chemotherapy_output, m_params, c("rate_longterm"), 
                            evsi_longterm, "gam")

#### STUDY 4: Retrospective Study of Hospital Data ####
# Data generation function
retrohosp_datagen_fn <- function(inputs, n = 500){
  # Load the data
  p_died <- inputs[, "p_died"]
  lambda_hosp <- inputs[, "lambda_hosp"]
  rate_recover_hosp <- -log(1 -lambda_hosp )
  
  X_dead <- N_recover_hospital <- vector("numeric", length = dim(inputs)[1])
  T_hosp <- matrix(NA, nrow = dim(inputs)[1], ncol = n)
  for(i in 1:dim(inputs)[1]){
    # Simulate the number of patients die
    X_dead[i] <- rbinom(1, n, p_died[i])
    
    ## Simulate recovery times for patients
    N_recover_hospital[i] <- n - X_dead[i]
    T_hosp[i, 1:N_recover_hospital[i]] <- rexp(N_recover_hospital[i], rate_recover_hosp[i])
  }
  
  data_save_dat <- data.frame(cbind(X_dead = X_dead, T_hosp = T_hosp))
  return(data_save_dat)
}

# Analysis function based on JAGS
retrohosp_analysis_fn <- function(data, args, pars){
  # Create the data list for JAGS
  data_jags <- list(X_dead = data$X_dead,
                    N_recover_hosp = args$n - data$X_dead,
                    T_hosp = as.vector(data[, 1 + (1:(args$n - data$X_dead))]),
                    n = args$n,
                    p_recovery_hosp_alpha = betaPar(args$p_recovery_hosp_mu, 
                                                    args$p_recovery_hosp_sd)$alpha,
                    p_recovery_hosp_beta = betaPar(args$p_recovery_hosp_mu,
                                                   args$p_recovery_hosp_sd)$beta,
                    n_died = args$n_died,
                    n_hospitalised = args$n_hospitalised)
  
  LogOR_addoutcomes_trial <- function(){
    
    ## Models for the data
    X_dead ~ dbin(p_died, n)
    
    rate_recover_hosp <- -log(1 - lambda_hosp)
    for(i in 1:N_recover_hosp){
      T_hosp[i] ~ dexp(rate_recover_hosp)
    }
    
    # Probability that a patient dies over the time horizon given they were 
    # hospitalised
    p_died ~ dbeta(1 + n_died, 1 + n_hospitalised - n_died)
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
    n.chains = 3, 
    n.iter = args$n.iter, 
    n.thin = 1, 
    n.burnin = 250, progress.bar = "none") 
  
  return(data.frame(p_died = bugs.data$BUGSoutput$sims.matrix[, "p_died"],
                    lambda_hosp = bugs.data$BUGSoutput$sims.matrix[, "lambda_hosp"]))
}

# EVSI calculation using the momemt matching method.
evsi_retrohosp <- evsi(outputs = chemotherapy_output,
                       inputs = m_params,
                       pars = c("p_died", "lambda_hosp"),
                       n = seq(500, 1500, by = 200),
                       method = "mm",
                       datagen_fn = retrohosp_datagen_fn,
                       model_fn = calculate_costs_effects,
                       analysis_args = list(n = 500,
                                            betaPar = betaPar,
                                            p_recovery_hosp_mu = p_recovery_hosp_mu,
                                            p_recovery_hosp_sd = p_recovery_hosp_sd,
                                            n.iter = 5000,
                                            n_died = n_died,
                                            n_hospitalised = n_hospitalised),
                       analysis_fn = retrohosp_analysis_fn, 
                       par_fn = generate_psa_parameters)

#### STUDY 5: Registry Study of Observational Data ####
# Data generation function
registry_datagen_fn <- function(inputs, n = 500){
  # Load the data
  p_hospitalised_total <- inputs[, "p_hospitalised_total"]
  p_died <- inputs[, "p_died"]
  lambda_home <- inputs[, "lambda_home"]
  lambda_hosp <- inputs[, "lambda_hosp"]
  rate_recover_hosp <- -log(1 -lambda_hosp )
  rate_recover_home <- -log(1 - lambda_home)
  
  X_hosp <- X_dead <- N_recover_home <-  N_recover_hospital <- vector("numeric", length = dim(inputs)[1])
  T_home <- T_hosp <- matrix(NA, nrow = dim(inputs)[1], ncol = 2 * n)
  for(i in 1:dim(inputs)[1]){
    # Simulate the number of patients hospitalised 
    X_hosp[i] <- rbinom(1, n, p_hospitalised_total[i])
    # Simulate the number of patients die
    X_dead[i] <- rbinom(1, X_hosp[i], p_died[i])
    
    ## Simulate recovery times for patients
    N_recover_home[i] <- n - X_hosp[i]
    T_home[i, 1:N_recover_home[i]] <- rexp(N_recover_home[i], rate_recover_home[i])
    N_recover_hospital[i] <- X_hosp[i] - X_dead[i]
    T_hosp[i, 1:N_recover_hospital[i]] <- rexp(N_recover_hospital[i], rate_recover_hosp[i])
    
  }
  
  data_save_dat <- data.frame(cbind(X_hosp = X_hosp, X_dead = X_dead,
                                    N_recover_home = N_recover_home,
                                    N_recover_hospital = N_recover_hospital,
                                    T_home = T_home, T_hosp = T_hosp))
  return(data_save_dat)
}

# Analysis function based on JAGS
registry_analysis_fn <- function(data, args, pars){
  # Create the data list for JAGS
  data_jags <- list(X_hosp = data$X_hosp,
                    X_dead = data$X_dead,
                    T_home = as.vector(data[, (1:data$N_recover_home) + 6]),
                    T_hosp = as.vector(data[, (6 + 2 * args$n) + (1:data$N_recover_hospital)]),
                    N_recover_home = data$N_recover_home,
                    N_recover_hosp = data$N_recover_hospital,
                    n = args$n,
                    n_side_effects = args$n_side_effects,
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
    X_hosp ~ dbinom(p_hospitalised_total, n)
    X_dead ~ dbin(p_died, X_hosp)
    
    rate_recover_home <- -log(1 - lambda_home)
    rate_recover_hosp <- -log(1 - lambda_hosp)
    
    for(i in 1:N_recover_home){
      T_home[i] ~ dexp(rate_recover_home)
    }
    
    for(i in 1:N_recover_hosp){
      T_hosp[i] ~ dexp(rate_recover_hosp)
    }
  
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
    n.chains = 3, 
    n.iter = args$n.iter, 
    n.thin = 1, 
    n.burnin = 250, progress.bar = "none") 
  
  return(data.frame(p_hospitalised_total= bugs.data$BUGSoutput$sims.matrix[, "p_hospitalised_total"],
                    p_died = bugs.data$BUGSoutput$sims.matrix[, "p_died"],
                    lambda_home = bugs.data$BUGSoutput$sims.matrix[, "lambda_home"],
                    lambda_hosp = bugs.data$BUGSoutput$sims.matrix[, "lambda_hosp"]))
}

# EVSI calculation using the momemt matching method.
evsi_registry <- evsi(outputs = chemotherapy_output,
                          inputs = m_params,
                          pars = c("p_hospitalised_total","p_died", "lambda_home",
                                   "lambda_hosp"),
                          n = seq(500, 1500, by = 200),
                          method = "mm",
                          datagen_fn = registry_datagen_fn,
                          model_fn = calculate_costs_effects,
                          analysis_args = list(n_side_effects = n_side_effects,
                                               n = 500,
                                               betaPar = betaPar,
                                               p_recovery_home_mu = p_recovery_home_mu,
                                               p_recovery_home_sd = p_recovery_home_sd,
                                               p_recovery_hosp_mu = p_recovery_hosp_mu,
                                               p_recovery_hosp_sd = p_recovery_hosp_sd,
                                               n.iter = 5000,
                                               n_died = n_died,
                                               n_hospitalised = n_hospitalised),
                          analysis_fn = registry_analysis_fn, 
                          par_fn = generate_psa_parameters)

#### STUDY 6: A Cost Analysis ####
# Data generation function
cost_datagen_fn <- function(inputs, n = 500, 
                            v_home_care_fun = function(){return(sqrt(rgamma(1, shape = 3, scale = 39)))},
                            v_hospital_fun = function(){return(sqrt(rgamma(1, shape = 10, scale = 45)))},
                            v_death_fun = function(){return(sqrt(rgamma(1, shape = 15, scale = 112.5)))}
                            ){
  lognormPar = function(m,s) {
    # m: Mean of Log-Normal distribution
    # s: Standard deiviation of Log-Normal distribution
    
    var <- s^2
    meanlog <- log(m) - 0.5 * log(1 + var/m^2)
    varlog <- log(1 + (var/m^2))
    sdlog <- sqrt(varlog)
    
    return(
      list(meanlog = meanlog, sdlog = sdlog)
    )
  }

  
  X_home_care <- X_hospital <- X_death <- matrix(NA, nrow = dim(inputs)[1], ncol = n[1])
  for(i in 1:dim(inputs)[1]){
    # Load the data
    m_home_care <- inputs[i, "c_home_care"]
    m_hospital <- inputs[i, "c_hospital"]
    m_death <- inputs[i, "c_death"]
    v_home_care <- v_home_care_fun()
    v_hospital <- v_hospital_fun()
    v_death <- v_death_fun()
    par_home_care <- lognormPar(m_home_care, sqrt(v_home_care))
    par_hospital <- lognormPar(m_hospital, sqrt(v_hospital))
    par_death <- lognormPar(m_death, sqrt(v_death))
    # Simulate the costs 
    X_home_care[i, ] <- rlnorm(n[1], par_home_care$meanlog, par_home_care$sdlog)
    X_hospital[i, ] <- rlnorm(n[1], par_hospital$meanlog, par_hospital$sdlog)
    X_death[i, ] <- rlnorm(n[1], par_death$meanlog, par_death$sdlog)
    
  }
  
  data_save_dat <- data.frame(cbind(X_home_care = X_home_care,
                                    X_hospital = X_hospital,
                                    X_death = X_death))
  return(data_save_dat)
}

dat_try <- cost_datagen_fn(m_params)



# Analysis function based on JAGS
cost_analysis_fn <- function(data, args, pars){
  # Create the data list for JAGS
  data_jags <- list(X_home_care = as.vector(data[, (1:args$n)]),
                    X_hospital = as.vector(data[, args$n + (1:args$n)]),
                    X_death = as.vector(data[, 2*args$n + (1:args$n)]),
                    n = args$n,
                    mu_home_care = args$lognPar(args$c_home_care_mu, 
                                                    args$c_home_care_sd)$meanlog,
                    t_home_care = 1/args$lognPar(args$c_home_care_mu, 
                                                   args$c_home_care_sd)$sdlog^2,
                    mu_hospital = args$lognPar(args$c_hospital_mu, 
                                                    args$c_hospital_sd)$meanlog,
                    t_hospital = 1 / args$lognPar(args$c_hospital_mu,
                                                   args$c_hospital_sd)$sdlog^2,
                    mu_death = args$lognPar(args$c_death_mu,
                                       args$c_death_sd)$meanlog,
                    t_death = 1 / args$lognPar(args$c_death_mu,
                                       args$c_death_sd)$sdlog^2,
                    a_home_care = 3,
                    b_home_care = 1 / 39,
                    a_hospital = 10,
                    b_hospital = 1 / 45,
                    a_death = 15,
                    b_death = 1 / 112.5)
  
  LogOR_addoutcomes_trial <- function(){
    for(i in 1:n){
      X_home_care[i] ~ dlnorm(m_home_care, tau_ind_home_care)
      X_hospital[i] ~ dlnorm(m_hospital, tau_ind_hospital)
      X_death[i] ~ dlnorm(m_death, tau_ind_death)
    }
    
    c_home_care ~ dlnorm(mu_home_care, t_home_care)
    c_hospital ~ dlnorm(mu_hospital, t_hospital)
    c_death ~ dlnorm(mu_death, t_death)
    
    v_ind_home_care ~ dgamma(a_home_care, b_home_care)
    v_ind_hospital ~ dgamma(a_hospital, b_hospital)
    v_ind_death ~ dgamma(a_death, b_death)
    
    m_home_care <- log(c_home_care) - 0.5 * 
      log(1 + v_ind_home_care/c_home_care^2)
    tau_ind_home_care <- 1 / log(1 + (v_ind_home_care/m_home_care^2))
    m_hospital <- log(c_hospital) - 0.5 * 
      log(1 + v_ind_hospital/c_hospital^2)
    tau_ind_hospital <- 1 / log(1 + (v_ind_hospital/m_hospital^2))
    m_death <- log(c_death) - 0.5 * 
      log(1 + v_ind_death/c_death^2)
    tau_ind_death <- 1 / log(1 + (v_ind_death/m_death^2))

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
  
  return(data.frame(c_home_care= bugs.data$BUGSoutput$sims.matrix[, "c_home_care"],
                    c_hospital = bugs.data$BUGSoutput$sims.matrix[, "c_hospital"],
                    c_death = bugs.data$BUGSoutput$sims.matrix[, "c_death"]))
}


# EVSI calculation using the momemt matching method.
evsi_costs <- evsi(outputs = chemotherapy_output,
                         inputs = m_params,
                         pars = c("c_home_care", "c_hospital", "c_death"),
                         n = seq(30, 1000, by = 200),
                         method = "mm",
                         datagen_fn = cost_datagen_fn,
                         model_fn = calculate_costs_effects,
                         analysis_args = list(n = 20,
                                              lognPar = lognPar,
                                              c_home_care_mu = c_home_care_mu,
                                              c_home_care_sd = c_home_care_sd,
                                              c_hospital_mu = c_hospital_mu,
                                              c_hospital_sd = c_hospital_sd,
                                              c_death_mu = c_death_mu,
                                              c_death_sd = c_death_sd,
                                              n.iter = 2000),
                         analysis_fn = cost_analysis_fn, 
                         par_fn = generate_psa_parameters,
                      Q = 50)

#### STUDY 7: A Utility Analysis ####
# Data generation function
utility_datagen_fn <- function(inputs, n = 500, 
                            sd_recovery_fun = function(){return(runif(1, 0.000001, 0.00005))},
                            sd_home_care_fun = function(){return(runif(1, 0.00001, 0.005))},
                            sd_hospital_fun = function(){return(runif(1, 0.00001, 0.01))}
){
  betaPar <- function(m, s) {
    # m:  Mean of the Beta distribution
    # m: Standard deviation of the Beta distribution
    
    var <- s ^ 2
    alpha <- ((1 - m) / var - 1 / m) * m ^ 2
    beta <- alpha * (1 / m - 1)
    
    return(
      list(alpha = alpha, beta = beta)
    )
  }
  # Load the data
  X_home_care <- X_hospital <- X_recovery <- matrix(NA, nrow = dim(inputs)[1], ncol = n[1])
  for(i in 1:dim(inputs)[1]){
    set.seed(123 + i)
    m_recovery <- inputs[i, "u_recovery"]
    m_home_care <- inputs[i, "u_home_care"]
    m_hospital <- inputs[i, "u_hospital"]
    sd_recovery <- sd_recovery_fun()
    sd_home_care <- sd_home_care_fun()
    sd_hospital <- sd_hospital_fun()
    
    par_recovery <- betaPar(m_recovery, sd_recovery)
    par_home_care <- betaPar(m_home_care, sd_home_care)
    par_hospital <- betaPar(m_hospital, sd_hospital)
    
    # Simulate the costs 
    X_recovery[i, ] <- rbeta(n[1], par_recovery$alpha, par_recovery$beta)
    X_home_care[i, ] <- rbeta(n[1], par_home_care$alpha, par_home_care$beta)
    X_hospital[i, ] <- rbeta(n[1], par_hospital$alpha, par_hospital$beta)

    
  }
  
  data_save_dat <- data.frame(cbind(X_recovery = X_recovery,
                                    X_home_care = X_home_care,
                                    X_hospital = X_hospital))
  return(data_save_dat)
}
dat_try <- utility_datagen_fn(m_params)

# Analysis function based on JAGS
utility_analysis_fn <- function(data, args, pars){
  # Create the data list for JAGS
  data_jags <- list(X_recovery = as.vector(data[, (1:args$n)]),
                    X_home_care = as.vector(data[, args$n + (1:args$n)]),
                    X_hospital = as.vector(data[, 2*args$n + (1:args$n)]),
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
      X_recovery[i] ~ dbeta(a_recovery, b_recovery)
      X_home_care[i] ~ dbeta(a_home_care, b_home_care)
      X_hospital[i] ~ dbeta(a_hospital, b_hospital)
    }
    
    u_recovery ~ dbeta(alpha_recovery, beta_recovery)
    u_home_care ~ dbeta(alpha_home_care, beta_home_care)
    u_hospital ~ dbeta(alpha_hospital, beta_hospital)
    
    sd_recovery ~ dunif(0.000001, 0.00005)
    sd_home_care ~ dunif(0.00001, 0.005)
    sd_hospital ~ dunif(0.00001, 0.01)
    
    v_recovery <- sd_recovery ^ 2
    a_recovery <- ((1 - u_recovery) / v_recovery - 1 / u_recovery) * u_recovery ^ 2
    b_recovery <- a_recovery * (1 / u_recovery - 1)
    
    v_home_care <- sd_home_care ^ 2
    a_home_care <- ((1 - u_home_care) / v_home_care - 1 / u_home_care) * u_home_care ^ 2
    b_home_care <- a_home_care * (1 / u_home_care - 1)
    
    v_hospital <- sd_hospital ^ 2
    a_hospital <- ((1 - u_hospital) / v_hospital - 1 / u_hospital) * u_hospital ^ 2
    b_hospital <- a_hospital * (1 / u_hospital - 1)
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
undebug(evsi)
evsi_utility <- evsi(outputs = chemotherapy_output,
                   inputs = m_params,
                   pars = c("u_recovery", "u_home_care", "u_hospital"),
                   n = seq(50, 1500, by = 50),
                   method = "mm",
                   datagen_fn = utility_datagen_fn,
                   model_fn = calculate_costs_effects,
                   analysis_args = list(n = 50,
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

plotting <- evsi.plot.adapt(chemotherapy_output, m_params, c("u_recovery", "u_home_care", "u_hospital"), 
                            evsi_utility, "earth")
evsi.wtp.plot(plotting)
pop.adjust <- 46000 * (1 / (1 + 0.035)^3)
evsi.enbs.plot(plotting, c(1260000, 1400000), 2 * c(1560.55, 1600), 
               k = 20000, Pop = pop.adjust, Time = 10)
optim.ss(plotting, c(1260000, 1400000), 2 * c(1560.55, 1600), 
         k = 20000, Pop = pop.adjust, Time = 10)
coss(plotting, c(1260000, 1400000), 2 * c(1560.55, 1600), 
     Pop = pop.adjust, Time = 10)

