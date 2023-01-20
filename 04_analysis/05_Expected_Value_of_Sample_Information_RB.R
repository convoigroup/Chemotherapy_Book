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
                     n = seq(500, 1500, by = 200),
                     method = "gam",
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
  p_side_effects_t2 <- inputs[, "p_side_effects_t2"]
  X1 <- rbinom(length(p_side_effects_t1), n, p_side_effects_t1)
  X2 <- rbinom(length(p_side_effects_t2), n, p_side_effects_t2)
  OR <- (n - X2) / X2 / ((n - X1) / X1)
  data_save <- data.frame(OR = OR)
  return(data_save)
}


# EVSI calculation using the momemt matching method.
evsi_OR <- evsi(outputs = chemotherapy_output,
                inputs = m_params,
                pars = c("p_side_effects_t1", "p_side_effects_t2"),
                n = seq(500, 1500, by = 200),
                method = "gam",
                datagen_fn = OR_datagen_fn,
                model_fn = calculate_costs_effects,
                analysis_fn = OR_analysis_fn, 
                par_fn = generate_psa_parameters)

#### STUDY 2: Randomised Trial with Multiple Outcomes ####
# Data generation function
full_datagen_fn <- function(inputs, n = 500){
  # Load the data
  p_side_effects_t1 <- inputs[, "p_side_effects_t1"]
  p_side_effects_t2 <- inputs[, "p_side_effects_t2"]
  p_hospitalised_total <- inputs[, "p_hospitalised_total"]
  p_died <- inputs[, "p_died"]
  lambda_home <- inputs[, "lambda_home"]
  lambda_hosp <- inputs[, "lambda_hosp"]
  rate_recover_hosp <- -log(1 -lambda_hosp )
  rate_recover_home <- -log(1 - lambda_home)
  
  X1 <- X2 <- X_hosp <- X_dead <- N_recover_home <- 
    N_recover_hospital <- T_home_sum <- T_hosp_sum <- vector("numeric", length = dim(inputs)[1])
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
    T_home <- rexp(N_recover_home[i], rate_recover_home[i])
    N_recover_hospital[i] <- X_hosp[i] - X_dead[i]
    T_hosp <- rexp(N_recover_hospital[i], rate_recover_hosp[i])
    
    T_home_sum[i] <- sum(T_home)
    T_hosp_sum[i] <- sum(T_hosp)
  }
  OR <- (n - X2) / X2 / ((n - X1) / X1)
  p_hosp <- X_hosp / (X1 + X2)
  p_dead <- X_dead / X_hosp

  
  data_save_dat <- data.frame(cbind(OR = OR, 
                                    p_hosp = p_hosp, p_dead = p_dead,
                                    T_home_sum = T_home_sum, 
                                    T_hosp_sum = T_hosp_sum))
  return(data_save_dat)
}

# EVSI calculation using the momemt matching method.
evsi_OR_allout <- evsi(outputs = chemotherapy_output,
                       inputs = m_params,
                       pars = c("p_side_effects_t1", "p_side_effects_t2",
                                "p_hospitalised_total", "p_died",
                                "lambda_home", "lambda_hosp"),
                       n = seq(500, 1500, by = 200),
                       method = "earth",
                       datagen_fn = full_datagen_fn,
                       par_fn = generate_psa_parameters,
                       npreg_method = "inla")

#### STUDY 3: Long-term Survival ####
longterm_datagen_fn <- function(inputs, n = 40000){
  rate_longterm <- inputs[, "rate_longterm"]
  sum_of_surv <- rgamma(dim(inputs)[1], shape = 2 * n, scale = n / rate_longterm)
  return(data.frame(surv_sum = sum_of_surv))
}

# EVSI calculation using the momemt matching method.
evsi_longterm <- evsi(outputs = chemotherapy_output,
                      inputs = m_params,
                      pars = c("rate_longterm"),
                      n = seq(40000, 50000, by = 500),
                      method = "gam",
                      datagen_fn = longterm_datagen_fn,
                      model_fn = calculate_costs_effects,
                      analysis_args = list(n = 40000,
                                           gammaPar = gammaPar,
                                           rate_longterm_mu = rate_longterm_mu,
                                           rate_longterm_sd = rate_longterm_sd,
                                           n.iter = 2000),
                      analysis_fn = longterm_analysis_fn, 
                      par_fn = generate_psa_parameters)

#### STUDY 4: Retrospective Study of Hospital Data ####
# Data generation function
retrohosp_datagen_fn <- function(inputs, n = 500){
  # Load the data
  p_died <- inputs[, "p_died"]
  lambda_hosp <- inputs[, "lambda_hosp"]
  rate_recover_hosp <- -log(1 -lambda_hosp )
  
  X_dead <-T_hosp_sum <- N_recover_hospital <- vector("numeric", length = dim(inputs)[1])
  for(i in 1:dim(inputs)[1]){
    # Simulate the number of patients die
    X_dead[i] <- rbinom(1, n, p_died[i])
    
    ## Simulate recovery times for patients
    N_recover_hospital[i] <- n - X_dead[i]
    T_hosp <- rexp(N_recover_hospital[i], rate_recover_hosp[i])
    T_hosp_sum[i] <- sum(T_hosp)
  }
  
  data_save_dat <- data.frame(cbind(p_dead = X_dead / n, T_hosp_sum = T_hosp_sum))
  return(data_save_dat)
}

# EVSI calculation using the momemt matching method.
evsi_retrohosp <- evsi(outputs = chemotherapy_output,
                       inputs = m_params,
                       pars = c("p_died", "lambda_hosp"),
                       n = seq(500, 1500, by = 200),
                       method = "gam",
                       datagen_fn = retrohosp_datagen_fn,
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
  
  X_hosp <- X_dead <- N_recover_home <- T_home_sum <- 
    T_hosp_sum <- N_recover_hospital <- vector("numeric", length = dim(inputs)[1])
  for(i in 1:dim(inputs)[1]){
    # Simulate the number of patients hospitalised 
    X_hosp[i] <- rbinom(1, n, p_hospitalised_total[i])
    # Simulate the number of patients die
    X_dead[i] <- rbinom(1, X_hosp[i], p_died[i])
    
    ## Simulate recovery times for patients
    N_recover_home[i] <- n - X_hosp[i]
    T_home <- rexp(N_recover_home[i], rate_recover_home[i])
    T_home_sum[i] <- sum(T_home)
    N_recover_hospital[i] <- X_hosp[i] - X_dead[i]
    T_hosp <- rexp(N_recover_hospital[i], rate_recover_hosp[i])
    T_hosp_sum[i] <- sum(T_hosp)
    
  }
  
  data_save_dat <- data.frame(cbind(p_hosp = X_hosp / n, p_dead = X_dead / X_hosp,
                                    T_home_sum = T_home_sum, T_hosp_sum = T_hosp_sum))
  return(data_save_dat)
}

# EVSI calculation using the momemt matching method.
evsi_registry <- evsi(outputs = chemotherapy_output,
                      inputs = m_params,
                      pars = c("p_hospitalised_total","p_died", "lambda_home",
                               "lambda_hosp"),
                      n = seq(500, 1500, by = 200),
                      method = "gam",
                      datagen_fn = registry_datagen_fn,
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
  X_home_care_suff1 <- X_home_care_suff2 <- X_hospital_suff1 <- X_hospital_suff2 <- 
    X_death_suff1 <- X_death_suff2 <- vector("numeric", dim(inputs)[1])
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
    
    X_home_care_suff1[i] <- sum(log(X_home_care[i, ]))
    X_home_care_suff2[i] <- sum(log(X_home_care[i, ])^2)
    X_hospital_suff1[i] <- sum(log(X_hospital[i, ]))
    X_hospital_suff2[i] <- sum(log(X_hospital[i, ])^2)
    X_death_suff1[i] <- sum(log(X_death[i, ]))
    X_death_suff2[i] <- sum(log(X_death[i, ])^2)
    
  }
  
  data_save_dat <- data.frame(cbind(X_home_care_suff1 = X_home_care_suff1,
                                    X_home_care_suff2 = X_home_care_suff2,
                                    X_hospital_suff1 = X_hospital_suff1,
                                    X_hospital_suff2 = X_hospital_suff2,
                                    X_death_suff1 = X_death_suff1,
                                    X_death_suff2 = X_death_suff2))
  return(data_save_dat)
}

# EVSI calculation using the momemt matching method.
evsi_costs <- evsi(outputs = chemotherapy_output,
                   inputs = m_params,
                   pars = c("c_home_care", "c_hospital", "c_death"),
                   n = seq(30, 1000, by = 200),
                   method = "earth",
                   datagen_fn = cost_datagen_fn,
                   par_fn = generate_psa_parameters)

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
  
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  
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
    
    par_recovery <- betaPar(m_recovery, sd_recovery)
    par_home_care <- betaPar(m_home_care, sd_home_care)
    par_hospital <- betaPar(m_hospital, sd_hospital)
    
    # Simulate the costs 
    X_recovery[i, ] <- rbeta(n[1], par_recovery$alpha, par_recovery$beta)
    X_home_care[i, ] <- rbeta(n[1], par_home_care$alpha, par_home_care$beta)
    X_hospital[i, ] <- rbeta(n[1], par_hospital$alpha, par_hospital$beta)
    
    X_recovery_mean1[i] <- gm_mean(X_recovery[i, ])
    X_recovery_mean2[i] <- gm_mean(1 - X_recovery[i, ])
    
    X_home_care_mean1[i] <- gm_mean(X_home_care[i, ])
    X_home_care_mean2[i] <- gm_mean(1 - X_home_care[i, ])
    
    X_hospital_mean1[i] <- gm_mean(X_hospital[i, ])
    X_hospital_mean2[i] <- gm_mean(1 - X_hospital[i, ])
  }
  
  
  
  data_save_dat <- data.frame(cbind(X_recovery_mean1 = X_recovery_mean1,
                                    X_recovery_mean2 = X_recovery_mean2,
                                    X_home_care_mean1 = X_home_care_mean1,
                                    X_home_care_mean2 = X_home_care_mean2,
                                    X_hospital_mean1 = X_hospital_mean1,
                                    X_hospital_mean2 = X_hospital_mean2))
  return(data_save_dat)
}

# EVSI calculation using the momemt matching method.
evsi_utility <- evsi(outputs = chemotherapy_output,
                     inputs = m_params,
                     pars = c("u_recovery", "u_home_care", "u_hospital"),
                     n = seq(30, 1000, by = 200),
                     method = "earth",
                     datagen_fn = utility_datagen_fn,
                     par_fn = generate_psa_parameters)


