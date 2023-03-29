########################################################################
### Expected Value of Sample Information Analysis - Regression Based ###
########################################################################
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
# EVSI calculation with GAM regression
evsi_default <- evsi(outputs = chemotherapy_output,
                     inputs = m_params,
                     study = "trial_binary",
                     pars = c("p_side_effects_t1", "p_side_effects_t2"),
                     n = seq(500, 1500, by = 200),
                     method = "gam")

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

#### STUDY 2: Randomised Trial with Multiple Outcomes ####
# Data generation function
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

evsi.wtp.plot(plotting)
pop.adjust <- 46000 * (1 / (1 + 0.035)^1)
optim.ss(plotting, c(60000, 60000), c(0, 0), 
         k = 30000, Pop = pop.adjust, Time = 12)
evsi.prob.plot(plotting, c(60000, 60000), c(0, 0),
               k = 20000, Pop = c(0, pop.adjust * 2), c(0, 15))

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

# EVSI calculation using GAM regression
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

# EVSI calculation using GAM regression
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
    
    ## Sufficient statistics for the log-normal distribution
    # Data provides information on mean and variance of log-normal
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

# EVSI calculation using MARS regression due to large number of parameters.
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
    
    ## Sufficient statistic for beta distribution is geometric mean of X and (1-X)
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

# EVSI calculation using MARS regression due to large number of parameters.
evsi_utility <- evsi(outputs = chemotherapy_output,
                     inputs = m_params,
                     pars = c("u_recovery", "u_home_care", "u_hospital"),
                     n = seq(30, 1000, by = 200),
                     method = "earth",
                     datagen_fn = utility_datagen_fn)


