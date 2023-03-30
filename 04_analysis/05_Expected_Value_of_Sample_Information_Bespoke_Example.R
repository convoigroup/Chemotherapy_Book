#######################################################################
### Expected Value of Sample Information Analysis - Example         ###
#######################################################################
### Packages
library(voi)
library(ggplot2)
library(dplyr)
library(R2jags)

## Run the model
source("04_analysis/02_baseline_model_output.R")
source("03_R/01_misc_functions.R")

## Baseline Cost-Effectiveness Formatting
# The output from the cost-effectiveness model should be formatted for 
# the voi package.
# Use BCEA package to create a cost-effectiveness object.
chemotherapy_output <- list(e = m_costs_effects[, "Effects", ],
                            c = m_costs_effects[, "Costs", ],
                            k = seq(0, 50000, length.out = 501))

# Data generation function - aggregate data (for regression method)
utility_datagen_fn_agg <- function(inputs, n = 500){
  dat_indiv <- utility_datagen_fn_indiv(inputs, n = n)
  X_hospital_mean <- rowMeans(dat_indiv)
  data_save_dat <- data.frame(X_hospital_mean = X_hospital_mean)
  return(data_save_dat)
}

# Data generation function - individual data (for other EVSI methods)
utility_datagen_fn_indiv <- function(inputs, n = 500){
  # Load the data
  X_hospital <- matrix(nrow = nrow(inputs), ncol = n[1])
  X_hospital_mean1 <- X_hospital_mean2 <- numeric(nrow(inputs))
  for(i in 1:nrow(inputs)){
    set.seed(123 + i)
    m_hospital <- inputs[i, "u_hospital"]
    sd_hospital <- inputs[i, "sd_iid_hospital"]
    X_hospital[i, ] <- truncnorm::rtruncnorm(n[1], 
                                             mean = m_hospital, 
                                             sd = sd_hospital,
                                             a = -Inf, b = 1)
  }
  data_save_dat <- data.frame(cbind(X_hospital = X_hospital))
  return(data_save_dat)
}

m_params$sd_iid_hospital <- runif(nrow(m_params), 0.00001, 0.4)

## Regression Based Method
evsi_utility <- evsi(outputs = chemotherapy_output,
                     inputs = m_params,
                     pars = c("u_hospital"),
                     n = seq(50, 1000, by = 200),
                     method = "gam",
                     datagen_fn = utility_datagen_fn_rb)

## Moment Matching Method
# Analysis function based on JAGS
utility_analysis_fn <- function(data, args, pars){
  # Create the data list for JAGS
  data_jags <- list(X_hospital = as.vector(data),
                    n = args$n,
                    alpha_hospital = betaPar(
                      args$u_hospital_mu,
                      args$u_hospital_sd
                    )$alpha,
                    beta_hospital = betaPar(
                      args$u_hospital_mu,
                      args$u_hospital_sd
                    )$beta)
  
  trial <- "
  model { 
    for(i in 1:n){
      X_hospital[i] ~ dnorm(u_hospital, tau_hospital)T(, 1)
    }
    u_hospital ~ dbeta(alpha_hospital, beta_hospital)
    sd_hospital ~ dunif(0.00001, 0.4)
    tau_hospital <- 1 / sd_hospital ^ 2
  }
  "
  
  filein <- file.path(tempdir(),fileext="datmodel.txt")
  cat(trial, file=filein)
  
  # Perform the MCMC simulation with JAGS.
  bugs.data <- jags(
    data =  data_jags,
    parameters.to.save = pars,
    model.file = filein, 
    n.chains = 1, 
    n.iter = args$n.iter, 
    n.thin = 1, 
    n.burnin = 250, 
    quiet=TRUE, progress.bar = "none") 
  
  u_hospital <- bugs.data$BUGSoutput$sims.matrix[, "u_hospital"]
  return(data.frame(u_hospital = u_hospital))
}

analysis_args <- list(n = 30,
                      u_hospital_mu = u_hospital_mu,
                      u_hospital_sd = u_hospital_sd,
                      n.iter = 2000)

evsi_utility <- evsi(outputs = chemotherapy_output,
                     inputs = m_params,
                     pars = c("u_hospital"),
                     pars_datagen = c("u_hospital","sd_iid_hospital"),
                     n = seq(30, 1000, by = 200),
                     method = "mm",
                     datagen_fn = utility_datagen_fn_indiv,
                     model_fn = calculate_costs_effects,
                     analysis_args = analysis_args,
                     analysis_fn = utility_analysis_fn, 
                     par_fn = generate_psa_parameters,
                     Q = 50)

## Importance Sampling
# Likelihood function
utility_likelihood <- function(data, inputs){
  # Load the data
  ll <- numeric(nrow(inputs)) 
  data_vec <- unlist(data)

  for(i in 1:nrow(inputs)){
    m_hospital <- inputs[i, "u_hospital"]
    sd_hospital <- inputs[i, "sd_iid_hospital"]
    ll[i] <- exp(
      sum(
        log(
          truncnorm::dtruncnorm(data_vec, 
                                mean = m_hospital, 
                                sd = sd_hospital,
                                a = -Inf, b = 1)
        )))
  }
  return(ll)
}

# Importance Sampling - EVSI.  (this is slow).
evsi_utility <- evsi(outputs = chemotherapy_output,
                     inputs = m_params,
                     pars = c("u_hospital"),
                     pars_datagen = c("u_hospital", "sd_iid_hospital"),
                     n = seq(50, 1000, by = 200),
                     method = "is",
                     nsim = 1000, 
                     datagen_fn = utility_datagen_fn_is,
                     likelihood = utility_likelihood)


### Using the "trial_binary" built-in study design

# EVSI calculation using GAM regression.
evsi_builtin_rb <- evsi(outputs = chemotherapy_output,
                        inputs = m_params,
                        study = "trial_binary",
                        pars = c("p_side_effects_t1", 
                                 "p_side_effects_t2"),
                        n = seq(50, 500, by = 50),
                        method = "gam")

# EVSI calculation using Importance Sampling
evsi_builtin_is <- evsi(outputs = chemotherapy_output,
                        inputs = m_params,
                        study = "trial_binary",
                        pars = c("p_side_effects_t1", 
                                 "p_side_effects_t2"),
                        n = seq(50, 500, by = 50),
                        method = "is")


# Beta prior for standard care is set using the number of events
beta_params_t1 <- c(1 + n_side_effects,
                    1 + n_patients - n_side_effects)
# Beta prior for the novel intervention is approximated from the 
# mean andstandard deviation of the PA distribution for the 
# probability of side effects.
beta_params_t2 <- betaPar(mean(m_params$p_side_effects_t2),
                          sd(m_params$p_side_effects_t2))
# EVSI calculation with moment matching method
evsi_builtin_mm <- evsi(outputs = chemotherapy_output,
                     inputs = m_params,
                     study = "trial_binary",
                     pars = c("p_side_effects_t1", "p_side_effects_t2"),
                     n = seq(50, 500, by = 50),
                     method = "mm",
                     model_fn = calculate_costs_effects,
                     analysis_args = list(a1 = beta_params_t1[1],
                                          b1 = beta_params_t1[2],
                                          a2 = beta_params_t2$alpha,
                                          b2 = beta_params_t2$beta),
                     par_fn = generate_psa_parameters)

## Plotting
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
                pars = c("logor_side_effects"),
                pars_datagen = c("p_side_effects_t1", 
                                 "p_side_effects_t2"),
                n = seq(50, 500, by = 10),
                method = "gam",
                datagen_fn = OR_datagen_fn,
                par_fn = generate_psa_parameters)


# Using a bespoke analysis function - trial only updates odds ratio.
# Data generation function
OR_datagen_fn <- function(inputs, n = 500){
  p_side_effects_t1 <- inputs[, "p_side_effects_t1"]
  p_side_effects_t2 <- inputs[, "p_side_effects_t2"]
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
  
  LogOR_trial <- "
  model { 
    # Probability of side effects under treatment 1
    p_side_effects_t1 ~ dbeta(1 + n_side_effects, 
                              1 + n_patients - n_side_effects)
    
    # Log odds of side effects on treatment 2
    logor_side_effects ~ dnorm(logor_side_effects_mu, 
                               logor_side_effects_sd)
    # Odds of side effects on treatment 1
    odds_side_effects_t1 <- p_side_effects_t1 / (1 - p_side_effects_t1)
    # Odds for side effects on treatment 2
    odds_side_effects_t2 <- odds_side_effects_t1 * exp(logor_side_effects)
    
    # Probability of side effects under treatment 2
    p_side_effects_t2 <- 
      odds_side_effects_t2 / (1 + odds_side_effects_t2)
    
    X1 ~ dbin(p_side_effects_t1, n)
    X2 ~ dbin(p_side_effects_t2, n)
  }
  "
  filein <- file.path(tempdir(), fileext="datmodel.txt")
  cat(LogOR_trial,file=filein)
  
  # Perform the MCMC simulation with JAGS
  bugs.data <- jags(
    data =  data_jags,
    parameters.to.save = pars,
    model.file = filein, 
    n.chains = 3, 
    n.iter = args$n.iter, 
    n.thin = 1, 
    n.burnin = 250, 
    quiet=TRUE, progress.bar = "none") 
  
  # Resample treatment 1 from the prior as study will not update t1
  nsam <- length(bugs.data$BUGSoutput$sims.matrix[, pars[2]])
  resample_t1 <- rbeta(nsam,
                       1 + args$n_side_effects,
                       1 + args$n_patients - args$n_side_effects)
  resample_t2 <- bugs.data$BUGSoutput$sims.matrix[, pars[2]]
  return(data.frame(p_side_effects_t1 = resample_t1,
                    p_side_effects_t2 = resample_t2))
}

# EVSI calculation using the moment matching method.
analysis_args <- list(n_side_effects = n_side_effects,
                      n_patients = n_patients,
                      n = 500,
                      logor_side_effects_mu = logor_side_effects_mu,
                      logor_side_effects_sd = logor_side_effects_sd,
                      n.iter = 7500)
evsi_OR <- evsi(outputs = chemotherapy_output,
                inputs = m_params,
                pars = c("p_side_effects_t1", "p_side_effects_t2"),
                n = seq(50, 500, by = 10),
                method = "mm",
                datagen_fn = OR_datagen_fn,
                model_fn = calculate_costs_effects,
                analysis_args = analysis_args, 
                analysis_fn = OR_analysis_fn, 
                par_fn = generate_psa_parameters)
