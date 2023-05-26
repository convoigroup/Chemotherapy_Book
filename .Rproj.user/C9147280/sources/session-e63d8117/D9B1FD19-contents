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

### Developing a bespoke study design
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

## Add the individual level standard deviation
m_params$sd_iid_hospital <- runif(nrow(m_params), 0.00001, 0.4)

## Regression Based Method
evsi_utility <- evsi(outputs = chemotherapy_output,
                     inputs = m_params,
                     pars = c("u_hospital"),
                     n = seq(50, 1000, by = 200),
                     method = "gam",
                     datagen_fn = utility_datagen_fn_agg)

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
                     datagen_fn = utility_datagen_fn_indiv,
                     likelihood = utility_likelihood)
