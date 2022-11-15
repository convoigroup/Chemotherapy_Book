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

#### Randomised Trial for Log-Odds ratio ####
beta_params_t1 <- c(1 + n_side_effects, 
                    1 + n_patients - n_side_effects)
beta_params_t2 <- betaPar(mean(m_params$p_side_effects_t2),
        sd(m_params$p_side_effects_t2))

evsi_default <- evsi(outputs = chemotherapy_output,
                     inputs = m_params,
                     study = "trial_binary",
                     pars = c("p_side_effects_t1", "p_side_effects_t2"),
                     n = seq(500, 1500, by = 20),
                     method = "mm",
                     model_fn = calculate_costs_effects,
                     analysis_args = list(a1 = beta_params_t1[1],
                                          b1 = beta_params_t1[2],
                                          a2 = beta_params_t2$alpha,
                                          b2 = beta_params_t2$beta),
                     par_fn = generate_psa_parameters)

debug(evsi)
### Written Out ###
OR_datagen_fn <- function(inputs, n = 500){
  p_side_effects_t1 <- inputs[, "p_side_effects_t1"]
  p_side_effects_t2 <- inputs[, "p_side_effects_t2"]
  X1 <- rbinom(length(p_side_effects_t1), n, p_side_effects_t1)
  X2 <- rbinom(length(p_side_effects_t2), n, p_side_effects_t2)
  data_save <- data.frame(X1 = X1, X2 = X2)
  return(data_save)
}

### HAD TO LEAVE OUT THE n ARGS - NEED TO BE ABLE TO HAVE n IN CHECK_ANALYSIS_FUN.
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
    n.chains = 3, 
    n.iter = args$n.iter, 
    n.thin = 1, 
    n.burnin = 250, progress.bar = "none") 
  
  resample_t1 <- rbeta(length(bugs.data$BUGSoutput$sims.matrix[, pars[2]]),
                       1 + args$n_side_effects,
                       1 + args$n_patients - args$n_side_effects)
  
  return(data.frame(p_side_effects_t1 = resample_t1,
         p_side_effects_t2 = bugs.data$BUGSoutput$sims.matrix[, pars[2]]))
}




evsi_OR <- evsi(outputs = chemotherapy_output,
                inputs = m_params,
                pars = c("p_side_effects_t1", "p_side_effects_t2"),
                n = seq(500, 1500, by = 20),
                method = "mm",
                datagen_fn = OR_datagen_fn,
                model_fn = calculate_costs_effects,
                analysis_args = list(n_side_effects = n_side_effects,
                                     n_patients = n_patients,
                                     logor_side_effects_mu = logor_side_effects_mu,
                                     logor_side_effects_sd = logor_side_effects_sd,
                                     n.iter = 5000),
                analysis_fn = OR_analysis_fn, 
                par_fn = generate_psa_parameters)


# Plot
EVPI <- evpi(chemotherapy_output)
plot(EVPI,
     xlab = "Willingness-to-Pay",
     ylab = "EVPI",
     main = "Expected Value of Perfect Information",
     type = "l")
points(chemotherapy_output$k, evsi_OR$evppi, type = "l", lty = 2)
x <- evsi_default %>% filter(n == 760) %>% select(evsi)
points(chemotherapy_output$k, x$evsi, type = "l", lty = 4, col = "red")
points(chemotherapy_output$k, evsi_OR$evsi, type = "l", lty = 3)
