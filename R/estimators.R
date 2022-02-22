require(tidyverse)
require(maxLik)



#### Estimators (in: data, out: est. of event association with vaccination z)

# The "naive" estimator, in which we pretend y_star == y (no misclassfication):
estimator_naive <- function(dat) {
  fit_glm <- glm(y_star ~ z, data = dat, family = binomial)
  
  relative_risk <- with(dat, {
    tapply(y_star, z, mean) %>% 
      log %>% diff %>% exp
  })
  
  tibble(est_alpha_glm_naive = coef(fit_glm)['(Intercept)'],
         est_beta_glm_naive = coef(fit_glm)['z'],
         est_rr_naive = relative_risk)
}

# A maximum-likelihood estimator based on the correct misclassification model:
# (The function get_* returns the estimator, based on a specific tau and lambda.)
get_estimator_misclassification_maxlik <- function(tau, lambda) {
  function(dat) {
    
    # Define a log-likelihood function to optimize later
    loglik <- function(params) {
      
      alpha <- params['alpha'] # Intercept
      beta <- params['beta'] # Slope of z -> y
      
      y_star <- dat$y_star # Observed event
      z <- dat$z # Observed vaccination
    
      # Probability that (unobserved) true event happened:
      pi_y <- plogis(alpha + beta * z)
      
      # Marginal prob. of observing y_star == 1:
      p_y_star <- pi_y * plogis(tau + lambda) + (1 - pi_y) * plogis(tau)
      
      # Casewise log-likelihoods:
      logp <- dbinom(y_star, 1, prob = p_y_star, log = TRUE)
      
      return(logp) # maxLik can work with casewise likelihoods
    }
    
    # Use the maxLik package to find some estimates
    res <- maxLik(loglik, 
                  start = c('alpha' = 0, 'beta' = 0), 
                  method = "BHHH")
    
    # Output results in a nice tibble
    tibble(est_alpha_MLE = coef(res)['alpha'],
           est_beta_MLE = coef(res)['beta'],
           se_alpha_MLE = stdEr(res)['alpha'],
           se_beta_MLE = stdEr(res)['beta'])
  }
}

estimator_misclassification_EM <- function(dat, condition,
                                           max_iterations = 200, 
                                           params = NULL,
                                           verbose = TRUE,
                                           standard_errors = FALSE) {
  
  if(is.null(params)) params <- c("alpha" = rnorm(1, sd = 0.7), 
                                  "beta" = rnorm(1, sd = 0.7))
  
  tau <- condition$tau 
  lambda <- condition$lambda
  
  # Work with frequency data (TODO: this can be much more 
  #  efficient by generating data this way in the first place)
  dat_freq <- xtabs(~., dat) %>% 
    as.data.frame() %>% mutate_all(function(x) as.numeric(x) - 1)
  dat_freq <- bind_rows(dat_freq, dat_freq)
  dat_freq$X <- rep(c(0,1), each = nrow(dat_freq)/2)

  # Define a log-likelihood function to optimize later
  E_step <- function(params) {

    alpha <- params['alpha'] # Intercept
    beta <- params['beta'] # Slope of z -> y

    with(dat_freq, {

      P_X.Z <- plogis(alpha + beta * z)
      P_X.Z <- ifelse(X == 0, (1 - P_X.Z), P_X.Z)

      P_Y.X <- plogis(tau + lambda * X)
      P_Y.X <- ifelse(y_star == 0, (1 - P_Y.X), P_Y.X)

      P_YX.Z <- P_Y.X * P_X.Z
      P_Y.Z <- tapply(P_YX.Z, interaction(y_star, z), sum)
      dim(P_Y.Z) <- NULL

      # Posterior
      P_X.YZ <- P_YX.Z / P_Y.Z

      return(P_X.YZ)
    })
  }

  M_step <- function(posterior) {

    fit_structural <- 
      glm(X ~ z,
          weights = round(Freq * posterior),
          family = "binomial", data = dat_freq)
    
    coef(fit_structural)
  }
  
  log_likelihood <- function(params) {
    dat_freq_sub <- dat_freq %>% filter(X == 0)
    
    alpha <- params['alpha'] # Intercept
    beta <- params['beta'] # Slope of z -> y
    
    y_star <- dat_freq_sub$y_star # Observed event
    z <- dat_freq_sub$z # Observed vaccination
    
    # Probability that (unobserved) true event happened:
    pi_y <- plogis(alpha + beta * z)
    
    # Marginal prob. of observing y_star == 1:
    p_y_star <- pi_y * plogis(tau + lambda) + (1 - pi_y) * plogis(tau)
    
    # Casewise log-likelihoods:
    logp <- dbinom(y_star, 1, prob = p_y_star, log = TRUE)
    
    sum(dat_freq_sub$Freq * logp)
  }

  loglik_previous <- log_likelihood(params)
  converged <- FALSE
  
  for(i in 1:max_iterations) {
    # E_step
    posterior <- E_step(params)

    # M_step
    params[] <- M_step(posterior)
    
    loglik_current <- log_likelihood(params)
    
    loglik_increase <- (loglik_current - loglik_previous)
    if(verbose)
      cat(sprintf("%d\tloglik increase:\t%5.6f\n", i, loglik_increase))
    
    if(abs(loglik_increase) < 1e-5) {
      if(verbose)
        cat(sprintf("Reached convergence at loglik %5.6f.\n", loglik_current))
      converged <- TRUE
      break
    }
    
    loglik_previous <- loglik_current
  }
  if(!converged) warning("EM did not reach convergence criterion.\n")
  
  # Output results in a nice tibble
  output <- tibble(est_alpha_EM = params['alpha'], 
         est_beta_EM = params['beta'])

  if(standard_errors) {
    H <- numDeriv::hessian(log_likelihood, params)
    ses <- sqrt(diag(solve(-H)))
    output$se_alpha_EM <- ses[1]
    output$se_beta_EM <- ses[2]
  }
  
  output
}

