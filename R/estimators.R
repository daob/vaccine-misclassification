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

