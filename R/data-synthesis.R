require(tidyverse)

#### Synthetic data generation 

simulate_dataset <- function(condition) {
  with(condition, {
    
    # z is an indicator for whether the person is vaccinated (1) or not (0)
    z <- rbinom(n, 1, prob = p_vaccinated)
    
    # Ignore covariates x for now
    # x <- rbinom(n, 1, prob = plogis(z - 0.5))
    
    # True linear predictor in GLM for event (e.g. myocarditis)
    eta_y_true <- alpha_true + beta_true * z
    
    # y is an indicator for the _unobserved_ true event (1) or not (0)
    y <- rbinom(n, 1, prob = plogis(eta_y_true))
    
    # summary(glm(y ~ z, family = binomial)) # CHECK
    
    # y_star is the _observed_ event (1) or not (0)
    # The misclassifcation process is parameterized as a GLM for y_star given y:
    eta_ystar_true <- tau + lambda * y
    y_star <- rbinom(n, 1, prob = plogis(eta_ystar_true))
    
    # You can check that as n increases, this gives back the TP and TN rates:
    # table(y_star, y) %>% prop.table(2)
    
    # Put the observed data together in a data_frame:
    dat_samp <- tibble(y_star = y_star, z = z)
    
    dat_samp
  })
}
