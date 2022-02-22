library(tidyverse)

# Load estimators and data synthesis functions:
source("R/data-synthesis.R")
source("R/estimators.R") # Requires maxLik package

# For reproducibility of the simulations
set.seed(342)

# The number of simulation replicates:
nsim <- 500 

### Defining simulation conditions
conditions <- expand.grid(
  # Sample size:
  n = c(200, 1000, 2000, 10000),
  # Proportion of vaccinated:
  p_vaccinated = 0.70,
  # True logistic regression coefficient of vaccination on event:
  beta_true = c(0, 0.05, 0.1, 0.2),
  # True logistic intercept of event (e.g. myocarditis):
  alpha_true = c(-4.5, -3, -2), #c(-9, -7, -4.5, -2),
  # Sensitivity Pr(observed = 1 | true = 1) of event registration:
  p_true_positive = c(0.99,  0.9, 0.8, 0.7, 0.6),
  # Specificity Pr(observed = 0 | true = 0) of event registration:
  p_true_negative = c(0.99,  0.8)
) %>% tibble

# For convenience, the relationship between true and observed event is
#    parameterized as a GLM with an intercept tau and slope ("loading") lambda:
#
# So, Pr(observed = 1 | true = 1) == plogis(tau + lambda)
#    and Pr(observed = 1 | true = 0) == plogis(tau)
conditions$tau <-
  with(conditions, log(1 - p_true_negative) - log(p_true_negative))

conditions$lambda <-
  with(conditions, log(p_true_positive) - log(1 - p_true_positive) - tau)


### Running the simulation

# Function that takes one condition, simulations nsim datasets, and runs all 
#   models on each
run_condition <- function(cond, nsim = 100) {
  
  # The estimator depends on the condition, because it depends 
  #   on the fixed classification error probabilities (through tau & lambda)
  estimator_misclassification_maxlik <-
    get_estimator_misclassification_maxlik(tau = cond$tau, lambda = cond$lambda)
  
  # Loop over simulation replicates
  map_df(1:nsim, function(isim) {
    
    dat_samp <- simulate_dataset(cond)
    
    # Produce a few estimates using various methods:
    res_mle <- estimator_misclassification_maxlik(dat_samp)
    res_naive <- estimator_naive(dat_samp)
    
    # Return the results together
    bind_cols(res_naive, res_mle)
  })
  
}


# Loop over all possible conditions, run simulation replicates 

# Define a dummy result that should be returned when one of the 
#   the conditions fails due to an error thrown
dummy_result <- run_condition(conditions[218,], nsim = 1)
dummy_result[1, ] <- NA

# For the moment I limit the simulation to look at true effect size = 0.2
#  Just comment out the next line to go through all conditions
conditions_subset <- conditions %>% filter(beta_true == 0.2)

# Run the simulation only if there are no previously saved results:
if(file.exists("data/res_overall.rdata")) {
  load(file = "data/res_overall.rdata")
} else { 
  system.time({ # Keep track of time
    # res_overall will contain the results of all replications & conditions
    res_overall <- 
      map_df(1:nrow(conditions_subset), 
             function(icond) {
               cat(sprintf("%d\n", icond)) # A little progress indicator
               current_condition <- conditions_subset[icond, ] 
               res <- tryCatch(run_condition(current_condition, nsim = nsim),
                               # If any error is thrown during estimation, 
                               #   the condition will return the dummy result:
                               error = function(e) { dummy_result })
               # Return the resulting estimates together with the condition:
               res %>% bind_cols(current_condition) 
             })
  })
  # For convenience in analyzing results, save the results for later
  save(res_overall, file = "res_overall.rdata")
}

# Make sure results are grouped by condition:
res_overall <- res_overall %>%
  group_by(n,
           p_vaccinated,
           beta_true,
           alpha_true,
           p_true_negative,
           p_true_positive) 

# Calculate MSE and absolute error of the estimates
absolute_error <- res_overall %>% summarize(
  mae_naive = median(abs(est_beta_glm_naive - beta_true)),
  mae_MLE = median(abs(est_beta_MLE - beta_true)),
  rmse_naive = sqrt(mean( (est_beta_glm_naive - beta_true)^2 )),
  rmse_MLE = sqrt(mean((est_beta_MLE - beta_true)^2))
)

# Plot mean absolute error (lower is better):
absolute_error %>% 
  pivot_longer(starts_with("mae")) %>%
  mutate(Estimator = gsub("mae_", "", name),
         `Mean absolute error in β` = value) %>%
  ggplot(aes(n, `Mean absolute error in β`, color = Estimator)) +
  scale_y_continuous(trans='log10')  +
  facet_grid(alpha_true ~ p_true_negative + p_true_positive, scales = "free_y") +
  ggthemes::theme_few() +
  ggthemes::scale_color_colorblind() +
  theme(axis.text.x = element_text(angle = 90, size = 8, vjust = 1, hjust=1),
        legend.position="top") + 
  geom_point() + geom_line() +
  ggtitle("Median absolute error in estimates of β (lower is better)")

# ggsave("images/mean_absolute_error_big.pdf", height = 5, width = 10)


# Mean of estimates over replications (per condition)
res_overall %>% 
  summarize(across(starts_with("est"), mean, na.rm = TRUE)) %>%
  ungroup %>% select(n, starts_with("est")) 

# Median of estimates over replications (per condition)
ests_median <- res_overall %>%
  summarize(across(starts_with("est"), median, na.rm = TRUE))

# Plot bias in estimates by condition (closer to zero is better):
ests_median %>% 
  pivot_longer(c("est_beta_glm_naive", "est_beta_MLE")) %>%
  mutate(Estimator = gsub("est_beta.*_", "", name)) %>%
  mutate(
    bias = value - beta_true,
    abs_bias = abs(bias),
    pct_relbias = 100 * abs_bias / beta_true) %>%
  ggplot(aes(n, bias, color = Estimator)) + 
  scale_x_continuous(breaks = unique(ests_median$n), trans = "sqrt") +
  facet_grid(alpha_true ~   p_true_negative + p_true_positive) +
  ggthemes::theme_few() +
  ggthemes::scale_color_colorblind() +
  theme(axis.text.x = element_text(angle = 90, size = 8),
        legend.position="top") + 
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() + geom_line() +
  ggtitle("Bias in estimates of β (closer to zero line is better)")


# ggsave("images/bias_more.pdf", height = 5, width = 10)

# Simulation standard deviations and average estimated standard errors
res_sd <- res_overall %>%
  summarize(across(starts_with("est"), mad, na.rm = TRUE)) %>%
  ungroup %>% 
  bind_cols(
    res_overall %>%
      summarize(across(starts_with("se"), median, na.rm = TRUE)) %>%
      ungroup %>% select(starts_with("se"))
  )

# Plot relative "efficiency" of estimators 
#   (but WARNING: they're not necessarily unbiased!)
res_sd %>% 
  mutate(`Relative efficiency` = est_beta_MLE/est_beta_glm_naive) %>%
  ggplot(aes(n, `Relative efficiency`, color = factor(alpha_true))) + 
  scale_x_continuous(breaks = unique(ests_median$n), trans = "sqrt") +
  scale_y_continuous(trans = "log10") +
  facet_grid(p_true_negative ~ p_true_positive, scales = "free_y") +
  ggthemes::theme_few() +
  ggthemes::scale_color_colorblind() +
  theme(axis.text.x = element_text(angle = 90, size = 8),
        legend.position="top") + 
  geom_hline(yintercept = 1, linetype = 2) +
  geom_point() + geom_line() + 
  ggtitle("How much more variable is the MLE?")

# ggsave("images/efficiency_more.pdf", height = 5, width = 10)


# Plot accuracy of MLE theoretical standard errors as measure of simulation sd's
res_sd %>% 
  mutate(`se:sd ratio` = se_beta_MLE / est_beta_MLE) %>%
  ggplot(aes(n, `se:sd ratio`, color = factor(alpha_true))) + 
  scale_x_continuous(breaks = unique(ests_median$n), trans = "sqrt") +
  scale_y_continuous(trans = "log10") +
  facet_grid(p_true_negative ~ p_true_positive, scales = "free_y") +
  ggthemes::theme_few() +
  ggthemes::scale_color_colorblind() +
  theme(axis.text.x = element_text(angle = 90, size = 8),
        legend.position="top") + 
  geom_hline(yintercept = 1, linetype = 2) +
  geom_point() + geom_line() + 
  ggtitle("How accurate are MLE-estimated standard errors?")

# ggsave("images/accuracy_ses_more.pdf", height = 5, width = 10)
