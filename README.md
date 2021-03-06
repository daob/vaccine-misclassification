# Misclassification models

### How to use

1. Clone the reposity
2. Open the Rproj file in Rstudio (for example by double-clicking it)
3. Open the `R/simulation.R` file
4. Run it

### Description of setup

* _Simulates_ data with the following setup: 
    - Outcome is "event" (0/1)
    - Outcome is misclassified. So there is an unobserved true outcome, `y`, and an observed, potentially error-prone, outcome `ystar`
    - Predictor `x` is "vaccine" (0/1)
    - Everything, including misclassfication table, is parameterized as a logistic regression. So parameters in these tables are log-linear.
* Defines _estimators_ you might apply to such data. Currently the following are implemented:
    - "Naive", pretending there is no misclassfication: a logistic regression of oobserved `ystar` on `x`. And the relative risk calculated on the observed `ystar`
    - "MLE": this formulates the (true) model, including misclassification. This model is estimated by maximizing the marginal log-likelihood directly. 
    - "EM": this formulates the same (true) model as MLE, but uses a different estimation method (EM), which is usually more stable.  
* Defines and runs a _simulation study_ (experiment) in which the following factors are varieed:
    - Sample size `n`
    - The sensitivity and specificity of the event registration (reworked into loglinear paramters `tau` and `lambda`)
    - The overall base event rate, parameterized with logit parameter `alpha`
    - The true effect size, with logistic regression coefficient `beta` (fixed at 0.2 for the moment)
    - The overall vaccination rate (fixed at 0.7 for the moment)

At the moment, only nondifferential error is examined.

### Results so far

"Obviously", the MLE is unbiased, especially as `n` increases and/or `alpha` gets closer to zero (more events observed). However, in terms of mean-squared-error, it is almost never worthwhile to use the MLE. This is because, in this simple nondifferential setup, only the specificity, which we can assume to be excellent, is relevant to bias in `beta`. So there is little to no payoff for trading unbiasedness for the considerably larger variance in the MSE relative to the naive estimator. 

![Mean absolute error illustration](https://github.com/daob/vaccine-misclassification/blob/dcbe1cff67edcddcac11420ae1a4f42bb9b2ae23/images/mean_absolute_error_big.pdf)

The same results may not hold for relative risk (not examined yet).

The same results will not hold for differential error.