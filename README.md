
<!-- README.md is generated from README.Rmd. Please edit that file -->

# survtrans

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/survtrans)](https://CRAN.R-project.org/package=survtrans)
[![R-CMD-check](https://github.com/SignorinoY/survtrans/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SignorinoY/survtrans/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/SignorinoY/survtrans/graph/badge.svg)](https://app.codecov.io/gh/SignorinoY/survtrans)
<!-- badges: end -->

The goal of **survtrans** is to provide a framework for transferring
survival information from source domain(s) to a target domain. The
package currently supports the Cox proportional hazards model with both
global and local transfer learning.

## Installation

You can install the development version of **survtrans** with:

``` r
# install.packages("pak")
pak::pak("SignorinoY/survtrans")
```

## Example

This is a basic example showing how to transfer survival information
from multiple source domains to a target domain using the Cox
proportional hazards model:

``` r
library(survtrans)
formula <- survival::Surv(time, status) ~ . - group - id
fit <- coxtrans(
  formula, sim2, sim2$group, 1,
  lambda1 = 0.075, lambda2 = 0.04, lambda3 = 0.04, penalty = "SCAD"
)
summary(fit)
#> Call:
#> coxtrans(formula = formula, data = sim2, group = sim2$group, 
#>     target = 1, lambda1 = 0.075, lambda2 = 0.04, lambda3 = 0.04, 
#>     penalty = "SCAD")
#> 
#>   n=500, number of events=422
#> 
#>       coef exp(coef) se(coef)     z Pr(>|z|)    
#> X1 0.35052   1.41981  0.05366 6.533 6.47e-11 ***
#> X2 0.35914   1.43209  0.05424 6.621 3.56e-11 ***
#> X3 0.34485   1.41178  0.05431 6.350 2.15e-10 ***
#> X4 0.32870   1.38916  0.05190 6.333 2.40e-10 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>    exp(coef) exp(-coef) lower .95 upper .95
#> X1 1.4198    0.7043     1.2781    1.5773   
#> X2 1.4321    0.6983     1.2877    1.5927   
#> X3 1.4118    0.7083     1.2692    1.5703   
#> X4 1.3892    0.7199     1.2548    1.5379
```

You can also plot the estimated cumulative baseline hazard function and
compare it to the true function:

``` r
library(ggplot2)
basehaz_pred <- basehaz(fit)
basehaz_pred$color <- ifelse(
  as.numeric(basehaz_pred$strata) %% 2 == 0, "Group 2", "Group 1"
)
ggplot(
  basehaz_pred,
  aes(
    x = time,
    y = basehaz,
    group = strata,
    color = factor(color),
    linetype = "Estimates"
  )
) +
  geom_line() +
  geom_line(
    aes(x = time, y = time^2 / 2, color = "Group 1", linetype = "True")
  ) +
  geom_line(
    aes(x = time, y = time^3 / 3, color = "Group 2", linetype = "True")
  ) +
  labs(
    title = "Cumulative Baseline Hazard Function (Estimated vs. True)",
    x = expression(t),
    y = expression(Lambda[0](t))
  ) +
  scale_linetype_manual(values = c("Estimates" = "dashed", "True" = "solid")) +
  guides(
    color = guide_legend(title = "Strata"),
    linetype = guide_legend(title = "Type")
  )
```

<img src="man/figures/README-basehaz-1.png" alt="Estimated vs. True Cumulative Baseline Hazard Function" width="100%" />
