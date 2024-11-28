
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

The goal of survtrans is to provide a framework for transferring
survival information from source domain(s) to target domain. The package
now only supports the Cox proportional hazards model with global and
local transfer learning.

## Installation

You can install the development version of survtrans like so:

``` r
# install.packages("pak")
pak::pak("SignorinoY/survtrans")
```

## Example

This is a basic example which shows you how to transfer survival
information from multiple source domains to a target domain using the
Cox proportional hazards model:

``` r
library(survtrans)
formula <- survival::Surv(time, status) ~ . - group - id
fit <- coxtrans(
  formula, sim2, sim2$group,
  lambda1 = 0.03, lambda2 = 0.01, lambda3 = 0.01, penalty = "SCAD"
)
summary(fit)
#> Call:
#> coxtrans(formula = formula, data = sim2, group = sim2$group, 
#>     lambda1 = 0.03, lambda2 = 0.01, lambda3 = 0.01, penalty = "SCAD")
#> 
#>   n=500, number of events=422
#> 
#>               coef exp(coef) se(coef)      z Pr(>|z|)    
#> X1 (1)     0.34791   1.41611  0.05130  6.781 1.19e-11 ***
#> X1 (2, 4)  0.94904   2.58322  0.08576 11.066  < 2e-16 ***
#> X1 (3, 5) -0.25321   0.77630  0.06888 -3.676 0.000237 ***
#> X2 (1)     0.36098   1.43474  0.05370  6.723 1.78e-11 ***
#> X2 (2, 4)  0.96454   2.62358  0.08449 11.415  < 2e-16 ***
#> X2 (3, 5) -0.24257   0.78461  0.08020 -3.025 0.002490 ** 
#> X3 (ALL)   0.34501   1.41200  0.05541  6.227 4.76e-10 ***
#> X4 (ALL)   0.32763   1.38767  0.05335  6.141 8.22e-10 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>           exp(coef) exp(-coef) lower .95 upper .95
#> X1 (1)    1.4161    0.7062     1.2806    1.5659   
#> X1 (2, 4) 2.5832    0.3871     2.1835    3.0561   
#> X1 (3, 5) 0.7763    1.2882     0.6783    0.8885   
#> X2 (1)    1.4347    0.6970     1.2914    1.5940   
#> X2 (2, 4) 2.6236    0.3812     2.2232    3.0961   
#> X2 (3, 5) 0.7846    1.2745     0.6705    0.9182   
#> X3 (ALL)  1.4120    0.7082     1.2667    1.5740   
#> X4 (ALL)  1.3877    0.7206     1.2499    1.5406
```

We can also give the estimated cumulative hazard function as follows:

``` r
library(ggplot2)
basehaz_pred <- basehaz(fit)
basehaz_pred$color <- as.numeric(basehaz_pred$strata) %% 2
ggplot(basehaz_pred, aes(x = time, y = basehaz, group = strata, color = factor(color))) +
  geom_line() +
  geom_line(aes(x = time, y = time^2 / 2, color = "True"), linetype = "dotted") +
  geom_line(aes(x = time, y = time^3 / 3, color = "True"), linetype = "dotted") +
  labs(
    title = "Cumulative Baseline Hazard Function (Esimtated vs. True)",
    x = expression(t),
    y = expression(Lambda[0](t))
  ) +
  guides(color = guide_legend(title = "Strata"))
```

<img src="man/figures/README-basehaz-1.png" alt="Estimated vs. True Cumulative Baseline Hazard Function" width="100%" />
