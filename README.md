
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
#>              coef exp(coef) se(coef)      z Pr(>|z|)    
#> X1 (1)    0.34638   1.41394  0.05332  6.497 8.20e-11 ***
#> X1 (2)    0.93371   2.54394  0.10765  8.673  < 2e-16 ***
#> X1 (3)   -0.28992   0.74832  0.10673 -2.716   0.0066 ** 
#> X1 (4)    0.95170   2.59011  0.15661  6.077 1.22e-09 ***
#> X1 (5)   -0.21037   0.81028  0.09827 -2.141   0.0323 *  
#> X2 (1)    0.36122   1.43508  0.05483  6.588 4.46e-11 ***
#> X2 (2)    0.97637   2.65481  0.11734  8.321  < 2e-16 ***
#> X2 (3)   -0.24155   0.78541  0.11991 -2.014   0.0440 *  
#> X2 (4)    0.97474   2.65048  0.13296  7.331 2.29e-13 ***
#> X2 (5)   -0.26555   0.76678  0.11850 -2.241   0.0250 *  
#> X3 (ALL)  0.35856   1.43127  0.05650  6.347 2.20e-10 ***
#> X4 (ALL)  0.35181   1.42163  0.05574  6.312 2.76e-10 ***
#> X13 (4)   0.01383   1.01393  0.13008  0.106   0.9153    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>          exp(coef) exp(-coef) lower .95 upper .95
#> X1 (1)   1.4139    0.7072     1.2736    1.5697   
#> X1 (2)   2.5439    0.3931     2.0600    3.1415   
#> X1 (3)   0.7483    1.3363     0.6071    0.9224   
#> X1 (4)   2.5901    0.3861     1.9055    3.5207   
#> X1 (5)   0.8103    1.2341     0.6683    0.9824   
#> X2 (1)   1.4351    0.6968     1.2889    1.5979   
#> X2 (2)   2.6548    0.3767     2.1094    3.3413   
#> X2 (3)   0.7854    1.2732     0.6209    0.9935   
#> X2 (4)   2.6505    0.3773     2.0424    3.4396   
#> X2 (5)   0.7668    1.3042     0.6079    0.9673   
#> X3 (ALL) 1.4313    0.6987     1.2812    1.5989   
#> X4 (ALL) 1.4216    0.7034     1.2745    1.5857   
#> X13 (4)  1.0139    0.9863     0.7857    1.3084
```

We can also give the estimated cumulative hazard function as follows:

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
