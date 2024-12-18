
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dmpredict

## Overview

dmpredict is a package written to implement the prediction equations
reported in `NCD Risk Factor Collaboration (2023)` to estimate the
probability that an individual had HbA1c ≥6.5% based on measured fasting
plasma glucose (FPG), and vice versa. The details of the statistical
methods are provided in the citation below. In short, the prediction
equations are generalised linear mixed models fitted with package
`lme4`. The model objects in this package provide draws of the model
coefficients to allow estimation of the uncertainty associated with the
prediction equations.

NCD Risk Factor Collaboration (NCD-RisC). Global variation in diabetes
diagnosis and prevalence based on fasting glucose and hemoglobin A1c.
*Nature Medicine* 2023; **29**: 2885-2901.
<https://doi.org/10.1038/s41591-023-02610-2>

See also <https://zenodo.org/records/8169146>.

## Usage

The users must provide at least one of the following:

- `fpg_clean` cleaned FPG (mmol/L)
- `hba1c_clean` cleaned HbA1c (%)

Plus the following required predictors:

- `sex` coded as “female” and “male”
- `age` in years (current models have only been tested in adults aged
  18+ years)
- `study_id` unique identifier for each study (or study round)
- `iso` [ISO alpha-3
  code](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3) for the
  country the study was conducted in
- `medication` whether using medication for diabetes, coded as 1 = yes,
  0 = no, NA = missing

Optional covariates for prediction:

- `bmi_clean` cleaned body-mass index (kg/m<sup>2</sup>)
- `handheld_fpg` how FPG was measured, coded as 1 = using a POC handheld
  device, 0 = lab based measurement, NA = unknown
- `handheld_a1c` how HbA1c was measured, coded as 1 = using a POC
  handheld device, 0 = lab based measurement, NA = unknown

If `iso` is not recognised, provide a string with one of the following
in a variable named `Superregion`:

- “Central and Eastern Europe”
- “Central Asia, Middle East and North Africa”
- “High-income western”
- “Latin America and Caribbean”
- “Oceania”
- “Southeast Asia, East Asia and the Pacific”
- “South Asia”
- “Sub-Saharan Africa”

## Installation

To install the package, use:

``` r
remotes::install_github("NCD-RisC/dmpredict")
```

## Example

``` r
library(dmpredict)

# Make sure 'data' object has the minimum columns as stated above
# data <- readRDS("input_data.RDS")
df <- transform_data(data)

# Predict_missing_fpg or hba1c functions return a probability that the predicted biomarker
# is over the diagnostic threshold for diabetes, i.e. FPG >=7.0 mmol/L and HbA1c >=6.5%

data$predicted_a1c <- predict_missing_hba1c(df, cw_model_coefs, return_draws = FALSE, verbose = FALSE)
data$predicted_fpg <- predict_missing_fpg(df, cw_model_coefs, return_draws = FALSE, verbose = FALSE)

# using 'return_draws = TRUE' will return the N draws (default = 2000) of predictions in a matrix

res1 <- predict_missing_hba1c(df, cw_model_coefs, return_draws = TRUE, verbose = FALSE)
res2 <- predict_missing_fpg(df, cw_model_coefs, return_draws = TRUE, verbose = FALSE)

data$predicted_a1c2 <- res1$predicted_a1c
data$predicted_fpg2 <- res2$predicted_fpg
print(dim(res1$mat_predicted_a1c))
print(dim(res2$mat_predicted_fpg))

# using a different number of simulation draws

res3 <- predict_missing_hba1c(df, cw_model_coefs, N_sim = 5000, return_draws = TRUE, verbose = FALSE)
res4 <- predict_missing_fpg(df, cw_model_coefs, N_sim = 5000, return_draws = TRUE, verbose = FALSE)

print(dim(res3$mat_predicted_a1c))
print(dim(res4$mat_predicted_fpg))
```
