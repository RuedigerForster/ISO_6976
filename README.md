# ISO6976.2016

R implementation of **ISO 6976:2016** — *"Natural Gas — Calculation of
calorific values, density, relative density and Wobbe indices from
composition"*.

ISO 6976:2016 defines the method for deriving the thermodynamic and
combustion properties of a natural gas mixture from its molar composition,
together with the propagation of measurement uncertainties.  It is the
reference standard for custody transfer, billing, and energy measurement
in the natural gas industry.

## Installation

```r
# From CRAN
install.packages("ISO6976.2016")

# Development version from GitHub
remotes::install_github("RuedigerForster/ISO_6976")
```

## Documentation

- `?calculateProperties` — function reference
- `?GasComponents` — R6 class reference

## Status

[![R CMD check](https://github.com/RuedigerForster/ISO_6976/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/RuedigerForster/ISO_6976/actions/workflows/R-CMD-check.yaml)
