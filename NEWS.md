# ISO_6976 News

## ISO_6976 0.1-0 (2026-03-28)

Initial release — clean-room implementation of ISO 6976:2016 for CRAN.

### Scope

Implements all properties and uncertainties defined in ISO 6976:2016:

* Compression factor Z (§5, Eq. 1)
* Molar mass M (Eq. 5)
* Ideal-gas and real-gas gross and net calorific values on molar, mass, and
  volumetric bases (Eqs. 2–12)
* Ideal-gas and real-gas density and relative density (Eqs. 13–19)
* Ideal-gas and real-gas gross and net Wobbe indices (Eqs. 15–16, 20–21)
* Standard uncertainties of all properties (Annex B, Eqs. B.4–B.14)

### Design decisions

* Single Rcpp back-end (`iso6976_calc`) computes all properties in one call.
* `GasComponents` R6 class provides named component access without
  normalization (normalization is outside the scope of ISO 6976:2016).
* Component names are in English following ISO 6976:2016 Table A.2.
* Reference data sets `example1`, `example2`, `example3`, `example3_ex` are
  the worked examples from ISO 6976:2016 Annex D.
