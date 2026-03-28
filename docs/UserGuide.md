# ISO6976.2016 — User Guide

This package implements ISO 6976:2016 *"Natural Gas — Calculation of calorific
values, density, relative density and Wobbe indices from composition"*,
including the uncertainty propagation of Annex B.

---

## Installation

```r
# From CRAN (once published)
install.packages("ISO6976.2016")

# From GitHub
remotes::install_github("RuedigerForster/ISO_6976")
```

---

## Component ordering

The composition vector has **60 elements**, one per component, in the order of
ISO 6976:2016 Table A.2.  Use `componentNames()` to inspect the full list and
`componentIndex()` / `componentName()` to convert between names and indices.

```r
library(ISO6976.2016)

componentNames()[1:10]
# [1] "methane"    "ethane"     "propane"    "n-butane"   "isobutane"
# [6] "n-pentane"  "isopentane" "neopentane" "n-hexane"   "2-methylpentane"

componentIndex("nitrogen")       # 52
componentIndex("carbon dioxide") # 54
componentIndex("water")          # 42
componentName(41L)               # "hydrogen"
```

---

## Providing the composition

### Method A — bare vectors

Construct a length-60 vector directly.  All values must be in **mol/mol**
(not mol%).

```r
x   <- numeric(60)
u_x <- numeric(60)
r_x <- diag(60)          # identity: no inter-component correlations

x[componentIndex("methane")]        <- 0.93321
x[componentIndex("ethane")]         <- 0.02566
x[componentIndex("propane")]        <- 0.01537
x[componentIndex("nitrogen")]       <- 0.01035
x[componentIndex("carbon dioxide")] <- 0.01541

u_x[componentIndex("methane")]        <- 0.000346
u_x[componentIndex("ethane")]         <- 0.000243
u_x[componentIndex("propane")]        <- 0.000148
u_x[componentIndex("nitrogen")]       <- 0.000195
u_x[componentIndex("carbon dioxide")] <- 0.000111
```

### Method B — GasComponents class

`GasComponents` is an R6 container that accepts components by name or index
and enforces basic consistency on the correlation matrix.  It is useful when
building a composition incrementally, e.g. while reading a chromatograph
result file.

```r
gc <- GasComponents$new()

gc$setFraction("methane",        0.93321)
gc$setFraction("ethane",         0.02566)
gc$setFraction("propane",        0.01537)
gc$setFraction("nitrogen",       0.01035)
gc$setFraction("carbon dioxide", 0.01541)

gc$setUncertainty("methane",        0.000346)
gc$setUncertainty("ethane",         0.000243)
gc$setUncertainty("propane",        0.000148)
gc$setUncertainty("nitrogen",       0.000195)
gc$setUncertainty("carbon dioxide", 0.000111)

# Pass the three fields to calculateProperties:
res <- calculateProperties(
  gc$fractions, gc$uncertainties, gc$correlations,
  combustionTemperature = 25,
  volumeTemperature     = 0
)
```

### Method C — built-in Annex D datasets

The package ships with the four reference mixtures from Annex D.  Each
dataset is a named list with elements `fractionArray`, `uncertaintyArray`,
and `correlationMatrix`.

| Dataset      | Mixture | Correlation matrix |
|--------------|---------|--------------------|
| `example1`   | 5-component (Annex D.2) | identity |
| `example2`   | 11-component with water vapour (Annex D.3) | identity |
| `example3`   | 11-component (Annex D) | identity |
| `example3_ex`| same as `example3` | full GC calibration matrix |

```r
data("example1")
res <- calculateProperties(
  example1$fractionArray,
  example1$uncertaintyArray,
  example1$correlationMatrix,
  combustionTemperature = 15,
  volumeTemperature     = 15
)
```

---

## Correlation matrix

When a GC calibration delivers a full covariance matrix, convert it to a
correlation matrix and pass it as `correlationMatrix`.  All off-diagonal
elements must lie in (−1, 1).  Correlations are typically negative between
the dominant component (methane) and the heavier fractions.

```r
gc$setCorrelation("methane", "ethane",   -0.65)
gc$setCorrelation("methane", "propane",  -0.48)
gc$setCorrelation("ethane",  "propane",   0.31)

# Or replace the whole matrix at once:
gc$setCorrelationMatrix(R)   # R is a 60x60 matrix
```

`setCorrelation()` automatically enforces symmetry.

When correlations are unknown, pass `diag(60)`.  Uncertainties will then be
slightly conservative (identity gives larger u values than a realistic
negative-correlation matrix).

---

## Reference conditions

`calculateProperties()` accepts any combination of combustion temperature
and volume reference temperature permitted by ISO 6976:2016 §5.

| `combustionTemperature` | `volumeTemperature` | Typical use |
|------------------------|---------------------|-------------|
| 25 °C | 0 °C | Germany / DVGW G 685 |
| 25 °C | 15 °C | — |
| 15 °C | 15 °C | UK / legacy European |
| 15.55 °C | 15.55 °C | 60 °F / North American |
| 0 °C | 0 °C | — |

```r
# German standard
r_de <- calculateProperties(
  x, u_x, r_x,
  combustionTemperature = 25,
  volumeTemperature     = 0
)

# UK standard
r_uk <- calculateProperties(
  x, u_x, r_x,
  combustionTemperature = 15,
  volumeTemperature     = 15
)
```

The reference pressure defaults to 101.325 kPa and must remain within
90–110 kPa.

---

## Output

`calculateProperties()` returns a named list.  All calorific values and
Wobbe indices are in **MJ/m³** (volumetric) or **MJ/kg** (mass basis) or
**kJ/mol** (molar basis).  Density is in **kg/m³**, molar mass in
**kg/kmol**.

| Name | Description |
|------|-------------|
| `M` | Molar mass [kg/kmol] |
| `Z` | Compression factor [—] |
| `G_o`, `D_o` | Ideal-gas relative density [—], density [kg/m³] |
| `G`, `D` | Real-gas relative density, density |
| `Hcg`, `Hcn` | Molar gross / net CV [kJ/mol] |
| `Hmg`, `Hmn` | Mass-basis gross / net CV [MJ/kg] |
| `Hvg_o`, `Hvn_o` | Ideal-gas volumetric gross / net CV [MJ/m³] |
| `Hvg`, `Hvn` | Real-gas volumetric gross / net CV [MJ/m³] |
| `Wg_o`, `Wn_o` | Ideal-gas gross / net Wobbe index [MJ/m³] |
| `Wg`, `Wn` | Real-gas gross / net Wobbe index [MJ/m³] |

Every quantity that has an associated uncertainty carries a `u_` companion
(`u_G`, `u_D`, `u_Hcg`, `u_Hcn`, `u_Hmg`, `u_Hmn`, `u_Hvg_o`, `u_Hvn_o`,
`u_Hvg`, `u_Hvn`, `u_Wg`, `u_Wn`).

```r
cat("Hvg  =", round(res$Hvg,   4), "MJ/m³\n")
cat("u_Hvg=", round(res$u_Hvg, 5), "MJ/m³  (k = 1)\n")
cat("Wg   =", round(res$Wg,    4), "MJ/m³\n")
cat("u_Wg =", round(res$u_Wg,  5), "MJ/m³  (k = 1)\n")
```

---

## Uncertainties and coverage factor

All `u_` outputs are **standard uncertainties** (*k* = 1) by default,
propagated from `uncertaintyArray` and `correlationMatrix` via the
variance-covariance method of ISO 6976:2016 Annex B.

For **expanded uncertainties** multiply by the appropriate coverage factor.
Pass `coverage = 2` for approximately 95 % confidence (normal distribution):

```r
res_k2 <- calculateProperties(
  x, u_x, r_x,
  combustionTemperature = 25,
  volumeTemperature     = 0,
  coverage              = 2
)
# res_k2$u_Wg  is the expanded uncertainty U(Wg) at k = 2
```

Non-uncertainty outputs (`M`, `Z`, `Hvg`, `Wg`, …) are unaffected by the
coverage factor.

---

## Verification against ISO 6976:2016 Annex D

The four reference calculations in Annex D can be reproduced exactly.  The
table below shows Example 1 (Annex D.2, 15/15 °C).

```r
data("example1")
res <- calculateProperties(
  example1$fractionArray,
  example1$uncertaintyArray,
  example1$correlationMatrix,
  combustionTemperature = 15,
  volumeTemperature     = 15,
  coverage              = 1
)
```

| Property | Computed | ISO 6976:2016 |
|----------|:--------:|:-------------:|
| M [kg/kmol] | 17.3884301 | 17.3884301 |
| Z | 0.99776224 | 0.99776224 |
| Hcg [kJ/mol] | 906.1799588 | 906.1799588 |
| u(Hcg) | 0.615609872 | 0.615609872 |
| Hmg [MJ/kg] | 52.113961 | 52.113961 |
| u(Hmg) | 0.024301 | 0.024301 |
| Hvg [MJ/m³] | 38.410611 | 38.410611 |
| u(Hvg) | 0.026267 | 0.026267 |

---

## Common pitfalls

**Fractions in mol% instead of mol/mol.**  The function expects mol/mol
(sum ≈ 1).  Divide by 100 before passing if your GC software reports mol%.

```r
x_pct <- c(93.3212, 2.5656, ...)   # mol%
x     <- x_pct / 100               # mol/mol
```

**Z ≤ 0.9.**  The calculation is outside the application range of the
standard.  Check that fractions sum to 1 and that no component fraction is
unrealistically large.

**Water vapour.**  Water is component 42.  If the analysis is on a dry basis,
leave index 42 at zero.  If the gas is water-saturated, add the vapour
pressure fraction explicitly before calling `calculateProperties()`.

**Non-physical correlation matrix.**  A covariance matrix from GC calibration
must be converted to a correlation matrix (divide each element by the product
of the corresponding standard deviations) before passing it to the function.
`setCorrelationMatrix()` rejects any element outside [−1, 1] but does not
check positive semi-definiteness.

---

## Quick reference

```r
library(ISO6976.2016)

# Build composition
gc <- GasComponents$new()
gc$setFraction("methane",  0.9332)
gc$setFraction("nitrogen", 0.0668)
gc$setUncertainty("methane",  3.3e-4)
gc$setUncertainty("nitrogen", 1.9e-4)
gc$setCorrelation("methane", "nitrogen", -0.4)

# Calculate at 25/0 °C, expanded uncertainty k = 2
res <- calculateProperties(
  gc$fractions, gc$uncertainties, gc$correlations,
  combustionTemperature = 25,
  volumeTemperature     = 0,
  coverage              = 2
)

res$Hvg    # real-gas volumetric GCV [MJ/m³]
res$u_Hvg  # expanded uncertainty U(Hvg), k = 2
res$Wg     # real-gas gross Wobbe index [MJ/m³]
res$u_Wg   # expanded uncertainty U(Wg), k = 2
```
