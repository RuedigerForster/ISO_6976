# Tests against the worked examples in Annex D of ISO 6976:2016
# All expected values are taken directly from the standard.

library(ISO6976.2016)

################################################################################
# Example 1 (Annex D.2): 5-component mixture, 15/15 °C
################################################################################
test_that("Example 1 — 15/15 °C: molar mass, Z, molar GCV, mass GCV, vol. GCV", {
  data("example1", envir = environment())
  res <- calculateProperties(example1$fractionArray, example1$uncertaintyArray, example1$correlationMatrix,
                             combustionTemperature = 15,
                             volumeTemperature = 15,
                             coverage = 1)

  expect_equal(res$M,       17.3884301, tolerance = 1e-7)   # molar mass
  expect_equal(res$Z,        0.99776224, tolerance = 1e-8)  # compression factor
  expect_equal(res$Hcg,    906.1799588, tolerance = 1e-7)   # molar gross CV
  expect_equal(res$u_Hcg,    0.615609872, tolerance = 1e-9) # u(molar GCV)
  expect_equal(res$Hmg,     52.113961,  tolerance = 1e-6)   # mass gross CV
  expect_equal(res$u_Hmg,    0.024301,  tolerance = 1e-6)   # u(mass GCV)
  expect_equal(res$Hvg,     38.410611,  tolerance = 1e-6)   # real-gas vol. GCV
  expect_equal(res$u_Hvg,    0.026267,  tolerance = 1e-6)   # u(vol. GCV)
})

################################################################################
# Example 2 (Annex D.3): mixture with water vapour, 15.55/15.55 °C
# Known mismatch vs. standard for Z: -1.78e-5
################################################################################
test_that("Example 2 — 15.55/15.55 °C (water vapour)", {
  data("example2", envir = environment())
  res <- calculateProperties(example2$fractionArray, example2$uncertaintyArray, example2$correlationMatrix,
                             combustionTemperature = 15.55,
                             volumeTemperature = 15.55,
                             coverage = 1)

  expect_equal(res$M,       16.9891697,   tolerance = 1e-7)
  expect_equal(res$Z,        0.9975690,   tolerance = 1e-7)   # mismatch -1.78e-5
  expect_equal(res$Hcg,    871.443916,    tolerance = 1e-7)
  expect_equal(res$u_Hcg,    0.522493911, tolerance = 1e-9)
  expect_equal(res$Hmg,     51.294085,    tolerance = 1e-6)
  expect_equal(res$u_Hmg,    0.025938,    tolerance = 1e-6)
  expect_equal(res$Hvg,     36.874304,    tolerance = 1e-3)
  expect_equal(res$u_Hvg,    0.022289,    tolerance = 1e-6)
})

################################################################################
# Example 3 (Annex D): 11-component mixture, identity matrix, 15/15 °C
################################################################################
test_that("Example 3 — identity matrix, 15/15 °C", {
  data("example3", envir = environment())
  res <- calculateProperties(example3$fractionArray, example3$uncertaintyArray, example3$correlationMatrix,
                             combustionTemperature = 15,
                             volumeTemperature = 15,
                             coverage = 1)

  expect_equal(res$Hvg,    39.73351,  tolerance = 1e-5)
  expect_equal(res$u_Hvg,   0.026917, tolerance = 1e-6)
  expect_equal(res$Hvn,    35.86811,  tolerance = 1e-5)
  expect_equal(res$u_Hvn,   0.024757, tolerance = 1e-6)
  expect_equal(res$D,       0.76462,  tolerance = 1e-5)
  expect_equal(res$u_D,     0.000586, tolerance = 1e-6)
  expect_equal(res$G,       0.62391,  tolerance = 1e-5)
  expect_equal(res$u_G,     0.000478, tolerance = 1e-6)
  expect_equal(res$Wg,     50.30318,  tolerance = 1e-5)
  expect_equal(res$u_Wg,    0.021588, tolerance = 1e-6)
  expect_equal(res$Wn,     45.40954,  tolerance = 1e-5)
  expect_equal(res$u_Wn,    0.020151, tolerance = 1e-6)
})

################################################################################
# Example 3, identity matrix, 25/0 °C
################################################################################
test_that("Example 3 — identity matrix, 25/0 °C", {
  data("example3", envir = environment())
  res <- calculateProperties(example3$fractionArray, example3$uncertaintyArray, example3$correlationMatrix,
                             combustionTemperature = 25,
                             volumeTemperature = 0,
                             coverage = 1)

  expect_equal(res$Hvg,    41.89360,  tolerance = 1e-5)
  expect_equal(res$u_Hvg,   0.028425, tolerance = 1e-6)
  expect_equal(res$Hvn,    37.85228,  tolerance = 1e-5)
  expect_equal(res$u_Hvn,   0.026164, tolerance = 1e-6)
  expect_equal(res$D,       0.80701,  tolerance = 1e-5)
  expect_equal(res$u_D,     0.000619, tolerance = 1e-6)
  expect_equal(res$G,       0.62411,  tolerance = 1e-5)
  expect_equal(res$u_G,     0.000479, tolerance = 1e-6)
  expect_equal(res$Wg,     53.02930,  tolerance = 1e-5)
  expect_equal(res$u_Wg,    0.022783, tolerance = 1e-6)
  expect_equal(res$Wn,     47.91376,  tolerance = 1e-5)
  expect_equal(res$u_Wn,    0.021278, tolerance = 1e-6)
})

################################################################################
# Example 3, full correlation matrix, 15/15 °C
################################################################################
test_that("Example 3 — full correlation matrix, 15/15 °C", {
  data("example3_ex", envir = environment())
  res <- calculateProperties(example3_ex$fractionArray, example3_ex$uncertaintyArray, example3_ex$correlationMatrix,
                             combustionTemperature = 15,
                             volumeTemperature = 15)

  expect_equal(res$Hvg,    39.73351, tolerance = 1e-5)
  expect_equal(res$u_Hvg,   0.016316, tolerance = 1e-6)
  expect_equal(res$Hvn,    35.86811, tolerance = 1e-5)
  expect_equal(res$u_Hvn,   0.015305, tolerance = 1e-6)
  expect_equal(res$D,       0.76462,  tolerance = 1e-5)
  expect_equal(res$u_D,     0.000277, tolerance = 1e-6)
  expect_equal(res$G,       0.62391,  tolerance = 1e-5)
  expect_equal(res$u_G,     0.000226, tolerance = 1e-6)
  expect_equal(res$Wg,     50.30318,  tolerance = 1e-5)
  expect_equal(res$u_Wg,    0.019823, tolerance = 1e-6)
  expect_equal(res$Wn,     45.40954,  tolerance = 1e-5)
  expect_equal(res$u_Wn,    0.018498, tolerance = 1e-6)
})

################################################################################
# Example 3, full correlation matrix, 25/0 °C
################################################################################
test_that("Example 3 — full correlation matrix, 25/0 °C", {
  data("example3_ex", envir = environment())
  res <- calculateProperties(example3_ex$fractionArray, example3_ex$uncertaintyArray, example3_ex$correlationMatrix,
                             combustionTemperature = 25,
                             volumeTemperature = 0,
                             coverage = 1)

  expect_equal(res$Hvg,    41.89360,  tolerance = 1e-5)
  expect_equal(res$u_Hvg,   0.017241, tolerance = 1e-6)
  expect_equal(res$Hvn,    37.85228,  tolerance = 1e-5)
  expect_equal(res$u_Hvn,   0.016181, tolerance = 1e-6)
  expect_equal(res$D,       0.80701,  tolerance = 1e-5)
  expect_equal(res$u_D,     0.000293, tolerance = 1e-6)
  expect_equal(res$G,       0.62411,  tolerance = 1e-5)
  expect_equal(res$u_G,     0.000227, tolerance = 1e-6)
  expect_equal(res$Wg,     53.02930,  tolerance = 1e-5)
  expect_equal(res$u_Wg,    0.020914, tolerance = 1e-6)
  expect_equal(res$Wn,     47.91376,  tolerance = 1e-5)
  expect_equal(res$u_Wn,    0.019528, tolerance = 1e-6)
})

################################################################################
# Ideal-gas relational checks (Example 3, identity matrix, 15/15 °C)
# Verifies internal consistency of the ideal-gas properties.
################################################################################
test_that("Example 3 — ideal-gas property relationships", {
  data("example3", envir = environment())
  res <- calculateProperties(example3$fractionArray, example3$uncertaintyArray, example3$correlationMatrix,
                             combustionTemperature = 15,
                             volumeTemperature = 15,
                             coverage = 1)

  expect_equal(res$G_o, res$M / 28.96546,           tolerance = 1e-9)  # G_o = M/M_air
  expect_equal(res$D_o, res$D * res$Z,               tolerance = 1e-6)  # D_o = D * Z
  expect_equal(res$Hvg_o, res$Hvg * res$Z,           tolerance = 1e-6)  # Hv_o = Hv * Z
  expect_equal(res$Hvn_o, res$Hvn * res$Z,           tolerance = 1e-6)
  expect_equal(res$Wg_o, res$Hvg_o / sqrt(res$G_o),  tolerance = 1e-6)
  expect_equal(res$Wn_o, res$Hvn_o / sqrt(res$G_o),  tolerance = 1e-6)
})

################################################################################
# Ideal-gas uncertainty checks (Example 1, 15/15 °C)
# u_Hvg_o and u_Hvn_o are not tabulated in ISO 6976:2016 Annex D.
# The standard states Hv = Hv_o / Z, so u(Hv) = u(Hv_o) / Z, giving
# u(Hv_o) = u(Hv) * Z.  We verify this structural identity and that the
# values are positive.
################################################################################
test_that("Example 1 — ideal-gas CV uncertainty structural identity", {
  data("example1", envir = environment())
  res <- calculateProperties(example1$fractionArray, example1$uncertaintyArray,
                             example1$correlationMatrix,
                             combustionTemperature = 15, volumeTemperature = 15,
                             coverage = 1)

  expect_equal(res$u_Hvg_o, res$u_Hvg * res$Z, tolerance = 1e-9)
  expect_equal(res$u_Hvn_o, res$u_Hvn * res$Z, tolerance = 1e-9)
  expect_true(res$u_Hvg_o > 0)
  expect_true(res$u_Hvn_o > 0)
})

################################################################################
# NCV sanity checks (Example 1, 15/15 °C)
################################################################################
test_that("Example 1 — NCV values finite, positive, and below GCV", {
  data("example1", envir = environment())
  res <- calculateProperties(example1$fractionArray, example1$uncertaintyArray, example1$correlationMatrix,
                             combustionTemperature = 15,
                             volumeTemperature = 15,
                             coverage = 1)

  expect_true(is.finite(res$Hcn)   && res$Hcn   > 0)
  expect_true(is.finite(res$Hmn)   && res$Hmn   > 0)
  expect_true(is.finite(res$u_Hcn) && res$u_Hcn > 0)
  expect_true(is.finite(res$u_Hmn) && res$u_Hmn > 0)
  expect_lt(res$Hcn, res$Hcg)
  expect_lt(res$Hmn, res$Hmg)
})

################################################################################
# componentIndex / componentName / componentNames
################################################################################
test_that("componentNames() returns length-60 character vector", {
  nm <- componentNames()
  expect_type(nm, "character")
  expect_length(nm, 60)
})

test_that("componentIndex returns correct indices for boundary components", {
  expect_equal(componentIndex("methane"),       1L)
  expect_equal(componentIndex("ethane"),        2L)
  expect_equal(componentIndex("nitrogen"),      52L)
  expect_equal(componentIndex("carbon dioxide"), 54L)
  expect_equal(componentIndex("n-pentadecane"), 60L)
})

test_that("componentName returns correct names for boundary indices", {
  expect_equal(componentName(1L),  "methane")
  expect_equal(componentName(52L), "nitrogen")
  expect_equal(componentName(54L), "carbon dioxide")
  expect_equal(componentName(60L), "n-pentadecane")
})

test_that("componentIndex / componentName round-trip", {
  for (i in c(1L, 10L, 30L, 54L, 60L)) {
    expect_equal(componentIndex(componentName(i)), i)
  }
})

test_that("componentIndex errors on unknown name", {
  expect_error(componentIndex("unobtainium"), "Unknown component")
})

test_that("componentName errors on out-of-range index", {
  expect_error(componentName(0L),  "index must be between 1 and 60")
  expect_error(componentName(61L), "index must be between 1 and 60")
})

################################################################################
# GasComponents R6 class — construction and defaults
################################################################################
test_that("GasComponents initialises with zeros and identity matrix", {
  gc <- GasComponents$new()
  expect_equal(gc$fractions,    numeric(60))
  expect_equal(gc$uncertainties, numeric(60))
  expect_equal(gc$correlations,  diag(60))
})

################################################################################
# GasComponents — single-component getters/setters by index and by name
################################################################################
test_that("setFraction / getFraction work by index and by name", {
  gc <- GasComponents$new()
  gc$setFraction(1L, 0.9)
  expect_equal(gc$getFraction(1L), 0.9)
  expect_equal(gc$getFraction("methane"), 0.9)

  gc$setFraction("nitrogen", 0.05)
  expect_equal(gc$getFraction(52L),       0.05)
  expect_equal(gc$getFraction("nitrogen"), 0.05)
})

test_that("setUncertainty / getUncertainty work by index and by name", {
  gc <- GasComponents$new()
  gc$setUncertainty(1L, 0.001)
  expect_equal(gc$getUncertainty(1L),       0.001)
  expect_equal(gc$getUncertainty("methane"), 0.001)

  gc$setUncertainty("carbon dioxide", 0.0005)
  expect_equal(gc$getUncertainty(54L), 0.0005)
})

test_that("setCorrelation / getCorrelation are symmetric", {
  gc <- GasComponents$new()
  gc$setCorrelation(1L, 2L, 0.3)
  expect_equal(gc$getCorrelation(1L, 2L), 0.3)
  expect_equal(gc$getCorrelation(2L, 1L), 0.3)   # symmetry
  expect_equal(gc$correlations[1, 2],     0.3)
  expect_equal(gc$correlations[2, 1],     0.3)
})

################################################################################
# GasComponents — array setters
################################################################################
test_that("setFractionArray accepts length-60 vector", {
  gc <- GasComponents$new()
  x <- rep(0, 60)
  x[1] <- 0.85
  x[52] <- 0.15
  gc$setFractionArray(x)
  expect_equal(gc$fractions, x)
})

test_that("setUncertaintyArray accepts length-60 vector", {
  gc <- GasComponents$new()
  u  <- rep(0.001, 60)
  gc$setUncertaintyArray(u)
  expect_equal(gc$uncertainties, u)
})

test_that("setCorrelationMatrix accepts valid 60x60 matrix", {
  gc <- GasComponents$new()
  r <- diag(60)
  r[1, 2] <- r[2, 1] <- -0.5
  gc$setCorrelationMatrix(r)
  expect_equal(gc$correlations[1, 2], -0.5)
  expect_equal(gc$correlations[2, 1], -0.5)
})

################################################################################
# GasComponents — input validation errors
################################################################################
test_that("setFractionArray rejects wrong-length input", {
  gc <- GasComponents$new()
  expect_error(gc$setFractionArray(numeric(59)), "length 60")
  expect_error(gc$setFractionArray(numeric(61)), "length 60")
})

test_that("setUncertaintyArray rejects wrong-length input", {
  gc <- GasComponents$new()
  expect_error(gc$setUncertaintyArray(numeric(1)), "length 60")
})

test_that("setCorrelationMatrix rejects non-matrix and wrong-size inputs", {
  gc <- GasComponents$new()
  expect_error(gc$setCorrelationMatrix(numeric(3600)), "60x60 matrix")
  expect_error(gc$setCorrelationMatrix(matrix(0, 59, 60)), "60x60 matrix")
})

test_that("setCorrelationMatrix rejects out-of-range coefficients", {
  gc  <- GasComponents$new()
  bad <- diag(60)
  bad[1, 2] <- 1.1
  expect_error(gc$setCorrelationMatrix(bad), "\\[-1, 1\\]")
})

test_that(".resolve errors on index out of range", {
  gc <- GasComponents$new()
  expect_error(gc$getFraction(0L),  "between 1 and 60")
  expect_error(gc$getFraction(61L), "between 1 and 60")
})

test_that(".resolve errors on unknown name", {
  gc <- GasComponents$new()
  expect_error(gc$getFraction("unobtainium"), "Unknown component")
})

################################################################################
# calculateProperties() — input validation
################################################################################
test_that("calculateProperties errors on wrong-length compositionArray", {
  expect_error(
    calculateProperties(numeric(59), numeric(60), diag(60)),
    "length 60"
  )
})

test_that("calculateProperties errors on wrong-length uncertaintyArray", {
  expect_error(
    calculateProperties(numeric(60), numeric(61), diag(60)),
    "length 60"
  )
})

test_that("calculateProperties errors on non-matrix correlationMatrix", {
  expect_error(
    calculateProperties(numeric(60), numeric(60), numeric(3600)),
    "60x60 matrix"
  )
})

test_that("calculateProperties errors on wrong-size correlationMatrix", {
  expect_error(
    calculateProperties(numeric(60), numeric(60), matrix(0, 59, 60)),
    "60x60 matrix"
  )
})

test_that("calculateProperties errors on invalid combustionTemperature", {
  expect_error(
    calculateProperties(numeric(60), numeric(60), diag(60),
                        combustionTemperature = 30),
    "combustionTemperature"
  )
})

test_that("calculateProperties errors on invalid volumeTemperature", {
  expect_error(
    calculateProperties(numeric(60), numeric(60), diag(60),
                        volumeTemperature = 10),
    "volumeTemperature"
  )
})

test_that("calculateProperties errors on out-of-range pressure", {
  expect_error(
    calculateProperties(numeric(60), numeric(60), diag(60), pressure = 50),
    "pressure"
  )
  expect_error(
    calculateProperties(numeric(60), numeric(60), diag(60), pressure = 200),
    "pressure"
  )
})

################################################################################
# calculateProperties() — coverage factor scaling
################################################################################
test_that("coverage factor k doubles all uncertainty outputs", {
  data("example1", envir = environment())
  r1 <- calculateProperties(example1$fractionArray, example1$uncertaintyArray, example1$correlationMatrix,
                             combustionTemperature = 15, volumeTemperature = 15,
                             coverage = 1)
  r2 <- calculateProperties(example1$fractionArray, example1$uncertaintyArray, example1$correlationMatrix,
                             combustionTemperature = 15, volumeTemperature = 15,
                             coverage = 2)

  u_names <- c("u_G", "u_D", "u_Hcg", "u_Hcn", "u_Hmg", "u_Hmn",
                "u_Hvg_o", "u_Hvn_o", "u_Hvg", "u_Hvn", "u_Wg", "u_Wn")
  for (nm in u_names) {
    expect_equal(r2[[nm]], 2 * r1[[nm]], tolerance = 1e-12,
                 label = paste0("2 * u (k=1) == u (k=2) for ", nm))
  }
})

test_that("coverage factor does not affect non-uncertainty outputs", {
  data("example1", envir = environment())
  r1 <- calculateProperties(example1$fractionArray, example1$uncertaintyArray, example1$correlationMatrix,
                             combustionTemperature = 15, volumeTemperature = 15,
                             coverage = 1)
  r5 <- calculateProperties(example1$fractionArray, example1$uncertaintyArray, example1$correlationMatrix,
                             combustionTemperature = 15, volumeTemperature = 15,
                             coverage = 5)

  det_names <- c("M", "Z", "G_o", "D_o", "G", "D",
                 "Hcg", "Hcn", "Hmg", "Hmn",
                 "Hvg_o", "Hvn_o", "Hvg", "Hvn",
                 "Wg_o", "Wn_o", "Wg", "Wn")
  for (nm in det_names) {
    expect_equal(r5[[nm]], r1[[nm]], tolerance = 1e-12,
                 label = paste0(nm, " unchanged by coverage factor"))
  }
})
