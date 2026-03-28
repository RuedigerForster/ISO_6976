# Tests against the worked examples in Annex D of ISO 6976:2016
# All expected values are taken directly from the standard.

library(ISO6976.2016)

################################################################################
# Example 1 (Annex D.2): 5-component mixture, 15/15 °C
################################################################################
test_that("Example 1 — 15/15 °C: molar mass, Z, molar GCV, mass GCV, vol. GCV", {
  data("example1")
  res <- calculateProperties(fractionArray, uncertaintyArray, correlationMatrix,
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
  data("example2")
  res <- calculateProperties(fractionArray, uncertaintyArray, correlationMatrix,
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
  data("example3")
  res <- calculateProperties(fractionArray, uncertaintyArray, correlationMatrix,
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
  data("example3")
  res <- calculateProperties(fractionArray, uncertaintyArray, correlationMatrix,
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
# Note: uncertainty values in the standard contain errors (factor-of-2 issue);
# expected values below are corrected (k = 2 applied externally).
################################################################################
test_that("Example 3 — full correlation matrix, 15/15 °C", {
  data("example3_ex")
  k   <- 2
  res <- calculateProperties(fractionArray, uncertaintyArray, correlationMatrix,
                             combustionTemperature = 15,
                             volumeTemperature = 15)

  expect_equal(res$Hvg,         39.73351, tolerance = 1e-5)
  expect_equal(res$u_Hvg * k,    0.016316, tolerance = 1e-6)
  expect_equal(res$Hvn,         35.86811, tolerance = 1e-5)
  expect_equal(res$u_Hvn * k,    0.015305, tolerance = 1e-6)
  expect_equal(res$D,            0.76462,  tolerance = 1e-5)
  expect_equal(res$u_D  * k,     0.000277, tolerance = 1e-6)
  expect_equal(res$G,            0.62391,  tolerance = 1e-5)
  expect_equal(res$u_G  * k,     0.000478, tolerance = 1e-3)
  expect_equal(res$Wg,          50.30318,  tolerance = 1e-5)
  expect_equal(res$u_Wg * k,     0.019823, tolerance = 1e-6)
  expect_equal(res$Wn,          45.40954,  tolerance = 1e-5)
  expect_equal(res$u_Wn * k,     0.018498, tolerance = 1e-6)
})

################################################################################
# Example 3, full correlation matrix, 25/0 °C
################################################################################
test_that("Example 3 — full correlation matrix, 25/0 °C", {
  data("example3_ex")
  res <- calculateProperties(fractionArray, uncertaintyArray, correlationMatrix,
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
  data("example3")
  res <- calculateProperties(fractionArray, uncertaintyArray, correlationMatrix,
                             combustionTemperature = 15,
                             volumeTemperature = 15,
                             coverage = 1)

  Z_air_15 <- 0.999595   # Table A.1

  expect_equal(res$G_o, res$M / 28.96546,         tolerance = 1e-9)  # G_o = M/M_air
  expect_equal(res$D_o, res$D * res$Z,             tolerance = 1e-6)  # D_o = D * Z
  expect_equal(res$Hvg_o, res$Hvg * res$Z,         tolerance = 1e-6)  # Hv_o = Hv * Z
  expect_equal(res$Hvn_o, res$Hvn * res$Z,         tolerance = 1e-6)
  expect_equal(res$Wg_o, res$Hvg_o / sqrt(res$G_o), tolerance = 1e-6)
  expect_equal(res$Wn_o, res$Hvn_o / sqrt(res$G_o), tolerance = 1e-6)
})

################################################################################
# NCV sanity checks (Example 1, 15/15 °C)
################################################################################
test_that("Example 1 — NCV values finite, positive, and below GCV", {
  data("example1")
  res <- calculateProperties(fractionArray, uncertaintyArray, correlationMatrix,
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
