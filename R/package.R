#' ISO6976.2016: Calorific Values and Properties of Natural Gas per ISO 6976:2016
#'
#' Calculates calorific values (gross and net, molar, mass and volumetric
#' bases), density, relative density, and Wobbe indices together with their
#' standard uncertainties from natural gas composition, following
#' ISO 6976:2016 "Natural Gas — Calculation of calorific values, density,
#' relative density and Wobbe indices from composition".
#'
#' Uncertainty propagation is implemented according to Annex B of that
#' standard (variance-covariance method).
#'
#' **Application restrictions (ISO 6976:2016 §5):**
#' * Combustion temperature: 0, 15, 15.55 (60 °F), 20, or 25 °C.
#' * Volume reference temperature: 0, 15, 15.55 (60 °F), or 20 °C.
#' * Reference pressure: 90–110 kPa.
#' * Compression factor Z must be > 0.9.
#'
#' @references
#' ISO 6976:2016 "Natural Gas — Calculation of calorific values, density,
#' relative density and Wobbe indices from composition".
#'
#' @author Rüdiger Forster \email{meticulous.measurements@@gmail.com}
#'
#' @name ISO6976.2016-package
#' @aliases ISO6976.2016-package
#' @useDynLib ISO6976.2016, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom R6 R6Class
"_PACKAGE"
