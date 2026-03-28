#' ISO 6976:2016 Annex D.2 — Example 1 (5-component mixture, 15/15 °C)
#'
#' A five-component natural gas mixture from Annex D.2 of ISO 6976:2016.
#' Non-zero components: methane, ethane, propane, nitrogen, carbon dioxide.
#'
#' @format A named list with three elements:
#' \describe{
#'   \item{fractionArray}{Numeric vector (length 60) of mole fractions [mol/mol].}
#'   \item{uncertaintyArray}{Numeric vector (length 60) of standard uncertainties [mol/mol].}
#'   \item{correlationMatrix}{60×60 identity correlation matrix.}
#' }
#' @source ISO 6976:2016 Annex D.2, Table D.2.
#' @name example1
NULL

#' ISO 6976:2016 Annex D.3 — Example 2 (mixture with water vapour, 15.55/15.55 °C)
#'
#' An eleven-component natural gas mixture including water vapour from
#' Annex D.3 of ISO 6976:2016.
#'
#' @format A named list with three elements (same structure as \code{\link{example1}}).
#' @source ISO 6976:2016 Annex D.3, Table D.3.
#' @name example2
NULL

#' ISO 6976:2016 Annex D — Example 3 (11-component mixture, identity matrix)
#'
#' An eleven-component natural gas mixture from Annex D of ISO 6976:2016
#' with an identity correlation matrix (no inter-component correlations).
#' Used to verify calculations at 15/15 °C and 25/0 °C.
#'
#' @format A named list with three elements (same structure as \code{\link{example1}}).
#' @source ISO 6976:2016 Annex D.
#' @name example3
NULL

#' ISO 6976:2016 Annex D — Example 3, full correlation matrix
#'
#' The same eleven-component mixture as \code{\link{example3}} but with a full
#' inter-component correlation matrix from a GC calibration covariance analysis,
#' as given in Annex D of ISO 6976:2016.
#'
#' @format A named list with three elements (same structure as \code{\link{example1}}).
#' @source ISO 6976:2016 Annex D (extended correlation case).
#' @name example3_ex
NULL
