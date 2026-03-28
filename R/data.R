#' Example 1 from ISO 6976:2016 Annex D.2
#'
#' A simple five-component mixture at 15/15 °C.  Used as the primary
#' reference case for testing the calculation engine.
#'
#' @format Three objects are loaded:
#' \describe{
#'   \item{fractionArray}{Numeric vector (length 60) of mole fractions [mol/mol].}
#'   \item{uncertaintyArray}{Numeric vector (length 60) of standard uncertainties.}
#'   \item{correlationMatrix}{60x60 identity correlation matrix.}
#' }
#' @source ISO 6976:2016, Annex D, Example 1 (Table D.2).
#' @examples
#' data("example1")
#' res <- calculateProperties(fractionArray, uncertaintyArray, correlationMatrix,
#'                            combustionTemperature = 15, volumeTemperature = 15)
"example1"

#' Example 2 from ISO 6976:2016 Annex D.3
#'
#' A mixture containing water vapour, at 15.55/15.55 °C.
#'
#' @format Three objects are loaded (same structure as \code{\link{example1}}).
#' @source ISO 6976:2016, Annex D, Example 2 (Table D.3).
#' @examples
#' data("example2")
#' res <- calculateProperties(fractionArray, uncertaintyArray, correlationMatrix,
#'                            combustionTemperature = 15.55,
#'                            volumeTemperature = 15.55)
"example2"

#' Example 3 from ISO 6976:2016 Annex D — identity correlation matrix
#'
#' An 11-component mixture with an identity (uncorrelated) correlation matrix.
#'
#' @format Three objects are loaded (same structure as \code{\link{example1}}).
#' @source ISO 6976:2016, Annex D, Example 3.
#' @examples
#' data("example3")
#' res <- calculateProperties(fractionArray, uncertaintyArray, correlationMatrix,
#'                            combustionTemperature = 15, volumeTemperature = 15)
"example3"

#' Example 3 from ISO 6976:2016 Annex D — full correlation matrix
#'
#' The same 11-component mixture as \code{\link{example3}} but with a full
#' (non-identity) correlation matrix between components.
#'
#' @format Three objects are loaded (same structure as \code{\link{example1}}).
#' @source ISO 6976:2016, Annex D, Example 3 (extended case).
#' @examples
#' data("example3_ex")
#' res <- calculateProperties(fractionArray, uncertaintyArray, correlationMatrix,
#'                            combustionTemperature = 15, volumeTemperature = 15)
"example3_ex"
