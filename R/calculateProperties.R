#' Calculate natural gas properties per ISO 6976:2016
#'
#' Computes all combustion and volumetric properties of a natural gas mixture
#' together with their standard (or expanded) uncertainties, following the
#' method and tables of ISO 6976:2016 and the uncertainty propagation formulae
#' of its Annex B.
#'
#' @param compositionArray Numeric vector of length 60: mole fractions
#'   \[mol/mol\] in the component order of ISO 6976:2016 Table A.2.
#'   See \code{\link{componentNames}} for the ordering.
#'   If values are given in mol\%, divide by 100 before passing them here.
#' @param uncertaintyArray Numeric vector of length 60: standard uncertainties
#'   of the mole fractions (same units as \code{compositionArray}).
#' @param correlationMatrix 60x60 numeric matrix of correlation coefficients
#'   between component mole fractions.  Use \code{diag(60)} when correlations
#'   are unknown or assumed to be zero.
#' @param combustionTemperature Combustion reference temperature in °C.
#'   Permitted values: 0, 15, 15.55, 20, 25.  Default: 25.
#' @param volumeTemperature Volume reference temperature in °C.
#'   Permitted values: 0, 15, 15.55, 20.  Default: 15.
#' @param pressure Reference pressure in kPa.  Must be in \[90, 110\].
#'   Default: 101.325.
#' @param coverage Coverage factor \eqn{k}; uncertainties are multiplied by
#'   \eqn{k} before being returned.  Default: 1 (standard uncertainty).
#'
#' @return A named list with the following elements (all numeric scalars):
#'
#' | Name      | Description                                         | Unit       |
#' |-----------|-----------------------------------------------------|------------|
#' | M         | Molar mass                                          | kg/kmol    |
#' | Z         | Compression factor                                  | —          |
#' | G_o       | Ideal-gas relative density                          | —          |
#' | D_o       | Ideal-gas density                                   | kg/m³      |
#' | G, u_G    | Real-gas relative density and uncertainty           | —          |
#' | D, u_D    | Real-gas density and uncertainty                    | kg/m³      |
#' | Hcg, u_Hcg| Molar gross calorific value and uncertainty         | kJ/mol     |
#' | Hcn, u_Hcn| Molar net calorific value and uncertainty           | kJ/mol     |
#' | Hmg, u_Hmg| Mass-basis gross calorific value and uncertainty    | MJ/kg      |
#' | Hmn, u_Hmn| Mass-basis net calorific value and uncertainty      | MJ/kg      |
#' | Hvg_o, u_Hvg_o | Ideal-gas vol. gross CV and uncertainty        | MJ/m³      |
#' | Hvn_o, u_Hvn_o | Ideal-gas vol. net CV and uncertainty          | MJ/m³      |
#' | Hvg, u_Hvg| Real-gas vol. gross CV and uncertainty              | MJ/m³      |
#' | Hvn, u_Hvn| Real-gas vol. net CV and uncertainty                | MJ/m³      |
#' | Wg_o      | Ideal-gas gross Wobbe index                         | MJ/m³      |
#' | Wn_o      | Ideal-gas net Wobbe index                           | MJ/m³      |
#' | Wg, u_Wg  | Real-gas gross Wobbe index and uncertainty          | MJ/m³      |
#' | Wn, u_Wn  | Real-gas net Wobbe index and uncertainty            | MJ/m³      |
#'
#' @examples
#' \dontrun{
#' data("example1")
#' res <- calculateProperties(fractionArray, uncertaintyArray, correlationMatrix,
#'                            combustionTemperature = 15, volumeTemperature = 15)
#' res$M      # molar mass [kg/kmol]
#' res$Hvg    # real-gas vol. gross CV [MJ/m^3]
#' res$u_Wg   # standard uncertainty of gross Wobbe index
#' }
#'
#' @seealso \code{\link{GasComponents}}, \code{\link{componentNames}}
#' @references ISO 6976:2016 "Natural Gas — Calculation of calorific values,
#'   density, relative density and Wobbe indices from composition"
#' @export
calculateProperties <- function(
    compositionArray,
    uncertaintyArray,
    correlationMatrix,
    combustionTemperature = 25,
    volumeTemperature     = 15,
    pressure              = 101.325,
    coverage              = 1
) {
  # Input validation
  if (length(compositionArray) != 60)
    stop("compositionArray must have length 60")
  if (length(uncertaintyArray) != 60)
    stop("uncertaintyArray must have length 60")
  if (!is.matrix(correlationMatrix) ||
      nrow(correlationMatrix) != 60 || ncol(correlationMatrix) != 60)
    stop("correlationMatrix must be a 60x60 matrix")
  if (!combustionTemperature %in% c(0, 15, 15.55, 20, 25))
    stop("combustionTemperature must be 0, 15, 15.55, 20, or 25 \u00b0C")
  if (!volumeTemperature %in% c(0, 15, 15.55, 20))
    stop("volumeTemperature must be 0, 15, 15.55, or 20 \u00b0C")
  if (pressure < 90 || pressure > 110)
    stop("pressure must be in the range 90\u2013110 kPa")

  iso6976_calc(
    x   = compositionArray,
    u_x = uncertaintyArray,
    r_x = correlationMatrix,
    t1  = combustionTemperature,
    t2  = volumeTemperature,
    p2  = pressure,
    k   = coverage
  )
}
