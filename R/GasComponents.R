#' R6 class for natural gas composition input
#'
#' \code{GasComponents} holds the mole fractions, their standard uncertainties,
#' and the inter-component correlation matrix for a natural gas mixture.
#' All three quantities are required as input to \code{\link{calculateProperties}}.
#'
#' Components are identified either by their integer index (1–60, matching
#' ISO 6976:2016 Table A.2) or by English name (see \code{\link{componentNames}}).
#'
#' @export
GasComponents <- R6::R6Class("GasComponents",
  public = list(
    #' @field fractions Numeric vector (length 60) of mole fractions [mol/mol].
    fractions = numeric(60),

    #' @field uncertainties Numeric vector (length 60) of standard uncertainties.
    uncertainties = numeric(60),

    #' @field correlations 60x60 correlation matrix (default: identity matrix).
    correlations = diag(60),

    # ------------------------------------------------------------------ getters

    #' @description Get the mole fraction of a single component.
    #' @param name Component name (English) or integer index 1–60.
    #' @return Numeric scalar.
    getFraction = function(name) {
      self$fractions[[.resolve(name)]]
    },

    #' @description Get the standard uncertainty of a single component.
    #' @param name Component name (English) or integer index 1–60.
    #' @return Numeric scalar.
    getUncertainty = function(name) {
      self$uncertainties[[.resolve(name)]]
    },

    #' @description Get the correlation between two components.
    #' @param name1 First component name or index.
    #' @param name2 Second component name or index.
    #' @return Numeric scalar in \[-1, 1\].
    getCorrelation = function(name1, name2) {
      self$correlations[.resolve(name1), .resolve(name2)]
    },

    # ------------------------------------------------------------------ setters

    #' @description Set the mole fraction of a single component.
    #' @param name Component name or index.
    #' @param value Mole fraction [mol/mol].
    setFraction = function(name, value) {
      self$fractions[[.resolve(name)]] <- value
    },

    #' @description Set the standard uncertainty of a single component.
    #' @param name Component name or index.
    #' @param value Standard uncertainty.
    setUncertainty = function(name, value) {
      self$uncertainties[[.resolve(name)]] <- value
    },

    #' @description Set the correlation between two components.
    #' @param name1 First component name or index.
    #' @param name2 Second component name or index.
    #' @param value Correlation coefficient in \[-1, 1\].
    setCorrelation = function(name1, name2, value) {
      i <- .resolve(name1); j <- .resolve(name2)
      self$correlations[i, j] <- value
      self$correlations[j, i] <- value
    },

    # ------------------------------------------------------------------ array setters

    #' @description Set all mole fractions at once.
    #' @param x Numeric vector of length 60.
    setFractionArray = function(x) {
      if (length(x) != 60) stop("x must have length 60")
      self$fractions <- x
    },

    #' @description Set all uncertainties at once.
    #' @param u Numeric vector of length 60.
    setUncertaintyArray = function(u) {
      if (length(u) != 60) stop("u must have length 60")
      self$uncertainties <- u
    },

    #' @description Set the full correlation matrix.
    #' @param r 60x60 numeric matrix; values must lie in \[-1, 1\].
    setCorrelationMatrix = function(r) {
      if (!is.matrix(r) || nrow(r) != 60 || ncol(r) != 60)
        stop("r must be a 60x60 matrix")
      if (any(r < -1) || any(r > 1))
        stop("All correlation coefficients must be in [-1, 1]")
      self$correlations <- r
    }
  )
)

# Internal helper: resolve name or integer to an index
.resolve <- function(name) {
  if (is.numeric(name)) {
    idx <- as.integer(name)
    if (idx < 1L || idx > 60L) stop("index must be between 1 and 60")
    return(idx)
  }
  componentIndex(name)
}
