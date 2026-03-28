# Component name list — 60 components per ISO 6976:2016 Table A.2
.component_names <- c(
  "methane",           #  1
  "ethane",            #  2
  "propane",           #  3
  "n-butane",          #  4
  "isobutane",         #  5
  "n-pentane",         #  6
  "isopentane",        #  7
  "neopentane",        #  8
  "n-hexane",          #  9
  "2-methylpentane",   # 10
  "3-methylpentane",   # 11
  "2,2-dimethylbutane", # 12
  "2,3-dimethylbutane", # 13
  "n-heptane",         # 14
  "n-octane",          # 15
  "n-nonane",          # 16
  "n-decane",          # 17
  "ethylene",          # 18
  "propylene",         # 19
  "1-butene",          # 20
  "cis-2-butene",      # 21
  "trans-2-butene",    # 22
  "isobutylene",       # 23
  "1-pentene",         # 24
  "propadiene",        # 25
  "1,2-butadiene",     # 26
  "1,3-butadiene",     # 27
  "acetylene",         # 28
  "cyclopentane",      # 29
  "methylcyclopentane", # 30
  "ethylcyclopentane", # 31
  "cyclohexane",       # 32
  "methylcyclohexane", # 33
  "ethylcyclohexane",  # 34
  "benzene",           # 35
  "toluene",           # 36
  "ethylbenzene",      # 37
  "o-xylene",          # 38
  "methanol",          # 39
  "methanethiol",      # 40
  "hydrogen",          # 41
  "water",             # 42
  "hydrogen sulphide", # 43
  "ammonia",           # 44
  "hydrogen cyanide",  # 45
  "carbon monoxide",   # 46
  "carbonyl sulphide", # 47
  "carbon disulphide", # 48
  "helium",            # 49
  "neon",              # 50
  "argon",             # 51
  "nitrogen",          # 52
  "oxygen",            # 53
  "carbon dioxide",    # 54
  "sulphur dioxide",   # 55
  "n-undecane",        # 56
  "n-dodecane",        # 57
  "n-tridecane",       # 58
  "n-tetradecane",     # 59
  "n-pentadecane"      # 60
)

#' Return the index of a gas component by name
#'
#' @param name Character string: component name (English). See
#'   \code{\link{componentNames}} for the full list.
#' @return Integer index (1–60) into the composition vector.
#' @export
componentIndex <- function(name) {
  idx <- match(name, .component_names)
  if (is.na(idx)) stop("Unknown component: ", name)
  idx
}

#' Return the name of a gas component by index
#'
#' @param index Integer 1–60.
#' @return Character string with the English component name.
#' @export
componentName <- function(index) {
  if (index < 1L || index > 60L) stop("index must be between 1 and 60")
  .component_names[[index]]
}

#' Return the names of all 60 gas components
#'
#' Returns the ordered character vector of all 60 natural gas components
#' recognised by ISO 6976:2016 (Table A.2). The position of each name
#' corresponds to the index used in \code{compositionArray},
#' \code{uncertaintyArray}, and \code{correlationMatrix}.
#'
#' @return Character vector of length 60.
#' @seealso \code{\link{componentIndex}}, \code{\link{componentName}}
#' @export
componentNames <- function() .component_names
