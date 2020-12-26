#' Species Distribution
#'
#' A Species Distribution is a (preferably named) vector containing species abundances or probabilities.
#' `SpeciesDistribution` objects include [AbdVector] and [ProbaVector] objects.
#'
#' `as.AbdVector` counts the number of individuals (points) per species (in `marks$PointType`).
#'
#' `as.ProbaVector` normalizes the vector so that it sums to 1. If `Correction` is not `"None"`, the observed abundance distribution is used to estimate the actual species distribution. The list of species will be changed: zero-abundance species will be cleared, and some unobserved species will be added. First, observed species probabilities are estimated folllowing \insertCite{Chao2003;textual}{SpatDiv}, *i.e.* input probabilities are multiplied by the sample coverage, or according to more sophisticated models: \insertCite{Chao2013;textual}{SpatDiv}, single-parameter model or \insertCite{Chao2015;textual}{SpatDiv}, two-parameter model. The total probability of observed species equals the sample coverage. Then, the distribution of unobserved species can be unveiled: their number is estimated according to `RCorrection` (if the Jackknife estimator is chosen, the `JackOver` argument allows using the order immediately over the optimal one). The coverage deficit (1 minus the sample coverage) is shared by the unobserved species equally: `Unveiling = "unif"`, \insertCite{Chao2013;textual}{SpatDiv} or according to a geometric distribution: `Unveiling = "geom"`, \insertCite{Chao2015;textual}{SpatDiv}.
#'
#' `SpeciesDistribution` objects can be plotted. The `plot` method returns the estimated parameters of the fitted distribution. The broken stick has no parameter, so the maximum abundance is returned.
#'
#' @param x A [wmppp.object] with `PointType` values as species names, or a vector of factors or characters containing species names of each individual.
#' @param ... Further arguments. Unsused.
#' @name SpeciesDistributions
#' @return A vector of species abundances ([AbdVector]) or probabilities ([ProbaVector]).
#' @references
#' \insertAllCited{}
NULL



#' @rdname SpeciesDistributions
#' @importFrom entropart as.SpeciesDistribution
#' @method as.SpeciesDistribution wmppp
#' @export
as.SpeciesDistribution.wmppp <-
function (x, ...)
{
  # Table counts the number of individuals per species. It returns an array (1d) that must be converted to a vector.
  spD <- as.numeric(table(x$marks$PointType))
  class(spD) <- c("SpeciesDistribution", class(spD))
  return(spD)
}



#' @rdname SpeciesDistributions
#' @importFrom entropart as.SpeciesDistribution
#' @method as.SpeciesDistribution factor
#' @export
as.SpeciesDistribution.factor <-
function (x, ...)
{
  # Table counts the number of individuals per species. It returns an array (1d) that must be converted to a vector.
  spD <- as.numeric(table(x))
  class(spD) <- c("SpeciesDistribution", class(spD))
  return(spD)
}



#' @rdname SpeciesDistributions
#' @importFrom entropart as.SpeciesDistribution
#' @method as.SpeciesDistribution character
#' @export
as.SpeciesDistribution.character <-
function (x, ...)
{
  # Table counts the number of individuals per species. It returns an array (1d) that must be converted to a vector.
  spD <- as.numeric(table(x))
  class(spD) <- c("SpeciesDistribution", class(spD))
  return(spD)
}



#' @rdname SpeciesDistributions
#' @importFrom entropart as.ProbaVector
#' @method as.ProbaVector wmppp
#' @export
as.ProbaVector.wmppp  <-
function (x, ...)
{
  spD <- as.SpeciesDistribution.wmppp(x)

  class(spD) <- c("ProbaVector", class(spD))
  return(spD)
}



#' @rdname SpeciesDistributions
#' @importFrom entropart as.ProbaVector
#' @method as.ProbaVector factor
#' @export
as.ProbaVector.factor  <-
function (x, ...)
{
  spD <- as.SpeciesDistribution.factor(x)

  class(spD) <- c("ProbaVector", class(spD))
  return(spD)
}



#' @rdname SpeciesDistributions
#' @importFrom entropart as.ProbaVector
#' @method as.ProbaVector character
#' @export
as.ProbaVector.character  <-
function (x, ...)
{
  spD <- as.SpeciesDistribution.character(x)

  class(spD) <- c("ProbaVector", class(spD))
  return(spD)
}



#' @rdname SpeciesDistributions
#' @importFrom entropart as.AbdVector
#' @method as.AbdVector wmppp
#' @export
as.AbdVector.wmppp  <-
function (x, ...)
{
  spD <- as.SpeciesDistribution.wmppp(x)

  class(spD) <- c("AbdVector", class(spD))
  return(spD)
}



#' @rdname SpeciesDistributions
#' @importFrom entropart as.AbdVector
#' @method as.AbdVector factor
#' @export
as.AbdVector.factor  <-
function (x, ...)
{
  spD <- as.SpeciesDistribution.factor(x)

  class(spD) <- c("AbdVector", class(spD))
  return(spD)
}



#' @rdname SpeciesDistributions
#' @importFrom entropart as.AbdVector
#' @method as.AbdVector character
#' @export
as.AbdVector.character  <-
function (x, ...)
{
  spD <- as.SpeciesDistribution.character(x)

  class(spD) <- c("AbdVector", class(spD))
  return(spD)
}
