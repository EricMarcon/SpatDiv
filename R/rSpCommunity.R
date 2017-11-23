#' Random Spatialized Community
#'
#' @param n The number of communities to draw.
#' @param size The number of individuals to draw in each community.
#' @param NorP A numeric vector or a two-column matrix. Contains either abundances or probabilities. Two-column matrices should contain the observed abundances (or probabilities) in the first column and the expected ones in the second column, to allow using beta diversity functions.
#' @param BootstrapMethod  The method used to obtain the probabilities to generate bootstrapped communities from observed abundances. If \code{"Marcon"}, the probabilities are simply the abundances divided by the total number of individuals (Marcon \emph{et al.}, 2012). If \code{"Chao2013"} or \code{"Chao2015"} (by default), a more sophisticated approach is used (see \code{\link{as.ProbaVector}}) following Chao \emph{et al.} (2013) or Chao \emph{et al.} (2015).
#' @param S The number of species.
#' @param Distribution The distribution of species frequencies. May be \code{"lnorm"} (log-normal), \code{"lseries"} (log-series), \code{"geom"} (geometric) or \code{"bstick"} (broken stick).
#' @param sd The simulated distribution standard deviation. For the log-normal distribution, this is the standard deviation on the log scale.
#' @param prob The probability of success in each trial.
#' @param alpha Fisher's alpha.
#' @param Spatial TODO
#' @param win TODO
#' @param Species TODO
#' @param Sizes TODO
#' @inheritParams as.ProbaVector.wmppp
#'
#' @return An object of class "wmppp" (\code{\link{wmppp.object}}), with \code{PointType} values as species names if \code{n}=1.
#'         An object of class "SpCommunities", which is a list of \code{\link{wmppp.object}}s, is returned if \code{n}>1.
#' @export
#'
#' @examples TODO
rSpCommunity <-
function(n, size = sum(NorP), NorP = 1, BootstrapMethod = "Chao2015",
         S = 300, Distribution = "lnorm", sd = 1, prob = 0.1, alpha=40,
         Spatial = "Binomial", win=spatstat::owin(),
         Species = NULL,
         Sizes = "Uniform",
         CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckSpatDivArguments()

  # Species
  if (is.null(Species)) {
    Format <- ceiling(log10(S+1))
    Species <- paste("sp", formatC(1:S, width=Format, flag="0"), sep="")
  }

  # Sizes
  if (Sizes == "Uniform") {
    Sizes <- rep(1, S)
  }

  # Spatial distribution
  if (Spatial == "Binomial") {
    # Binomial point process: call rCommunity and runifpoint.
    Community <- entropart::rCommunity(n=n, size=size, NorP=NorP, BootstrapMethod=BootstrapMethod, S=S, Distribution=Distribution, sd=sd, prob=prob, alpha=alpha)
    # Community is an AbdVector if n=1 or a MetaCommunity if n>1. Make an array in all cases.
    if (n == 1) {
      Community = matrix(Community, ncol=1)
    } else {
      Community <- Community$Nsi
    }
    Binomial <- function(i) {
      X <- dbmss::as.wmppp(spatstat::runifpoint(sum(Community[, i]), win=win))
      # Associate species and points
      X$marks$PointType <- as.factor(rep(Species, Community[, i]))
      # Associate sizes and points
      X$marks$PointWeight <- Sizes
      return(X)
    }
    # Loop to simulate several point processes
    listX <- lapply(1:n, Binomial)
  }

  if (n == 1) {
    # Return a wmppp
    return(listX[[1]])
  } else {
    # Return a list of wmppp
    class(listX) <- c("SpCommunities", class(listX))
    return(listX)
  }
}
