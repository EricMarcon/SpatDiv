#' Random Spatialized Community
#' 
#' This function extends [rCommunity] by spatializing the randomized community.
#'
#'
#' @param n The number of communities to draw.
#' @param size The number of individuals to draw in each community.
#' @param NorP A numeric vector or a two-column matrix. 
#' Contains either abundances or probabilities. 
#' Two-column matrices should contain the observed abundances (or probabilities) in the first column 
#' and the expected ones in the second column, to allow using beta diversity functions.
#' @param BootstrapMethod  The method used to obtain the probabilities to generate bootstrapped communities 
#' from observed abundances. 
#' If `"Marcon"`, the probabilities are simply the abundances divided by the total number of individuals \insertCite{Marcon2012a}{SpatDiv}. 
#' If `"Chao2013"` or `"Chao2015"` (by default), a more sophisticated approach is used (see [as.ProbaVector]) 
#' following \insertCite{Chao2013;textual}{SpatDiv} or \insertCite{Chao2015;textual}{SpatDiv}.
#' @param S The number of species.
#' @param Distribution The distribution of species frequencies. 
#' May be `"lnorm"` (log-normal), `"lseries"` (log-series), `"geom"` (geometric) or `"bstick"` (broken stick).
#' @param sd The simulated distribution standard deviation. For the log-normal distribution, this is the standard deviation on the log scale.
#' @param prob The probability of success in each trial.
#' @param alpha Fisher's alpha.
#' @param Spatial The spatial distribution of points. 
#' May be "Binomial" (a completely random point pattern except for its fixed number of points) or 
#' "Thomas" for a clustered point pattern with parameters `scale` and `mu`
#' @param scale In Thomas point patterns, the standard deviation of random displacement (along each coordinate axis) of a point from its cluster center.
#' @param mu In Thomas point patterns, the mean number of points per cluster.
#' The intensity of the Poisson process of cluster centers is calculated as the number of points (`size`) per area divided by `mu`.
#' @param win The window containing the point pattern. It is an [owin] object.
#' @param Species A vector of characters or of factors containing the possible species.
#' @param Sizes The distribution of point sizes.
#' May be "Uniform" for a uniform distribution between `MinSize` and `MaxSize`.
#' By default, all sizes are 1.
#' @param MinSize The minimum size in a uniform distribution.
#' @param MaxSize The maximum size in a uniform distribution.
#' @inheritParams EntAccum
#'
#' @return A [wmppp.object], with `PointType` values as species names if `n`=1.
#'         An object of class "SpCommunities", which is a list of [wmppp.object]s, is returned if `n`>1.
#' @export
#'
#' @references
#' \insertAllCited{}
#' @examples
#' spCommunity <- rSpCommunity(1, size=30, S=5)
#' autoplot(spCommunity, which.marks = "PointType")
#' 
rSpCommunity <-
function(n, size = sum(NorP), NorP = 1, BootstrapMethod = "Chao2015",
         S = 300, Distribution = "lnorm", sd = 1, prob = 0.1, alpha=40,
         Spatial = "Binomial", 
         scale = 0.2, mu = 10,
         win=spatstat.geom::owin(),
         Species = NULL,
         Sizes = "Uniform", MinSize = 1, MaxSize = 1,
         CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckSpatDivArguments()

  # Species
  if (is.null(Species)) {
    Format <- ceiling(log10(S+1))
    Species <- paste("sp", formatC(1:S, width=Format, flag="0"), sep="")
  }
  
  # Species abundance: call rCommunity
  Community <- entropart::rCommunity(n=n, size=size, NorP=NorP, BootstrapMethod=BootstrapMethod, S=S, Distribution=Distribution, sd=sd, prob=prob, alpha=alpha, CheckArguments=FALSE)
  # Community is an AbdVector if n=1 or a MetaCommunity if n>1. Make an array in all cases.
  if (n == 1) {
    Community = matrix(Community, ncol=1)
  } else {
    Community <- Community$Nsi
  }
  
  # Spatial distribution
  if (Spatial == "Binomial") {
    # Sizes
    if (Sizes == "Uniform") {
      Sizes <- stats::runif(size, min=MinSize, max=MaxSize)
    }
    Binomial <- function(i) {
      X <- dbmss::as.wmppp(spatstat.random::runifpoint(sum(Community[, i]), win=win))
      # Associate species and points
      X$marks$PointType <- as.factor(rep(Species, Community[, i]))
      # Associate sizes and points
      X$marks$PointWeight <- Sizes
      return(X)
    }
    # Loop to simulate several point processes
    listX <- lapply(1:n, Binomial)
  }
  
  if (Spatial == "Thomas") {
    area <- spatstat.geom::area.owin(win)
    Thomas <- function(i) {
      # Prepare an empty point pattern
      X <- NULL
      # Draw each species
      for (s in 1:S) {
        Xs <- spatstat.random::rThomas(kappa=Community[s, i]/area/mu, 
                                    scale=scale,
                                    mu=mu, win=win)
        # Associate species and points
        PointType <- as.factor(rep(Species[s], Xs$n))
        # Associate sizes and points
        if (Sizes == "Uniform") {
          PointWeight <- stats::runif(Xs$n, min=MinSize, max=MaxSize)
        }
        # Add the marks
        spatstat.geom::marks(Xs) <- data.frame(PointType, PointWeight)
        # Add the species to the point pattern
        if (is.null(X)) {
          X <- Xs
        } else {
          X <- spatstat.geom::superimpose(X, Xs)
        }
      }
      return(dbmss::as.wmppp(X))
    }
    # Loop to simulate several point processes
    listX <- lapply(1:n, Thomas)
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
