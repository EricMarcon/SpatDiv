% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpCommunity.R
\name{rSpCommunity}
\alias{rSpCommunity}
\title{Random Spatialized Community}
\usage{
rSpCommunity(
  n,
  size = sum(NorP),
  NorP = 1,
  BootstrapMethod = "Chao2015",
  S = 300,
  Distribution = "lnorm",
  sd = 1,
  prob = 0.1,
  alpha = 40,
  Spatial = "Binomial",
  win = spatstat::owin(),
  Species = NULL,
  Sizes = "Uniform",
  CheckArguments = TRUE
)
}
\arguments{
\item{n}{The number of communities to draw.}

\item{size}{The number of individuals to draw in each community.}

\item{NorP}{A numeric vector or a two-column matrix. Contains either abundances or probabilities. Two-column matrices should contain the observed abundances (or probabilities) in the first column and the expected ones in the second column, to allow using beta diversity functions.}

\item{BootstrapMethod}{The method used to obtain the probabilities to generate bootstrapped communities from observed abundances. If \code{"Marcon"}, the probabilities are simply the abundances divided by the total number of individuals \insertCite{Marcon2012a}{SpatDiv}. If \code{"Chao2013"} or \code{"Chao2015"} (by default), a more sophisticated approach is used (see \link{as.ProbaVector}) following \insertCite{Chao2013;textual}{SpatDiv} or \insertCite{Chao2015;textual}{SpatDiv}.}

\item{S}{The number of species.}

\item{Distribution}{The distribution of species frequencies. May be \code{"lnorm"} (log-normal), \code{"lseries"} (log-series), \code{"geom"} (geometric) or \code{"bstick"} (broken stick).}

\item{sd}{The simulated distribution standard deviation. For the log-normal distribution, this is the standard deviation on the log scale.}

\item{prob}{The probability of success in each trial.}

\item{alpha}{Fisher's alpha.}

\item{Spatial}{TODO}

\item{win}{The window containing the point pattern. It is an \link{owin} object.}

\item{Species}{A vector of characters or of factors containing the possible species.}

\item{Sizes}{TODO}

\item{CheckArguments}{If \code{TRUE} (default), the function arguments are verified. Should be set to \code{FALSE} to save time in simulations for example, when the arguments have been checked elsewhere.}
}
\value{
A \link{wmppp.object}, with \code{PointType} values as species names if \code{n}=1.
An object of class "SpCommunities", which is a list of \link{wmppp.object}s, is returned if \code{n}>1.
}
\description{
This function extends \link{rCommunity} by spatializing the randomized community.
}
\examples{
spCommunity <- rSpCommunity(1, size=30, S=5)
plot(spCommunity, which.marks = "PointType")

}
\references{
\insertAllCited{}
}
