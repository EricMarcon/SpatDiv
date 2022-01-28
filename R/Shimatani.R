#' Spatially Explicit Simpson's Entropy
#' 
#' Simpson's entropy of the neighborhood of individuals, up to a distance \insertCite{Shimatani2001}{SpatDiv}.
#'
#' @param spCommunity A spatialized community (A [wmppp.object] with `PointType` values as species names.)
#' @param r A vector of distances. If `NULL` accumulation is along `n`, else neighbors are accumulated in circles of radius `r`.
#' @param spCorrection The edge-effect correction to apply when estimating the K function with [Kest].
#'        Default is "isotropic".
#' @param CheckArguments If `TRUE` (default), the function arguments are verified. Should be set to `FALSE` to save time in simulations for example, when the arguments have been checked elsewhere.
#' 
#' @name Simpson_r
#' @references
#' \insertAllCited{}
NULL


#' @rdname Simpson_r
#' @return `Simpson_r` returns an object of class `fv`, see [fv.object].
#' There are methods for print and plot for this class.
#' It contains the value of the spatially explicit Simpson's entropy for each distance in `r`.
#' @export
#'
#' @examples
#' # Generate a random community
#' spCommunity <- rSpCommunity(1, size=1000, S=3)
#' # Calculate the entropy and plot it
#' autoplot(Simpson_r(spCommunity))
#' 
Simpson_r <- function(spCommunity, r = NULL, spCorrection = "isotropic", 
                      CheckArguments = TRUE) {
  
  if (CheckArguments)
    CheckSpatDivArguments()
  
  fvCorrection <- function(x) {
    switch(spCorrection ,
           "isotropic" = x$iso,
           "translate" = x$trans,
           "none" = x$un
    )
  }
  
  # Summary
  Ns <- tapply(spCommunity$marks$PointType, spCommunity$marks$PointType, length)
  Ns <- Ns[!is.na(Ns)]
  Nall <- sum(Ns)
  Ps <- Ns/Nall
  
  # r
  if (is.null(r)) {
    rMax <- spatstat.geom::diameter(spCommunity$window)
    r <- rMax * c(0, 1:20, seq(22, 40, 2), seq(45, 100, 5), 
                  seq(110, 200, 10), seq(220, 400, 20))/800
  }
  # K all points
  Kall <- fvCorrection(spatstat.core::Kest(spCommunity, r=r, correction=spCorrection))
  
  # The point pattern is separated into a list of ppp for each mark
  pppList <- split(spCommunity, as.factor(spCommunity$marks$PointType))
  # K for each ppp
  KList <- lapply(pppList, spatstat.core::Kest, r=r, correction=spCorrection)
  Ks <- as.data.frame(lapply(KList, fvCorrection))
  # Ks is NA for species with a a single point. Should be 0
  Ks[is.na(Ks)] <- 0
  
  # Result
  Shr <- (1 - rowSums((Ks*rep(Ns*(Ns-1), each=dim(Ks)[1]))/
                        (Kall*Nall*(Nall-1)))) * (Nall-1)/Nall
  
  # Build a dataframe with r, theoretical value and S(r)
  ShiEstimate <- data.frame(r, as.numeric(entropart::Simpson(Ns)), Shr)
  ColNames <- c("r", "Simpson", "S_r")
  colnames(ShiEstimate) <- ColNames
  
  # Return the values of Shimatani(r)
  Labl <- c("r", "hat(%s)", "hat(%s)(r)")
  Desc <- c("Distance argument r", "Asymptotic %s", "Estimated %s")
  S <- spatstat.core::fv(
    ShiEstimate, 
    argu="r", 
    ylab=quote(Shimatani(r)), 
    valu="S_r", 
    fmla= ". ~ r", 
    alim=c(0, max(r)), 
    labl=Labl, 
    desc=Desc, 
    unitname=spCommunity$window$unit, 
    fname="Simpson's Entropy")
  spatstat.core::fvnames(S, ".") <- ColNames[-1]
  return (S)
}



#' @rdname Simpson_r
#' @param NumberOfSimulations The number of simulations to run, 100 by default.
#' @param Alpha The risk level, 5% by default.
#' @param SimulationType A string describing the null hypothesis to simulate. 
#' The null hypothesis may be "RandomPosition": points are drawn in a Poisson process (default); 
#' "RandomLabeling": randomizes point types, keeping locations unchanged.
#' @param Global Logical; if `TRUE`, a global envelope sensu \insertCite{Duranton2005}{SpatDiv} is calculated.
#'
#' @return `Simpson_rEnvelope` returns an envelope object [envelope]. 
#' There are methods for print and plot for this class.
#' It contains the observed value of the function, its average simulated value and the confidence envelope.
#' @export
#'
#' @examples
#' # Generate a random community
#' spCommunity <- rSpCommunity(1, size=1000, S=3)
#' # Calculate the entropy and plot it
#' autoplot(Simpson_rEnvelope(spCommunity, NumberOfSimulations=10))
#' 
Simpson_rEnvelope <- function(spCommunity, r = NULL, NumberOfSimulations = 100, 
                              Alpha = 0.05, spCorrection = "isotropic", 
                              SimulationType = "RandomLabeling", Global = FALSE, 
                              CheckArguments = TRUE) {
    
  if (CheckArguments)
    CheckSpatDivArguments()
  
  # Choose the null hypothesis
  SimulatedPP <- switch (SimulationType,
                         RandomPosition = expression(dbmss::rRandomPositionK(spCommunity, CheckArguments = FALSE)),
                         RandomLabeling = expression(dbmss::rRandomLabeling(spCommunity, CheckArguments = FALSE))
  )
  if (is.null(SimulatedPP))
    stop(paste("The null hypothesis", sQuote(SimulationType), "has not been recognized."))
  # local envelope, keep extreme values for lo and hi (nrank=1)
  Envelope <- spatstat.core::envelope(spCommunity, fun=Simpson_r, nsim=NumberOfSimulations, nrank=1,
                                      r=r, spCorrection=spCorrection,
                                      CheckArguments = FALSE,
                                      simulate=SimulatedPP, savefuns=TRUE
                                      )
  attr(Envelope, "einfo")$H0 <- switch (SimulationType,
                                        RandomPosition = "Random Position",
                                        RandomLabeling = "Random Labeling"
                                        )
  # Calculate confidence intervals
  Envelope <- dbmss::FillEnvelope(Envelope, Alpha, Global)
  # Return the envelope
  return (Envelope)
}
