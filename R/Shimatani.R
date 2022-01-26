#' Title
#'
#' @param X 
#' @param r 
#' @param correction 
#' @param Biased 
#' @param CheckArguments 
#'
#' @return
#' @export
#'
#' @examples
Simpson_r <- function(X, r = NULL, correction = "isotropic", Biased = TRUE, CheckArguments = TRUE) {
  
  if (CheckArguments)
    CheckSpatDivArguments()
  
  fvCorrection <- function(x) {
    switch(correction,
           "isotropic" = x$iso,
           "translate" = x$trans,
           "none" = x$un
    )
  }
  
  # Summary
  Ns <- tapply(X$marks$PointType, X$marks$PointType, length)
  Ns <- Ns[!is.na(Ns)]
  Nall <- sum(Ns)
  Ps <- Ns/Nall
  
  # r
  if (is.null(r)) {
    rMax <- spatstat.geom::diameter(X$window)
    r <- rMax * c(0, 1:20, seq(22, 40, 2), seq(45, 100, 5), 
                  seq(110, 200, 10), seq(220, 400, 20))/800
  }
  # K all points
  Kall <- fvCorrection(spatstat.core::Kest(X, r=r, correction=correction))
  
  # The point pattern is separated into a list of ppp for each mark
  pppList <- split(X, as.factor(X$marks$PointType))
  # K for each ppp
  KList <- lapply(pppList, spatstat.core::Kest, r=r, correction=correction)
  Ks <- as.data.frame(lapply(KList, fvCorrection))
  # Ks is NA for species with a a single point. Should be 0
  Ks[is.na(Ks)] <- 0
  
  # Result
  Shr <- (1 - rowSums((Ks*rep(Ns*(Ns-1), each=dim(Ks)[1]))/(Kall*Nall*(Nall-1)))) * ifelse(Biased, 1, (Nall-1)/Nall)
  
  # Build a dataframe with r, theoretical value and S(r)
  ShiEstimate <- data.frame(r, entropart::Simpson(Ps), Shr)
  colnames(ShiEstimate) <- c("r", "Simpson", "S")
  
  # Return the values of Shimatani(r)
  return (spatstat.core::fv(ShiEstimate, argu="r", ylab=quote(Shimatani(r)), valu="S", fmla= ". ~ r", alim=c(0, max(r)), labl=c("r", "Simpson", paste("hat(%s)(r)", sep="")), desc=c("distance argument r", "Simpson", "Estimated S(r)"), unitname=X$window$unit, fname="S"))
}
