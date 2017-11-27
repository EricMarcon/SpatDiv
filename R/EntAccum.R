#' Entropy Accumulation
#'
#' @param spCommunity A spatialized community (An object of class "wmppp" (\code{\link{wmppp.object}}), with \code{PointType} values as species names.)
#' @param divCorrection A string containing one of the possible corrections to calculate diversity, see \code{\link{Tsallis}} for all possible values. \code{"None"} uses the plugin estimator. \code{"Best"} is the default value.
#' @param q.seq A numeric vector: the sequence of diversity orders to address. Default is from 0 to 2.
#' @param n.seq A vector of integers. Entropy will be accumulated along this number of neighbors around each individual. Default is 10\% of the individuals.
#' @param r.seq A vector of distances. If \code{NULL} accumulation is along \code{n}, else neighbors are accumulated in circles of radius \code{r}.
#' @param spCorrection A string containing one of the possible corrections for edge effects. When a point is close to the window's border, its neighborhood is underestimated. "None" is for no correction.
#' @param correction The edge-effect correction to apply when estimating the entropy of a neighborhood community that does not fit in the window. 
#'        Does not apply if neihborhoods are defined by the number of neighbors. Default is "None".
#' @param ShowProgressBar If `TRUE` (default), a progress bar is shown.
#' @inheritParams as.ProbaVector.wmppp
#'
#' @return A matrix containing average entropy.
#'
#' @importFrom spatstat nnwhich
#' @importFrom parallel mclapply
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @export
#'
#' @examples #TODO
EntAccum <-
function(spCommunity, q.seq = seq(0,2,by=0.1), divCorrection = "None", n.seq = 1:ceiling(spCommunity$n/10), r.seq = NULL, spCorrection = "None", 
         ShowProgressBar = TRUE, CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckSpatDivArguments()

  if (is.null(r.seq)) {
    # n nearest neighbors. Find them.
    nNeighbors <- spatstat::nnwhich(spCommunity, k=n.seq)
    # Add the reference point to get a table: center points in line, neoghbors in columns, including the point itself in the first column
    nNeighbors <- cbind(Reference=1:spCommunity$n, nNeighbors)

    # Compute the average entropy of nNeighbors
    nNeighborhoodEntropies <- function (nNeighbors) {
      # Neighbor communities: 1 community per column
      NeighborCommunities <- apply(nNeighbors, 1, function(Neighbors) spCommunity$marks$PointType[Neighbors])
      # Calculate entropy of each community and all q values
      qCommunityEntropies <-  apply(NeighborCommunities, 2, function(Community) sapply(q.seq, function(q) Tsallis(Community, q=q, Correction=divCorrection)))
      # Mean entropy. If qCommunityEntropies is a vector (i.e. a single value of q is provided), transpose it to get a 1-row matrix.
      if (is.null(dim(qCommunityEntropies))) qCommunityEntropies <- t(qCommunityEntropies)
      return(apply(t(t(qCommunityEntropies)), 1, mean))
    }

    # Prepare a progress bar and the result matrix
    if (ShowProgressBar) ProgressBar <- utils::txtProgressBar(min=0, max=length(n.seq))
    qEntropies <- matrix(nrow=length(q.seq), ncol=1+length(n.seq))
    for (k in n.seq) {
      # At each number of neighbors, calculate the entropy of all points' neighborhood for each q
      qEntropies[, k+1] <- nNeighborhoodEntropies(nNeighbors[, 1:(k+1)])
      if (ShowProgressBar) utils::setTxtProgressBar(ProgressBar, k)
    }
    # Entropy of a single individual is 0
    qEntropies[, 1] <- 0
    # Name lines and columns
    rownames(qEntropies) <- q.seq
    colnames(qEntropies) <- c(1, 1+n.seq)

  } else {
    # neighbors up to distance r. Distances are in r.seq. The first distance is 0 (verified by CheckArguments)

    # The max value of the factors is needed
    NbSpecies <- max(as.integer(spCommunity$marks$PointType))
    # Run C++ routine to fill a 3D array. Rows are points, columns are r, the 3rd dimension has a z-value per species. Values are the number (weights) of neighbors of each point, up to ditance r, of species z.
    rNeighbors <- parallelCountNbd(r=r.seq, NbSpecies,
                         x=spCommunity$x, y=spCommunity$y,
                         Type=spCommunity$marks$PointType, Weight=spCommunity$marks$PointWeight)
    # The array of neighbor communities is built from the vector returned.
    dim(rNeighbors) <- c(spCommunity$n, length(r.seq), NbSpecies)

    # Compute the average entropy of neighborhoods.
    rNeighborhoodEntropies <- function (rNeighbors) {
      # Neighbor communities: 1 community per column
      NeighborCommunities <- sapply(1:NbSpecies, function(i) colSums(rNeighbors[i, , ]))
      # Calculate entropy of each community and all q values
      if (spCorrection == "None") {
        # No edge-effect correction
        qCommunityEntropies <-  apply(NeighborCommunities, 2, function(Community) sapply(q.seq, function(q) Tsallis(Community, q=q, Correction=divCorrection)))
      } else {
        # TODO
      }
      # Mean entropy. If qCommunityEntropies is a vector (i.e. a single value of q is provided), transpose it to get a 1-row matrix.
      if (is.null(dim(qCommunityEntropies))) qCommunityEntropies <- t(qCommunityEntropies)
      return(apply(t(t(qCommunityEntropies)), 1, mean))
    }

    # Prepare a progress bar and the result matrix
    if (ShowProgressBar) ProgressBar <- utils::txtProgressBar(min=1, max=length(r.seq))
    qEntropies <- matrix(nrow=length(q.seq), ncol=length(r.seq))
    for (k in 2:length(r.seq)) {
      # At each distance, calculate the entropy of all points' neighborhood for each q
      qEntropies[, k] <- rNeighborhoodEntropies(rNeighbors[, 1:k, ])
      if (ShowProgressBar) utils::setTxtProgressBar(ProgressBar, k)
    }
    # Entropy at r=0 is 0
    qEntropies[, 1] <- 0
    # Name lines and columns
    rownames(qEntropies) <- q.seq
    colnames(qEntropies) <- r.seq

  }

  class(qEntropies) <- c("EntAccum", "Accumulation", class(qEntropies))
  return(qEntropies)
}


#' Diversity Accumulation
#'
#' @inheritParams EntAccum
#' @param H0 The null hypothesis to compare the distribution of `spCommunity` to. If "none", the default value, no null hypothesis is tested.
#'        "Binomial" means the community will be rarefied down to the number of neighbors of `n.seq`.
#'        "RandomLocation" means the points will we randomly permuted accross their actual locations.
#' @param Alpha The risk level of the envelope of the null hypothesis. Default is 5\%.
#' @param Simulations The number of bootstraps to build confidence intervals. Default is 50.
#'
#' @return A 3-dimensional array containing diversity accumulation curves.

#' @importFrom dbmss rRandomLocation
#' @importFrom iNEXT iNEXT
#' @importFrom spatstat area
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @export
#'
#' @examples #TODO
DivAccum <-
function(spCommunity, q.seq = seq(0,2,by=0.1), divCorrection = "None", n.seq = 1:ceiling(spCommunity$n/10), r.seq = NULL, spCorrection = "None",
         H0 = "None", Alpha = 0.05, Simulations = 50,
         ShowProgressBar = TRUE, CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckSpatDivArguments()

  # Prepare an array to store data
  if (is.null(r.seq)) {
    # Neighborhoods defined as the number of neighbors + a comlumn for no neighbor.
    nCols <- 1+length(n.seq)
    seq <- c(1, 1+n.seq)
  } else {
    # Neighborhoods defined by distances. The first distance in r.seq is 0.
    nCols <- length(r.seq)
    seq <- r.seq
  }
  if (H0 == "none") {
    # No H0: qDiversity contains q rows, r in columns and a single z-value: the observed data.
    qDiversities <- array(0.0, dim=c(length(q.seq), nCols, 1), dimnames = list(q=q.seq, nr=seq, Values="Observed"))
  } else {
    # 4 z-values: observed, expected under H0, lower and upper bounds of H0.
    qDiversities <- array(0.0, dim=c(length(q.seq), nCols, 4), dimnames = list(q=q.seq, nr=seq, Values=c("Observed", "Theoretical", "Lower bound", "Upper bound")))
  }
  # Get entropy
  qDiversities[, , 1] <- EntAccum(spCommunity=spCommunity, q.seq=q.seq, divCorrection=divCorrection, n.seq=n.seq, r.seq=r.seq, spCorrection=spCorrection, 
                                  ShowProgressBar=(ShowProgressBar & (H0 == "None" | H0 == "Binomial")), CheckArguments=FALSE)

  # Calculate Hill Numbers, by row
  for (i in 1:length(q.seq)) {
    # Transform entropy to diversity, by row (where q does not change)
    qDiversities[i, , 1] <- entropart::expq(qDiversities[i, , 1], q.seq[i])
  }
  
  # Null distribution
  if (H0 == "Binomial") {
    # Rarefy the community with iNEXT
    if (!is.null(r.seq)) stop("The 'Binomial' null hypothesis only applies to accumulation by number of neighbors.")
    # Prepare a progress bar 
    if (ShowProgressBar) ProgressBar <- utils::txtProgressBar(min=0, max=length(q.seq))
    # Prepare the distribution of the abundances of species.
    Ns <- as.numeric(table(spCommunity$marks$PointType))
    # Sample size of the neighborhoods
    Sizes <- c(1, 1+n.seq)
    for (i in 1:length(q.seq)) {
      # Rarefy the who community to the sizes of neighborhoods
      H0Values <- suppressWarnings(iNEXT::iNEXT(Ns, q=as.numeric(q.seq[i]), size=Sizes, conf=1-Alpha, nboot=Simulations)$iNextEst)
      # Extract the results from the object returnes by iNEXT.
      qDiversities[i, , 2] <- H0Values$qD
      qDiversities[i, , 3] <- H0Values$qD.LCL
      qDiversities[i, , 4] <- H0Values$qD.UCL
      if (ShowProgressBar) utils::setTxtProgressBar(ProgressBar, i)
    }
  }
  if (H0 == "RandomLocation") {
    if (is.null(r.seq)) stop("The 'RandomLocation' null hypothesis only applies to accumulation by distance.")
    # Prepare a progress bar 
    if (ShowProgressBar) ProgressBar <- utils::txtProgressBar(min=0, max=Simulations)
    # Prepare a 3-D array to store results. Rows are q, columns are r, z-values are for each simulation.
    H0qDiversities <- array(0.0, dim=c(length(q.seq), nCols, Simulations))
    # Simulate communities according to H0
    for (i in (1:Simulations)) {
      # Random community
      H0spCommunity <- dbmss::rRandomLocation(spCommunity, CheckArguments=FALSE)
      # Calculate its accumulated diversity
      H0qDiversities[, , i] <- DivAccum(H0spCommunity, q.seq=q.seq, divCorrection=divCorrection, n.seq=NULL, r.seq=r.seq, spCorrection=spCorrection, H0="None", ShowProgressBar=FALSE, CheckArguments=FALSE)[, , 1]
      if (ShowProgressBar) utils::setTxtProgressBar(ProgressBar, i)
    }
    # Calculate quantiles
    for (q in 1:length(q.seq)) {
      for (r in 1:length(r.seq)) {
        qDiversities[q, r, 3:4] <- quantile(H0qDiversities[q, r, ], c(Alpha, 1-Alpha))
        qDiversities[q, r, 2] <- mean(H0qDiversities[q, r, ])
      }
    }
  }

  class(qDiversities) <- c("DivAccum", "Accumulation", class(qDiversities))
  return(qDiversities)
}



#' Mixing index
#'
#' The mixing index is the ratio of observed diversity (effective number of species) to its theoretical, null-hypothesis value.
#'
#' @inheritParams DivAccum
#'
#' @return A 3-dimensional array containing mixing index curves.
#' @export
#'
#' @examples #TODO
Mixing <-
  function(spCommunity, q.seq = seq(0,2,by=0.1), divCorrection = "None", n.seq = 1:ceiling(spCommunity$n/10), r.seq = NULL, spCorrection = "None",
           H0 = "None", Alpha = 0.05, Simulations = 50,
           ShowProgressBar = TRUE, CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckSpatDivArguments()

  # Get the diversity accumulation
  qMixing <- DivAccum(spCommunity=spCommunity, q.seq=q.seq, divCorrection=divCorrection, n.seq=n.seq, r.seq=r.seq, spCorrection=spCorrection, 
                      H0=H0, Alpha=Alpha, Simulations=Simulations, ShowProgressBar=ShowProgressBar, CheckArguments=FALSE)

  # Normalize it
  qMixing[, , 1] <- qMixing[, , 1] / qMixing[, , 2]
  qMixing[, , 3] <- qMixing[, , 3] / qMixing[, , 2]
  qMixing[, , 4] <- qMixing[, , 4] / qMixing[, , 2]
  qMixing[, , 2] <- 1

  class(qMixing) <- c("Mixing", "Accumulation", class(qMixing))
  return(qMixing)

}



#' Plot Diversity Accumulation
#'
#' @param An \code{\link{Accumulation}} object that cat be accumulation of diversity (\code{\link{DivAccum}}), entropy (\code{\link{EntAccum}}) or the Mixing index (\code{\link{Mixing}}).
#' @param ... Further plotting arguments.
#' @param q The order of Diversity
#' @param type Plotting parameter. Default is "l".
#' @param main Main title of the plot.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param ylim Limits of the Y-axis, as a vector of two numeric values.
#' @param LineWidth Width of the Diversity Accumulation Curve line.
#' @param LineWidth Width of the Diversity Accumulation Curve line.
#' @param ShadeColor The color of the shaded confidence envelope.
#' @param BorderColor The color of the borders of the confidence envelope.
#'
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom graphics lines
#' @export
#' @method plot Accumulation
#'
#' @examples TODO
plot.Accumulation <-
function(x, ..., q = 0,
         type = "l",  main = "Accumulation of ...", xlab = "Sample size", ylab = "Diversity...", ylim = NULL,
         LineWidth = 2, ShadeColor = "grey75", BorderColor = "red")
{
  # Find the row in the accumulation table
  Whichq <- which(dimnames(x)$q==q)
  if (length(Whichq) != 1)
    stop("The value of q does not correspond to any accumulation curve.")

  if (is.null(ylim)) {
    # Evaluate ylim if not set by an argument
    ymin <- min(x[Whichq, , ])
    ymax <- max(x[Whichq, , ])
  } else {
    ymin <- ylim[1]
    ymax <- ylim[2]
  }

  if (main == "Accumulation of ...") {
    # Prepare the main title
    if (inherits(x, "EntAccum")) main <- paste("Accumulation of Entropy of order", q)
    if (inherits(x, "DivAccum")) {
      if (q == 0)
        main <- "Species Accumulation Curve"
      else
        main <- paste("Accumulation of Diversity of order", q)
    }
    if (inherits(x, "Mixing")) main <- paste("Mixing index of order", q)
  }

  if (ylab == "Diversity...") {
    # Prepare Y-axis
    if (inherits(x, "EntAccum")) ylab <-"Diversity"
    if (inherits(x, "DivAccum")) {
      if (q == 0)
        ylab <- "Richness"
      else
        ylab <- "Diversity"
    }
    if (inherits(x, "Mixing")) ylab <- "Mixing index"
  }

  # Prepare the plot
  graphics::plot(x=dimnames(x)$n, y=x[Whichq, , 1], ylim=c(ymin, ymax),
                 type=type, main=main, xlab=xlab, ylab=ylab)

  # Confidence envelope
  graphics::polygon(c(rev(dimnames(x)$n), dimnames(x)$n), c(rev(x[Whichq, , 4]), x[Whichq, , 3]), col=ShadeColor, border=FALSE)
  # Add red lines on borders of polygon
  graphics::lines(dimnames(x)$n, x[Whichq, , 4], col=BorderColor, lty=2)
  graphics::lines(dimnames(x)$n, x[Whichq, , 3], col=BorderColor, lty=2)
  # Redraw the SAC
  graphics::lines(x=dimnames(x)$n, y=x[Whichq, , 1], lwd=LineWidth, ...)

}
