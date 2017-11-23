#' Entropy Accumulation
#'
#' @param spCommunity A spatialized community (An object of class "wmppp" (\code{\link{wmppp.object}}), with \code{PointType} values as species names.)
#' @param q.seq A numeric vector: the sequence of diversity orders to address. Default is from 0 to 2.
#' @param n A vector of integers. Entropy will be accumulated along this number of neighbors around each individual. Default is 10\% of the individuals.
#' @param r A vector of distances. If \code{NULL} accumulation is along \code{n}, else neighbors are accumulated in circles of radius \code{r}.
#' @inheritParams as.ProbaVector.wmppp
#'
#' @return A matrix containing average entropy. Columns
#'
#' @importFrom spatstat nnwhich
#' @importFrom parallel mclapply
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @export
#'
#' @examples TODO
EntAccum <-
function(spCommunity, q.seq = seq(0,2,by=0.1), Correction = "None", n.seq = 1:ceiling(spCommunity$n/10), r = NULL, CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckSpatDivArguments()

  if (is.null(r)) {
    # n nearest neighbors. Find them
    nNeighbors <- spatstat::nnwhich(spCommunity, k=n.seq)
    # Add the reference point to get a table: center points in line, neoghbors in columns, including the point itself in the first column
    nNeighbors <- cbind(Reference=1:spCommunity$n, nNeighbors)

    # Compute the average entropy of nNeighbors
    NeighborhoodEntropies <- function (nNeighbors) {
      # Neighbor communities: 1 community per column
      NeighborCommunities <- apply(nNeighbors, 1, function(Neighbors) spCommunity$marks$PointType[Neighbors])
      # Calculate entropy of each community and all q values
      qCommunityEntropies <-  apply(NeighborCommunities, 2, function(Community) sapply(q.seq, function(q) Tsallis(Community, q=q, Correction=Correction)))
      # Mean entropy. If qCommunityEntropies is a vector (i.e. a single value of q is provided), transpose it to get a 1-row matrix.
      if (is.null(dim(qCommunityEntropies))) qCommunityEntropies <- t(qCommunityEntropies)
      return(apply(t(t(qCommunityEntropies)), 1, mean))
    }

    # Apply NeighborhoodEntropies to increasing number of neighbors, starting from 1 (i.e. 2 individuals including the reference point)
#    qEntropies <- parallel::mclapply(n.seq, function(k) NeighborhoodEntropies(nNeighbors[, 1:(k+1)]))
#    qEntropies <- simplify2array(qEntropies)
    ProgressBar <- utils::txtProgressBar(min=0, max=length(n.seq))
    qEntropies <- matrix(nrow=length(q.seq), ncol=1+length(n.seq))
    for (k in n.seq) {
      qEntropies[, k+1] <- NeighborhoodEntropies(nNeighbors[, 1:(k+1)])
      utils::setTxtProgressBar(ProgressBar, k)
    }
    # Entropy of a single individual is 0
    qEntropies[, 1] <- 0
    # Name lines and columns
    rownames(qEntropies) <- q.seq
    colnames(qEntropies) <- c(1, 1+n.seq)

    class(qEntropies) <- c("EntAccum", class(qEntropies))
    return(qEntropies)

  } else {
    # neighbors up to distance r
  }
}


#' Diversity Accumulation
#'
#' @inheritParams EntAccum
#' @param H0 A spatialized community (An object of class "wmppp" (\code{\link{wmppp.object}}), with \code{PointType} values as species names.)
#' @param Alpha The risk level of the envelope of the null hypothesis. Default is 5\%.
#' @param Simulations The number of bootstraps to build confidence intervals. Default is 50.
#'
#' @return A 3-dimensional array containing diversity accumulation curves.

#' @importFrom iNEXT iNEXT
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @export
#'
#' @examples TODO
DivAccum <-
function(spCommunity, q.seq = seq(0,2,by=0.1), Correction = "None", n.seq = 1:ceiling(spCommunity$n/10), r = NULL,
         H0 = FALSE, Alpha = 0.05, Simulations = 50,
         CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckSpatDivArguments()

  # Prepare an array to store data
  if (H0) {
    qDiversities <- array(NA, dim=c(length(q.seq), 1+length(n.seq), 4), dimnames = list(q=q.seq, n=c(1, 1+n.seq), Values=c("Observed", "Theoretical", "Lower bound", "Upper bound")))
  } else {
    qDiversities <- array(NA, dim=c(length(q.seq), 1+length(n.seq), 1), dimnames = list(q=q.seq, n=c(1, 1+n.seq), Values="Observed"))
  }
  # Get entropy
  qDiversities[, , 1] <- EntAccum(spCommunity=spCommunity, q.seq=q.seq, Correction=Correction, n.seq=n.seq, r=r, CheckArguments=FALSE)

  # Calculate Hill Numbers, by line, and the null distribution
  if (H0) {
    ProgressBar <- utils::txtProgressBar(min=0, max=length(q.seq))
    Ns <- as.numeric(table(spCommunity$marks$PointType))
  }
  for (i in 1:length(q.seq)) {
    qDiversities[i, , 1] <- entropart::expq(qDiversities[i, , 1], q.seq[i])
    if (H0) {
      H0Values <- suppressWarnings(iNEXT::iNEXT(Ns, q=as.numeric(q.seq[i]), size=c(1, 1+n.seq), conf=1-Alpha, nboot=Simulations)$iNextEst)
      qDiversities[i, , 2] <- H0Values$qD
      qDiversities[i, , 3] <- H0Values$qD.LCL
      qDiversities[i, , 4] <- H0Values$qD.UCL
      utils::setTxtProgressBar(ProgressBar, i)
    }
  }

  class(qDiversities) <- c("DivAccum", class(qDiversities))
  return(qDiversities)
}



#' Plot Diversity Accumulation
#'
#' @param A \code{\link{DivAccum}} object
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
#'
#' @examples TODO
plot.DivAccum <-
function(x, ..., q = 0,
         type = "l",  main = paste("Accumulation of Diversity of Order", q), xlab = "Sample size", ylab = "Diversity", ylim = NULL,
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
