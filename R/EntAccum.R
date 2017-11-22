#' Entropy Accumulation
#'
#' @param spCommunity A spatialized community (An object of class "wmppp" (\code{\link{wmppp.object}}), with \code{PointType} values as species names.)
#' @param q.seq A numeric vector: the sequence of diversity orders to address. Default is from 0 to 2.
#' @param n A vector of integers. Entropy will be accumulated along this number of neighbors around each individual. Default is 10\% of the individuals.
#' @param r A vector of distances. If \code{NULL} accumulation is along \code{n}, else neighbors are accumulated in circles of radius \code{r}.
#' @inheritParams as.ProbaVector.wmppp
#'
#' @return A matrix containing average entropy. Columns
#' @export
#' @importFrom spatstat nnwhich
#' @importFrom parallel mclapply
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
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

    return(qEntropies)

  } else {
    # neighbors up to distance r
  }
}


#' Diversity Accumulation
#'
#' @inheritParams EntAccum
#'
#' @return A matrix containing diversity accumulation curves
#' @export
#'
#' @examples TODO
DivAccum <-
  function(spCommunity, q.seq = seq(0,2,by=0.1), Correction = "None", n.seq = 1:ceiling(spCommunity$n/10), r = NULL, CheckArguments = TRUE)
  {
    if (CheckArguments)
      CheckSpatDivArguments()

    # Get entropy
    qDiversities <- EntAccum(spCommunity=spCommunity, q.seq=q.seq, Correction=Correction, n.seq=n.seq, r=r, CheckArguments=FALSE)

    # Calculate Hill Numbers, by line
    for (i in 1:length(q.seq)) {
      qDiversities[i, ] <- entropart::expq(qDiversities[i, ], q.seq[i])
    }

    return(qDiversities)
  }
