#' alpha-shape calculation
#'
#' Calculate a window containing all points of a point pattern.
#' The window is not convex but as close as possible to the points.
#' 
#' The typical use of this function is to define a narrow window around a point pattern that has been created with a default, rectangle window.
#' 
#' The window is built by the [alphahull::ashape()] function and then transformed into a [spatstat.geom::owin.object].
#' The `alpha` parameter determines the smallest size of zones excluded from the window.
#' If it is not specified, a first attempt is 1/256 of the diameter of the existing window of `X`.
#' If the shape cannot be calculated, `alpha` is doubled and a new attempt is made.
#' 
#' @param X A planar point pattern ([spatstat.geom::ppp.object]).
#' @param alpha A smoothing parameter to delimit concave polygons.
#' @param CheckArguments If `TRUE` (default), the function arguments are verified. 
#' Should be set to `FALSE` to save time in simulations for example, when the arguments have been checked elsewhere.
#'
#' @return A window, i.e. a [spatstat.geom::owin.object].
#' @export
#' 
#' @examples
#' # Simulate a point pattern
#' if (require(spatstat.random)) {
#'   X <- rpoispp(10)
#'   plot(X)
#'   # Calculate its border
#'   X$window <- alphahull(X)
#'   plot(X)
#' }
alphahull <- function(X, alpha = NULL, CheckArguments = TRUE) {
  if (CheckArguments)
    CheckSpatDivArguments()
  
  # At least 3 points are needed
  if (X$n < 3) 
    return(X$window)
  if (X$n == 3)
    # Window is a convex hull
    return(spatstat.geom::convexhull.xy(X$x, X$y))
  
  # Prepare 
  is_validated_alpha <- FALSE
  # Unique points are needed
  xy <- unique(data.frame(x=X$x, y=X$y))
  # Size of the window
  Diameter <- spatstat.geom::diameter(X$window)
  # Use alphahull::ashape to obtain a concave hull. 
  # Parameter alpha must be small for accuracy (argument alpha)
  if (is.null(alpha)) {
    alpha <- Diameter/256 
  }
  # but large enough to be able to build a correct graph from the points: 
  # alpha is multiplied by 2 until success.
  while (!is_validated_alpha){
    AlphaShape <- alphahull::ashape(xy, alpha=alpha)
    # Convert alpha shape into polygon (https://rpubs.com/geospacedman/alphasimple)
    # Make a graph with edges, library igraph
    AlphaShape_graph <- igraph::graph.edgelist(cbind(
      as.character(AlphaShape$edges[, "ind1"]), 
      as.character(AlphaShape$edges[, "ind2"])), 
      directed = FALSE)
    # Tests: the graph must be connected and circular. If it is not, increase alpha.
    Error <- ""
    if (AlphaShape$length == 0) {
      Error <- "No edges in alpha shape"
    } else if (!igraph::is.connected(AlphaShape_graph)) {
      Error <- "Graph not connected"
    } else if (any(igraph::degree(AlphaShape_graph) != 2)) {
      Error <- "Graph not circular"
    } else if (igraph::clusters(AlphaShape_graph)$no > 1) {
      Error <- "Graph composed of more than one circle"
    }
    if (Error == "") {
      is_validated_alpha <- TRUE
    } else {
      if (alpha > Diameter) {
        # Unable to make a circular graph: give up.
        warning(paste("Unable to build an alpha hull:", Error))
      }
      else # Try to double alpha
        alpha <- 2*alpha
    }
  }
  
  # Eliminate the first node to destroy circularity
  Cut_graph <- AlphaShape_graph - igraph::E(AlphaShape_graph)[1]
  # Find chain end points
  ends <- names(which(igraph::degree(Cut_graph) == 1))
  path <- igraph::get.shortest.paths(Cut_graph, ends[1], ends[2])[[1]]
  # This is an index into the points
  pathX <- as.numeric(igraph::V(Cut_graph)[unlist(path)]$name)
  # Join the ends to restore circularity
  pathX = c(pathX, pathX[1])
  
  # Get the points from the ashape object, make an owin. Manage reverse by tryCatch
  the_window <- tryCatch(
    spatstat.geom::owin(poly=list(x=AlphaShape$x[pathX, ][, 1],y=AlphaShape$x[pathX, ][, 2])),
    error=function(e) # Error if the polygon is traversed clockwise
      spatstat.geom::owin(poly=list(x=AlphaShape$x[rev(pathX), ][, 1],y=AlphaShape$x[rev(pathX), ][, 2]))
  )
  
  return(the_window)
}
