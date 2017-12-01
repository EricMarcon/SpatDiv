#' MapPlot
#'
#' Map Spatialized Communities
#' 
#' Map an "Accum" object.
#' 
#' @param x An object to map.
#' @param ... Further parameters.
#' @export
MapPlot <-
function (x, ...) 
{
  UseMethod("MapPlot")
}



#' @rdname MapPlot
#' 
#' @inheritParams DivAccum
#' @param Order The order of diversity. It can be an integer, interpreted as the index of the rows of the Accumulation array, or a string, interpreted as the value of $q$.
#' @param NeighborHood The neignborhood size. It can be an integer, interpreted as the index of the column of the Accumulation array, or a string, interpreted as the value of the number of neighobors or the distance.
#' @param Nbx The number of columns (pixels) of the resulting map, 128 by default.
#' @param Nby The number of rows (pixels) of the resulting map, 128 by default.
#' @param Contour If \code{TRUE}, contours are added to the map.
#' @param Palette The color palette of the map.
#' @param ContourColor The color of the contour lines.
#'
#' @return \code{MapPlot.Accumulation} returns an \code{\link{autoKrige}} object that can be used to produce alternative maps.
#' 
#' @importFrom automap autoKrige
#' @importFrom graphics contour
#' @importFrom graphics image
#' @importFrom grDevices topo.colors
#' @importFrom sp gridded
#' @importFrom sp SpatialPoints
#' @importFrom sp SpatialPointsDataFrame
#' @importFrom spatstat inside.owin
#' @importFrom spatstat inside.owin
#' @method MapPlot Accumulation
#' @export
#' 
#' @examples
#' # Generate a random community
#' spCommunity <- rSpCommunity(1, size=50, S=10)
#' # Calculate the species accumulation curve
#' accum <- DivAccum(spCommunity, q.seq=0, n.seq=4, Individual=TRUE)
#' # Plot the local richness, accumulated up to 5 individuals.
#' MapPlot(accum, Order="0", NeighborHood="5")
#' 
MapPlot.Accumulation <-
function (x, Order, NeighborHood, Nbx = 128, Nby = 128, Contour = TRUE, 
          Palette = grDevices::topo.colors(128, alpha=1), ContourColor = "dark red",
          ..., CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckSpatDivArguments()
  
  if (is.null(dim(x$Neighborhoods)))
    stop("The Accumulation object does not contain individual data to plot.")
  
  # Convert the SpCommunity to a SpatialPointsDataFrame
  sdfCommunity <- sp::SpatialPointsDataFrame(coords=data.frame(x=x$SpCommunity$x, y=x$SpCommunity$y), 
                                         data=data.frame(Accumulation=x$Neighborhoods[Order, NeighborHood,]))
  # Prepare a grid
  xy <- spatstat::gridcentres(x$SpCommunity, Nbx, Nby)
  ok <- spatstat::inside.owin(xy$x, xy$y, x$SpCommunity$window)
  xygrid <- sp::SpatialPoints(cbind(xy$x[ok], xy$y[ok]))
  sp::gridded(xygrid) <- TRUE
  # Proceed to krigeing
  krigedCommunity <- automap::autoKrige(Accumulation~1, sdfCommunity, new_data=xygrid)
  # Map
  graphics::image(krigedCommunity$krige_output, col=Palette, asp=1)
  if (Contour)
    graphics::contour(krigedCommunity$krige_output, add=TRUE, col=ContourColor)

  # Return the kriged community to allow further processing
  return(invisible(krigedCommunity))
}

