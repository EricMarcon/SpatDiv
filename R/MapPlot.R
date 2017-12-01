#' MapPlot
#'
#' Methods to map spatialized communities.
#'
#' @title Map Spatialized Communities
#' @param x An object to map.
#' @param ... Further parameters.
#' @export
MapPlot <-
function (x, ...) 
{
  UseMethod("MapPlot")
}



#' @rdname MapPlot
#' @return \code{MapPlot.Accumulation} returns an \code{\link{autoKrige}} object that can be used to produce alternative maps.
#' @importFrom automap autoKrige
#' @importFrom sp SpatialPoints
#' @importFrom spatstat inside.owin
#' @importFrom spatstat inside.owin
#' @method MapPlot Accumulation
#' @export
MapPlot.Accumulation <-
function (x, Order, NeighborHood, Nbx = 128, Nby = 128, Contour = TRUE, 
          Palette = topo.colors(128, alpha=1), ContourColor = "dark red",
          ..., CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckSpatDivArguments()
  
  if (is.null(dim(x$Neighborhoods)))
    stop("The Accumulation object does not contain individual data to plot.")
  
  # Convert the SpCommunity to a SpatialPointsDataFrame
  sdfCommunity <- SpatialPointsDataFrame(coords=data.frame(x=x$SpCommunity$x, y=x$SpCommunity$y), 
                                         data=data.frame(Accumulation=x$Neighborhoods[Order, NeighborHood,]))
  # Prepare a grid
  xy <- spatstat::gridcentres(x$SpCommunity, Nbx, Nby)
  ok <- spatstat::inside.owin(xy$x, xy$y, x$SpCommunity$window)
  xygrid <- sp::SpatialPoints(cbind(xy$x[ok], xy$y[ok]))
  gridded(xygrid) <- TRUE
  # Proceed to krigeing
  krigedCommunity <- automap::autoKrige(Accumulation~1, sdfCommunity, new_data=xygrid)
  # Map
  image(krigedCommunity$krige_output, col=Palette, asp=1)
  if (Contour)
    contour(krigedCommunity$krige_output, add=TRUE, col=ContourColor)

  # Return the kriged community to allow further processing
  return(invisible(krigedCommunity))
}

