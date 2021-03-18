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
#' @param Order The order of diversity. It can be a number or a character string, interpreted as the value of $q$.
#' @param NeighborHood The neignborhood size. It can be a number or a character string, interpreted as the value of the number of neighobors or the distance.
#' @param AllowJitter If `TRUE`, duplicated points are jittered to avoid their elimination by the krigeing procedure.
#' @param Nbx The number of columns (pixels) of the resulting map, 128 by default.
#' @param Nby The number of rows (pixels) of the resulting map, 128 by default.
#' @param Palette The color palette of the map.
#' @param Contour If `TRUE`, contours are added to the map.
#' @param Contournlevels The number of levels of contours.
#' @param Contourcol The color of the contour lines.
#'
#' @return `MapPlot.Accumulation` returns an [autoKrige] object that can be used to produce alternative maps.
#' 
#' @method MapPlot Accumulation
#' @export
#' 
#' @examples
#' # Generate a random community
#' spCommunity <- rSpCommunity(1, size=50, S=10)
#' # Calculate the species accumulation curve
#' accum <- DivAccum(spCommunity, q.seq=c(0,1), n.seq=4, Individual=TRUE)
#' # Plot the local richness, accumulated up to 5 individuals.
#' MapPlot(accum, Order=0, NeighborHood=5)
#' 
#' @section TODO: MapPlot fails if q.seq is a scalar
MapPlot.Accumulation <-
function (x, Order, NeighborHood, AllowJitter = TRUE,
          Nbx = 128, Nby = 128, Contour = TRUE, 
          Palette = grDevices::topo.colors(128, alpha=1), 
          Contournlevels = 10, Contourcol = "dark red",
          ..., CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckSpatDivArguments()
  
  if (is.null(dim(x$Neighborhoods)))
    stop("The Accumulation object does not contain individual data to plot.")
  
  # Jitter
  if (AllowJitter) {
    # Find duplicates
    Dups <- spatstat.geom::duplicated.ppp(x$SpCommunity, rule="unmark")
    if (sum(Dups)>0) {
      # Extract the duplicates and jitter them
      Dupswmppp <- spatstat.geom::rjitter(x$SpCommunity[Dups])
      # Put the coordinates back into the original wmppp
      x$SpCommunity$x[Dups] <- Dupswmppp$x
      x$SpCommunity$y[Dups] <- Dupswmppp$y
    }
  }
  
  # Convert numeric values of Order and Neighborhood into their index
  if (is.numeric(Order)) Order <- which(as.numeric(rownames(x$Neighborhoods)) == Order)
  if (is.numeric(NeighborHood)) NeighborHood <- which(as.numeric(colnames(x$Neighborhoods)) == NeighborHood)
  # Verify that values exist
  if (length(dim(x$Neighborhoods[Order, , ])) != 2) stop("Incorrect Order.") 
  if (length(dim(x$Neighborhoods[, NeighborHood, ])) != 2) stop("Incorrect Neighborhood.") 

  # Convert the SpCommunity to a SpatialPointsDataFrame
  sdfCommunity <- sp::SpatialPointsDataFrame(coords=data.frame(x=x$SpCommunity$x, y=x$SpCommunity$y), 
                                         data=data.frame(Accumulation=x$Neighborhoods[Order, NeighborHood,]))
  # Prepare a grid
  xy <- spatstat.geom::gridcentres(x$SpCommunity, Nbx, Nby)
  ok <- spatstat.geom::inside.owin(xy$x, xy$y, x$SpCommunity$window)
  xygrid <- sp::SpatialPoints(cbind(xy$x[ok], xy$y[ok]))
  sp::gridded(xygrid) <- TRUE
  # Proceed to krigeing
  krigedCommunity <- automap::autoKrige(Accumulation~1, sdfCommunity, new_data=xygrid)
  # Map
  graphics::image(krigedCommunity$krige_output, col=Palette, asp=1)
  if (Contour)
    graphics::contour(krigedCommunity$krige_output, add=TRUE, nlevels=Contournlevels, col=Contourcol)

  # Return the kriged community to allow further processing
  return(invisible(krigedCommunity))
}

