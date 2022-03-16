#' MapPlot
#'
#' Map Spatialized Communities
#' 
#' Map an "Accum" object.
#' 
#' @param x An object to map.
#' @param ... Further parameters.
#' @param AllowJitter If `TRUE`, duplicated points are jittered to avoid their elimination by the kriging procedure.
#' @param Nbx The number of columns (pixels) of the resulting map, 128 by default.
#' @param Nby The number of rows (pixels) of the resulting map, 128 by default.
#' @param Palette The color palette of the map.
#' @param Contour If `TRUE`, contours are added to the map.
#' @param Contournlevels The number of levels of contours.
#' @param Contourcol The color of the contour lines.
#' 
#' @return An [autoKrige] object that can be used to produce alternative maps.
#' 
#' @export
MapPlot <-
function (x, ...) 
{
  UseMethod("MapPlot")
}



#' @rdname MapPlot
#' 
#' @inheritParams DivAccum
#' @param Order The order of diversity, i.e. the value of $q$.
#' @param NeighborHood The neighborhood size, i.e. the number of neighbors or the distance to consider.
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
  Order <- which(as.numeric(rownames(x$Neighborhoods)) == Order)
  NeighborHood <- which(as.numeric(colnames(x$Neighborhoods)) == NeighborHood)
  # Verify that values exist: if which() did not match, we get integer(0) for Order or NeighborHood
  # then data is of length 0.
  if (length(x$Neighborhoods[Order, , ]) == 0) stop("Incorrect Order.") 
  if (length(x$Neighborhoods[, NeighborHood, ]) == 0) stop("Incorrect Neighborhood.") 

  # Convert the SpCommunity to a SpatialPointsDataFrame
  sdfCommunity <- sp::SpatialPointsDataFrame(
    coords=data.frame(x=x$SpCommunity$x, y=x$SpCommunity$y), 
    data=data.frame(Accumulation=x$Neighborhoods[Order, NeighborHood,])
    )
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



#' @rdname MapPlot
#' 
#' @param r The distance at which the function value must be considered.
#' 
#' @method MapPlot fv
#' @export
#' 
#' @examples
#' # Generate a random community
#' spCommunity <- rSpCommunity(1, size=50, S=9, Species = paste0("sp",1:9))
#' library("dbmss")
#' M_i <- Mhat(spCommunity, ReferenceType = "sp1", Individual = TRUE)
#' MapPlot(M_i, spCommunity[spCommunity$marks$PointType=="sp1"], r=0.5)
#' 
MapPlot.fv <-
function (x, spCommunity, r, AllowJitter = TRUE,
          Nbx = 128, Nby = 128, Contour = TRUE, 
          Palette = grDevices::topo.colors(128, alpha=1), 
          Contournlevels = 10, Contourcol = "dark red",
          ..., CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckSpatDivArguments()
      
  # TODO
  # Check consistency between x and spCommunity
  # spCommunity points filtered to fit x
  # Test x has value_xxx columns
  
  # Jitter
  if (AllowJitter) {
    # Find duplicates
    Dups <- spatstat.geom::duplicated.ppp(spCommunity, rule="unmark")
    if (sum(Dups)>0) {
      # Extract the duplicates and jitter them
      Dupswmppp <- spatstat.geom::rjitter(spCommunity[Dups])
      # Put the coordinates back into the original wmppp
      spCommunity$x[Dups] <- Dupswmppp$x
      spCommunity$y[Dups] <- Dupswmppp$y
    }
  }
  
  # Find the max r value of x lower than or equal to argument r
  r_to_plot <- max(x$r[x$r<=r])

  # Convert the data to a SpatialPointsDataFrame
  df <- as.data.frame(x)
  # Pivot longer. Columns are named M_1, etc. for function M
  df <- tidyr::pivot_longer(
    df,
    cols = tidyselect::starts_with(paste0(attr(x, "valu"),"_")),
    names_to = "point",
    values_to = "dbmss"
  )
  # Filter the appropriate distance
  df <- dplyr::filter(
    df,
    r == r_to_plot
  )
  # Select the function value column only
  df <- dplyr::select(
    df,
    dbmss
  )
  # Make a SpatialPointsDataFrame
  sdfCommunity <- sp::SpatialPointsDataFrame(
    coords=data.frame(x=spCommunity$x, y=spCommunity$y),
    data=df
    )
  
  # Prepare a grid
  xy <- spatstat.geom::gridcentres(spCommunity, Nbx, Nby)
  ok <- spatstat.geom::inside.owin(xy$x, xy$y, spCommunity$window)
  xygrid <- sp::SpatialPoints(cbind(xy$x[ok], xy$y[ok]))
  sp::gridded(xygrid) <- TRUE
  # Proceed to krigeing
  krigedCommunity <- automap::autoKrige(dbmss~1, sdfCommunity, new_data=xygrid)
  # Map
  graphics::image(krigedCommunity$krige_output, col=Palette, asp=1)
  if (Contour)
    graphics::contour(krigedCommunity$krige_output, add=TRUE, nlevels=Contournlevels, col=Contourcol)
  
  # Return the kriged community to allow further processing
  return(invisible(krigedCommunity))
}