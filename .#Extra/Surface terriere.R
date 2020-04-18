t<- dirichlet(rjitter(Paracou6))
plot(t)
plot(Paracou6, which=1, add = TRUE, col="red")
a <- sapply(t$tiles, area)
sum(a)
g <- Paracou6$marks$PointWeight * a /sum(a) * sum(Paracou6$marks$PointWeight)

Nbx = 128; Nby = 128 ;
Palette = grDevices::topo.colors(128, alpha=1) ;
Contournlevels = 10 ; Contourcol = "dark red"
# Convert the SpCommunity to a SpatialPointsDataFrame
sdfCommunity <- sp::SpatialPointsDataFrame(coords=data.frame(x=Paracou6$x, y=Paracou6$y), 
                                           data=data.frame(Accumulation=g))
# Prepare a grid
xy <- spatstat::gridcentres(Paracou6, Nbx, Nby)
ok <- spatstat::inside.owin(xy$x, xy$y, Paracou6$window)
xygrid <- sp::SpatialPoints(cbind(xy$x[ok], xy$y[ok]))
sp::gridded(xygrid) <- TRUE
# Proceed to krigeing
krigedCommunity <- automap::autoKrige(Accumulation~1, sdfCommunity, new_data=xygrid)
# Map
graphics::image(krigedCommunity$krige_output, col=Palette, asp=1)
graphics::contour(krigedCommunity$krige_output, add=TRUE, nlevels=Contournlevels, col=Contourcol)

# Alternative
d <- density(Paracou6, weights = Paracou6$marks$PointWeight, leaveoneout=FALSE)
plot(d)
