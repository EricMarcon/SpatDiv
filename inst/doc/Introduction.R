## ----install, eval=FALSE---------------------------------------------------
#  library("devtools")
#  devtools::install_gitub("EricMarcon/SpatDiv")

## ----rSpCommunity----------------------------------------------------------
library("SpatDiv")
rSpCommunity(n=1, size=100, S=10) -> spCommunity
plot(spCommunity, which.marks = "PointType")

## ----Whittaker-------------------------------------------------------------
plot(as.AbdVector(spCommunity))

## ----plotDivAccum----------------------------------------------------------
divAccum <- DivAccum(spCommunity, n.seq = 1:50, q.seq=0:2, H0="Multinomial")
plot(divAccum, q = 1)

## ----Mixing----------------------------------------------------------------
mixing <- Mixing(spCommunity, n.seq = 1:50, q.seq=0:2, H0="Multinomial", Individual = TRUE)
plot(mixing, q = 1)

## ----MapPlot---------------------------------------------------------------
MapPlot(mixing, Order = 0, NeighborHood = 10)

