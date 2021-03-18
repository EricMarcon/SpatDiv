#### Test de rSpCommunity
library("SpatDiv")
data(Paracou618)

# 100 arbres, 10 espèces
# rCommunity(n=3, size = 100, S=10) -> Community

# Fabrication de la liste des espèces
Tree <- Paracou618.Taxonomy
# build a phylo object
phyTree <- ape::as.phylo(Tree)
# Liste des espèces dans tip.label
spList <- sample(phyTree$tip.label, size = 10)

# Fabrication d'une communauté spatiale
plot(rSpCommunity(n=1, size = 100, S=10, Species = spList) -> spCommunity, which.marks = "PointType")
plot(as.AbdVector(spCommunity))

# Plusieurs communautés
rSpCommunity(n=3, size = 100, S=10, Species = spList) -> spCommunityList
plot(spCommunityList[[2]])

# Distances: inutile
pd <- spatstat.geom::pairdist(spCommunity)
colnames(pd) <- rownames(pd) <- spCommunity$marks$PointType

# Plus proches voisins > création d'une liste de communautés.
nb <- spatstat.geom::nnwhich(spCommunity, k=1:10)
# ajouter le point de référence
nb <- cbind(Reference=1:spCommunity$n, nb)
# Extraction
nbX <- apply(nb, 1, function(Neighbors) spCommunity[Neighbors])
plot(nbX[[1]], which.marks = "PointType")
# Affichage du centre
points(nbX[[1]]$x[1], nbX[[1]]$y[1], col="red")


# Entropie
Diversity(spCommunity)

# Accumulation of entropy n
system.time(entAccum <- EntAccum(spCommunity, n.seq = 1:50, q.seq=0:2, Individual=T))

# Accumulation n
system.time(divAccum <- DivAccum(spCommunity, n.seq = 1:50, q.seq=0:2, H0="Binomial", Individual=T))
plot(divAccum, q = 0)
MapPlot(divAccum, Order = "0", NeighborHood = "10")

mixing <- Mixing(spCommunity, n.seq = 1:50, q.seq=0:2, H0="Binomial", Individual=T)
plot(mixing, q = 0)

# Accumulation r
system.time(divAccum <- DivAccum(spCommunity, r.seq=seq(0, 1, 0.1), H0="RandomLocation"))
plot(divAccum, q = 0)

#Indice de mélange
mixing <- Mixing(spCommunity, r.seq = seq(0, 1, .1), q.seq=0:2, H0="RandomLocation")
plot(mixing, q = 0)

