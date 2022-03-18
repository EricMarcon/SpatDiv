# SpatDiv 0.8-3

## New features

* `MapPlot.fv()`: maps of local concentration.


# SpatDiv 0.7-6

## New features

* `rSpSpecies()`: random spatialized species.


# SpatDiv 0.6-2

## New features

* `Simpson_r()`: Spatially Explicit Simpson's Entropy.


# SpatDiv 0.5-3

## New features

* `alphahull()` to build a concave hull around a set of points.


# SpatDiv 0.4-2

## New features

* Diversity accumulation with respect to distance, with edge-effect corrections.
* Multinomial rarefaction by _entropart_.

## Improvements

* Documentation by _pkgdown_.

## External changes
 
- Updates in the _spatstat_ package: _SpatDiv_ has been updated to address the creation of _spatstat.core_ et al.


# SpatDiv 0.3-3

## New features

* Vignettes by _pkgdown_.
* `Paracou6` datasets.
* Examples.
* `Richness`, `Shannon` and `Simpson` methods.
* _ggplot_ support with `autoplot` methods.

## Improvements

* Tuned imports.
* Documentation by _packagedoc_.


# SpatDiv 0.2-1

## New features

* `MapPlot()` to make maps.
* CI by Travis.
* `Mixing()` calculates the mixing index
* `plot.Accumumation()` replaces `plot.DivAccum()` to allow printing `EntAccum`, `DivAccum` and `Mixing`objects.
* accumulation by distance (sample area) works. No edge-effect correction yet. Null hypothesis is "RandomLocation".


# SpatDiv 0.1-1

## New features

* Species distribution S3 methods : `as.AbdVector()`, `as.ProbaVector()` for `wmppp` objects, factors and character vectors.
* `Tsallis()` and `Diversity()` for `wmppp` objects, factors and character vectors.
* `rSpCommunity()`generates a spatialized random community (in a `wmppp` object). In progress: limited to a binomial point process.
* `EntAccum()` and `DivAccum()` calculate entropy and diversity accumulations curves. DAC's can be plotted by `plot.DivAccum()`
