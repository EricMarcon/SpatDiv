# SpatDiv 0.3-2

## New features

* `Richness`, `Shannon` and `Simpson` methods.

# SpatDiv 0.3-1

## New features

* _ggplot_ support with `autoplot` methods.

## Improvements

* Tuned imports.


# SpatDiv 0.3-0

## New features

* Documentation by _packagedoc_.


# SpatDiv 0.2-1

## New features

* `MapPlot()` to make maps.

* CI by Travis.


# SpatDiv 0.1-2

## New features

* `Mixing()` calculates the mixing index

* `plot.Accumumation()` replaces `plot.DivAccum()` to allow printing `EntAccum`, `DivAccum` and `Mixing`objects.

* accumulation by distance (sample area) works. No edge-effect correction yet. Null hypothesis is "RandomLocation".


# SpatDiv 0.1-1

## New features

* Species distribution S3 methods : `as.AbdVector()`, `as.ProbaVector()` for `wmppp` objects, factors and character vectors.

* `Tsallis()` and `Diversity()` for `wmppp` objects, factors and character vectors.

* `rSpCommunity()`generates a spatialized random community (in a `wmppp` object). In progress: limited to a binomial point process.

* `EntAccum()` and `DivAccum()` calculate entropy and diversity accumulations curves. DAC's can be plotted by `plot.DivAccum()`
