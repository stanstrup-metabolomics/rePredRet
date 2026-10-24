# Local Scale of Residuals

Models the absolute residuals as a smooth function of x.

## Usage

``` r
local_residual_scale(x, e, w, k, floor_frac = 0.3)
```

## Arguments

- x:

  Source RTs.

- e:

  Residuals.

- w:

  Observation weights (sigmoid outlier weights).

- k:

  Basis dimension of the mean fit; the scale smooth uses at most 5.

- floor_frac:

  The scale is floored at this fraction of the global weighted scale
  (default 0.3).

## Value

A function returning the residual standard deviation at new x.

## Details

A Gamma GAM with log link is fitted to \|e\| with the weights `w`. Under
normal errors E\|e\| = sigma \* sqrt(2 / pi), so the fitted values are
scaled by sqrt(pi / 2). Absolute residuals are floored at 5% of their
mean, since exact zeros (e.g. RTs reported with few decimals) are not
allowed by the Gamma family. If the GAM fails, the global weighted scale
is used.

The GAM is only evaluated at the compounds that carry weight (`w` at
least 5% of its maximum) and interpolated linearly between them; beyond
them the scale is held constant. In gaps between compounds and next to
down-weighted outliers the log-link smooth is unconstrained and can grow
without bound (scales above 1e30 min were seen). The scale is capped at
sqrt(pi / 2) times the largest absolute residual.
