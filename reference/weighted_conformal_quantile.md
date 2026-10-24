# Weighted Conformal Quantile

The `1 - alpha` quantile of the scores, where each compound counts by
its outlier weight, so that compounds the fit treats as outliers (wrong
annotations) do not set the width of the interval.

## Usage

``` r
weighted_conformal_quantile(scores, w, alpha)
```

## Arguments

- scores:

  Normalized residuals \|e\| / s(x).

- w:

  Sigmoid outlier weights of the compounds.

- alpha:

  Significance level.

## Value

The multiplier q.

## Details

The smallest score at which the cumulative weight of the sorted scores
reaches `(1 - alpha) * (sum(w) + max(w))`; the extra `max(w)` stands for
the new compound, assumed to be well annotated. With equal weights this
is the split-conformal `ceiling((n + 1) * (1 - alpha))`-th smallest
score. If the total weight is too small to reach that level (fewer than
about `1 / alpha` well-fitting compounds), the largest score of the
compounds with at least half the maximum weight is used.

## References

Tibshirani, R. J., Barber, R. F., Candès, E. J. & Ramdas, A. (2019).
Conformal prediction under covariate shift. *Advances in Neural
Information Processing Systems* 32.
