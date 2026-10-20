# Get Model Details

Retrieves a pre-built model with its prediction grid and calibration
data.

## Usage

``` r
rePredRet_model(from_id, to_id, refresh = FALSE)
```

## Arguments

- from_id:

  Source system identifier.

- to_id:

  Target system identifier.

- refresh:

  If TRUE, download a fresh copy instead of using the cache.

## Value

A list with:

- sys1_id, sys2_id:

  Source and target system

- n_compounds:

  Number of calibration compounds

- method:

  Interval method

- predictions:

  Data frame with x, pred, lower, upper

- calibration:

  Data frame with x, y, compound

- stats:

  List with median_ci_width, median_error, r_squared, coverage (fraction
  of calibration compounds inside the PI), rel_mean_ci_width (mean PI
  width / target RT range) and build_time. Models built before these
  fields were added lack them.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- rePredRet_model("0001", "0002")
} # }
```
