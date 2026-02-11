# Check Input Data Quality

Validates InChI strings and retention times before model building.

## Usage

``` r
check_input_data(inchis, rts, fix = TRUE)
```

## Arguments

- inchis:

  Character vector of InChI strings.

- rts:

  Numeric vector of retention times.

- fix:

  If TRUE (default), attempts to fix common issues and returns corrected
  data. If FALSE, only reports issues.

## Value

If fix = TRUE, returns a list with corrected `inchis` and `rts`. If fix
= FALSE, returns a list of issues found.

## Details

Checks performed:

- Length mismatch between inchis and rts

- Invalid InChI format (must start with "InChI=")

- Missing or NA values

- Non-positive retention times

- Duplicate compounds (same InChI after stripping stereochemistry)

- Outlier RTs (\> 3 SD from mean)

## Examples

``` r
if (FALSE) { # \dontrun{
# Check data quality
check_input_data(my_inchis, my_rts, fix = FALSE)

# Fix issues automatically
fixed <- check_input_data(my_inchis, my_rts, fix = TRUE)
models <- build_user_models(fixed$inchis, fixed$rts, "RP")
} # }
```
