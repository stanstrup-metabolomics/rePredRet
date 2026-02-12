# Export Predictions to File

Exports predictions to CSV or Excel format.

## Usage

``` r
export_predictions(predictions, file, format = c("auto", "csv", "xlsx"))
```

## Arguments

- predictions:

  A UserRTPredictions object or data.frame of predictions.

- file:

  Output file path. Extension determines format if format is "auto".

- format:

  Output format: "auto" (from extension), "csv", or "xlsx".

## Value

Invisibly returns the file path.

## Details

For Excel export, the `writexl` package is required. Install with
`install.packages("writexl")` if needed.

## Examples

``` r
if (FALSE) { # \dontrun{
# Export to CSV
export_predictions(predictions, "my_predictions.csv")

# Export to Excel
export_predictions(predictions, "my_predictions.xlsx")

# Explicit format
export_predictions(predictions, "output.txt", format = "csv")
} # }
```
