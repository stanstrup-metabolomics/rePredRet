# Export Predictions to Files

Saves predictions to CSV and JSON files with full provenance.

Exports predictions to CSV or Excel format.

## Usage

``` r
export_predictions(predictions, file, format = c("auto", "csv", "xlsx"))

export_predictions(predictions, file, format = c("auto", "csv", "xlsx"))
```

## Arguments

- predictions:

  A UserRTPredictions object or data.frame of predictions.

- file:

  Output file path. Extension determines format if format is "auto".

- format:

  Output format: "auto" (from extension), "csv", or "xlsx".

- output_dir:

  Directory to save predictions.

- studies:

  Studies tibble for system metadata.

## Value

Invisibly returns the output directory path.

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
