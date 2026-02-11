# Export Statistics to Files

Saves computed statistics to CSV and JSON files for the website.

## Usage

``` r
export_statistics(stats_list, output_dir = "statistics", metadata = NULL)
```

## Arguments

- stats_list:

  A list containing model_stats, database_stats, and prediction_stats
  from the calculate\_\* functions.

- output_dir:

  Directory to save statistics.

- metadata:

  Additional metadata (e.g., RepoRT version, run date).

## Value

Invisibly returns the output directory path.
