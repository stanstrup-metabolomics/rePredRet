# Read Model Index

Reads model_index.csv in the current (sys1_id/sys2_id) or the older
(from_id/to_id) format.

## Usage

``` r
read_model_index(index_file)
```

## Arguments

- index_file:

  Path to model_index.csv, or NULL

## Value

data.frame with model_key and the index columns
