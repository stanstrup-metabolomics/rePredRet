# Search for a Compound Across All Systems

Searches for a compound by InChI and returns all systems where it has
been measured. Stereochemistry is automatically stripped for matching.

## Usage

``` r
search_compound(inchi, system_type = NULL, data_path = NULL)
```

## Arguments

- inchi:

  InChI string to search for. Stereochemistry layers will be
  automatically removed for matching.

- system_type:

  Optional filter by system type ("RP" or "HILIC").

- data_path:

  Optional path to local RepoRT data directory.

## Value

A data.frame with columns: system_id, system_name, system_type,
compound_name, rt, inchi.

## Details

The search automatically sanitizes InChI strings by removing
stereochemistry layers (/t, /b, /m, /s) before matching, ensuring
consistent results regardless of whether the query or database entries
include stereochemistry information.

## Examples

``` r
if (FALSE) { # \dontrun{
# Search by InChI (caffeine)
search_compound("InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3")

# Partial InChI also works
search_compound("InChI=1S/C8H10N4O2")

# Filter to RP systems only
search_compound("InChI=1S/C15H14O6", system_type = "RP")
} # }
```
