# Batch Search for Multiple Compounds

Searches for multiple compounds by InChI and returns all systems where
they have been measured.

## Usage

``` r
search_compounds(inchis, system_type = NULL, data_path = NULL)
```

## Arguments

- inchis:

  Character vector of InChI strings to search for.

- system_type:

  Optional filter by system type ("RP" or "HILIC").

- data_path:

  Optional path to local RepoRT data directory.

## Value

A data.frame with columns: query_inchi, system_id, system_name,
system_type, compound_name, rt, inchi.

## See also

[`search_compound`](https://stanstrup.github.io/rePredRet/reference/search_compound.md)
for single compound search

## Examples

``` r
if (FALSE) { # \dontrun{
# Search for multiple compounds
results <- search_compounds(c(
  "InChI=1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3",
  "InChI=1S/C15H14O6/c16-8-4-11(18)9-6-13(20)15(21-14(9)5-8)7-1-2-10(17)12(19)3-7"
))
} # }
```
