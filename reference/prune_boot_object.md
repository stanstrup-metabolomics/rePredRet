# Prune Bootstrap Object for Storage

Removes large components from a boot object to reduce storage size while
preserving essential information.

## Usage

``` r
prune_boot_object(boot_object)
```

## Arguments

- boot_object:

  A boot object from
  [`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html).

## Value

Pruned boot object with reduced size.

## Details

The function removes:

- `statistic`: References the calling environment (huge)

- `t`: All bootstrap iteration results (we only need final CI)

This typically reduces object size by ~95%.
