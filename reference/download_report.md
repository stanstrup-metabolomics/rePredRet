# Download RepoRT Repository

Downloads the RepoRT repository from GitHub.

## Usage

``` r
download_report(dest_dir = tempdir(), release = NULL, overwrite = FALSE)
```

## Arguments

- dest_dir:

  Directory to save the downloaded data.

- release:

  Optional release tag. If NULL (default), downloads the master branch.

- overwrite:

  If TRUE, overwrite existing files.

## Value

Path to the extracted RepoRT directory.

## Examples

``` r
if (FALSE) { # \dontrun{
report_path <- download_report(dest_dir = tempdir())
} # }
```
