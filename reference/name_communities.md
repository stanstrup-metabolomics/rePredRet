# Name Groups of Systems

Gives each group (community) of systems the letter and name of the named
reference group it matches. The reference groups
(`inst/extdata/system_groups.csv` and `system_group_names.csv`) were
named by hand from the studies that contributed their systems and the
classes of their compounds (RepoRT data of October 2026, see the Network
article). Since the groups are recomputed from the current data, a group
gets a reference name only when they match by mutual majority: at least
half of the group's systems are in the reference group and at least half
of the reference group's current systems are in the group. Each name is
used at most once.

## Usage

``` r
name_communities(
  ids,
  community,
  reference = utils::read.csv(system.file("extdata", "system_groups.csv", package =
    "rePredRet"), colClasses = "character"),
  reference_names = utils::read.csv(system.file("extdata", "system_group_names.csv",
    package = "rePredRet"), colClasses = "character")
)
```

## Arguments

- ids:

  System ids.

- community:

  Group of each system (factor).

- reference:

  data.frame with sys_id and group (letter).

- reference_names:

  data.frame with group and name.

## Value

data.frame with one row per community: community, group (letter or NA),
name (or NA) and label (letter, or "unnamed" with a number).
