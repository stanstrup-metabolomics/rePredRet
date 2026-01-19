# Tests for data loading functions

test_that("get_common_compounds handles empty data", {
  # Create mock datasets with no overlap
  dataset1 <- list(
    rtdata = data.frame(
      inchi.std = c("InChI=1S/A", "InChI=1S/B"),
      rt = c(1.0, 2.0)
    )
  )

  dataset2 <- list(
    rtdata = data.frame(
      inchi.std = c("InChI=1S/C", "InChI=1S/D"),
      rt = c(3.0, 4.0)
    )
  )

  result <- get_common_compounds(dataset1, dataset2)

  expect_equal(nrow(result), 0)
  expect_equal(ncol(result), 2)
})


test_that("get_common_compounds finds matching compounds", {
  # Create mock datasets with overlap
  dataset1 <- list(
    rtdata = data.frame(
      inchi.std = c("InChI=1S/A", "InChI=1S/B", "InChI=1S/C"),
      rt = c(1.0, 2.0, 3.0)
    )
  )

  dataset2 <- list(
    rtdata = data.frame(
      inchi.std = c("InChI=1S/B", "InChI=1S/C", "InChI=1S/D"),
      rt = c(5.0, 6.0, 7.0)
    )
  )

  result <- get_common_compounds(dataset1, dataset2)

  expect_equal(nrow(result), 2)  # B and C are common
  expect_equal(ncol(result), 2)
  expect_equal(colnames(result), c("rt_source", "rt_target"))
})


test_that("get_common_compounds handles duplicates by taking median", {
  # Create mock datasets with duplicates
  dataset1 <- list(
    rtdata = data.frame(
      inchi.std = c("InChI=1S/A", "InChI=1S/A", "InChI=1S/B"),
      rt = c(1.0, 1.2, 2.0)  # Two measurements of A
    )
  )

  dataset2 <- list(
    rtdata = data.frame(
      inchi.std = c("InChI=1S/A", "InChI=1S/B"),
      rt = c(5.0, 6.0)
    )
  )

  result <- get_common_compounds(dataset1, dataset2)

  # Should have 2 compounds, with median RT for A
  expect_equal(nrow(result), 2)

  # Find row for compound A (should have median RT = 1.1)
  a_row <- which(attr(result, "compounds") == "InChI=1S/A")
  expect_equal(as.numeric(result[a_row, "rt_source"]), 1.1, tolerance = 0.01)
})


test_that("get_common_compounds strips stereochemistry from InChI", {
  # InChI with stereochemistry layers
  dataset1 <- list(
    rtdata = data.frame(
      inchi.std = c("InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1"),
      rt = c(1.0)
    )
  )

  # Same compound without stereochemistry
  dataset2 <- list(
    rtdata = data.frame(
      inchi.std = c("InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2"),
      rt = c(2.0)
    )
  )

  result <- get_common_compounds(dataset1, dataset2)

  # Should match after stripping stereochemistry
  expect_equal(nrow(result), 1)
})
