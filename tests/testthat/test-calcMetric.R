# load an example
spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
rSeq <- seq(0, 50, length.out = 50)

## test function calcMetricPerFov
test_that("Output contains correction for discrete single mark", {
  rSeq <- seq(0, 50, length.out = 50)
  metricRes <- calcMetricPerFov(spe, "alpha",
                                subsetby = "image_number", fun = "Gest",
                                marks = "cell_type",
                                rSeq = rSeq, by = c(
                                  "patient_stage", "patient_id",
                                  "image_number"),
                                correction = "rs",
                                ncores = 1
  )
  expect_contains(colnames(metricRes), "rs")
  expect_contains(colnames(metricRes), "r")
  expect_contains(colnames(metricRes), "theo")
})

test_that("Output contains correction for continuous single mark", {
  # add continuous mark to colData
  protein <- "CD31"
  expr <- assay(spe, "exprs")[protein, ] %>%
    as.matrix() %>%
    data.frame() %>%
    rename("CD31" = ".")
  colData(spe) <- colData(spe) %>% cbind(expr)

  rSeq <- seq(0, 50, length.out = 50)
  metricRes <- calcMetricPerFov(spe, "alpha",
                                subsetby = "image_number", fun = "markcorr",
                                marks = protein,
                                rSeq = rSeq, by = c(
                                  "patient_stage", "patient_id",
                                  "image_number"),
                                correction = "iso",
                                ncores = 1,
                                continuous = TRUE
  )
  expect_contains(colnames(metricRes), "iso")
  expect_contains(colnames(metricRes), "r")
  expect_contains(colnames(metricRes), "theo")
})

test_that("Output contains correction for two marks", {
  metricRes <- calcMetricPerFov(spe, c("alpha", "beta"),
                                subsetby = "image_number", fun = "Gcross",
                                marks = "cell_type",
                                rSeq = rSeq, by = c(
                                  "patient_stage", "patient_id",
                                  "image_number"),
                                correction = "rs",
                                ncores = 1
  )
  expect_contains(colnames(metricRes), "rs")
  expect_contains(colnames(metricRes), "r")
  expect_contains(colnames(metricRes), "theo")
})

test_that("Output has correct dimensions", {
  metricRes <- calcMetricPerFov(spe, "alpha",
                                subsetby = "image_number", fun = "Gest",
                                marks = "cell_type",
                                rSeq = rSeq, by = c(
                                  "patient_stage", "patient_id",
                                  "image_number"),
                                correction = "rs",
                                ncores = 1
  )
  expect_length(metricRes$rs, length(rSeq) * length(unique(spe$image_name)))
})

test_that("Function fails if marks not in ColData", {
  expect_error(calcMetricPerFov(spe, c("alpha", "epsilon"),
                                subsetby = "image_number", fun = "Gcross",
                                marks = "cell_type",
                                rSeq = rSeq, by = c(
                                  "patient_stage", "patient_id",
                                  "image_number"
                                ),
                                ncores = 1
  ))
})

test_that("Function fails if fun not in spatstat.explore", {
  expect_error(calcMetricPerFov(spe, c("alpha", "beta"),
                                subsetby = "image_number", fun = "Mcross",
                                marks = "cell_type",
                                rSeq = rSeq, by = c(
                                  "patient_stage", "patient_id",
                                  "image_number"
                                ),
                                ncores = 1
  ))
})

# Test function calcCrossMetricPerFov
test_that("Cross function output has correct dimensions", {
  selection <- c("alpha", "beta", "delta")
  metricRes <- calcCrossMetricPerFov(spe, selection,
                                     subsetby = "image_number", fun = "Gcross",
                                     marks = "cell_type",
                                     rSeq = seq(0, 50, length.out = 50), by = c(
                                       "patient_stage", "patient_id",
                                       "image_number"
                                     ),
                                     ncores = 1
  )
  expect_length(metricRes$rs,
                length(rSeq) * length(unique(spe$image_name))
                * (length(selection)^2))
})

test_that("Cross function output has correct dimensions for Kdot", {
  selection <- c("alpha", "beta", "delta")
  metricRes <- calcCrossMetricPerFov(spe, selection,
                                     subsetby = "image_number", fun = "Kdot",
                                     marks = "cell_type",
                                     rSeq = seq(0, 50, length.out = 50), by = c(
                                       "patient_stage", "patient_id",
                                       "image_number"
                                     ),
                                     correction = "border",
                                     ncores = 1
  )
  expect_length(metricRes$border,
                length(rSeq) * length(unique(spe$image_name))
                * (length(selection)))
})
