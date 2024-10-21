# load an example
spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
rSeq <- seq(0, 50, length.out = 50)

## test function calcMetricPerFov
test_that("Output contains correction for single mark", {
  rSeq <- seq(0, 50, length.out = 50)
  metricRes <- calcMetricPerFov(spe, "alpha",
                                subsetby = "image_number", fun = "Gest",
                                marks = "cell_type",
                                rSeq = rSeq, by = c(
                                  "patient_stage", "patient_id",
                                  "image_number"),
                                correction = 'rs',
                                ncores = 1
  )
  expect_contains(colnames(metricRes), 'rs')
})

test_that("Output contains correction for two marks", {
  metricRes <- calcMetricPerFov(spe, c("alpha", "beta"),
                                subsetby = "image_number", fun = "Gcross",
                                marks = "cell_type",
                                rSeq = rSeq, by = c(
                                  "patient_stage", "patient_id",
                                  "image_number"),
                                correction = 'rs',
                                ncores = 1
  )
  expect_contains(colnames(metricRes), 'rs')
})

test_that("Output has correct dimensions", {
  metricRes <- calcMetricPerFov(spe, "alpha",
                                subsetby = "image_number", fun = "Gest",
                                marks = "cell_type",
                                rSeq = rSeq, by = c(
                                  "patient_stage", "patient_id",
                                  "image_number"),
                                correction = 'rs',
                                ncores = 1
  )
  expect_length(metricRes$rs, length(rSeq)*length(unique(spe$image_name)))
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

# Test function calcCrossMetricPerFov
test_that("Cross function output has correct dimensions", {
  metricRes <- calcCrossMetricPerFov(spe, c("alpha", "beta", "delta"),
                                     subsetby = "image_number", fun = "Gcross",
                                     marks = "cell_type",
                                     rSeq = seq(0, 50, length.out = 50), by = c(
                                       "patient_stage", "patient_id",
                                       "image_number"
                                     ),
                                     ncores = 1
  )
  expect_length(metricRes$rs,
                length(rSeq)*length(unique(spe$image_name))
                *(length(selection)^2))
})

test_that("Cross function output has correct dimensions for Kdot", {
  metricRes <- calcCrossMetricPerFov(spe, c("alpha", "beta", "delta"),
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
                length(rSeq)*length(unique(spe$image_name))
                * (length(selection)))
})
