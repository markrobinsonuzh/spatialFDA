library("tidyr")
library("stringr")
library("dplyr")
spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
# calculate the Gcross metric for alpha and beta cells
metricRes <- calcMetricPerFov(spe, c("alpha", "beta"),
                              subsetby = "image_number", fun = "Gcross",
                              marks = "cell_type", rSeq = seq(0, 50, length.out = 50),
                              c("patient_stage", "patient_id", "image_number"), ncores = 1
)
metricRes$ID <- paste0(
  metricRes$patient_stage, "x", metricRes$patient_id,
  "x", metricRes$image_number
)
# prepare data for FDA
dat <- prepData(metricRes, "r", "rs")

# drop rows with NA
dat <- dat |> drop_na()
# create meta info of the IDs
splitData <- str_split(dat$ID, "x")
dat$condition <- factor(sapply(splitData, function(x) x[1]))
dat$patient_id <- factor(sapply(splitData, function(x) x[2]))
dat$image_id <- factor(sapply(splitData, function(x) x[3]))
# calculate fPCA
mdl <- functionalPCA(
  data = dat, r = metricRes$r |> unique(),
  knots = 30, pve = 0.99
)

test_that("Correct output type", {
  # of the correct type
  expect_equal(is(mdl), "fpca")
  # contains the required output
  expect_true(!is.null(mdl$Yhat))
  expect_true(!is.null(mdl$scores))
  expect_true(!is.null(mdl$efunctions))
  expect_true(!is.null(mdl$evalues))
  expect_true(!is.null(mdl$npc))
  expect_true(!is.null(mdl$pve))
})

test_that("Can handle NA in response", {
  dat[1, 2][9] <- NA
  mdl <- functionalPCA(
    data = dat, r = metricRes$r |> unique(),
    knots = 30, pve = 0.99
  )
  expect_equal(is(mdl), "fpca")
  # contains the required output
  expect_true(!is.null(mdl$Yhat))
  expect_true(!is.null(mdl$scores))
  expect_true(!is.null(mdl$efunctions))
  expect_true(!is.null(mdl$evalues))
  expect_true(!is.null(mdl$npc))
  expect_true(!is.null(mdl$pve))
})
