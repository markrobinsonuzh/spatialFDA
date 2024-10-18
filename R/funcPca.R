#' Functional Principal Component Analysis
#'
#' A function that takes as input the output of `calcMetricPerFov` which has to
#' be converted into the correct format by `prepData`. The output is a list with
#' the `fpca.face` output from refund.
#'
#' @param dat a data object for functional data analysis containing at least the
#' the functional
#' @param r the functional domain
#' @param knots the number of knots
#' @param pve the proportion of variance explained
#'
#' @return a list with components of fpca.face
#' @export
#'
#' @examples
#' # load the pancreas dataset
#' library("tidyr")
#' library("stringr")
#' library("dplyr")
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' # calculate the Gcross metric for alpha and beta cells
#' metricRes <- calcMetricPerFov(spe, c("alpha", "beta"),
#'     subsetby = "image_number", fun = "Gcross",
#'     marks = "cell_type", rSeq = seq(0, 50, length.out = 50),
#'     c("patient_stage", "patient_id", "image_number"), ncores = 1
#' )
#' metricRes$ID <- paste0(
#'     metricRes$patient_stage, "x", metricRes$patient_id,
#'     "x", metricRes$image_number
#' )
#' # prepare data for FDA
#' dat <- prepData(metricRes, "r", "rs")
#'
#' # drop rows with NA
#' dat <- dat |> drop_na()
#' # create meta info of the IDs
#' splitData <- str_split(dat$ID, "x")
#' dat$condition <- factor(sapply(splitData, function(x) x[1]))
#' dat$patient_id <- factor(sapply(splitData, function(x) x[2]))
#' dat$image_id <- factor(sapply(splitData, function(x) x[3]))
#' # calculate fPCA
#' mdl <- functionalPCA(
#'     dat = dat, r = metricRes$r |> unique(),
#'     knots = 30, pve = 0.99
#' )
#' @import dplyr
functionalPCA <- function(dat, r, knots, pve = 0.95) {
    stopifnot(is(dat, 'data.frame'))
    stopifnot(is(r, 'vector'))
    stopifnot(is(knots, 'numeric'))
    stopifnot(is(pve, 'numeric'))
    # calculate the fPCA - this is a bit a pointless wrapper until now
    res <- refund::fpca.face(
        Y = dat$Y, center = TRUE, argvals = r,
        knots = knots, pve = pve
    )
    return(res)
}

#' Plot a biplot from an fPCA analysis
#'
#' A function that takes the output from the `functionalPCA` function and returns
#' a `ggplot` object of the first two dimensions of the PCA as biplot.
#'
#' @param dat a data object for functional data analysis containing at least the
#' the functional
#' @param res the output from the fPCA calculation
#' @param colourby the variable by which to colour the PCA plot by
#' @param labelby the variable by which to label the PCA plot by
#'
#' @return a list with components of fpca.face
#' @export
#'
#' @examples
#' # load the pancreas dataset
#' library("tidyr")
#' library("stringr")
#' library("dplyr")
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' # calculate the Gcross metric for alpha and beta cells
#' metricRes <- calcMetricPerFov(spe, c("alpha", "beta"),
#'     subsetby = "image_number", fun = "Gcross",
#'     marks = "cell_type", rSeq = seq(0, 50, length.out = 50),
#'     c("patient_stage", "patient_id", "image_number"), ncores = 1
#' )
#' metricRes$ID <- paste0(
#'     metricRes$patient_stage, "x", metricRes$patient_id,
#'     "x", metricRes$image_number
#' )
#' # prepare data for FDA
#' dat <- prepData(metricRes, "r", "rs")
#'
#' # drop rows with NA
#' dat <- dat |> drop_na()
#' # create meta info of the IDs
#' splitData <- str_split(dat$ID, "x")
#' dat$condition <- factor(sapply(splitData, function(x) x[1]))
#' dat$patient_id <- factor(sapply(splitData, function(x) x[2]))
#' dat$image_id <- factor(sapply(splitData, function(x) x[3]))
#' # calculate fPCA
#' mdl <- functionalPCA(
#'     dat = dat, r = metricRes$r |> unique(),
#'     knots = 30, pve = 0.99
#' )
#' p <- plotFpca(
#'     dat = dat, res = mdl, colourby = "condition",
#'     labelby = "patient_id"
#' )
#' print(p)
#' @import dplyr
plotFpca <- function(dat, res, colourby = NULL, labelby = NULL) {
    stopifnot(is(dat, 'data.frame'))
    stopifnot(is(res, 'fpca'))
    scoresDf <- res$scores %>% as.data.frame()
    # plot fCPA results - assumes same order of fPCA results and input data
    p <- ggplot(scoresDf, aes(scoresDf[, 1], scoresDf[, 2],
        colour = factor(dat[[colourby]])
    )) +
        geom_point() +
        coord_equal() +
        theme_light() +
        xlab("functional PC1") +
        ylab("functional PC2")
    if (!is.null(labelby)) {
        p <- p +
            geom_text(hjust = 0, vjust = 0, aes(label = factor(dat[[labelby]])))
    }
    return(p)
}
