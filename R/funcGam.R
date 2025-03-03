#' General additive model with functional response
#'
#' A function that takes the output of a metric calculation as done by
#' `calcMetricPerFov`. The data has to be prepared into the correct format for
#' the functional analysis by the `prepData` function. The output is a `pffr`
#' object as implemented by `refund`.
#'
#' @param data a dataframe with the following columns: Y = functional response;
#' sample_id = sample ID; image_id = image ID;
#' @param x the x-axis values of the functional response
#' @param designmat a design matrix as defined by model.matrix()
#' @param weights weights as the number of points per image. These weights are
#' normalised by the mean of the weights in the fitting process
#' @param formula the formula for the model. The colnames of the designmatrix
#' have to correspond to the variables in the formula.
#' @param family the distributional family as implemented in `family.mgcv`. For
#' fast computation the default is set to `gaussian`. If the covariance scales
#' e.g. as a function of the domain, this estimation can be improved with
#' `gaulss` - for more information see `family.mgcv`.
#' @param ... Other parameters passed to `pffr`
#'
#' @return a fitted pffr object which inherits from gam
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
#'
#' # create meta info of the IDs
#' splitData <- str_split(dat$ID, "x")
#' dat$condition <- factor(sapply(splitData, function(x) x[1]))
#' dat$patient_id <- factor(sapply(splitData, function(x) x[2]))
#' dat$image_id <- factor(sapply(splitData, function(x) x[3]))
#' # create a designmatrix
#' condition <- dat$condition
#' # relevel the condition - can set explicit contrasts here
#' condition <- relevel(condition, "Non-diabetic")
#' designmat <- model.matrix(~condition)
#' # colnames don't work with the '-' sign
#' colnames(designmat) <- c(
#'     "Intercept", "conditionLong_duration",
#'     "conditionOnset"
#' )
#' # fit the model
#' mdl <- functionalGam(
#'     data = dat, x = metricRes$r |> unique(),
#'     designmat = designmat, weights = dat$npoints,
#'     formula = formula(Y ~ conditionLong_duration +
#'         conditionOnset + s(patient_id, bs = "re"))
#' )
#' summary(mdl)

#' @import dplyr
#' @importFrom methods is
#' @importFrom stats terms
functionalGam <- function(data, x, designmat, weights, formula,
                          family = "gaussian", ...) {
    # type checking
    stopifnot(is(data, "data.frame"))
    stopifnot(is(x, "vector"))
    stopifnot(is(designmat, "matrix"))
    stopifnot(is(weights, "integer"))
    stopifnot(is(formula, "formula"))
    stopifnot(is(family, "character"))

    data <- cbind(data, designmat)
    # Test that length of number of points and nrow of designmat correspond
    stopifnot(length(weights) == nrow(designmat))
    # normalise the weights
    weights <- weights / mean(weights)
    # Test that the colnames of the designmat correspond to formula
    # with the exception of the intercept as this is inferred - there are
    # functional and constant intercepts
    # deselect "Intercept"
    colNames <- colnames(designmat)[!colnames(designmat) == "Intercept"]
    # stop if the rest terms of design mat are not in the formula
    stopifnot(colNames %in% attr(terms(formula), "term.labels"))
    mdl <- refund::pffr(formula,
        yind = x,
        data = data,
        weights = weights,
        family = family,
        ...
    )
    return(mdl)
}
