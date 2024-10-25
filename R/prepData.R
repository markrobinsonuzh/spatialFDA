#' Prepare data from calcMetricRes to be in the right format for FDA
#'
#' @param metricRes a dataframe as calculated by calcMetricRes - requires
#' the columns ID (unique identifier of each row)
#' @param x the name of the x-axis of the spatial metric
#' @param y the name of the y-axis of the spatial metric
#'
#' @return returns a list with three entries, the unique ID, the functional
#' response Y and the weights
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' metricRes <- calcMetricPerFov(spe, c("alpha", "beta"),
#'     subsetby = "image_number", fun = "Gcross", marks = "cell_type",
#'     rSeq = seq(0, 50, length.out = 50), by = c(
#'         "patient_stage", "patient_id",
#'         "image_number"
#'     ),
#'     ncores = 1
#' )
#'
#' # create a unique ID for each row
#' metricRes$ID <- paste0(
#'     metricRes$patient_stage, "x", metricRes$patient_id,
#'     "x", metricRes$image_number
#' )
#' dat <- prepData(metricRes, "r", "rs")
#' @import tidyr
#' @importFrom methods is
prepData <- function(metricRes, x, y) {
    # type checking
    stopifnot(is(metricRes, "data.frame"))
    stopifnot(is(x, "character"))
    stopifnot(is(y, "character"))
    # extract the functional response matrix
    mat <- metricRes %>%
        select("ID", x, y) %>%
        spread("ID", y) %>%
        select(!x)
    # extract the number of points as weights - are ordered differently,
    # thus order according to image ID
    weights <- metricRes %>%
        group_by("ID") %>%
        select("ID", "npoints") %>%
        unique() %>%
        arrange("ID")
    # create a dataframe as required by pffr
    dat <- data.frame(ID = colnames(mat))
    dat$Y <- t(mat)
    dat$weights <- weights
    return(dat)
}
