#' Compute a spatial metric on a SpatialExperiment object
#'
#' A function that takes a `SpatialExperiment` object and computes a spatial
#' statistics function as implemented in `spatstat`. The output is a `spatstat`
#' object.
#'
#' @param df A `dataframe` with the x and y coordinates from the corresponding
#' `SpatialExperiment` and the `colData`
#' @param selection the mark(s) you want to compare
#' @param fun the `spatstat` function to compute on the point pattern object
#' @param marks the marks to consider e.g. cell types
#' @param rSeq the range of r values to compute the function over
#' @param by the spe `colData` variable(s) to add to the meta data
#' @param continuous A boolean indicating whether the marks are continuous
#' defaults to FALSE
#' @param window a observation window for the point pattern of class `owin`.
#' @param ... Other parameters passed to `spatstat.explore` functions
#'
#' @return a `spatstat` metric object with the fov number, the number of
#' points and the centroid of the image
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' speSub <- subset(spe, , image_number == "138")
#' dfSub <- .speToDf(speSub)
#' metricRes <- .extractMetric(dfSub, c("alpha", "beta"),
#'     fun = "Gcross",
#'     marks = "cell_type", rSeq = seq(0, 1000, length.out = 100),
#'     by = c("patient_stage", "patient_id", "image_number")
#' )
#' @import spatstat.explore
.extractMetric <- function(df,
    selection,
    fun,
    marks = NULL,
    rSeq = NULL,
    by = NULL,
    continuous = FALSE,
    window = NULL,
    ...) {
    stopifnot(is(df, 'data.frame'))
    pp <- .dfToppp(df, marks = marks, continuous = continuous, window = window)
    if (!continuous) {
        ppSub <- subset(pp, marks %in% selection, drop = TRUE)
        metaData <- df[, by] %>% unique()
    } else {
        ppSub <- pp
        metaData <- df[, by] %>% unique()
        metaData$gene <- names(df)[names(df) %in% marks]
    }
    # small quality control to only consider pp that have more than 2 points per
    # fov and more than one unique mark and that each mark has more than one point
    if (spatstat.geom::npoints(ppSub) > 2 &&
        ((length(unique(
            spatstat.geom::marks(ppSub)
        )) > 1 &&
            sum(table(ppSub$marks) > 0) > 1) ||
            length(selection) == 1)) {
        metricRes <- tryCatch(
            {
                metricRes <- do.call(fun,
                    args = list(X = ppSub, r = rSeq, ...)
                )
            },
            warning = function(w) {
                print(w)
                metricRes <- do.call(fun,
                    args = list(X = ppSub, r = rSeq, ...)
                )
            },
            error = function(e) {
                print(e)
                metricRes <- data.frame(
                    r = rSeq,
                    fun = fun,
                    row.names = seq_along(rSeq)
                )
            }
        )
    }
    # This handles the case when we do cross functions for the same type
    else if (spatstat.geom::npoints(ppSub) > 2 &&
        length(unique(selection)) == 1 &&
        length(selection) > 1) {
        metricRes <- tryCatch(
            {
                # here we use pp, otherwise there are problems with the
                # mark connection function
                metricRes <- do.call(fun, args = list(
                    X = pp,
                    i = selection[1],
                    j = selection[2],
                    r = rSeq,
                    ...
                ))
            },
            warning = function(w) {
                print(w)
                metricRes <- do.call(fun, args = list(
                    X = pp,
                    i = selection[1],
                    j = selection[2],
                    r = rSeq,
                    ...
                ))
            },
            error = function(e) {
                print(e)
                metricRes <- data.frame(
                    r = rSeq,
                    fun = fun,
                    row.names = seq_along(rSeq)
                )
            }
        )
    } else {
        metricRes <- data.frame(
            r = rSeq,
            fun = fun,
            row.names = seq_along(rSeq)
        )
    }
    metricRes <- cbind(metricRes, metaData)
    metricRes$npoints <- spatstat.geom::npoints(ppSub)
    centroid <- spatstat.geom::centroid.owin(ppSub$window)
    metricRes$centroidx <- centroid$x
    metricRes$centroidy <- centroid$y
    return(metricRes)
}

#' Calculate a spatial metric on a `SpatialExperiment` object per field of view
#'
#' A function that takes a `SpatialExperiment` object as input and calculates a
#' spatial metric as implemented by `spatstat` per field of view.
#'
#' @param spe a `SpatialExperiment` object
#' @param selection the mark(s) you want to compare
#' @param subsetby the spe `colData` variable to subset the data by
#' @param fun the `spatstat` function to compute on the point pattern object
#' @param marks the marks to consider e.g. cell types
#' @param rSeq the range of r values to compute the function over
#' @param by the spe `colData` variable(s) to add to the meta data
#' @param continuous A boolean indicating whether the marks are continuous
#' defaults to FALSE
#' @param ncores the number of cores to use for parallel processing, default = 1
#' @param ... Other parameters passed to `spatstat.explore` functions
#'
#' @return a `dataframe` of the `spatstat` metric objects with the radius r, the
#' theoretical value of a Poisson process, the different border corrections
#' the fov number, the number of points and the centroid of the image
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
#' @import dplyr parallel
calcMetricPerFov <- function(spe, selection, subsetby = NULL, fun, marks = NULL,
    rSeq = NULL, by = NULL, continuous = FALSE, ncores = 1, ...) {
    # type checking of input
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot(is(fun, "character"))
    stopifnot(is(marks, "character"))
    stopifnot(is(ncores, "numeric"))

    df <- .speToDf(spe)
    # we have one case for discrete cell types where we have one column to subset
    if (length(subsetby) == 1) {
        dfLs <- split(df, df[[subsetby]])
    } else {
        dfLs <- purrr::map(subsetby, ~ df %>%
            select(all_of(setdiff(names(df), subsetby)), .x))
    }
    metricDf <- parallel::mclapply(dfLs, function(dfSub) {
        metricRes <- .extractMetric(
            df = dfSub,
            selection = selection,
            fun = fun,
            marks = marks,
            rSeq = rSeq,
            by = by,
            continuous = continuous,
            ...
        ) %>% as.data.frame()
        return(metricRes)
    }, mc.cores = ncores) %>% bind_rows()
    # store metadata of the calculation in the dataframe
    metricDf$ID <- paste0(metricDf[[by[1]]], "|", metricDf[[by[2]]])
    metricDf$fun <- fun
    metricDf$selection <- paste(selection, collapse = " and ")
    return(metricDf)
}


#' Calculate cross spatial metrics for all combinations per FOV
#'
#' A function that takes a `SpatialExperiment` object as input and calculates a
#' cross spatial metric as implemented by `spatstat` per field of view for all
#' combinations provided by the user.
#'
#' @param spe a `SpatialExperiment` object
#' @param selection the mark(s) you want to compare
#' @param subsetby the spe `colData` variable to subset the data by
#' @param fun the `spatstat` function to compute on the point pattern object
#' @param marks the marks to consider e.g. cell types
#' @param rSeq the range of r values to compute the function over
#' @param by the spe `colData` variable(s) to add to the meta data
#' @param ncores the number of cores to use for parallel processing, default = 1
#' @param continuous A boolean indicating whether the marks are continuous
#' defaults to FALSE
#' @param ... Other parameters passed to spatstat.explore functions
#'
#' @return a dataframe of the `spatstat` metric objects with the radius r, the
#' theoretical value of a Poisson process, the different border corrections
#' the fov number, the number of points and the centroid of the image
#' @export
#'
#' @examples
#' spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
#' metricRes <- calcCrossMetricPerFov(spe, c("alpha", "beta", "delta"),
#'     subsetby = "image_number", fun = "Gcross", marks = "cell_type",
#'     rSeq = seq(0, 50, length.out = 50), by = c(
#'         "patient_stage", "patient_id",
#'         "image_number"
#'     ),
#'     ncores = 1
#' )
calcCrossMetricPerFov <- function(
        spe, selection, subsetby = NULL, fun,
        marks = NULL, rSeq = NULL, by = NULL,
        ncores = 1, continuous = FALSE, ...) {
    # type checking of input
    stopifnot(is(spe, "SpatialExperiment"))
    stopifnot(is(fun, "character"))
    stopifnot(is(marks, "character"))
    stopifnot(is(ncores, "numeric"))

    # Special case of dot functions
    if (grepl("dot", fun)) {
        # one vs all other
        ls <- unique(selection)
        # calculate the metric per FOV
        resLs <- lapply(ls, function(x) {
            print(x)
            calcMetricPerFov(
                spe = spe, selection = x, subsetby = subsetby, fun = fun,
                marks = marks, rSeq = rSeq, by = by, ncores = ncores,
                continuous = continuous, ...
            )
        })
        # Bind the data and return
        return(bind_rows(resLs))
    } else {
        # This creates a grid with all possible 2 way combinations
        ls <- apply(expand.grid(selection, selection), 1, function(x) {
            return(c(x[1], x[2]))
        }) |> t()

        # calculate the metric per FOV
        resLs <- apply(ls, 1, function(x) {
            calcMetricPerFov(
                spe = spe, selection = x, subsetby = subsetby, fun = fun,
                marks = marks, rSeq = rSeq, by = by, ncores = ncores,
                continuous = continuous
            )
        })

        # Bind the data and return
        return(bind_rows(resLs))
    }
}
