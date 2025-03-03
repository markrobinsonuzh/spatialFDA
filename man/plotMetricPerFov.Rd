% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMetric.R
\name{plotMetricPerFov}
\alias{plotMetricPerFov}
\title{Plot a spatial metric per field of view}
\usage{
plotMetricPerFov(
  metricDf,
  theo = FALSE,
  correction = NULL,
  x = NULL,
  imageId = NULL,
  ID = NULL
)
}
\arguments{
\item{metricDf}{the metric \code{dataframe} as calculated by \code{calcMetricPerFov}}

\item{theo}{logical; if the theoretical line should be plotted}

\item{correction}{the border correction to plot}

\item{x}{the x-axis variable to plot}

\item{imageId}{the ID of the image/fov}

\item{ID}{the (optional) ID for plotting combinations}
}
\value{
a \code{ggplot} object
}
\description{
A function that plots the output of the function \code{calcMetricPerFov}. The plot
contains one curve per FOV and makes subplots by samples.
}
\examples{
spe <- imcdatasets::Damond_2019_Pancreas("spe", full_dataset = FALSE)
metricRes <- calcMetricPerFov(spe, c("alpha", "beta"),
    subsetby = "image_number", fun = "Gcross", marks = "cell_type",
    rSeq = seq(0, 50, length.out = 50), by = c(
        "patient_stage", "patient_id",
        "image_number"
    ),
    ncores = 1
)
p <- plotMetricPerFov(metricRes,
    correction = "rs", x = "r",
    imageId = "image_number", ID = "ID"
)
print(p)
}
