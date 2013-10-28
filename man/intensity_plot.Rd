\name{intensity_plot}
\alias{intensity_plot}
\title{Generates an intensity plot from x, y datasets}
\description{
  This function is used by \code{\link{match_alleles}} to
  generate an intensity plot. This function is currently only
  partially implemented, so we recommend that users do not
  bother with it.
}
\usage{
intensity_plot(x, y, strata, nbin = 20,
               xmax = max(x), xmin = min(x),
               ymax = max(y), ymin = min(y),
               strata_colours = c("black", "red", "turquoise3"),
               verbose = TRUE, xlab = "x", ylab = "y", ...)
}
\arguments{
  \item{x, y}{numerical vectors; the x and y coordiantes of the
    datapoints.}
  \item{strata}{logical vector; indicates whether the datapoint
    belongs to strata 1 or 2. If missing, all datapoints are
    assumed to belong to strata 1.}
  \item{nbin}{integer: the number of bins (categories) on the x
    and y axis.}
  \item{xmax, xmin, ymax, ymin}{numeric; the range of x and y
    values shown in the plot.}
  \item{strata_colours}{character vector of length \code{3},
    indicating the colours to be used for entries of strata 1, 2
    or mixed, respectively.}
  \item{verbose}{logical; determines whether a warning is printed
    in the console when datapoints are removed.}
  \item{xlab, ylab, \dots}{arguments passed to \code{\link{plot}}.}
}
\details{
  This function is intended as an alternative to the standard
  scatter plot for the allele-frequency correlation graph.
}
\value{
  An invisible object of type \code{list} with the following two
  components:
  \item{NA_removed}{The number of entries removed due to missing
    values}
  \item{outliers_removed}{The number of entries removed because
    they exceeded the thresholds specified by the min/max
    arguments.}
}
\seealso{\code{\link{match_alleles}}}
\keyword{distribution}
\keyword{univar}
