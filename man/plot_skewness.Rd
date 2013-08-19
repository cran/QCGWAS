\name{plot_skewness}
\alias{plot_skewness}
\title{GWAS Skewness vs. Kurtosis Plot}
\description{
  This function generates the skewness vs. kurtosis plot created
  by \code{\link{QC_series}}.
}
\usage{
plot_skewness(skewness,
              kurtosis,
              labels = paste("Study", 1:length(skewness)),
              plot_labels = "outliers",
              save_name = "Graph_skewness_kurtosis",
              save_dir = getwd(), ...)
}
\arguments{
  \item{skewness, kurtosis}{
    Vectors containing the skewness and kurtosis values of the
    datasets}
  \item{labels}{vector containing names or other identifiers for
    the datapoints, to be plotted in the graph. Note: it's best
    to keep these very short.}
  \item{plot_labels}{character string or logical determining
    whether the values in \code{labels} are plotted next to the
    data points. The possible settings are: \code{"none"} (or 
    \code{FALSE}); \code{"all"} (or \code{TRUE}); and
    \code{"outliers"} for outliers only.}
  \item{save_name}{character string; the filename, \emph{without}
    extension, for the graph file.}
  \item{save_dir}{character string; the directory where the
    graph is saved. Note that R uses \emph{forward}
    slash (/) where Windows uses backslash (\\).}
  \item{\dots}{arguments passed to \code{\link{plot}}.}
}
\details{
  When running a QC over multiple files, \code{\link{QC_series}}
  collects the values of the \code{skewness_HQ} and \code{kurtosis_HQ}
  output of \code{\link{QC_GWAS}} in a table, which is then
  passed to this function to convert it into a plot. Note that
  this values are calculated over high-quality SNPs only.
}
\value{An invisible \code{NULL}.}
\seealso{
  For calculating skewness and kurtosis: \code{\link{calc_kurtosis}}.
  
  For other plots comparing GWAS results files:
  \code{\link{plot_precision}} and \code{\link{plot_distribution}}.
}
\examples{
  value_S <- c(0.05, -0.27, 0.10, 0.11)
  value_K <- c( 6.7,  10.0, 10.1,  6.6)
  value_labels <- paste("cohort", 1:4)

  plot_skewness(skewness = value_S,
                kurtosis = value_K,
                labels = value_labels,
                plot_labels = "outliers",
                save_name = "sample_skewness_kurtosis")
}
\keyword{distribution}
