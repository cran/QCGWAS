\name{plot_precision}
\alias{plot_precision}
\title{GWAS Precision Plot}
\description{
  This function generates the precision plot created by
  \code{\link{QC_series}}. Precision is defined as:
    
  \code{1 / median standard-error}.
}
\usage{
plot_precision(SE, N,
               labels = NULL,
               save_name = "Graph_precision",
               save_dir = getwd(), ...)
}
\arguments{
  \item{SE, N}{numeric vectors containing the median standard-error
    and sample size of the datasets, respectively.}
  \item{labels}{vector containing names or other identifiers for
    the datapoints, to be plotted in the graph. Note: it's best
    to keep these very short. To disable labeling, set to
    \code{NULL} (default).}
  \item{save_name}{character string; the filename, \emph{without}
    extension, for the graph file.}
  \item{save_dir}{character string; the directory where the
    graph is saved. Note that R uses \emph{forward}
    slash (/) where Windows uses backslash (\\).}
  \item{\dots}{arguments passed to \code{\link{plot}}.}
}
\details{
  When running a QC over multiple files, \code{\link{QC_series}}
  collects the values of the \code{SE_HQ_median} and \code{sample_size_HQ}
  outputs of \code{\link{QC_GWAS}} in a table, which is then
  passed to this function to convert it into a plot.
}
\value{An invisible \code{NULL}.}
\seealso{
  \code{\link{plot_distribution}} and \code{\link{plot_skewness}}.
}
\examples{
  value_SE <- c(0.078, 0.189, 0.077, 0.040, 0.021, 0.072)
  value_N <- c(870, 830, 970, 690, 2200, 870)
  value_labels <- paste("cohort", 1:6)
  
  plot_precision(SE = value_SE, N = value_N,
                 labels = value_labels,
                 save_name = "sample_precision")
}
\keyword{distribution}
