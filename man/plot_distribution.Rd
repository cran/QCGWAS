\name{plot_distribution}
\alias{plot_distribution}
\title{GWAS effect-Size distribution plot}
\description{
  This function generates the effect-size distribution boxplot
  created by \code{\link{QC_series}}.
}
\usage{
plot_distribution(data_table,
                  names = 1:ncol(data_table),
                  include = TRUE,
                  plot_order = 1:ncol(data_table),
                  quantile_lines = FALSE,
                  save_name = "Graph_distribution",
                  save_dir = getwd(), ...)
}
\arguments{
  \item{data_table}{table with a column of effect
   sizes for every dataset.}
  \item{names}{vector; the names for the datasets, for use in the
    graph. Note: it's best to keep these very short, as long
    labels won't be plotted. The default is
    the column \emph{numbers} of \code{data_table}.}
  \item{include}{logical vector indicating which columns of
    \code{data_table} are included in the plot. The default
    setting is to include all.}
  \item{plot_order}{numeric vector determining the left-to-right
    order of plotting the datasets (columns). \code{\link{QC_series}}
    uses the sample size for this.}
  \item{quantile_lines}{logical; should lines representing the
    median and quartile values be included?}
  \item{save_name}{character string; the filename, \emph{without}
    extension, for the graph file.}
  \item{save_dir}{character string; the directory where the
    graph is saved. Note that R uses \emph{forward}
    slash (/) where Windows uses backslash (\\).}
  \item{\dots}{arguments passed to \code{\link{boxplot}}.}
}
\details{
  When running a QC over multiple files, \code{\link{QC_series}}
  collects the values of the \code{effectsize_HQ} output of
  \code{\link{QC_GWAS}} in a table, which is then passed to this
  function. If there are significant
  differences in the distribution of effect sizes, it usually
  indicates that the datasets did not use the same model or unit.
}
\value{An invisible \code{NULL}.}
\seealso{
  For comparing reported to expected effect-size distribution:
  \code{\link{QC_histogram}}.
  
  For other plots comparing GWAS results files:
  \code{\link{plot_precision}} and \code{\link{plot_skewness}}.
}
\examples{
  data("gwa_sample")

  chunk1 <- gwa_sample$EFFECT[1:1000]
  chunk2 <- gwa_sample$EFFECT[1001:2000]
  chunk3 <- gwa_sample$EFFECT[2001:3000]

  plot_distribution(
    data_table = data.frame(chunk1, chunk2, chunk3),
    names = c("chunk 1", "chunk 2", "chunk 3"),
    quantile_lines = TRUE,
    save_name = "sample_distribution")
}
\keyword{distribution}
