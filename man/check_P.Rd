\name{check_P}
\alias{check_P}
\title{Checking GWAS p-values}
\description{
  A simple test to check if the reported p-values in a GWAS
  results file match the other statistics. This function
  calculates an expected p-value (from the effect size and
  standard error) and then correlates it with the actual,
  reported p-value.
}
\usage{
check_P(dataset, HQ_subset,
        plot_correlation = FALSE, plot_if_threshold = FALSE,
        threshold_r = 0.99,
        save_name = "dataset", save_dir = getwd(),
        header_translations,
        use_log = FALSE, dataN = nrow(dataset), ...)
}
\arguments{
  \item{dataset}{table with at least three
    columns: p-value, effect size and standard error.}
  \item{HQ_subset}{an \emph{optional} logical or numeric vector
    indicating the rows in \code{dataset} that contain high
    quality SNPs.}
  \item{plot_correlation}{logical; should a scatterplot of
    the reported vs. calculated p-values be made? If \code{TRUE},
    the plot is saved as a .png file.}
  \item{plot_if_threshold}{logical; if \code{TRUE}, the
    scatterplot is only saved when the correlation between
    reported and calculated p-values is lower than
    \code{threshold_r}.}
  \item{threshold_r}{numeric; the correlation threshold for the
    scatterplot.}
  \item{save_name}{character string; the filename, \emph{without}
    extension, for the scatterplot.}
  \item{save_dir}{character string; the directory where the
    output files are saved. Note that R uses \emph{forward}
    slash (/) where Windows uses backslash (\\).}
  \item{header_translations}{translation table for column names
    See \code{\link{translate_header}} for more information. If
    the argument is left empty, \code{dataset} is assumed to use
    the standard column names used by \code{\link{QC_GWAS}}.}
  \item{use_log, dataN}{arguments used by \code{\link{QC_GWAS}};
    redundant when \code{check_P} is used separately.}
  \item{\dots}{arguments passed to \code{\link{plot}}.}
}
\details{
  \code{check_P} calculates the expected p-value by taking the
  chi-square (1 degree of freedom) of the effect size divided by
  the standard error squared.
  
  In a typical GWAS dataset, the expected and observed p-values
  should correlate perfectly. If this isn't the case, the
  problem either lies in a misidentified column, or the wrong
  values were used when generating the dataset.
}
\value{
  The correlation between expected and reported p-values.
}
\examples{
  data("gwa_sample")

  selected_SNPs <- HQ_filter(data = gwa_sample,
                             FRQ_val = 0.05,
                             cal_val = 0.95,
                             filter_NA = FALSE)
  # To calculate a correlation between predicted and actual p-values:
  check_P(gwa_sample, HQ_subset = selected_SNPs,
          plot_correlation = FALSE)
  
  # To plot the correlation:
  \dontrun{
    check_P(gwa_sample, HQ_subset = selected_SNPs,
            plot_correlation = TRUE, plot_if_threshold = FALSE,
            save_name = "sample")
  }
}
\keyword{univar}
