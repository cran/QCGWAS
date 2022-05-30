\name{QQ_plot}
\alias{QQ_plot}
\alias{QQ_series}
\title{QQ plot(s) of expected vs. reported p-values}
\description{
  \code{QQ_plot} generates a simple QQ plot of the expected and
  reported p-value distribution. It includes the option to
  filter the data with the high-quality filter. \code{QQ_series}
  generates a series of such QQ plots for multiple filter
  settings.
}
\usage{
QQ_plot(dataset, save_name = "dataset", save_dir = getwd(),
        filter_FRQ = NULL, filter_cal = NULL,
        filter_HWE = NULL, filter_imp = NULL,
        filter_NA = TRUE,
        filter_NA_FRQ = filter_NA, filter_NA_cal = filter_NA,
        filter_NA_HWE = filter_NA, filter_NA_imp = filter_NA,
        p_cutoff = 0.05, plot_QQ_bands = FALSE,
        header_translations,
        check_impstatus = FALSE, ignore_impstatus = FALSE,
        T_strings = c("1", "TRUE", "yes", "YES", "y", "Y"),
        F_strings = c("0", "FALSE", "no", "NO", "n", "N"),
        NA_strings = c(NA, "NA", ".", "-"), ...)
QQ_series(dataset, save_name = "dataset", save_dir = getwd(),
          filter_FRQ = NULL, filter_cal = NULL,
          filter_HWE = NULL, filter_imp = NULL,
          filter_NA = TRUE,
          filter_NA_FRQ = filter_NA, filter_NA_cal = filter_NA,
          filter_NA_HWE = filter_NA, filter_NA_imp = filter_NA,
          p_cutoff = 0.05, plot_QQ_bands = FALSE,
          header_translations,
          check_impstatus = FALSE, ignore_impstatus = FALSE,
          T_strings = c("1", "TRUE", "yes", "YES", "y", "Y"),
          F_strings = c("0", "FALSE", "no", "NO", "n", "N"),
          NA_strings = c(NA, "NA", ".", "-"), ...)
}
\arguments{
  \item{dataset}{a data frame containing
    the p-value column and (depending on the settings) columns
    for chromosome number, position, the quality parameters,
    sample size and imputation status.}
  \item{save_name}{for \code{QQ_plot}, a character string;
    for \code{QQ_series}, a vector of character strings;
    specifying the filename(s) of the graph, \emph{without}
    extension.}
  \item{save_dir}{character string; the directory where the
    output files are saved. Note that R uses \emph{forward}
    slash (/) where Windows uses the backslash (\\).}
  \item{filter_FRQ, filter_cal, filter_HWE, filter_imp}{
    Filter threshold-values for allele-frequency, callrate,
    HWE p-value and imputation quality, respectively. Passed to
    \code{\link{HQ_filter}}. \code{QQ_plot} takes only
    single values, but  \code{QQ_series} accepts vectors
    as well (see 'details').}
  \item{filter_NA}{logical; if \code{TRUE}, then missing filter
    variables will be excluded; if \code{FALSE}, they will be
    ignored. \code{QQ_plot} takes only single values, but 
    \code{QQ_series} accepts vectors as well (see
    'Details'). \code{filter_NA} is the default setting for all
    variables; variable-specific settings can be specified with
    the following arguments.}
  \item{filter_NA_FRQ, filter_NA_cal,
        filter_NA_HWE, filter_NA_imp}{
    logical; variable-specific settings for \code{filter_NA}.
    These arguments are passed to \code{\link{HQ_filter}}.}
  \item{p_cutoff}{numeric; the threshold of p-values to be
    shown in the QQ plot(s). Higher (less significant) p-values
    are excluded from the plot. The default setting is \code{0.05},
    which excludes 95\% of data-points. It's not recommended to
    increase the value above \code{0.05}, as this may
    dramatically increase running time and memory usage.}
  \item{plot_QQ_bands}{logical; should probability bands be
    added to the QQ plot?}
  \item{header_translations}{translation table for column names.
    See \code{\link{translate_header}} for more information. If
    the argument is left empty, \code{dataset} is assumed to use
    the standard column-names of \code{\link{QC_GWAS}}.}
  \item{check_impstatus}{logical; should the imputation-status
    column be passed to \code{\link{convert_impstatus}}?}
  \item{ignore_impstatus}{logical; if \code{FALSE}, HWE p-value
    and callrate filters are applied only to genotyped SNPs, and
    imputation quality filters only to imputed SNPs. If
    \code{TRUE}, the filters are applied to all SNPs regardless
    of the imputation status.}
  \item{T_strings, F_strings, NA_strings}{arguments passed
    to \code{\link{convert_impstatus}}.}
  \item{\dots}{arguments passed to \code{\link{plot}}.}
}
\details{
  \code{QQ_series} accepts multiple filter-values, and
  passes these one by one to \code{QQ_plot} to generate a
  series of plots. For example, specifying:
  
  \code{filter_FRQ = c(0.05, 0.10), filter_cal = c(0.90, 0.95)}
  
  will generate two plots. The first excludes SNPs with
  allele frequency < 0.05 or callrate < 0.90; the second allele 
  frequency < 0.10 or callrate < 0.95. The same principle
  applies to the \code{NA_filter} settings. If the vectors
  submitted to the filter arguments are of unequal length, the
  shorter vector will be recycled until it equals the length of
  the longer (if possible). To filter missing values only, set
  the filter to \code{NA} and the corresponding NA-filter
  argument to \code{TRUE}. Setting the filter argument to
  \code{NULL} will disable the filter entirely, regardless of
  the NA-filter setting.
}
\value{Both functions return an invisible value \code{NULL}.}
\seealso{
  \code{\link{QC_plots}} for generating more complex QQ plots
  as well as Manhattan plots.
  
  \code{\link{QC_histogram}} for creating histograms.
  
  \code{\link{check_P}} for comparing the reported p-values to
  the p expected from the effect size and standard error.
}
\examples{
  \dontrun{
    data("gwa_sample")
  
    QQ_plot(dataset = gwa_sample,
            save_name = "sample_QQ",
            filter_FRQ = 0.01, filter_cal = 0.95,
            filter_NA = FALSE)
  
    QQ_series(dataset = gwa_sample,
              save_name = "sample_QQ",
              filter_FRQ = c(NA, 0.01, 0.01),
              filter_cal = c(NA, 0.95, 0.95),
              filter_NA = c(FALSE, FALSE, TRUE))
  }
}
\keyword{distribution}
