\name{QC_histogram}
\alias{QC_histogram}
\alias{histogram_series}
\title{Histogram(s) of expected and observed data distribution}
\description{
  \code{QC_histogram} creates two histograms: one showing the
  observed data distribution of a numeric variable, and one
  showing the expected distribution.
  It includes the option to filter the data with the
  high-quality filter. \code{histogram_series} generates a
  series of such histograms for multiple filter settings.
}
\usage{
QC_histogram(dataset, data_col = 1,
             save_name = "dataset", save_dir = getwd(),
             export_outliers = FALSE,
             filter_FRQ = NULL, filter_cal = NULL,
             filter_HWE = NULL, filter_imp = NULL,
             filter_NA = TRUE,
             filter_NA_FRQ = filter_NA, filter_NA_cal = filter_NA,
             filter_NA_HWE = filter_NA, filter_NA_imp = filter_NA,
             breaks = "Sturges",
             graph_name = colnames(dataset)[data_col],
             header_translations, check_impstatus = FALSE,
             ignore_impstatus = FALSE,
             T_strings = c("1", "TRUE", "yes", "YES", "y", "Y"),
             F_strings = c("0", "FALSE", "no", "NO", "n", "N"),
             NA_strings = c(NA, "NA", ".", "-"), ...)
histogram_series(dataset, data_col = 1,
   save_name = paste0("dataset_F", 1:nrow(plot_table)),
   save_dir = getwd(), export_outliers = FALSE,
   filter_FRQ = NULL, filter_cal = NULL,
   filter_HWE = NULL, filter_imp = NULL,
   filter_NA = TRUE,
   filter_NA_FRQ = filter_NA, filter_NA_cal = filter_NA,
   filter_NA_HWE = filter_NA, filter_NA_imp = filter_NA,
   breaks = "Sturges",
   header_translations, ignore_impstatus = FALSE,
   check_impstatus = FALSE,
   T_strings = c("1", "TRUE", "yes", "YES", "y", "Y"),
   F_strings = c("0", "FALSE", "no", "NO", "n", "N"),
   NA_strings = c(NA, "NA", ".", "-"),
   ...)
}
\arguments{
  \item{dataset}{vector or table containing the variable of interest.}
  \item{data_col}{name or number of the column of \code{dataset}
    containing the variable of interest.}
  \item{save_name}{for \code{QC_histogram}, a character string;
    for \code{histogram_series}, a vector of character strings;
    specifying the filename(s) of the graph, \emph{without}
    extension.}
  \item{save_dir}{character string; the directory where the
    output files are saved. Note that R uses \emph{forward}
    slash (/) where Windows uses the backslash (\\).}
  \item{export_outliers}{logical or numeric value; should outlying
    entries (which are excluded from the plot) be exported to an
    output file? If numeric, the number specifies the max.
    number of entries that is exported.}
  \item{filter_FRQ, filter_cal, filter_HWE, filter_imp}{
    Filter threshold-values for allele-frequency, callrate,
    HWE p-value and imputation quality, respectively. Passed to
    \code{\link{HQ_filter}}. \code{QC_histogram} takes only
    single values, but  \code{histogram_series} accepts vectors
    as well (see 'details').}
  \item{filter_NA}{logical; if \code{TRUE}, then missing filter
    variables will be excluded; if \code{FALSE}, they will be
    ignored. \code{QC_histogram} takes only single values, but 
    \code{histogram_series} accepts vectors as well (see
    'Details'). \code{filter_NA} is the default setting for all
    variables; variable-specific settings can be specified with
    the following arguments.}
  \item{filter_NA_FRQ, filter_NA_cal, filter_NA_HWE, filter_NA_imp}{
    logical; variable-specific settings for \code{filter_NA}.
    These arguments are passed to \code{\link{HQ_filter}}.}
  \item{breaks}{argument passed to \code{\link{hist}}; determines
    the cell-borders in the histogram.}
  \item{graph_name}{character string; used in the title of the plot.}
  \item{header_translations}{translation table for column names.
    See \code{\link{translate_header}} for more information. If
    the argument is left empty, \code{dataset} is assumed to use
    the standard column-names used by \code{\link{QC_GWAS}}.}
  \item{check_impstatus}{logical; should
    \code{\link{convert_impstatus}} be called to convert the
    imputation-status column into standard values?}
  \item{ignore_impstatus}{logical; if \code{FALSE}, HWE p-value
    and callrate filters are applied only to genotyped SNPs, and
    imputation quality filters only to imputed SNPs. If
    \code{TRUE}, the filters are applied to all SNPs regardless
    of the imputation status.}
  \item{T_strings, F_strings, NA_strings}{arguments passed to
  \code{\link{convert_impstatus}}.}
  \item{\dots}{in \code{histogram_series}: arguments passed to
    \code{QC_histogram}; in \code{QC_histogram}, arguments passed
    to \code{\link{hist}}.}
}
\details{
  \code{histogram_series} accepts multiple filter-values, and
  passes these one by one to \code{QC_histogram} to generate a
  series of histograms. For example, specifying:
  
  \code{filter_FRQ = c(0.05, 0.10), filter_cal = c(0.90, 0.95)}
  
  will generate two histograms. The first excludes SNPs with
  allele frequency < 0.05 or callrate < 0.90; the second allele 
  frequency < 0.10 or callrate < 0.95. The same principle
  applies to the \code{NA_filter} settings. If the vectors
  submitted to the filter arguments are of unequal length, the
  shorter vector will be recycled until it equals the length of
  the longer (if possible). To filter missing values only, set
  the filter to \code{NA} and the corresponding NA-filter
  argument to \code{TRUE}. Setting the filter argument to
  \code{NULL} will disable the filter entirely, regardless of
  the NA filter setting.
}
\value{Both functions return an invisible value \code{NULL}.}
\seealso{For creating QQ plots: \code{\link{QQ_plot}}.}
\examples{
\dontrun{
  data("gwa_sample")

  QC_histogram(dataset = gwa_sample, data_col = "EFFECT",
             save_name = "sample_histogram",
             filter_FRQ = 0.01, filter_cal = 0.95,
             filter_NA = FALSE,
             graph_name = "Effect size histogram")

  histogram_series(dataset = gwa_sample, data_col = "EFFECT",
             save_name = "sample_histogram",
             filter_FRQ = c(NA, 0.01, 0.01),
             filter_cal = c(NA, 0.95, 0.95),
             filter_NA = c(FALSE, FALSE, TRUE))
}}
\keyword{distribution}
