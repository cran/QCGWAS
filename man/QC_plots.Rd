\name{QC_plots}
\alias{QC_plots}
\title{QQ and Manhattan plots}
\description{
  This function creates the most important graphs of the QC: the
  QQ plots and the Manhattan plot. It also calculates lambda,
  and determines the effect of the filters.
}
\usage{
QC_plots(dataset,
         plot_QQ = TRUE, plot_Man = TRUE,
         FRQfilter_values = NULL, FRQfilter_NA = filter_NA,
         HWEfilter_values = NULL, HWEfilter_NA = filter_NA,
         calfilter_values = NULL, calfilter_NA = filter_NA,
         impfilter_values = NULL, impfilter_NA = filter_NA,
         impfilter_min = min(dataset$IMP_QUALITY, na.rm = TRUE),
         manfilter_FRQ = NULL, manfilter_HWE = NULL,
         manfilter_cal = NULL, manfilter_imp = NULL,
         filter_NA = TRUE,
         plot_cutoff_p = 0.05, plot_names = FALSE,
         QQ_colors = c("red", "blue", "orange", "green3", "yellow"),
         plot_QQ_bands = FALSE,
         save_name = "dataset", save_dir = getwd(),
         header_translations, use_log = FALSE,
         check_impstatus = FALSE, ignore_impstatus = FALSE,
         T_strings = c("1", "TRUE", "yes", "YES", "y", "Y"),
         F_strings = c("0", "FALSE", "no", "NO", "n", "N"),
         NA_strings = c(NA, "NA", ".", "-"))
}
\arguments{
  \item{dataset}{vector of p-values or a data frame containing
    the p-value column and (depending on the settings) columns
    for chromosome number, position, the quality parameters,
    sample size and imputation status.}
  \item{plot_QQ, plot_Man}{logical; should QQ and Manhattan plots
    be saved?}
  \item{FRQfilter_values, HWEfilter_values,
        calfilter_values, impfilter_values}{
    numeric vectors; the threshold values for the QQ plot filters.
    The filters are for allele-frequency, HWE p-values, callrate
    and imputation-quality parameters, respectively. A maximum of
    five values can be specified per parameter.
    \itemize{
      \item Set to \code{NULL} to disable the QQ filter for that
        parameter.
      \item To filter missing values only, set the filter value
        to \code{NA} and the corresponding \code{filter_NA}
        argument to \code{TRUE}.
      \item The allele-frequency filter is two-sided: for a
        filter-value of \code{x}, it will exclude entries
        with \code{freq < x} and \code{freq > 1 - x}.
      \item Values >= 1 will be divided by the SNP's sample size.
        This allows sample-size dependent filtering of allele
        frequencies. Note that this uses the sample size reported
        in the sample-size column for that specific SNP. SNPs
        without sample size will be excluded if the corresponding
        \code{filter_NA} argument is \code{TRUE} and ignored if
        it is \code{FALSE}.
    }}
  \item{FRQfilter_NA, HWEfilter_NA, calfilter_NA,
        impfilter_NA, filter_NA}{
    logical; should the filters exclude (\code{TRUE}) or
    ignore (\code{FALSE}) missing values? \code{filter_NA}
    is the default setting, the others allow
    variable-specific settings.}
  \item{impfilter_min}{numeric; the lowest possible value for
    imputation-quality. This argument is currently redundant, as
    it is calculated automatically.}
  \item{manfilter_FRQ, manfilter_HWE, manfilter_cal, manfilter_imp}{
    single, numeric values; the filter-settings for allele-frequency,
    HWE p-values, callrate and imputation quality respectively,
    for the Manhattan plot. The arguments are passed to
    \code{\link{HQ_filter}}.
    To filter missing values only, set to \code{NA} and the
    corresponding \code{filter_NA} argument to \code{TRUE}. To
    disable filtering entirely, set to \code{NULL}.}
  \item{plot_cutoff_p}{numeric; the threshold of p-values to be
    shown in the QQ & Manhattan plots. Higher (less
    significant) p-values are excluded from the plot. The default
    setting is \code{0.05}, which excludes 95\% of data-points.
    It's not recommended to increase the value above \code{0.05},
    as this may dramatically increase running time and memory usage.}
  \item{plot_names}{argument currently redundant.}
  \item{QQ_colors}{vector of R color-values; the color of the
    QQ filter-plots. The unfiltered data is black by default.
    This argument sets the colors of the least (first value) to
    most (last value) stringent filters. (For this setting,
    filter values \code{>= 1} (i.e. sample-size based filtering)
    are considered less stringent than values \code{< 1}.)}
  \item{plot_QQ_bands}{logical; should probability bands be
    added to the QQ plot?}
  \item{save_name}{character string; the filename, \emph{without}
    extension, for the graphs.}
  \item{save_dir}{character string; the directory where the
    graphs are saved. Note that R uses \emph{forward}
    slash (/) where Windows uses the backslash (\\).}
  \item{header_translations}{translation table for column names.
    See \code{\link{translate_header}} for more information. If
    the argument is left empty, \code{dataset} is assumed to use
    the standard column-names of \code{\link{QC_GWAS}}.}
  \item{use_log}{argument used by \code{\link{QC_GWAS}};
    redundant when \code{QC_plots} is used separately.}
  \item{check_impstatus}{logical; should the imputation-status
    column be passed to \code{\link{convert_impstatus}}?}
  \item{ignore_impstatus}{logical; if \code{FALSE}, HWE p-value
    and callrate filters are applied only to genotyped SNPs, and
    imputation quality filters only to imputed SNPs. If
    \code{TRUE}, the filters are applied to all SNPs regardless
    of the imputation status.}
  \item{T_strings, F_strings, NA_strings}{arguments passed
    to \code{\link{convert_impstatus}}.}
}
\details{
  The function \code{QC_plots} grew out of phase 4 of
  \code{\link{QC_GWAS}}. It carries out three functions, hence
  the vague name: it calculates  lambda, it applies the
  QQ filters, and it creates the QQ and Manhattan plots (a
  separate function is available for creating
  regional-association plots: see below). The function schematic
  is as follows:
  \itemize{
    \item Preparing the dataset: this step involves translating
      the dataset header to the standard column-names (by
      \code{\link{identify_column}}) and converting imputation
      status (by \code{\link{convert_impstatus}}). Both steps
      are optional, and are disabled by default. If the function
      cannot identify the imputation status column, it will
      generate a warning message and disable the
      imputation-status dependent filters.
    \item Calculating the QC stats: here it generates the
      filters an calculates how many SNPs are removed. Lambda
      is also calculated at this point.
    \item Creating a QQ graph of every variable for which
      filters have been specified. Every graph contains an
      unfiltered plot, plus plots for every effective filter.
      ("Effective" means "excludes more SNPs than the
      previous, less-stringent filter".)
    \item Creating the Manhattan plot. The default Manhattan
      plot covers chromosomes 1 to 23 (X). Fields for XY, Y and
      M are added when such SNPs are present.
  }
}
\value{
  An object of class 'list' with the following components:
  \item{lambda}{vector of the lambda values of all SNPs,
    genotyped SNPs and imputed SNPs, respectively.}
  \item{ignore_impstatus}{logical value indicating whether
    imputation status was used when applying the filters.}
  \item{FRQfilter_names, HWEfilter_names, calfilter_names, impfilter_names}{
    character vectors naming the specified QQ filters.}
  \item{FRQfilter_N, HWEfilter_N, calfilter_N, impfilter_N}{
    numeric vectors; the number of SNPs removed by the specified
    filters. Note that the filters are sorted before being
    applied, so the order may not match that of the input.
    Check the \code{filter_names} output to see the order
    that was used inside \code{QC_plots}.}
  \item{Manfilter_N}{numeric; the number of SNPs removed by the
    Manhattan filter. This does not include those SNPs removed
    because they lacked p or chromosome/position-values, or failed
    the p-cutoff threshold.}
}
\note{
  By default, \code{QC_plots} expects \code{dataset} to use the
  standard column-names used by \code{\link{QC_GWAS}}. A
  translation table can be specified in \code{header_translations}
  to allow non-standard names. See \code{\link{translate_header}}
  for more information.
  
  The function accepts both integer and character chromosome
  values. Character values of \code{"X"}, \code{"Y"}, \code{"XY"}
  and \code{"M"} are automatically converted to integers. By
  default, the Manhattan plot shows all autosomal
  chromosomes and chromosome X. Fields for Y, XY and M are
  added only when such SNPs are present.
  
  There must be more than 10 p-values at or below the
  \code{plot_cutoff_p} threshold for the QQ and Manhattan
  plots to be created.  
}
\seealso{
  \code{\link{plot_regional}} for creating a regional association
  plot.
  
  \code{\link{check_P}} for comparing the reported p-values to
  the p expected from the effect size and standard error.
  
  \code{\link{QQ_plot}} for generating simpler QQ plots.
}
\examples{
  \dontrun{
  data("gwa_sample")

  QC_plots(dataset = gwa_sample,
    plot_QQ = TRUE, plot_QQ_bands = TRUE, plot_Man = TRUE,
    FRQfilter_values = c(NA, 0.01, 0.05, 3),
    calfilter_values = c(NA, 0.95, 0.99),
    manfilter_FRQ = 0.05, manfilter_cal = 0.95,
    filter_NA = TRUE, save_name = "sample_plots")
    }
}
\keyword{distribution}
\keyword{univar}
