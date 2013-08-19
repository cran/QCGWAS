\name{HQ_filter}
\alias{HQ_filter}
\title{Select high-quality data in GWAS datasets}
\description{
  This function accepts a \code{\link{QC_GWAS}} dataset and
  returns a vector of logical values indicating which entries
  meet the quality criteria.
}
\usage{
HQ_filter(data,
          ignore_impstatus = FALSE,
          FRQ_val = NULL, HWE_val = NULL,
          cal_val = NULL, imp_val = NULL,
          filter_NA = TRUE,
          FRQ_NA = filter_NA, HWE_NA = filter_NA,
          cal_NA = filter_NA, imp_NA = filter_NA)
}
\arguments{
  \item{data}{table to be filtered. \code{HQ_filter}
    assumes the dataset uses the standard \code{\link{QC_GWAS}}
    column names.}
  \item{ignore_impstatus}{logical; if \code{FALSE}, HWE p-value
    and callrate filters are applied only to genotyped SNPs, and
    imputation quality filters only to imputed SNPs. If
    \code{TRUE}, the filters are applied to all SNPs regardless
    of the imputation status.}
  \item{FRQ_val, HWE_val, cal_val, imp_val}{
    numeric; the \emph{minimal} required value for allele frequency,
    HWE p-value, callrate and imputation quality respectively.
    Note that the allele-frequency filter is two-sided: for a
    filter-value of \code{x}, it will exclude entries with
    \code{freq < x} and \code{freq > 1 - x}.}
  \item{filter_NA}{
    logical; if \code{TRUE}, then missing filter variables
    will be excluded; if \code{FALSE}, they will be ignored.
    \code{filter_NA} is the default setting for all variables.
    Variable-specific settings can be specified with the
    following arguments.}
  \item{FRQ_NA, HWE_NA, cal_NA, imp_NA}{logical; variable-specific 
    settings for \code{filter_NA}.}
}
\details{
  A SNP is considered high-quality if it meets all quality
  criteria. The thresholds are inclusive; i.e. SNPs that have a
  value equal or higher than the threshold will be considered
  high-quality.

  To filter missing values only, set the filter argument to
  \code{NA}, and the corresponding NA-filter to \code{TRUE}.

  To disable filtering entirely, set to \code{NULL}. This
  disables the filtering of missing values as well.

  When imputation status is missing or invalid (and
  \code{ignore_impstatus} is \code{FALSE}), only the
  allele-frequency filter will be applied.
}
\value{
  A vector of logical values, indicating which values in
  \code{data} meet (\code{TRUE}) or fail (\code{FALSE}) the
  quality criteria.
}
\note{
  The table entered in the \code{data} argument must use the
  standard column names of \code{\link{QC_GWAS}}. Functions
  using \code{HQ_filter} usually allow the user to specify a
  translation table. If not, \code{\link{translate_header}} can
  be used to translate the header manually.
}
\examples{
  data("gwa_sample")

  selected_SNPs <- HQ_filter(data = gwa_sample,
                             FRQ_val = 0.01,
                             cal_val = 0.95,
                             filter_NA = FALSE)
  summary(gwa_sample[selected_SNPs, ])
  
  selected_SNPs <- HQ_filter(data = gwa_sample,
                             FRQ_val = 0.01,
                             cal_val = 0.95,
                             filter_NA = FALSE,
                             ignore_impstatus = TRUE)
  summary(gwa_sample[selected_SNPs, ])
}
