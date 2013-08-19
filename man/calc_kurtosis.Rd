\name{calc_kurtosis}
\alias{calc_kurtosis}
\alias{calc_skewness}
\title{Skewness and Kurtosis}
\description{Functions for calculating skewness and kurtosis}
\usage{
calc_kurtosis(input,
              FRQ_val = NULL, HWE_val = NULL,
              cal_val = NULL, imp_val = NULL, ...)
calc_skewness(input,
              FRQ_val = NULL, HWE_val = NULL,
              cal_val = NULL, imp_val = NULL, ...)
}
\arguments{
  \item{input}{either a vector of effect sizes or a data frame
    using the standard column names.}
  \item{FRQ_val, HWE_val, cal_val, imp_val, \dots}{arguments
    passed to \code{\link{HQ_filter}}.}
}
\details{
Kurtosis is calculated as:
    
\code{sum( (ES - mean(ES))^4) / ((length(ES)-1) * sd(ES)^4 )}
  
Skewness is calculated as:
  
\code{sum( (ES - mean(ES))^3) / ((length(ES)-1) * sd(ES)^3 )}
}
\value{
  Respectively the kurtosis and skewness of the input
  effect-size distribution.
}
\note{
  Both functions accept vectors as \code{input}. If \code{input}
  is a data frame, the column names must match the standard
  names used by \code{\link{QC_GWAS}} (\code{"EFFECT"} for
  effect sizes, \code{"EFF_ALL_FREQ"} for allele frequency, etc.)
}
\seealso{
  For plotting skewness and kurtosis:
  \code{\link{plot_skewness}}.
}
\examples{
  data("gwa_sample")  
  
  calc_kurtosis(gwa_sample$EFFECT)
  calc_kurtosis(gwa_sample)
  calc_kurtosis(gwa_sample$EFF_ALL_FREQ)
  calc_kurtosis(gwa_sample,
                FRQ_val = 0.05, cal_val = 0.95,
                filter_NA = FALSE)
}
\keyword{univar}
