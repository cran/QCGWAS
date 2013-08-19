\name{convert_impstatus}
\alias{convert_impstatus}
\title{Convert imputation-status values to the QCGWAS standard}
\description{
  Converts imputation-status data to the standard values used
  by \code{\link{QC_GWAS}}
}
\usage{
convert_impstatus(impstatus,
                  T_strings = c("1", "TRUE", "T"),
                  F_strings = c("0", "FALSE", "F"),
                  NA_strings = c(NA, "NA", ".", "-"),
                  use_log = FALSE, ...)
}
\arguments{
  \item{impstatus}{vector of imputation-status data.}
  \item{T_strings}{character-string(s) indicating imputed data.}
  \item{F_strings}{character-string(s) indicating genotyped
    data.}
  \item{NA_strings}{character-string(s) indicating missing
    data.}
  \item{use_log, \dots}{arguments used by \code{\link{QC_GWAS}};
    redundant when \code{convert_impstatus} is used separately.}
}
\details{
  This function is used to convert the imputation-status column
  into the standard format, where \code{0} is genotyped and
  \code{1} is imputed. Untranslated values (i.e. strings that do
  not appear in any of the string arguments) will trigger a
  warning message and are set to \code{NA}.
}
\value{
  Numeric vector with values \code{0} for genotyped, \code{1}
  for imputed and \code{NA} for unknown data.
}
\note{
  The implementation of this function has changed. Previously,
  only character data was translated using the string
  arguments. Since version 1.0-4, the string arguments are used
  for all data types, so the user can determine the conversion
  of logical and numeric values as well.

  However, this does mean that values \code{1} and \code{0} are
  no longer converted automatically. They must be specified in
  the string arguments, or else the function will report them as
  untranslated and converts them to \code{NA}. The same applies
  to \code{TRUE}, \code{FALSE} and \code{NA}. Note the
  difference between the character string \code{"NA"} and the
  value \code{NA}.

  Finally, if the imputation-status column contains character
  strings, the main QC function \code{\link{QC_GWAS}} requires
  that all values are translated. If not, \code{QC_GWAS} will
  abort the QC.
}
\examples{
  status1 <- c(0,0,0,1,1,1,2,NA)
  status2 <- c(TRUE, FALSE, TRUE, FALSE, NA)
  status3 <- c("imputed", "genotyped", "NA", NA)
  status4 <- c(status1, status2, status3)
  
  ( outcome1 <- convert_impstatus(status1,
                                  T_strings = 1,
                                  F_strings = 0,
                                  NA_strings = NA) )
  # status 1 contains an untranslated value "2", which is
  # converted to NA. To avoid the warning message:
  
  ( outcome1 <- convert_impstatus(status1,
                                  T_strings = 1,
                                  F_strings = 0,
                                  NA_strings = c(NA, 2)) )
  
  
  ( outcome2 <- convert_impstatus(status2,
                                  T_strings = TRUE,
                                  F_strings = FALSE,
                                  NA_strings = NA) )
  
  ( outcome3 <- convert_impstatus(status3,
                                  T_strings = "imputed",
                                  F_strings = "genotyped",
                                  NA_strings = c("NA", NA)) )
  # Note that NA_strings includes both the character-string "NA"
  # and the value NA. Otherwise, one of the two would go 
  # "untranslated" and trigger a warning message.
  
  
  # And to check them all together
  ( outcome4 <- convert_impstatus(status4,
                          T_strings = c(1, TRUE, "imputed"),
                          F_strings = c(0, FALSE, "genotyped"),
                          NA_strings = c("NA", NA, 2)) )
}
\keyword{classes}
