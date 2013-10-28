\name{header_translations}
\alias{header_translations}
\docType{data}
\title{
  Translation table for GWAS dataset headers
}
\description{
  This is a sample translation table, as used by
  \code{\link{translate_header}} to translate dataset column
  names to the \code{QCGWAS} standard. An editable .txt version
  can be found in "R\\library\\QCGWAS\\doc".
}
\usage{data(header_translations)}
\format{
  A data frame with 104 observations on the following 2 variables.
  \describe{
    \item{\code{STANDARD}}{a character vector; standard column
      names as used by QCGWAS.}
    \item{\code{ALTERNATIVE}}{a character vector; alternative
      column names, as used in the datasets.}
  }
}
\keyword{datasets}
