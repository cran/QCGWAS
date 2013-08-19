\name{switch_strand}
\alias{switch_strand}
\title{Convert alleles to the opposing DNA strand}
\description{
  This function is a subroutine of \code{\link{QC_GWAS}} and
  \code{\link{match_alleles}}. It converts allele-pairs to
  the configuration on the opposing DNA strand.
}
\usage{switch_strand(input, strand_col = FALSE)}
\arguments{
  \item{input}{table with the alleles in column 1 and 2, and
    (optionally) the strand-information (coded as \code{"+"} or
    \code{"-"}) in column 3. Note: the alleles must be uppercase
    characters.}
  \item{strand_col}{logical; if strand-information is present,
    this switches the sign in column 3 as well.}
}
\value{
  A table with two or three columns, depending on \code{strand_col}.
}
\seealso{\code{\link{match_alleles}}}
\examples{
  data("gwa_sample")

  switched_data <- gwa_sample[ , c("MARKER", "EFFECT_ALL",
                                   "OTHER_ALL", "STRAND")]
  switched_data[1:10, ]

  switched_data[ , 2:4] <- switch_strand(input =
            switched_data[ , 2:4], strand_col = TRUE)
  switched_data[1:10, ]
}
