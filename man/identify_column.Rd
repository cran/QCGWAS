\name{identify_column}
\alias{identify_column}
\title{Identify non-standard column names}
\description{
  This function is a subroutine of \code{\link{translate_header}}
  and several other functions. It is used to translate non-standard
  column names into standard ones.
}
\usage{
identify_column(std_name, alt_names, header)
}
\arguments{
  \item{std_name}{character string; the standard name for
    the data-column.}
  \item{alt_names}{translation table; with in the left column
    the standard name(s) and in the right column possible
    alternatives. See \code{\link{translate_header}} for more
    details.}
  \item{header}{the column names of the dataset.
    The names should be entirely in uppercase.}
}
\details{
  The purpose of \code{identify_column} is essentially to look
  up in the translation table (\code{alt_names}) which of the
  names in \code{header} can be translated into \code{std_name}.
}
\value{
  An integer vector of the entry(s) in \code{header} (i.e. the
  column-numbers) that can be translated into \code{std_name}.
}
\seealso{
  \code{\link{translate_header}}
}
\examples{
  sample_data <-
    data.frame(SNP = paste("rs", 1:10, sep = ""),
               chrom = 2,
               effect = 1:10/10,
               misc = NA)
  sample_header <- toupper(names(sample_data))


  alt_headers <-
    data.frame(
      standard = c("MARKER", "MARKER", "CHR", "CHR"),
      alternative = c("MARKER", "SNP", "CHR", "CHROM"),
      stringsAsFactors = FALSE)
  
  identify_column(std_name = "EFFECT", alt_names = alt_headers,
                  header = sample_header)
  identify_column(std_name = "MARKER", alt_names = alt_headers,
                  header = sample_header)
  identify_column(std_name = "CHR", alt_names = alt_headers,
                  header = sample_header)
  identify_column(std_name = "MISC", alt_names = alt_headers,
                  header = sample_header)
}
