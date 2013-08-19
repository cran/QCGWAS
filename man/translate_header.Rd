\name{translate_header}
\alias{translate_header}
\title{Translate column names into standard names}
\description{
  This function is used to translate non-standard column names
  into the standard ones used by \code{\link{QC_GWAS}} and other
  functions. It can also translate the standard names into
  user-specified names (via the \code{out_header} argument of
  \code{QC_GWAS}).
}
\usage{
translate_header(header,
    standard = c("MARKER", "CHR", "POSITION", "EFFECT_ALL",
                 "OTHER_ALL", "STRAND", "EFFECT", "STDERR",
                 "PVALUE", "EFF_ALL_FREQ", "HWE_PVAL",
                 "CALLRATE", "N_TOTAL", "IMPUTED",
                 "USED_FOR_IMP", "IMP_QUALITY"),
    alternative)
}
\arguments{
  \item{header}{character vector; the header to be translated.}
  \item{standard}{character vector; the names \code{header}
    should be translated into.}
  \item{alternative}{translation table; see 'Details' for more information.}
}
\details{
  In a nutshell: the \code{header} argument contains the names you
  have; the \code{standard} argument contains the names you want; and
  \code{alternative} is the conversion table.
  
  The table in \code{alternative} should have two columns. The
  left column contains the standard names; the right column
  possible alternatives. Only one alternative name should be
  listed per row. \code{translate_header} automatically changes
  the contents of \code{header} to uppercase, so \code{standard}
  and the right column of \code{alternative} should be uppercase
  as well.
  
  A sample translation table can be found in
  "R\\library\\QCGWAS\\doc". (It not recommended to load it into
  R with the \code{\link{data}} function, as this automatically
  converts the characters into factors. Use
  \code{\link{read.table}} instead, or put it in the reference
  directory and pass the filename to \code{\link{QC_GWAS}} by
  the \code{header_translations} argument. \code{QC_GWAS} will
  then load the table automatically.)
}
\value{
  \code{translate_header} returns an object of class 'list' with 6 components:
  \item{header_h}{character vector; the translated header.
    Unknown columns are included under their old names.}
  \item{missing_h}{character vector; the standard column names
    that were not found. If none, this returns \code{NULL}.}
  \item{unknown_h}{character vector; column names that could not
    be converted to a standard name. Note that these columns
    are also included in \code{header_h}. If none, this returns
    \code{NULL}.}
  \item{header_N, missing_N, unknown_N}{integer; the lengths of
    the above three vectors}
}
\seealso{\code{\link{identify_column}}}
\examples{
  sample_data <-
    data.frame(SNP = paste("rs", 1:10, sep = ""),
               chrom = 2,
               effect = 1:10/10,
               misc = NA,
               stringsAsFactors = FALSE)
  # Creates a table with four columns:
  #   SNP, chrom, effect and misc.

  ( alt_headers <-
      data.frame(
        standard = c("MARKER", "MARKER", "CHR", "CHR"),
        alternative = c("MARKER", "SNP", "CHR", "CHROM"),
        stringsAsFactors = FALSE) )
  # Creates the translation table
  #  with the standard names in column 1 and the alternatives
  #  in column 2.
  
  ( header_info <- 
      translate_header(header = names(sample_data),
        standard = c("MARKER", "CHR", "EFFECT"),
        alternative = alt_headers)  )
        
  # Despite not being in the translation table, EFFECT is
  #  changed to uppercase because it is present in standard.
  #  misc is neither in standard or the translation table, so
  #  it is marked as unknown and left unchanged.

  names(sample_data) <- header_info$header_h
}
