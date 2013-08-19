\name{gwa_sample}
\alias{gwa_sample}
\docType{data}
\title{
  Sample dataset for the QCGWAS package
}
\description{
  A fake GWAS results dataset for use in the examples.
}
\usage{data(gwa_sample)}
\format{
  A data frame with 10000 observations on the following 15 variables.
  \describe{
    \item{\code{MARKER}}{a character vector; the SNP IDs.}
    \item{\code{STRAND}}{a character vector; the DNA strand of
      the listed alleles.}
    \item{\code{CHR}}{a character vector; the chromosome of the
      listed SNP.}
    \item{\code{POSITION}}{a numeric vector; the basepair
      position of the listed SNP.}
    \item{\code{EFFECT_ALL}}{a character vector; the effect
      allele.}
    \item{\code{OTHER_ALL}}{a character vector; the non-effect
      allele.}
    \item{\code{N_TOTAL}}{a numeric vector; the sample size.}
    \item{\code{EFF_ALL_FREQ}}{a numeric vector; frequency of
      the effect-allele in the dataset.}
    \item{\code{HWE_PVAL}}{a numeric vector; Hardy-Weinberg
      (HWE) p-value of the listed alleles.}
    \item{\code{CALLRATE}}{a numeric vector; SNP callrate.}
    \item{\code{EFFECT}}{a numeric vector; effect size of the
      listed effect allele.}
    \item{\code{STDERR}}{a numeric vector; standard error of the
      effect.}
    \item{\code{PVALUE}}{a numeric vector; p-value of the
      effect.}
    \item{\code{IMPUTED}}{a numeric vector; imputation status -
      i.e. whether the SNP is genotyped (\code{0}) or imputed
      (\code{1}).}
    \item{\code{IMP_QUALITY}}{a numeric vector; imputation
      quality.}
  }
}
\keyword{datasets}
