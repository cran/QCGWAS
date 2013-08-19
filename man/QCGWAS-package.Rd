\name{QCGWAS-package}
\alias{QCGWAS-package}
\alias{QCGWAS}
\docType{package}
\title{Quality Control of Genome Wide Association Study results}
\description{
  Functions for automated and manual quality control of Genome
  Wide Association Study results.
}
\details{
\tabular{ll}{
Package: \tab QCGWAS\cr
Type: \tab Package\cr
Version: \tab 1.0-6\cr
Date: \tab 2013-08-19\cr
License: \tab GPL (>= 3)\cr
}
  The core of the package is the function \code{\link{QC_GWAS}}.
  This function carries out an automated quality-control (QC) of
  a Genome Wide Association Study (GWAS) results file, reporting
  on the data-distribution, checking the SNPs for missing and
  invalid data, comparing the alleles and allele-frequency to a
  reference, and creating QQ and Manhattan plots.
  
  Although the number of arguments in \code{QC_GWAS} may seem
  overwhelming, only three of them are required to run a basic
  QC. The name of the file to be QC'ed should be passed to the
  \code{filename} argument; the directory of said file to the
  \code{dir_data} argument; and a header-translation table to
  the \code{header_translations} argument. The results
  will be saved in a number of files and graphs in the data
  directory. For a quick introduction to QCGWAS, read the quick
  reference guide that can be found in "R\\library\\QCGWAS\\doc".
}
\author{
Peter J. van der Most and Ilja M. Nolte

Maintainer: Peter J. van der Most <p.j.van.der.most@umcg.nl>
}
\references{
  Van der Most, P.J., Vaez, A., Prins, B., Loretto Munoz, M.,
  Snieder, H., Alizadeh, B.Z. & Nolte, I.M. (In preparation).
  QCGWAS: A flexible R package for automated quality control of
  Genome-Wide Association results.
}
\keyword{ package }
\seealso{
  \code{\link{QC_GWAS}}
}
