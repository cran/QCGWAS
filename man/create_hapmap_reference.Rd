\name{create_hapmap_reference}
\alias{create_hapmap_reference}
\title{Create an allele-reference file from HapMap data}
\description{
  This function creates the standard allele reference file, as
  used by \code{\link{QC_GWAS}} and \code{\link{match_alleles}},
  from data publicly available at the website of the
  international HapMap project (see 'References').
}
\usage{
create_hapmap_reference(dir = getwd(),
   download_hapmap = FALSE, download_subset,
   hapmap_files = list.files(path = dir, pattern = "freqs_chr"),
   filename = "allele_reference_HapMap",
   save_txt = TRUE, save_rdata = !save_txt,
   return_reference = FALSE)
}
\arguments{
  \item{dir}{character string; the directory of
    the input and output files. Note that R uses \emph{forward}
    slash (/) where Windows uses the backslash (\\).}
  \item{download_hapmap}{logical; if \code{TRUE}, the required
    allele-frequency files are downloaded from the HapMap website
    into \code{dir}, and then turned into a reference. If
    \code{FALSE}, the files specified in \code{hapmap_files} are
    used.}
  \item{download_subset}{character-string; indicates the
    population to download for creating the reference. Options
    are: ASW, CEU, CHB, CHD, GIH, JPT, LWK, MEX, MKK, TSI, YRI.}
\item{hapmap_files}{character vector of the filenames of
    HapMap frequency-files to be included in the reference. The
    default option includes all files with the string
    "freqs_chr" in the filename. (This argument is only
    used when \code{download_hapmap} is \code{FALSE}.)}
  \item{filename}{character string; the name of the output file,
    \emph{without} file-extension.}
  \item{save_txt, save_rdata}{logical; should the reference be
    saved as a tab-delimitated text file and/or an RData file?
    If saved as RData, the object name \code{allele_ref_std} is
    used for the reference table.}
  \item{return_reference}{logical; should the function return
    the reference as it output value?}
}
\details{
  The function removes SNPs with invalid alleles and with allele
  frequencies that do not add up to \code{1}. It also removes
  all instances of duplicate SNPids. If such entries are
  encountered, a warning is printed in the R console and the
  entries are saved in a .txt file in the output directory.

  Like the \code{QC_GWAS}, \code{create_hapmap_reference} codes
  the X chromosome as \code{23}, Y as \code{24}, XY (not
  available on HapMap website) as \code{25} and M as \code{26}.

  Both the .RData export and  the function return store the
  alleles as factors rather than character strings.
}
\value{
  If \code{return_reference} is \code{TRUE}, the function
  returns the generated reference table. If \code{FALSE}, it
  returns an invisible \code{NULL}.
}
\references{
  The required data is available at the Website of the
  International HapMap project, under bulk data downloads >
  bulk data > frequencies

  http://hapmap.ncbi.nlm.nih.gov
  
  The HapMap files downloaded by this function are subject to
  the HapMap terms and policies. See:
  http://hapmap.ncbi.nlm.nih.gov/datareleasepolicy.html
}
\seealso{
  \code{\link{match_alleles}}
}
\examples{
  # This command will download the CEU HapMap dataset and use
  # it to generate an allele-reference. Create a folder
  # "new_hapmap" to store the data and make sure there is
  # sufficient disk space and a reasonably fast internet
  # connection.

  #new_hapmap <- create_hapmap_reference(dir = "C:/new_hapmap",
  # download_hapmap = TRUE, download_subset = "CEU",
  # filename = "new_hapmap", save_txt = TRUE,
  # return_reference = TRUE)
}
