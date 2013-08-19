\name{load_GWAS}
\alias{load_GWAS}
\alias{load_test}
\title{Easy loading of GWAS results files}
\description{
  \code{load_GWAS} is wrapper-function of \code{\link{read.table}}
  that makes loading large GWAS results files less of a hassle. It
  automatically unpacks .zip and .gz files and uses
  \code{load_test} to determine which column separator the file
  uses.
}
\usage{
load_GWAS(filename, dir = getwd(),
          column_separators = c("\\t", " ", "", ",", ";"),
          test_nrows = 1000,
          header = TRUE, nrows = -1,
          comment.char = "", na.strings = c("NA", "."),
          stringsAsFactors = FALSE, ...)

load_test(filename, dir = getwd(),
          column_separators = c("\\t", " ", "", ",", ";"),
          test_nrows = 1000, ...)
}
\arguments{
  \item{filename}{character string; the complete filename of the
    file to be loaded. Note that compressed files (.gz or .zip
    files) can only be unpacked if the filename of the archive
    contains the extension of the archived file. For example, if
    the archived file is named \code{"data1.csv"}, the archive
    should be \code{"data1.csv.zip"}.}
  \item{dir}{character string; the directory containing the file.
    Note that R uses \emph{forward} slash (/) where Windows uses
    backslash (\\).}
  \item{column_separators}{character string or vector of
    the column-separators to be tried by \code{load_test}.
    White-space can be specified by \code{""}, but it is
    recommended you try tab (\code{"\\t"}) and space
    (\code{" "}) first.}
  \item{test_nrows}{integer; the number of lines that
    \code{load_test} checks in the trail-load. A smaller number
    means faster loading, but also makes it more likely that
    errors slip through. To check the entire dataset, set to
    \code{-1}.}
  \item{header, nrows, comment.char, na.strings,
        stringsAsFactors, \dots}{
    Arguments passed to \code{\link{read.table}}.}
}
\details{
  \code{load_test} determines the correct column separator
  simply by trying them individually until it finds one that
  works (that is: one that results in a dataset with an equal
  number of cells in every row AND at least five or more
  columns). If none work, it reports the error-message generated
  by the last column separator tried.
  
  The column separators are tried in the order specified by the
  \code{column_separators} argument.
  
  By default, \code{load_test} only checks the first 1000 lines
  (adjustable by the \code{test_nrows} argument); if the problem
  lies further down in the dataset, it will not catch it. In such
  a case, \code{load_GWAS} and \code{\link{QC_GWAS}} will crash
  when attempting to load the dataset.
  
  A common problem is employing white-space (\code{""}) as
  column separator for a file that uses empty fields to indicate
  missing values. The separators surrounding an empty field are
  adjacent, so R parses them as a single column separator. In
  this particular example, specifying a single space
  (\code{" "}) or tab (\code{"\\t"}) as column separator solves
  the problem (this is why the default setting of
  \code{column_separators} puts these values before white-space).
}
\value{
  \code{load_GWAS} returns the table imported from the specified file.
  
  \code{load_test} returns a list with 4 components:
  \item{success}{logical; whether \code{load_test} was able to
    load a dataset with five or more columns.}
  \item{error}{character string; if unable to load the file, this
    returns the error-message of the last column separator to be
    tried.}
  \item{file_type}{character string; the last three characters
    of \code{filename}.}
  \item{sep}{the first column-separator that succeeded in loading
    a dataset with five or more columns.}
}
\note{
  \code{load_GWAS} uses the same default loading-settings as
  \code{QC_GWAS}. \code{load_test}, on the other hand, has no
  default values for \code{header}, \code{comment.char},
  \code{na.strings} and \code{stringsAsFactors}, and uses the
  \code{\link{read.table}} defaults instead.
}
\examples{
  ## As the function requires a GWAS file to work,
  ## the following code should be adjusted before execution.
  ## Because this is a demonstration, the nrows argument is used
  ## to read only the first 100 rows.
  
  # data_GWAS <-
  #    load_GWAS("GWA_results1.txt.zip",
  #              dir = "C:/GWAS_results",
  #              nrows = 100)
}
\keyword{files}
