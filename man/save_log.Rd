\name{save_log}
\alias{save_log}
\title{Log entries for QC_GWAS}
\description{
A subroutine for reporting problems encountered while running \code{\link{QC_GWAS}}. It's not designed for use outside of this context, so we recommend that users do not bother with it. It is described here for the sake of completeness only.
}
\usage{
save_log(phaseL, checkL, typeL, SNPL = allSNPs,
         allSNPs = 1L, actionL, noteL = "", fileL)
}
\arguments{
  \item{phaseL, checkL, typeL}{character-strings indicating (very briefly) what occurred to trigger a log entry.}
  \item{SNPL}{integer value indicating the number of SNPs affected.}
  \item{allSNPs}{integer value indicating the total number of SNPs in the dataset.}
  \item{actionL}{character-string describing very briefly the response to situation.}
  \item{noteL}{character-string for longer explanations of what happened.}
  \item{fileL}{character-string of the directory path and the file name,
    \emph{without} file extension, of the log file.}
}
\details{
\code{save_log} is used to add entries to the log file the moment they are triggered (as opposed to waiting until the QC concludes and then saving the entire log). In case of a fatal crash, the log file should therefor give some indication of where in the QC it occurred.

\code{save_log} does not create the log file (this is done by \code{\link{QC_GWAS}}): it merely appends entries to the bottom of the specified file, regardless of what that file is. Hence it is not recommended to use the function outside of \code{QC_GWAS}.
}
\value{\code{save_log} returns an invisible \code{NULL}. The relevant output is the entry added to the log file.}
\keyword{IO}
\keyword{print}
