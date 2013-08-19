\name{plot_regional}
\alias{plot_regional}
\title{Regional Association Plot}
\description{
  A regional association plot is essentially a zoomed-in
  Manhattan plot, allowing the researcher to look at associations
  in a small, pre-defined area of the genome.
}
\usage{
plot_regional(dataset,
              chr, start_pos, end_pos,
              plot_cutoff_p = 1,
              name_cutoff_p,
              data_name = NULL,
              save_name = "regional_association_plot",
              save_dir = getwd(),
              header_translations,
              main = "Regional association plot", ...)
}
\arguments{
  \item{dataset}{data frame containing the SNPs chromosome number,
    position, p-value and (if a \code{name_cutoff_p} is specified)
    their SNP IDs.}
  \item{chr}{character or numeric value indicating the chromosome
    of interest.}
  \item{start_pos, end_pos}{Numeric; the position values of the
    beginning and end of the region of interest.}
  \item{plot_cutoff_p}{numeric - the threshold of p-values to be
    shown in the QQ & Manhattan plots. Higher (less
    significant) p-values are not included in the plot. Note
    that, unlike in \code{\link{QC_GWAS}} or \code{\link{QC_plots}},
    the default value is \code{1} (i.e. includes everything),
    rather than \code{0.05}. Setting it to \code{0.05} excludes
    95\% of data-points, which significantly reduces running time
    and memory usage in larger datasets.}
  \item{name_cutoff_p}{numeric; SNPs with p-values lower than or
    equal to this value will have their names plotted in the
    graph. If set to \code{NULL} (default), no names will be
    plotted.}
  \item{data_name}{character string; the subtitle for the plot.}
  \item{save_name}{character string; the filename, \emph{without}
    extension, for the graph file.}
  \item{save_dir}{character string; the directory where the
    graph is saved. Note that R uses \emph{forward}
    slash (/) where Windows uses backslash (\\).}
  \item{header_translations}{translation table for column names.
    See \code{\link{translate_header}} for more information. If
    the argument is left empty, \code{dataset} is assumed to use
    the standard column-names of \code{\link{QC_GWAS}}.}
  \item{main, \dots}{arguments passed to \code{\link{plot}}.}
}
\details{
  By default, \code{plot_regional} expects \code{dataset} to use the
  standard column-names used by \code{\link{QC_GWAS}}. A
  translation table can be specified in \code{header_translations}
  to allow non-standard names. See \code{\link{translate_header}}
  for more information.
}
\value{An invisible \code{NULL}.}
\seealso{
  For creating a Manhattan plot: \code{\link{QC_plots}}.
}
\examples{
  data("gwa_sample")

  plot_regional(dataset = gwa_sample,
          chr = 2, start_pos = 55000000, end_pos = 75000000,
          data_name = "QC GWAS sample data",
          save_name = "sample_regional_association")
}
\keyword{distribution}
\keyword{univar}
