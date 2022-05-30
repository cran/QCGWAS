\name{match_alleles}
\alias{match_alleles}
\title{Check and correct alleles in GWAS result files}
\description{
  This function checks the reported alleles and allele frequencies
  in GWAS results data by comparing them to a reference table.
  It will also uniformize the dataset by switching all SNPs to
  the positive strand and flipping the alleles so that the
  effect allele matches the reference minor allele.
}
\usage{
match_alleles(dataset, ref_set, HQ_subset,
              dataname = "dataset", ref_name = "reference",
       unmatched_data = !all(dataset$MARKER \%in\% ref_set$SNP),
              check_strand = FALSE,
              save_mismatches = TRUE, delete_mismatches = FALSE,
              delete_diffEAF = FALSE, threshold_diffEAF = 0.15,
              check_FRQ = TRUE, check_ambiguous = FALSE,
              plot_FRQ = FALSE, plot_intensity = FALSE, 
              plot_if_threshold = FALSE,
              threshold_r = 0.95,
              return_SNPs = FALSE, return_ref_values = FALSE,
              header_translations, header_reference,
              save_name = dataname, save_dir = getwd(),
              use_log = FALSE, log_SNPall = nrow(dataset))
}
\arguments{
  \item{dataset}{
    table containing the allele data. \code{dataset}
    should always contain columns for the SNPID, the effect
    allele and the other allele. Strand and allele frequency
    may be required, depending upon the settings, while
    effect size is optional. \code{match_alleles} accepts
    non-standard column names, provided a translation table is
    specified in \code{header_translations}. The order of columns
    does not matter; nor does the presence of other columns.}
  \item{ref_set}{
    table containing the reference data. \code{ref_set}
    should always contain columns for the SNPID, the minor
    allele and the major allele. Minor allele frequency is only
    required if \code{check_FRQ} or \code{check_ambiguous} are
    \code{TRUE}. The standard column-names are \code{"SNP"},
    \code{"MINOR"}, \code{"MAJOR"} and \code{"MAF"}. Non-standard
    column names are accepted if a translation table is
    specified in \code{header_reference}. All SNPs \emph{must}
    be aligned to the positive strand.}
  \item{HQ_subset}{an \emph{optional} logical or numeric vector
    indicating the rows in \code{dataset} that contain high
    quality SNPs.}
  \item{dataname, ref_name}{character strings; the names of the
    dataset and reference, respectively. Used as identifiers
    in the output files.}
  \item{unmatched_data}{logical; are there SNPs in the dataset
    that do not appear in the reference? This argument is
    currently redundant: the function will determine it
    automatically.}
  \item{check_strand}{logical; should the function check for
    negative-strand SNPs? If \code{FALSE}, all SNPs are assumed
    to be on the positive strand.}
  \item{save_mismatches}{logical; should mismatching entries be
    exported to a .txt file before they are corrected?}
  \item{delete_mismatches}{logical; should mismatching SNPs
    (that could not be corrected by strand-switching) have
    their effect allele set to missing?}
  \item{delete_diffEAF}{logical; should SNPs that exceed the
    \code{threshold_diffEAF} have their effect allele set to
    missing?}
  \item{threshold_diffEAF}{numeric; the max. allowed difference
    between reported and reference allele frequency.}
  \item{check_FRQ}{logical; should the function correlate
    the reported allele-frequency with that of the reference?}
  \item{check_ambiguous}{logical; should the function do separate
    frequency correlations and create separate plots for SNPs
    with a strand-independent allele-pair (i.e. an A/T or C/G
    configuration)?}
  \item{plot_FRQ}{logical; should a scatterplot of reported vs.
    reference allele-frequency be made?}
  \item{plot_intensity}{logical; if \code{TRUE}, instead of a
    scatterplot an intensity plot is generated. This option is
    currently only partially implemented. Leave to \code{FALSE}
    for now.}
  \item{plot_if_threshold}{logical; if \code{TRUE}, the
    scatterplot is only made when frequency correlation is below
    the threshold specified by \code{threshold_r}.}
  \item{threshold_r}{numeric; the correlation threshold value.}
  \item{return_SNPs}{logical; should the return value include
    the relevant columns of \code{dataset}?}
  \item{return_ref_values}{logical; should the return-value
    include the matching entries in \code{ref_set}?}
  \item{header_translations, header_reference}{
    translation tables for converting the column names of
    \code{dataset} and \code{ref_set} to standard names,
    respectively. See \code{\link{translate_header}} for more
    information.}
  \item{save_name}{
    character string; the filename, \emph{without}
    extension, for the various output files.}
  \item{save_dir}{character string; the directory where the
    output files are saved. Note that R uses \emph{forward}
    slash (/) where Windows uses backslash (\\).}
  \item{use_log, log_SNPall}{
    arguments used by \code{\link{QC_GWAS}};
    redundant when \code{match_alleles} is used separately.}
}
\details{
  \code{match_alleles} is one of the more complicated functions
  of \code{QCGWAS}. However, what it does is quite simple:
  \itemize{
    \item Check for incorrect allele pairs
    \item Uniformize the output so that all SNPs are on the
      positive strand, and identical SNPs will have the same
      effect allele and other allele in all datasets
    \item Check the allele frequency
  }
  The complexity stems from the fact that these three tasks have
  to be carried out together and often overlap.
  So the actual function schematic looks like this:
  \itemize{
    \item Switch negative-strand SNPs (i.e. SNPs with \code{"-"}
      in the strand-column) to the positive strand. This step can be
      disabled by setting \code{check-strand} to \code{FALSE}.
    \item Correct (if possible) mismatching alleles. A mismatch
      is when the reported allele-pair does not match that in
      the reference. \code{match_alleles} will attempt to fix
      the mismatch by "strand-switching" the alleles. The
      assumption is that dataset merely reported the wrong
      strand, so converting them to the opposing strand should
      solve the mismatch. If it does not, the SNPs are truly
      mismatches. If \code{delete_mismatches} is \code{TRUE},
      the effect alleles are set to \code{NA}; if \code{FALSE},
      they are restored to their original configuration and
      excluded from the allele-frequency test. If
      \code{save_mismatches} is \code{TRUE}, the entries are
      exported in a .txt before being changed. Note that
      \code{save_mismatches} only exports true mismatches (i.e.
      not those that were fixed after strand-switching).
    \item Align the allele-pairs with the reference. In order to
      have the same effect allele with the same SNP in every
      dataset, SNPs are "flipped" so that the effect allele
      matches the reference minor allele. Flipped alleles will
      also have their allele frequency and effect size inverted.
    \item Check for undetected strand-mismatch by counting the
      number of ambiguous and (if \code{check_FRQ} is \code{TRUE})
      suspect SNPs. Ambiguous SNPs are SNPs with an allele-pair
      that is identical on both strands (i.e. A/T or C/G).
      Suspect SNPs are ambiguous SNPs whose allele-frequency
      differs strongly from that of the reference.
    \item Check allele-frequencies by correlating and/or
      plotting them against the reference. If
      \code{check_ambiguous} is \code{TRUE}, additional
      scatterplots will be made for the subsets of ambiguous and
      non-ambiguous SNPs. If \code{delete_diffEAF} is \code{TRUE},
      SNPs whose allele-frequency differs from the reference by
      more than \code{threshold_diffEAF} have their effect alleles
      set to \code{NA} as well. This entire step can be disabled
      by setting \code{check_FRQ} to \code{FALSE}.
  }  
}
\value{
  An object of class 'list' with the following components:
  
  \item{FRQ_cor, FRQ_cor_ambiguous, FRQ_cor_nonambi}{
    Allele-frequency correlations for all, ambiguous, and
    non-ambiguous SNPs respectively.}
  \item{n_SNPs}{Total number of SNPs in \code{dataset}}
  \item{n_missing, n_missing_data, n_missing_ref}{
    Number of SNPs with missing allele-data in either
    \code{dataset} or \code{reference}, dataset only, and
    reference only, respectively.}
  \item{n_negative_strand, n_negative_switch, n_negative_mismatch}{
    Number of negative-strand SNPs, the subset of negative-strand
    SNPs that were strand-switched twice because they did not match
    the reference, and the subset of double-switched SNPs that
    were still mismatching after the second strand-switch.}
  \item{n_strandswitch, n_mismatch}{
    Number of SNPs that was strand-switched because they did not
    match the reference, and the subset of those that still did
    not match after the strand-switch.}
  \item{n_flipped}{Number of SNPs whose alleles were flipped to
    align them with the reference.}
  \item{n_ambiguous, n_suspect}{Number of ambiguous SNPs, and
    the subset of those that had a large allele-frequency aberration.}
  \item{n_diffEAF}{Number of SNPs whose allele-frequency differs from
    the reference by more than \code{threshold_diffEAF}.}
  \item{MARKER}{When \code{return_SNPs} and/or \code{return_ref_values}
    is \code{TRUE}, this returns the column of \code{dataset}
    containing the SNP IDs. If not, this returns \code{NULL}.}
  \item{EFFECT_ALL, OTHER_ALL, STRAND, EFFECT, EFF_ALL_FREQ}{
    If \code{return_SNPs} is \code{TRUE}, these elements return
    the corrected columns of data-set. If \code{FALSE}, these
    return \code{NULL}. Note: \code{match_alleles} only returns
    those columns that were checked; if \code{check_FRQ} is
    \code{FALSE}, \code{EFF_ALL_FREQ} return \code{NULL}. The
    same goes for \code{check_strand} and \code{STAND}.
    \code{EFFECT} is only returned if present in \code{dataset}.
  }
  \item{ref_MINOR, ref_MAJOR, ref_MAF}{
    If \code{return_ref_values} is \code{TRUE}, these elements
    return the reference minor and major alleles and allele
    frequency column for the SNPs in \code{MARKER}.If
    \code{FALSE}, these return \code{NULL}. \code{ref_MAF} is
    only returned when \code{check_FRQ} is \code{TRUE}.
  }
}
\section{Interpreting the output}{
  The output of \code{match_alleles} may seem a bit overwhelming
  at first, so here is a short explanation of what it means and
  what you should pay attention to.
  
  The columns included in the return value when \code{return_SNPs}
  is \code{TRUE} are the  post-matching dataset. This is only
  relevant if you want to  continue working with the corrected
  dataset. Similarly, the output of \code{return_ref_values} is
  only used for comparing the post-matching dataset to the
  reference.
  
  \code{n_missing}, \code{n_missing_data} and
  \code{n_missing_ref} report the prevalence of missing allele
  data, but are otherwise irrelevant.
  
  The majority of return values serve to check whether
  strand-switching was performed correctly. \code{n_strandswitch}
  indicates how many SNPs were converted to the other strand
  because of a mismatch with the reference. In our experience,
  many cohorts do not include stand data, or simply set all SNPs
  to \code{"+"}, so the presence of strand-switched SNPs
  isn't an indicator of problems by itself. However, if the
  strand-switching did not fix the mismatch, there may a problem.
  The subset of strand-switched SNPs that could not be fixed is
  reported as \code{n_mismatch}, and indicates incorrect allele
  data or, possibly, triallelic SNPs.
  
  Depending on the argument \code{save_mismatches}, mismatching
  entries are exported
  as a .txt file, together with the reference data. This allows
  the user to see which SNPs are affected.
  
  Another sign of trouble is when negative-strand SNPs
  (\code{n_negative_strand}) are present (i.e. the cohort included
  real strand data, rather just setting it to \code{"+"}), but
  strand-switching still occurred. Negative-strand SNPs
  are converted to the positive strand before their alleles are
  compared to the reference, so they should not appear here.
  If they do, it means that either the strand-column data is
  incorrect, or it is an ordinary mismatch (see above).
  
  Negative-strand SNPs that are "strand-switched" will revert
  to their original allele configuration (but the strand column
  now reports them as being on positive strand). The output of
  \code{QC_GWAS} calls them double strand-switches, but here
  they are reported as \code{n_negative_switch}. The subset of
  those that could not be fixed is \code{n_negative_mismatch}.
  
  Just to be clear: the relevant output is still
  \code{n_strandswitch} and \code{n_negative_strand}, not 
  \code{n_negative_switch} and \code{n_negative_mismatch}.
  Whether it was the negative-strand SNPs that were switched or
  not is not important: the important thing is that there were
  negative-strand SNPs (i.e. the cohort included real strand
  data rather setting everything to \code{"+"}); yet
  strand-switches were still necessary and cannot be attributed
  to mismatch (i.e. the strand data is incorrect). 
  
  \code{n_flipped} counts how many SNPs had their alleles
  reversed to match the effect allele with the reference
  minor-allele. This is merely recorded for "administrative"
  purposes, and shouldn't concern the user.
  
  \code{n_ambiguous} and \code{n_suspect} are another test of
  the strand information. Ambiguous SNPs are SNPs that have
  the same allele pair on the positive and negative strands (i.e.
  A/T or C/G). Matching them with the allele-reference therefor
  won't detected incorrect strand-information. In a normal-sized
  GWAS results file, about 15\% of SNPs will be ambiguous.
  
  Suspect SNPs are the subset of ambiguous SNPs whose allele
  frequency is significantly different from that in the reference
  ( < 0.35. vs > 0.65 or visa versa). In our experience,
  a GWAS results file with 2.5M SNPs will have only a few dozen
  suspect SNPs. However, if it's a sizable proportion of all
  ambiguous SNPs, it indicates that the ambiguous SNPs are
  listed for the wrong strand. This will also have resulted
  in the wrong SNPs being flipped in the previous step, so
  it should be visible in the allele-frequency correlation test
  as well.
  
  \code{n_diffEAF} counts SNPs with significantly different
  allele-frequencies. A large number here indicates either
  that the allele-frequencies are incorrect or listed for the
  wrong allele (see below), or that the population used in the
  dataset does not match that of the reference.
  
  The \code{FRQ_cor} value is the correlation between the
  reported and reference allele-frequencies. If allele frequency
  is correct, the correlation should be near \code{1}. If it's
  close to \code{-1}, the listed frequency is that of the other
  (i.e. non-effect) allele.
  
  the \code{FRQ_cor_ambiguous} and \code{FRQ_cor_nonambi} values
  are the same test for the subsets of ambiguous and
  non-ambiguous SNPs. If ambiguous SNPs are listed on the wrong
  strand, then they will have been flipped incorrectly, so
  allele-frequency correlation should also move towards \code{-1}.
}
\note{
  The function does not delete SNPs, regardless of the
  \code{delete_mismatches} \code{delete_diffEAF} arguments.
  Setting these to \code{TRUE} means that any such
  SNPs are marked by having their effect allele set to \code{NA}.
  The actual deletion takes place inside \code{\link{QC_GWAS}}.
}
\seealso{
  \code{\link{create_hapmap_reference}} for creating an allele
  reference from publicly-available HapMap data
  
  \code{\link{switch_strand}}
}
\examples{
# In order to keep the QCGWAS package small, no allele reference
# is included. Use the create_hapmap_reference function (see
# above) to create one.

\dontrun{
  data("gwa_sample")
  hapmap_ref <- read.table("C:/new_hapmap/new_hapmap.txt",
                  header = TRUE, stringsAsFactors = FALSE)

  match_alleles(gwa_sample, hapmap_ref,
    dataname = "sample data", ref_name = "HapMap",
    save_name = "test_allele1", save_dir = "C:/new_hapmap",
    check_strand = TRUE, plot_FRQ = TRUE)

 HQ_SNPs <- HQ_filter(data = gwa_sample, filter_NA = TRUE,
                          filter_FRQ = 0.01, filter_cal = 0.95)
  match_alleles(gwa_sample, hapmap_ref,
    HQ_subset = HQ_SNPs,
    dataname = "sample data", ref_name = "HapMap",
    save_name = "test_allele2", save_dir = "C:/new_hapmap",
    check_strand = TRUE, plot_FRQ = TRUE)

  match_output <-
    match_alleles(gwa_sample, hapmap_ref,
      HQ_subset = HQ_SNPs,
      delete_mismatches = TRUE, return_SNPs = TRUE,
      delete_diffEAF = TRUE, threshold_diffEAF = 0.15,
      dataname = "sample data", ref_name = "HapMap",
      save_name = "test_allele3", save_dir = "C:/new_hapmap",
      check_strand = TRUE, plot_FRQ = TRUE)
  
  if(sum(match_output$n_negative_strand,
         match_output$n_strandswitch, match_output$n_mismatch,
         match_output$n_flipped, match_output$n_diffEAF) > 0){
    gwa_sample$EFFECT_ALL   <- match_output$EFFECT_ALL
    gwa_sample$OTHER_ALL    <- match_output$OTHER_ALL
    gwa_sample$STRAND       <- match_output$STRAND
    gwa_sample$EFFECT       <- match_output$EFFECT
    gwa_sample$EFF_ALL_FREQ <- match_output$EFF_ALL_FREQ
  }
  }
}
\keyword{distribution}
\keyword{univar}
