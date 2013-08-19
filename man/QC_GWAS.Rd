\name{QC_GWAS}
\alias{QC_GWAS}
\alias{QC_series}
\title{Automated Quality Control of GWAS Results files}
\description{
  \code{QC_GWAS} runs a full quality control (QC) over a single
  GWAS results file. It removes missing and invalid data, checks
  the alleles and allele frequency with a reference, tests the
  reported p-value against both calculated and expected values,
  creates QQ and Manhattan plots and reports the distribution of
  the quality-parameters within the dataset, as well as various
  QC statistics.
  
  \code{QC_series} does the same thing for multiple GWAS results
  files. It's mainly a wrapper that passes individual files to
  \code{QC_GWAS}, but it has a few extra features, such as
  making a checklist of important QC stats
  and creating several graphs to compare the QC'ed files.
  
  Although the number of arguments in \code{QC_GWAS} may seem
  overwhelming, only three of them are required to run a basic
  QC. The name of the file to be QC'ed should be passed to the
  \code{filename} argument; the directory of said file to the
  \code{dir_data} argument; and a header-translation table to
  the \code{header_translations} argument.
  For a quick introduction to QCGWAS, read the quick
  reference guide that can be found in "R\\library\\QCGWAS\\doc".
}
\usage{
QC_GWAS(filename,
        filename_output = paste("QC_", filename, sep = ""),
        dir_data = getwd(),
        dir_output = paste(dir_data, "QCGWASed", sep = "/"),
        dir_references = dir_data,
        header_translations,
        column_separators = c("\t", " ", "", ",", ";"),
        nrows = -1, nrows_test = 1000,
        header = TRUE, comment.char = "",
        na.strings = c("NA", "nan", "NaN", "."),
        imputed_T = c("1", "TRUE", "T"),
        imputed_F = c("0", "FALSE", "F"),
        imputed_NA = c(NA, "-"),
        save_final_dataset = TRUE, order_columns = FALSE,
        out_header = "standard",
        out_quote = FALSE, out_sep = "\t", out_eol = "\n",
        out_na = "NA", out_dec = ".", out_qmethod = "escape",
        out_rownames = FALSE, out_colnames = TRUE,
        return_HQ_effectsizes = FALSE,
        remove_X = FALSE, remove_Y = FALSE,
        remove_XY = remove_Y, remove_M = FALSE,
        calculate_missing_p = FALSE,
        make_plots = TRUE, only_plot_if_threshold = TRUE,
        threshold_allele_freq_correlation = 0.95,
        threshold_p_correlation = 0.99,
        plot_histograms = make_plots, plot_QQ = make_plots,
        plot_QQ_bands = TRUE, plot_Manhattan = make_plots,
        plot_cutoff_p = 0.05,
        allele_ref_std, allele_name_std,
        allele_ref_alt, allele_name_alt,
        update_alt = FALSE, update_savename,
        update_as_rdata = FALSE, backup_alt = FALSE,
        remove_mismatches = FALSE,
        remove_mismatches_std = remove_mismatches,
        remove_mismatches_alt = remove_mismatches,
        threshold_diffEAF = 0.15, remove_diffEAF = FALSE,
        remove_diffEAF_std = remove_diffEAF,
        remove_diffEAF_alt = remove_diffEAF,
        check_ambiguous_alleles = TRUE,
        use_threshold = 0.1,
        useFRQ_threshold = use_threshold,
        useHWE_threshold = use_threshold,
        useCal_threshold = use_threshold,
        useImp_threshold = use_threshold,
        useMan_threshold = use_threshold,
        HQfilter_FRQ = NULL, HQfilter_HWE = NULL,
        HQfilter_cal = NULL, HQfilter_imp = NULL,
        QQfilter_FRQ = NULL, QQfilter_HWE = NULL,
        QQfilter_cal = NULL, QQfilter_imp = NULL,
        NAfilter = TRUE,
        NAfilter_FRQ = NAfilter, NAfilter_HWE = NAfilter,
        NAfilter_cal = NAfilter, NAfilter_imp = NAfilter,
        ignore_impstatus = FALSE,
        minimal_impQ_value = -0.5, maximal_impQ_value = 1.5,
        logI = 1L, logN = 1L, ...)
QC_series(data_files, datafile_tag, output_filenames,
          dir_data = getwd(),
          dir_output = paste(dir_data, "QCGWASed", sep = "/"),
          dir_references = dir_data,
          header_translations, out_header = "standard",
          allele_ref_std, allele_name_std,
          allele_ref_alt, allele_name_alt,
          update_alt = FALSE, update_savename,
          update_as_rdata = FALSE, backup_alt = FALSE,
          plot_effectsizes = TRUE, lim_effectsizes = NULL,
          plot_SE = TRUE, label_SE = TRUE,
          plot_SK = TRUE, label_SK = "outliers",
          save_filtersettings = FALSE, ...)
}
\arguments{
  \item{filename, data_files, datafile_tag}{
    \code{filename} and \code{data_files} are, respectively, the
    name and names of the GWAS results file(s) to be QC'ed. If
    no \code{data_files} are specified, \code{QC_series} will
    process all filed in \code{dir_data} containing the string
    passed to \code{datafiles_tag} in their filename. See
    below for more information on the input requirements.}
  \item{filename_output, output_filenames}{
    respectively the filename or names of the output of the QC.
    This should not include an extension, since the QC will
    automatically add one. The default is to use the input
    filename with \code{"QC_"} prefixed.}
  \item{dir_data, dir_output, dir_references}{character strings
    specifying the directory dress of the folders for, respectively,
    the input file(s), the output file(s) and the auxillary files
    (header-translation tables and allele references). Note that
    R uses \emph{forward} slash (/) where Windows uses backslash (\\).
    If \code{dir_output} does not exist, it will be created.
    If no \code{dir_output} is specified, a folder named
    \code{"QCGWASed"} will be created in \code{dir_data}.}
  \item{header_translations}{Translation table for the column 
    names of the \emph{input} file. Alternatively, the name of a file
    in \code{dir_references} containing such a table. See
    \code{\link{translate_header}} for details.}
  \item{column_separators}{character string or vector; specifies
    the values used as column delimitator in the GWAS file. The
    argument is passed to \code{\link{load_test}}; see the
    description of that function for more information.}
  \item{nrows_test}{integer; the number of rows used for
    "trial-loading". Before loading the entire dataset, the
    function \code{\link{load_test}} is called to determine the
    dataset's file-format by reading the top \code{x} lines, where
    \code{x} is \code{nrows_test}. Setting \code{nrows_test} to a low number
    (e.g. \code{150}) means quick testing, but runs the risk of
    missing problems in lower rows. To test the entire dataset,
    set it to \code{-1}.}
  \item{nrows, header, comment.char}{arguments passed to
    \code{\link{read.table}} when importing the dataset.}
  \item{na.strings}{character vector describing the values that
    represent missing data in the dataset. Passed to
    \code{\link{read.table}}.}
  \item{imputed_T, imputed_F, imputed_NA}{character vectors;
    passed to \code{\link{convert_impstatus}} (as \code{T_strings},
    \code{F_strings} and \code{NA_strings}, respectively) to
    translate the imputation-status column. Note that the
    current version of \code{QC_GWAS} \emph{always} translates
    the imputation status. Even when the dataset already has
    the correct format, the user still needs to specify \code{1},
    \code{0} and \code{NA} for these arguments, respectively.
    Also note that R distinguishes between the value \code{NA}
    and the character string \code{"NA"}.}
  \item{save_final_dataset}{logical; should the post-QC dataset
    be saved?}
  \item{order_columns}{logical; should the post-QC dataset use
    the default column order?}
  \item{out_header}{Translation table for the column names of
    the \emph{output} file. This argument is the opposite of
    \code{header_translations}: it translates the standard
    column-names of \code{QC_GWAS} to user-defined ones.
    \code{output_header} can be one of three things:
    \itemize{
      \item A user specified table similar to the one used by
        \code{\link{translate_header}}. However, as this
        translates standard names into non-standard ones, the
        standard names should be in the right column, and the
        desired ones in the left. There is also no requirement
        for the names in the \emph{left} column to be uppercase.
      \item The name of a file in \code{dir_references} containing
        such a column.
      \item Character string specifying a standard form. See
        the section 'Output header' below for the options.
    }}
  \item{out_quote, out_sep, out_eol, out_na, out_dec, out_qmethod,
    out_rownames, out_colnames}{arguments passed to
    \code{\link{write.table}} when saving the final dataset.}
  \item{return_HQ_effectsizes}{logical; return a vector of (max.
    1000) high-quality effect-sizes? (In \code{QC_series}, this
    is set by the \code{plot_effectsizes} argument.)}
  \item{remove_X, remove_Y, remove_XY, remove_M}{logical;
    respectively whether X-chromosome, Y-chromosome,
    pseudo-autosomal and mitochondrial SNPs are removed from the
    dataset.}
  \item{calculate_missing_p}{logical; should the QC calculate
    missing/invalid p-values in the dataset?}
  \item{make_plots}{logical; should the QC generate and save
    QQ plots, a Manhattan plot, histograms of data distribution
    and scatter plots of correlation in \code{dir_reference}?}
  \item{only_plot_if_threshold}{logical; should the scatterplots
    only be made if the correlation is below a threshold value?}
  \item{threshold_allele_freq_correlation, threshold_p_correlation}{
    numeric; thresholds for reporting and plotting
    the correlation between respectively
    the allele frequency of the dataset and the reference, and
    the calculated and reported p-values.}
  \item{plot_histograms}{logical; should histograms of the effect
    sizes, standard errors, allele frequencies, HWE p-values,
    callrates and imputation quality be made?}
  \item{plot_QQ, plot_Manhattan}{logical; should QQ and Manhattan
    plots be made?}
  \item{plot_QQ_bands}{logical; include probability bands in the
    QQ plot?}
  \item{plot_cutoff_p}{numeric; significance threshold for
    inclusion in the QQ and Manhattan plots. The default value
    (\code{0.05}) excludes 95\% of SNPs, significantly reducing
    running-time and memory usage. For this reason it is not
    recommend to set a higher value when QC'ing a normal-sized
    GWAS dataset.}
  \item{allele_ref_std, allele_ref_alt}{the standard and alternative
    allele-reference tables. Alternatively, the name of a file
    in \code{dir_references} containing said table. Files in
    \code{.RData} format are accepted, but the table's object
    name must match the argument name. See \code{\link{match_alleles}}
    for more information on the input requirements.}
  \item{allele_name_std, allele_name_alt}{character strings;
    these name the standard and alternative allele reference in
    the output. If no values are given, the function will use the
    reference's filename (if specified) or a default name.}
  \item{update_alt}{logical; if the function encounters SNPs not
    included in either the standard or alternative reference,
    should these be added to the alternative reference file? If
    no alternative reference was specified, this creates one.}
  \item{update_savename}{character string; the filename for
    saving the updated alternative reference, \emph{without}
    extension. If \code{allele_ref_alt} is a filename, it is
    not necessary to specify this argument.}
  \item{update_as_rdata}{logical; should the updated alternative
    allele reference be saved as \code{.RData} (\code{TRUE}) or
    a tab-delimitated .txt file (\code{FALSE}).}
  \item{backup_alt}{logical; if the alternative allele reference is updated,
    should a back-up be made of the original reference file?}
  \item{remove_mismatches, remove_mismatches_std,
    remove_mismatches_alt}{
    logical; should SNPs with mismatching alleles be removed
    from the dataset? \code{remove_mismatches} serves as the
    default value; the other two arguments determine this
    setting for the standard and alternative references,
    respectively.}
  \item{threshold_diffEAF}{
    Numeric; the threshold for the difference between reported
    and reference allele-frequency. SNPs for which the
    difference exceeds the threshold are counted and
    (optionally) removed.}
  \item{remove_diffEAF, remove_diffEAF_std, remove_diffEAF_alt}{
    Logical; should SNPs that exceed the \code{threshold_diffEAF}
    be removed from the dataset? \code{remove_diffEAF} serves as
    the default value; the other two arguments determine this
    setting for the standard and alternative references,
    respectively.}
  \item{check_ambiguous_alleles}{logical; check for SNPs with
    strand-independent allele-configurations (i.e. A/T and C/G
    SNPs)?}
  \item{use_threshold}{
    numeric; threshold value. The relative or absolute number of
    usable values required for a variable to be used in the QC.
    These arguments prevent the QC from applying filters to variables
    with no data. If a variable has less non-missing, non-invalid
    values than specified in the threshold, it will be ignored;
    i.e. no filter will be applied to it and no plots will be made.
    Values \code{> 1} specify the absolute threshold, while
    values of \code{1} or lower specify the fraction of SNPs
    remaining in the dataset. This argument is the
    default threshold for all variables; variable-specific
    thresholds can be set with the following arguments.}
  \item{useFRQ_threshold, useHWE_threshold, useCal_threshold,
    useImp_threshold, useMan_threshold}{
    numeric; variable-specific thresholds for allele frequency,
    HWE p-value, callrate, imputation quality and Manhattan plot
    (i.e. chromosome & position values) respectively.}
  \item{HQfilter_FRQ, HQfilter_HWE, HQfilter_cal, HQfilter_imp}{
    numeric; threshold values for the high-quality SNP filter.
    SNPs that do not meet or exceed all four thresholds will be
    excluded from several QC tests. The filters are for allele
    frequency, HWE p-value, callrate & imputation quality,
    respectively, and are processed by \code{\link{HQ_filter}}.
    See 'Filter arguments' for more information. Note: the high-quality
    filter does not remove SNPs; it merely excludes them from
    several QC tests.}
  \item{QQfilter_FRQ, QQfilter_HWE, QQfilter_cal, QQfilter_imp}{
    numeric vector; threshold values for the QQ plot filters.
    SNPs that do not meet or exceed the value will be excluded
    from the QQ plot. Up to five values can be specified per
    filter. The filters are for allele-frequency, HWE p-value,
    callrate & imputation quality respectively, and are
    processed by \code{\link{QC_plots}}. See 'Filter arguments' for more
    information.}
  \item{NAfilter, NAfilter_FRQ, NAfilter_HWE, NAfilter_cal, NAfilter_imp}{
    logical; should the high-quality and QQ filters exclude
    (\code{TRUE}) or ignore (\code{FALSE}) missing values?
    \code{NAfilter} is the default setting; the others allow
    allow variable specific settings.}
  \item{ignore_impstatus}{logical; if \code{FALSE}, HWE p-value
    and callrate filters are applied only to genotyped SNPs, and
    imputation quality filters only to imputed SNPs. If
    \code{TRUE}, the filters are applied to all SNPs regardless
    of the imputation status.}
  \item{minimal_impQ_value, maximal_impQ_value}{
    numeric; the minimal and maximal possible (i.e. non-invalid)
    imputation quality values.}
  \item{logI, logN}{progress indicators used by \code{QC_series}:
    irrelevant for users.}
  \item{plot_effectsizes, plot_SE, plot_SK}{logical; additional
    plot options for \code{QC_series}. The arguments toggle the
    creation of plots comparing the effect-size distribution,
    precision and skewness vs. kurtosis of all successfully QC'ed
    datasets, respectively. See \code{\link{plot_distribution}},
    \code{\link{plot_precision}} and \code{\link{plot_skewness}}
    for more information.}
  \item{lim_effectsizes}{specifies the y-axis range of the
    effect-size distribution plot.}
  \item{label_SE}{logical; should the datapoints in the precision
    plot be labeled?}
  \item{label_SK}{character string; determines whether the
    datapoints in the skewness vs. kurtosis plot are labeled.
    Options are \code{"none"}, \code{"all"}, or \code{"outliers"}
    (outliers only).}
  \item{save_filtersettings}{logical; saves the filtersettings
    used by the high-quality filter to a file
    'Check_filtersettings.txt' in the output directory. If
    a file of that name already exists, the settings are added
    to the end (i.e. it updates rather than overwrites existing
    files). The file can be used as ini file by
    \code{\link{filter_GWAS}}.}
  \item{\dots}{in \code{QC_series}: arguments passed to
    \code{QC_GWAS}; in \code{QC_GWAS}: arguments passed to
    \code{\link{read.table}} when importing the dataset.}
}
\details{
  The full quality-control carried out by \code{QC_GWAS} consists
  of 5 phases. The function takes a single dataset (or, rather,
  the location and filename of a single dataset) and runs it
  through the following phases:
  \itemize{
    \item 1: Importing the dataset
    \item 2: Checking data integrity
    \item 3: Checking alleles
    \item 4: Generating QC statistics and graphs
    \item 5: Saving the output
  }
  
  \emph{Phase 1: importing the dataset}
  
  GWAS results files come in a variety of formats, so
  \code{QC_GWAS} is flexible about loading data. It
  uses an autoloader function (\code{\link{load_GWAS}}) to
  unpack \code{.zip} or \code{.gz} files and determine the
  column-separator used in the file. See the section
  'Requirement for the input dataset' for more information.
  
  Next, the function attempts to translate the dataset's column
  names (the header) to standard names, so that it knows what
  type of data a column contains. This is done by comparing
  the column names to a translation table (specified in the
  \code{header_translations} argument). See
  \code{\link{translate_header}} for more information.
  
  Note that only the SNP ID, alleles, effect-size and standard
  error columns are required. The absence of other standard
  columns (chromosome, position, strand, allele frequency, HWE
  p-value, callrate, imputation quality, imputation status and
  used for imputation) will not cause the QC to abort.
  Instead, a warning is printed on screen and in the log file,
  and a dummy column filled with \code{NA} values is added to
  the dataset.
  
  It is therefor important to check the log file: if a standard
  column is present but not identified (because it is missing or
  misspelled in the translation table) the QC will continue,
  but is unable to check/use the data inside. The unidentified
  column will be reported in the \code{columns_unidentified}
  value of the \code{QC_GWAS} return or in the
  "QC_checklist.txt" file generated by \code{QC_series}.
  
  \emph{Phase 2: checking data integrity}
  
  The purpose of phase 2 is to ensure that the dataset \emph{can}
  be QC'ed: that that all SNPs have the required data and that
  all columns contain only valid (or missing) values.
  
  The first step is to remove SNPs that won't be used: monomorphic
  SNPs and (if specified by the arguments) allosomal, 
  pseudo-autosomal and mitochondrial SNPs. The function considers
  SNPs monomorphic if they have a missing or invalid other
  (non-effect) allele, identical alleles or an allele-frequency
  of \code{1} or \code{0}.
  
  The second step is to check the imputation status column with
  the function \code{\link{convert_impstatus}}. See the section
  'Requirement for the input dataset' for more information.  
  Imputation status is one of the most important variables in
  the dataset: if unknown, the HWE p-value, callrate and imputation
  quality won't be used (unless \code{ignore_impstatus} is
  \code{TRUE}), as the function cannot determine which
  SNPs are imputed and which are not. For this reason, if
  \code{convert_impstatus} is unable to translate any character 
  string in the column, the QC will abort.
  
  The third step carries out three tests for all other standard
  variables:
  \itemize{
    \item Does the column contain the correct type of data
      (integer, numeric or character)?
    \item How many values are missing (\code{NA})?
    \item How many values are invalid (= impossible)?
  }
  The exact nature of the three tests differs per variable: see
  the documentation file in "R\\library\\QCGWAS\\doc" for more
  detail.
  
  The presence of the wrong data-type will cause the QC to abort.
  Wrong data-type indicates either a problem in the file
  itself, or with the way it was imported (in which case it is
  most likely due to a mistranslated header).
  
  The final step is the removal of the invalid values and of
  unusable SNPs. The variables MARKER, EFFECT_ALLELE, OTHER_ALLELE,
  EFFECT and STDERR are considered crucial. SNPs with missing or
  invalid values in any of these variables are removed the dataset.
  Missing values in the other variables are ignored, while
  invalid values are set to \code{NA}.
  
  \emph{Phase 3: checking alleles}
  
  This phase has three functions:
  \itemize{
    \item To check if the correct alleles are reported for each SNP
    \item To check if the allele-frequency is reported for the
      correct (effect) allele
    \item To ensure that SNPs are aligned to the positive strand
      and use the same effect-allele in all post-QC datasets
  }
  
  This is achieved by comparing the data to a reference, using
  the function \code{\link{match_alleles}}. First, all SNPs are
  switched to the positive strand (the alleles are converted to
  their opposing base and the strand-value is set to \code{"+"}).
  If there are SNPs whose allele pair doesn't match the
  reference, \code{match_alleles} assumes the information in
  the strand column is absent or incorrect, and will also
  switch those SNPs to the other (presumably positive) strand.
  This step is referred to as strand-switching in QC output, and
  is independent from the negative-strand SNP conversion. It is
  therefor possible that a SNP is switched twice: once because
  the strand-column indicates it is on the negative strand, and
  twice because of a mismatch. This is referred to as double
  strand-switch in the output, and indicates either the wrong
  value in the strand column, or a mismatch with the reference.
  In the latter case, it will most likely be picked up in the
  next step.
  
  If the strand-switch does not fix the mismatch, the SNPs
  are counted in the QC output as mismatches. Depending on the
  \code{remove_mismatches} arguments, the SNPs will either be
  removed from the dataset, or left in but excluded from the
  further tests of the allele-matching.
  
  Next, \code{match_alleles} "flips" SNPs so that their effect
  allele matches the reference minor allele. This ensures
  that a SNP will have the same effect allele in all post-QC
  datasets.
  
  \code{match_alleles} also counts the number of SNPs with a
  strand-independent allele configuration (A/T or C/G; these are
  designated as "ambiguous SNPs"), and the subset of those with
  an allele-frequency that is substantially different from the
  reference ("suspect SNPs"). If a substantial proportion of
  ambiguous SNPs is suspect, it indicates that the strand
  information is incorrect. In our experience, a regular, 2.5M
  SNP dataset usually consists of 15\% ambiguous SNPs, of which
  a few dozen will be suspect.
  
  The function also counts the number of SNPs whose allele
  frequency differs from the reference by more than a set amount
  (\code{threshold_difEAF}). If the relevant \code{remove_diffEAF}
  argument is \code{TRUE}, these SNP will be excluded from the
  dataset after the allele-matching is finished.
  
  The final step is to correlate the reported allele-frequency
  against the reference. If allele-frequency is reported for
  the correct (effect) allele, the correlation should be close
  to \code{1}. If the outcome is close to \code{-1}, the reported
  frequency is that of the other allele. Depending on the plot settings,
  a scatter plot of reported vs. reference frequency is made
  for all SNPs, and for the subsets of ambiguous and non-ambiguous
  SNPs.
  
  \emph{The standard and alternative allele reference}
  
  The above steps describe what happens when the dataset is
  compared to a single reference. However, we found that many
  GWAS datasets contain SNPs not present in our standard HapMap
  reference, so we added a second, flexible reference that can be
  updated with any unknown SNPs the QC encounters.
  
  SNPs that are not found in either reference are converted
  to the positive strand, and "flipped" if their allele frequency
  is > 0.50. If \code{update_alt} is \code{TRUE}, these SNPs are
  then added to the alternative
  reference and saved under the name \code{update_savename}.
  There are a few caveats to this system: see the section
  'Updating the alternative reference' for details.
  
  \emph{Phase 4: generating QC statistics and graphs}
  
  At this stage, no further changes will be made to the dataset
  (except, optionally, to recalculate missing p-values).
  The function will now start to calculate the QC statistics and
  generate the important graphs. These are:
  \itemize{
    \item Create histograms of variable distribution (optional)
    \item Check p-values by correlating them to a p calculated
      from the effect-size and standard-error (via the
      \code{\link{check_P}} function).
    \item Recalculate missing/invalid p-values (optional)
    \item Calculate QC statistics
    \item Create QQ and Manhattan plots (optional, see
      \code{\link{QC_plots}} function for more information).
  }
  
  \emph{Phase 5: saving the output}

  A series of tables is added to the bottom of
  the log file, reporting the QC statistics and the data
  distribution. If \code{save_final_dataset} is \code{TRUE}, the
  post-QC data will be exported as a .txt file. The column names
  and format of that file can be specified by the out arguments.
}
\value{
  The most important output of the QC is the log file. See the
  section 'QC output files' for more details. This section only
  describes the function return within R.
  
  \code{QC_series} returns a single, invisible, logical value,
  indicating whether the alternative allele-reference has
  been updated.
  
  \code{QC_GWAS} returns an object of class 'list'. If the QC
  was not successful, this list contains only five of the following
  components (\code{QC_successful}, \code{filename_input},
  \code{filename_output}, \code{all_ref_changed}, \code{effectsize_return}).
  If it was, it will contain all of them:
  
  \item{QC_successful}{logical; indicates whether the QC was
    completed. If \code{FALSE}, the function was either unable
    to load the dataset, encountered an unexpected datatype, or
    removed all SNPs during the QC. The log file and screen
    output will indicate what triggered the abort.}
  \item{filename_input, filename_output}{the filenames of the
    dataset pre- and post-QC respectively.}
  \item{sample_size, sample_size_HQ}{the highest reported sample
    size for all SNPs and high-quality SNPs only, respectively.}
  \item{lambda, lambda_geno, lambda_imp}{the lambda values of
    all, genotyped and imputed SNPs, respectively.}
  \item{SNP_N_input}{the number of SNPs in the original dataset.}
  \item{SNP_N_input_monomorphic}{the number of SNPs removed
    because they are monomorphic.}
  \item{SNP_N_input_monomorphic_identic_alleles}{the subset of
    above that had identical alleles, but allele-frequencies
    that were not \code{0} or \code{1}.}
  \item{SNP_N_input_chr}{the number of SNPs removed because they
    were X-chromosomal, Y-chromosomal, pseudo-autosomal or
    mitochondrial (depends on the remove-arguments). If all
    remove arguments were set to \code{FALSE}, this returns \code{NA}.}
  \item{SNP_N_preQC}{the number of SNPs that entered phase 2b
    (i.e. after removal of the monomorphic and excluded-chromosome
    SNPs).}
  \item{SNP_N_preQC_unusable}{the number of SNPs removed in phase
    2d, due to missing or invalid crucial variables.}
  \item{SNP_N_preQC_invalid}{the number of SNPs with invalid,
    non-crucial values in phase 2d.}
  \item{SNP_N_preQC_min}{the number of negative-strand SNPs in
    phase 2d.}
  \item{SNP_N_midQC}{the number of SNPs in the dataset during
    allele-matching (phase 3).}
  \item{SNP_N_midQC_min}{the number of negative strand SNPs in
    phase 3.}
  \item{SNP_N_midQC_min_std, SNP_N_midQC_min_alt,SNP_N_midQC_min_new}{
    the number of negative strand SNPs matched against,
    respectively, the standard allele reference, the alternative
    allele reference or neither.}
  \item{SNP_N_midQC_strandswitch_std, SNP_N_midQC_strandswitch_alt}{
    the number of SNPs that were strand-switched because of
    a mismatch with the standard and alternative allele reference,
    respectively.}
  \item{SNP_N_midQC_strandswitch_std_min, SNP_N_midQC_strandswitch_alt_min}{
    the subset of previous that were negative-strand SNPs. NOTE:
    at this point in the QC, negative-strand SNPs have already
    been converted to the positive strand, i.e. they should
    \emph{not} appear in this category. If they do, there is a
    problem with the reported strand, or with the reference table.}
  \item{SNP_N_midQC_mismatch}{the number of SNPs that were still
    mismatching after the strand-switch.}
  \item{SNP_N_midQC_mismatch_std, SNP_N_midQC_mismatch_alt}{
    subset of previous that were matched with the standard and
    alternative allele reference, respectively.}
  \item{SNP_N_midQC_mismatch_std_min, SNP_N_midQC_mismatch_alt_min}{
    subset of previous that are negative-strand SNPs.}
  \item{SNP_N_midQC_flip_std, SNP_N_midQC_flip_alt, SNP_N_midQC_flip_new}{
    Number of SNPs that were flipped (had their alleles reversed)
    when matched against, respectively, the standard allele reference,
    the alternative allele reference or neither.}
  \item{SNP_N_midQC_ambiguous}{the number of ambiguous SNPs}
  \item{SNP_N_midQC_ambiguous_std, SNP_N_midQC_ambiguous_alt,
    SNP_N_midQC_ambiguous_new}{subset of ambiguous SNPs matched
    against, respectively, the standard allele reference,
    the alternative allele reference or neither.}
  \item{SNP_N_midQC_suspect}{the subset of ambiguous SNPs whose
    allele frequencies differ strongly from those in the reference.}
  \item{SNP_N_midQC_suspect_std, SNP_N_midQC_suspect_alt}{
    the subsets of previous matched against the standard and
    alternative allele reference, respectively.}
  \item{SNP_N_midQC_diffEAF}{the number of SNPs whose allele
    frequency differs strongly from the reference.}
  \item{SNP_N_midQC_diffEAF_std, SNP_N_midQC_diffEAF_alt}{
    subset of previous that were matched with the standard and
    alternative allele reference, respectively.}
  \item{SNP_N_postQC}{the number of SNPs in the final dataset.}
  \item{SNP_N_postQC_geno, SNP_N_postQC_imp}{the number of
    genotyped and imputed SNPs in the final dataset.}
  \item{SNP_N_postQC_invalid}{the number of SNPs with invalid
    values remaining in the final dataset. Note: any invalid
    values have already been changed to \code{NA} at this point.
    This merely counts how many of those SNPs are still in the
    dataset.}
  \item{SNP_N_postQC_min}{the number of negative-strand SNPs
    in the final dataset. Note: all SNPs have been switched to
    the positive strand at this point. This merely counts how
    many of those SNPs are still in the dataset.}
  \item{SNP_N_postQC_HQ}{the number of high-quality SNPs in the
    final dataset.}
  \item{fixed_HWE, fixed_callrate, fixed_sampleN, fixed_impQ}{
    logical or character string; are HWE p-values, callrates,
    sample-size and imputation quality values identical for all
    relevant SNPs? If \code{TRUE}, it indicates that the parameters
    have not been calculated and are dummy values. If a parameter
    fails the threshold test, this returns \code{"insuf. data"}
    (\code{"no data"} for sample size).}
  \item{effect_25, effect_median, effect_75}{the quartile values
    of the effect-size distribution.}
  \item{effect_mean}{the mean of the effect-size distribution.}
  \item{SE_median, SE_median_HQ}{the median standard error of
    all SNPs and high-quality SNPs only, respectively.}
  \item{skewness, kurtosis}{the skewness and kurtosis value of
    the dataset.}
  \item{skewness_HQ, kurtosis_HQ}{the skewness and kurtosis value
    for high-quality SNPs only.}
  \item{all_ref_std_name, all_ref_alt_name}{the names used for
    the standard and alternative allele-references in the output.}
  \item{all_MAF_std_r, all_MAF_alt_r}{allele-frequency
    correlation with the standard and alternative allele references.}
  \item{all_ambiguous_MAF_std_r, all_ambiguous_MAF_alt_r}{
    allele-frequency correlations for the subset of ambiguous
    SNPs in the standard and alternative allele references,
    respectively.}
  \item{all_non_ambig_MAF_std_r, all_non_ambig_MAF_alt_r}{
    allele-frequency correlations for the subset of non-ambiguous
    SNPs in the standard and alternative allele references,
    respectively.}
  \item{all_ref_changed}{logical; has an updated alternative
    allele reference been saved?}
  \item{effectsize_return}{logical; is a vector of high-quality
    effect-sizes returned in \code{effectsizes_HQ}?}
  \item{effectsizes_HQ}{if \code{effectsize_return} is \code{TRUE},
    a vector of 1000 high-quality effect-sizes; if not, \code{NULL}.
    If a dataset contains less than 1000 high-quality SNPs, the
    vector is padded with NA's to bring it to 1000 values.}
  \item{pvalue_r}{the correlation between reported and calculated
    p-values.}
  \item{visschers_stat, visschers_stat_HQ}{the Visscher's statistic
    for all SNPs and high-quality SNPs only, respectively.}
  \item{columns_std_missing}{the names of any missing, standard
    columns: if no columns are missing, this returns \code{0}.}
  \item{columns_std_empty}{the names of any empty, standard
    columns: if no columns are empty, this returns \code{0}.}
  \item{columns_unidentified}{the names of any unidentified
    columns in the input dataset. If none, this returns \code{0}.}
  \item{outcome_useFRQ, outcome_useHWE, outcome_useCal,
  outcome_useImp, outcome_useMan}{logical; indicates whether
    the variable passed the threshold test.}
  \item{\dots}{the remaining 'setting' components return
    the the actual filter settings used in the QC (i.e.
    taking into account whether the variables passed the
    threshold test).}
}
\section{QC output files}{
  The results of the QC are reported in a variety of .txt and
  .jpeg files saved in \code{dir_output}. The files use the same
  output name as the dataset, with an extension to indicate
  their contents (i.e. '_log.txt', '_graph_QQ.jpg'). The .txt
  files are tab-delimited and are best viewed in a spreadsheet
  program. The most important one is the log file. This file
  summarizes the results of the QC and the data inside the file.
  
  \emph{The log file}
  
  The top of the file is table of log entries, reporting any
  (potential) problems encountered during the QC. Some of these
  are just routine updates; the removal of SNPs with missing
  data, for example. However, do check the other entries.
  These report important but non-fatal problems, relating to
  crucial missing data or invalid data. In such a case, and
  provided the QC did not abort, the affected SNPs will have
  been exported to a .txt file before being excluded, so the
  user can inspect them without having to reload the entire
  dataset. The .txt's are:
  
  \code{[filename_output]_SNPs_invalid_allele2.txt}
  
  \code{[filename_output]_SNPs_duplicates.txt}
  
  \code{[filename_output]_SNPs_removed.txt}
  
	\code{[filename_output]_SNPs_improbable_values.txt}
	
  (Note: the names of the files are slightly confusing: the
  "SNPs_removed" file contains all SNPs removed in phase 2d.
  This does not include monomorphic SNPs, or SNPs from excluded
  chromosomes, as these are removed in phase 2a. Also, the
  "SNPs_improbable_values" file does not include SNPs with
  invalid values for crucial
  variables, as these are already in the "SNPs_removed" file.)
  
  Another important but non-fatal problem is missing columns.
  \code{QC_GWAS} uses a translation table to determine the
  contents of a column. If the translation table is incomplete,
  or contains a typo, the function will be unable to translate
  (and therefor use) a column. If this involves, say, callrate,
  it merely means the function cannot apply the callrate filters,
  but the absence of p-values or imputation status will disable
  many features of the QC. If you know that a data column is
  present, yet the log reports it missing, check the translation
  table. The \code{QC_series} checklist output and the
  \code{columns_unidentified} value of the \code{QC_GWAS} return
  report the names of any unidentified columns in the dataset.
    
  If the QC aborts, the log file should give some indication why
  this occurred. However, if it was successful, there will be
  several other tables in the log file.
  
  The second table reports the number of SNPs in the dataset at
  various stages of the QC; as well as how many (and for what
  reason) SNPs were removed.
  
  The third table gives an overview of the data itself. It
  reports how many values were missing and invalid per variable
  in both the pre- and post-QC datasets, as well as the quartile
  and mean values of the post-QC data. A few notes on the nomenclature:
  invalid values will have been removed (for crucial variables) or set
  to missing (for non-crucial variables) in stage 2d. The post-QC
  'invalid' column merely records how many of these SNPs remain in
  the dataset. 'Unusable' values are the missing and invalid values
  combined (shown as percentage of the total data). Finally,
  pre-QC refers to the dataset during stage 2b-c, but
  after monomorphic SNPs and SNPs from excluded chromosomes (stage
  2a) have been removed.
  
  The fourth table reports on the allele-matching in phase 3.
  The concepts are explained in 'details' and the
  \code{\link{match_alleles}} function; here we just mention
  what the user needs to pay attention to. Strand-switching
  counts SNPs whose strand was switched because it mismatched
  with the reference. As many cohorts do not add strand data
  (or set every SNP to \code{"+"}), the presence of such SNPs is
  not a red flag by itself. However, if there are mismatching
  SNPs (the subset of strand-switched SNPs that could not be
  fixed), there is a problem with the allele data (or, possibly,
  a triallelic SNP). Check the
  \code{[filename_output]_SNPs_mismatches-[ref-name].txt} file
  to see the affected SNPs.
  
  Another red flag is if there are strand-switched SNPs in a
  dataset that \emph{also} contains negative-strand SNPs (i.e.
  the cohort included real strand data, rather just setting it
  to \code{"+"} for all SNPs). Negative-strand SNPs are converted
  to the positive strand beforehand, so they should not appear in
  this step (if they do, they are counted in the "double
  strand-switch" entry, but that is of minor importance). The real
  problem is that if a cohort includes negative-strand SNPs (i.e.
  real strand data), and there are still strand-switches, the strand
  data must be incorrect. Whether the strand-switches and the
  negative-strand SNPs overlap is unimportant.
  
  The third possible problem is when a large proportion of
  ambiguous SNPs is suspect: it indicates that they are on the
  wrong strand.
  
  Finally, a large number of SNPs with a deviating allele
  frequency indicates either that the frequency is reported
  for the wrong allele (see below) or that the dataset population
  does not match that of the reference.
  
  The fifth table reports the QC outcome statistics. The p-value
  and allele-frequency correlations should be close to \code{1}.
  An allele-frequency correlation of \code{-1} means that the
  frequency was reported for the wrong (non-effect) allele.
  
  The sixth table reports how many SNPs were removed by the
  various QQ plot filters.
  
  The seventh table reports the chromosomes and alleles present
  in the final dataset.
  
  The eighth table counts invalid values in the pre-
  and post-QC files for several variables. 'Extreme p' is a
  value that is only used when \code{calculate_missing_p} is
  \code{TRUE}. Any newly-calculated p-values that are
  \code{< 1E-300} will be set to \code{1E-300}, in order to
  exclude any values of \code{0} (\code{1E-300} is close to the
  smallest numeric value that R can handle safely).
  
  The final four tables contain the settings of the QC.
  
  
  \emph{The QC_series output}

  \code{QC_series} saves a checklist, showing the most
  important QC stats of the various files side by side, and
  (depending on the plot-arguments) several graphs comparing
  effect-size distribution, precision and skewness vs. kurtosis
  of the QC'ed files.
}
\section{Output header}{
  The standard-format values used by \code{out_header} are:
  \itemize{
    \item \code{"standard"} retains the default column names
      used by \code{QC_GWAS}.
    \item \code{"original"} restores the column names used in
      the input file.
    \item \code{"old"} uses the default column names of pre-v1.0b
      versions of \code{QCGWAS}.
    \item \code{"GWAMA"}, \code{"PLINK"}, \code{"META"} and
      \code{"GenABEL"} set the
      column names to those use by the respective programs. Note
      that META's alleleB corresponds to EFFECT_ALL in \code{QC_GWAS}.
  }
}
\section{Requirement for the input dataset}{
  \code{QC_GWAS} will automatically
  unpack \code{.gz} and \code{.zip} files, provided the filename
  includes the extension of the packed file. For example, if the
  data file is named \code{"data1.csv"}, the zip file should be
  \code{"data1.csv.zip"}. If it is named \code{"data1.zip"},
  \code{QC_GWAS} won't be able to "call" the file inside.
  
  \code{QC_GWAS} is flexible when it comes to file-format. By
  default, it can open datasets with a variety of column separators
  and NA values (the user can specify these via the
  \code{column_separators} and \code{na.strings} arguments).
  Read the documentation of the auto-loader function
  \code{\link{load_GWAS}} for more information.
  
  Chromosome values can be coded as numeric or character: 
  values of \code{"X"}, \code{"Y"}, \code{"XY"} and
  \code{"M"} will automatically be converted to \code{23},
  \code{24}, \code{25} and \code{26}, respectively.
  
  By default, imputation status is coded as integers \code{0}
  (genotyped) and \code{1} (imputed). As of version 1.0-4,
  imputation status will always be translated using the
  \code{imputed_T}, \code{imputed_F} and \code{imputed_NA}
  arguments. This means that the user must specify values for
  these arguments, even when the dataset already uses the
  standard format.
  Because of the importance of imputation status, if the function
  is unable to translate character values, the QC will abort.
  
  The minimal and maximal valid imputation quality values are
  determined by the \code{minimal_impQ_value} and
  \code{maximal_impQ_value} arguments.
  
  Standard errors and p-values of \code{0} are considered invalid and
  removed in phase 2d, while values of \code{-1} will be set to
  \code{NA}. Effect sizes of \code{-1} are accepted, unless the
  standard error and/or the p value are also \code{-1}, in which
  case it is also set to \code{NA}.  
}
\section{Filter arguments}{
  \code{QC_GWAS} has three sets arguments relating to filters:
  the arguments for the HQ (high-quality), QQ (plot) and NA
  (missing values) filters. The HQ and QQ arguments work
  mostly in the same way, but there are a few key differences.
  
  The high-quality arguments accept single, numeric values that
  determine the minimal values of allele-frequency (FRQ), HWE
  p-value (HWE), callrate (cal) and imputation quality (imp) for
  a SNP to be of high-quality. The high-quality filter is used
  for the effect-size boxplot and the Manhattan plot, as well
  as several QC stats.
  
  The QQ arguments accept a vector of max. 5 numeric values that
  are applied sequentially as filters in the QQ plots.
  
  Both of these use the NA filter argument(s) to determine
  whether to exclude or ignore missing values.
  
  Neither filter is used to \emph{remove} SNPs; they merely
  exclude them from several QC tests. Both HQ and QQ filter
  criteria are only applied if the variable passed the threshold
  test, i.e. if there are sufficient non-missing, non-invalid
  values for the filter to be applied (see the \code{use_threshold}
  argument for details). It's pointless to filter an empty column.
  
  If \code{ignore_impstatus} is \code{FALSE} (default), the 
  imputation-quality criterion is only applied to
  imputed SNPs, and the HWE p-value and callrate criteria only
  to genotyped SNPs. If \code{TRUE}, the filters are applied to
  all SNPs, regardless of the imputation status.
  
  The allele-frequency filters are two-sided:
  when set to value x, SNPs with frequency < x or > 1 - x are
  excluded.
  
  To filter missing values only, set the filter to \code{NA} and
  the corresponding NA filter argument to \code{TRUE}. To disable
  entirely, set to \code{NULL} (this means the NA-filter setting is
  ignored as well).
  
  The differences between the HQ and QQ filters are:
  
  The HQ filter arguments accept a single value, the QQ filters
  can accept up to 5.
  
  The HQ filter is a single filter: a SNP needs to meet all
  relevant criteria to be considered high-quality. The QQ filter
  values are applied separately.
  
  The QQ filter has an additional feature: if passed a value
  \code{x > 1}, it will calculate a filter value of
  \code{x / sample size}. This is to allow size-based filtering
  of allele frequencies. Note that this filter uses the sample
  size listed for that specific SNP. If the sample size is
  missing, the relevant NA-filter setting is used to determine
  whether it should be excluded.
}
\section{Updating the alternative reference}{
  There are two drawbacks to the way the function updates the
  alternative reference file. One is a technical issue, but the
  other can affect the QC of subsequent files.

  Firstly, the argument \code{update_alt} has a slightly misleading
  name: the alternative-reference \emph{file} is updated, but
  the reference inside the R memory is not. If the user wants to
  do further QCs with the updated reference, (s)he will have to
  manually reload the updated file into R.
  
  This is caused by the way R handles data-alterations that occur
  inside a function. Any changes made to data last only for the
  duration of the function. Once the function terminates, the memory
  reverts to its original state. In other words: the allele reference
  is updated inside \code{QC_GWAS}, but goes back to the pre-QC
  state once the QC is finished. \code{QC_series} deals with this by
  automatically reloading the reference file whenever it is
  updated, but, again, once the function terminates it will
  revert to its original state.
  
  The second drawback is that the content of the alternative
  reference is arbitrary, depending on which file
  an unknown SNP is encountered in first. For example, suppose
  that SNP rs31 has alleles A and G, an allele frequency that
  varies around 0.5, and does not appear in the standard reference.
  When it is added to the alternative reference, the allele
  listed as the minor one depends entirely on the allele
  frequency in the first file it is encountered in.
  
  More seriously, if the information in the first file is
  incorrect, the SNP may be strand-switched or excluded
  in subsequent files because it does not match the (incorrect)
  reference. This is another reason why it is important to check
  the log files: if there is a problem with a datafile's
  strand, alleles or allele-frequency, and the alternative
  reference was updated, the incorrect data may have been added
  to the reference. If so, one should go back to a previous
  reference file. The argument \code{backup_alt} is useful for
  this, though note that \code{QC_series} only does this the
  first time the reference is updated.
  
  Also, if one wants to QC a large number of files for a meta-GWAS,
  one should use the same alternative allele reference file (and
  let \code{QC_GWAS} update it) for every file, otherwise it is
  possible  that rs13 may have a different effect alleles in some
  post-QC files.
}
\seealso{
  For the plots created by \code{QC_series}:
  \code{\link{plot_distribution}},
  \code{\link{plot_precision}} and \code{\link{plot_skewness}}.
  
  For loading and preparing a GWAS dataset:
  \code{\link{load_GWAS}}, \code{\link{translate_header}},
  \code{\link{convert_impstatus}}.
  
  For carrying out separate steps of \code{QC_GWAS}:
  \code{\link{match_alleles}}, \code{\link{check_P}},
  \code{\link{QC_plots}}.
}
\examples{
## For instructions on how to run QC_GWAS and QC_series
## check the quick start guide in /R/library/QCGWAS/doc
}
