\name{filter_GWAS}
\alias{filter_GWAS}
\title{Automated filtering and reformatting of GWAS results files}
\description{
  This function was created as a convenient way to automate the
  removal of low-quality and non-autosomal SNPs. It
  includes the same formatting options as \code{\link{QC_GWAS}}.
}
\usage{
filter_GWAS(ini_file,
            GWAS_files, output_names,
            gzip_output = TRUE,
            dir_GWAS = getwd(), dir_output = dir_GWAS,
            FRQ_HQ = NULL, HWE_HQ = NULL,
            cal_HQ = NULL, imp_HQ = NULL,
            FRQ_NA = TRUE, HWE_NA = TRUE,
            cal_NA = TRUE, imp_NA = TRUE,
            ignore_impstatus = FALSE,
            remove_X = FALSE, remove_Y = FALSE,
            remove_XY = FALSE, remove_M = FALSE,
            header_translations,
            check_impstatus = FALSE,
            imputed_T = c("1", "TRUE", "yes", "YES", "y", "Y"),
            imputed_F = c("0", "FALSE", "no", "NO", "n", "N"),
            imputed_NA = NULL,
            column_separators = c("\t", " ", "", ",", ";"),
            header = TRUE, nrows = -1, nrows_test = 1000,
            comment.char = "", na.strings = c("NA", "."),
            out_header = "original", out_quote = FALSE,
            out_sep = "\t", out_eol = "\n", out_na = "NA",
            out_dec = ".", out_qmethod = "escape",
            out_rownames = FALSE, out_colnames = TRUE, ...)
}
\arguments{
  \item{ini_file}{
    (the filename of) a table listing the files to be
    processed and the filters to be applied. See 'Details'.}
  \item{GWAS_files}{character vector: when no \code{ini_file} is
    provided, this identifies the files to be processed. See
    'Details'.}
  \item{output_names}{character vector: the filenames for the
    output files. The default option is to use the input
    filenames. Note that, unlike with other \code{QCGWAS}
    functions, the file extensions should be included (However,
    the function will automatically add \code{".gz"} when the
    files are compressed.}
  \item{gzip_output}{logical; should the output files be
    compressed?}    
  \item{dir_GWAS, dir_output}{character-strings specifying the
    directory address of the folders for the input files and the
    output, respectively. Note that R uses \emph{forward} slash
    (/) where Windows uses backslash (\\).}
  \item{FRQ_HQ, HWE_HQ, cal_HQ, imp_HQ}{
    Numeric vectors. When no \code{ini_file} is provided, these
    arguments specify the filter threshold-values for allele
    frequency, HWE p-value, callrate and imputation quality,
    respectively. Passed to \code{\link{HQ_filter}}.}
  \item{FRQ_NA, HWE_NA, cal_NA, imp_NA}{
    Logical vectors. When no \code{ini_file} is provided, these
    arguments specify whether missing values (of allele frequency,
    HWE p-value, callrates and imputation quality, respectively)
    are excluded (\code{TRUE}) or ignored (\code{FALSE}).
    Passed to \code{\link{HQ_filter}}.}
  \item{ignore_impstatus}{
    Logical vector. When no \code{ini_file} is provided, this
    argument specifies whether imputation status is taken into
    account when applying the filters. If \code{FALSE}, HWE p-value
    and callrate filters are applied only to genotyped SNPs, and
    imputation quality filters only to imputed SNPs. If
    \code{TRUE}, the filters are applied to all SNPs regardless
    of the imputation status.}
  \item{remove_X, remove_Y, remove_XY, remove_M}{
    logical; respectively whether X-chromosome, Y-chromosome,
    pseudo-autosomal and mitochondrial SNPs are removed. Note:
    these arguments accept only a single \code{TRUE} or
    \code{FALSE} value.
    Unlike the above settings, it's not possible to specify
    them independently for every dataset.}
  \item{header_translations}{translation table for column names.
    See \code{\link{translate_header}} for more information. If
    the argument is left empty, \code{dataset} is assumed to use
    the standard column-names used by \code{\link{QC_GWAS}}.}
  \item{check_impstatus}{logical; should
    \code{\link{convert_impstatus}} be called to convert the
    imputation-status column into standard values?}
  \item{imputed_T, imputed_F, imputed_NA}{arguments passed to
  \code{\link{convert_impstatus}}.}
  \item{column_separators}{character string or vector; specifies
    the values used as column delimitator in the GWAS file(s). The
    argument is passed to \code{\link{load_test}}; see the
    description of that function for more information.}
  \item{nrows_test}{integer; the number of rows used for
    "trial-loading". Before loading the entire dataset, the
    function \code{\link{load_test}} is called to determine the
    dataset's file-format by reading the top \code{x} lines, where
    \code{x} is \code{nrows_test}. Setting \code{nrows_test} to
    a low number
    (e.g. \code{150}) means quick testing, but runs the risk of
    missing problems in lower rows. To test the entire dataset,
    set it to \code{-1}.}
  \item{header, nrows, comment.char, na.strings, \dots}{
    arguments passed to \code{\link{read.table}} when importing
    the dataset.}
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
      \item The name of a file in \code{dir_GWAS} containing
        such a table.
      \item Character string specifying a standard form. See
        \code{\link{QC_GWAS}}, section 'Output header' for the
        options.
    }}
  \item{out_quote, out_sep, out_eol, out_na, out_dec,
    out_qmethod, out_rownames, out_colnames}{
    arguments passed to
    \code{\link{write.table}} when saving the final dataset.}
}
\details{
  The easiest way to use \code{filter_GWAS} is by passing an ini
  file to the \code{ini_file} argument.
  The ini file can be generated by running \code{\link{QC_series}}
  with the \code{save_filtersettings} argument set to \code{TRUE}.
  The output will include a file 'Check_filtersettings.txt',
  describing the (high-quality) filter settings used for each
  file (taking into account whether there was enough data, i.e.
  whether the \code{use_threshold} was met, to apply the filters).
  
  The \code{ini_file} argument accepts both a table
  or the name of a file in \code{dir_GWAS} or the
  current R working directory.
  
  If no \code{ini_file} is specified, the function will use the
  \code{GWAS_files}, x_HQ, x_NA and \code{ignore_impstatus}
  arguments to construct such a table.
  \code{GWAS_files} can either be a character vector or a single
  value. If a single string, all filenames containing the string
  will be processed. The other arguments can also be a vector or
  a single value; if the latter, they will be recycled to create
  a vector of the correct length.  
  
  If neither \code{ini_file} nor \code{GWAS_files} are specified,
  the function will look for a file
  \code{Check_filtersettings.txt}
  in \code{dir_GWAS} and the current R working directory.
  
  Note that \code{ini_file} overrules the other filter settings,
  i.e. one cannot adjust \code{ini_file} through the other
  arguments.
}
\value{
  An invisible logical vector, indicating which files were
  successfully filtered.}
\note{
  R is not the optimal platform for filtering GWAS files. This
  function was added at the request of a user, but an UNIX script
  is likely to be faster.}
