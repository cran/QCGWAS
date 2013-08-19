match_alleles <-
function(dataset, ref_set, HQ_subset,
                          dataname = "dataset", ref_name = "reference",
                          unmatched_data = !all(dataset$MARKER %in% ref_set$SNP),
                          check_strand = FALSE, save_mismatches = TRUE, delete_mismatches = FALSE,
                          delete_diffEAF = FALSE, threshold_diffEAF = 0.15,
                          check_FRQ = TRUE, check_ambiguous = TRUE, plot_FRQ = FALSE, plot_if_threshold = FALSE, threshold_r = 0.95,
                          return_SNPs = FALSE, return_ref_values = FALSE, 
                          header_translations, header_reference,
                          save_name = dataname, save_dir = getwd(), use_log = FALSE, log_SNPall = nrow(dataset)
) {
  # Note that log_SNPall is meant exclusively for the log files. It
  # won't necessarily equal the number of SNPs passed to the function.
  
  # Part 1: checking headers & merging data	
  #	note that data_col/ref_col have different contents,
  #	depending on whether check_header/check_ref is TRUE or FALSE
  if(delete_diffEAF & !check_FRQ) stop("cannot remove abberant allele-frequencies when check_FRQ is FALSE")
  check_header <- !missing(header_translations)
  check_ref <- !missing( header_reference)
  col_effect <- TRUE
  
  
  data_std <- c("MARKER", "EFFECT_ALL", "OTHER_ALL", "STRAND", "EFFECT", "EFF_ALL_FREQ")[c(TRUE, TRUE, TRUE, check_strand, TRUE, check_FRQ)]
  if(check_header) {
    if(any(duplicated(header_translations[ ,2]))) stop("duplicated elements in header_translations, column 2")
    data_h <- toupper(colnames(dataset))
    data_col <- vector(mode = "integer", length = length(data_std))
    for(forI in 1:length(data_std)) {
      data_cur <- identify_column(data_std[forI], header_translations, data_h)
      if(length(data_cur) == 1L) {
        data_col[forI] <- data_cur
      } else {
        if(length(data_cur) == 0L) { 
          if(data_std[forI] == "EFFECT") { col_effect <- FALSE
          } else { stop(paste("Cannot identify data column:", data_std[forI])) }
        } else { stop(paste("Multiple data columns identified as:", data_std[forI])) }
      }
    }
    if(!col_effect) {
      data_std <- c("MARKER", "EFFECT_ALL", "OTHER_ALL", "STRAND", "EFFECT", "EFF_ALL_FREQ")[c(TRUE, TRUE, TRUE, check_strand, FALSE, check_FRQ)]
      data_col <- data_col[data_col != 0L]
    }
  } else {
    if(!("EFFECT" %in% colnames(dataset))) {
      col_effect <- FALSE
      data_std <- c("MARKER", "EFFECT_ALL", "OTHER_ALL", "STRAND", "EFFECT", "EFF_ALL_FREQ")[c(TRUE, TRUE, TRUE, check_strand, FALSE, check_FRQ)]
    }
    data_col <- data_std
  }
  
  ref_std <- c("SNP", "MINOR", "MAJOR", "MAF")[c(TRUE, TRUE, TRUE, check_FRQ)]
  if(check_ref) {
    if(any(duplicated(header_reference[ ,2]))) stop("duplicated elements in header_reference, column 2")
    ref_h	 <- toupper(colnames(ref_set))
    ref_col <- vector(mode = "integer", length = length(ref_std))
    for(forI in 1:length(ref_std)) {
      ref_cur <- identify_column(ref_std[forI], header_reference, ref_h)
      if(length(ref_cur) == 1L) {
        ref_col[forI] <- ref_cur
      } else {
        if(length(ref_cur) == 0L) { stop(paste("Cannot identify reference column:", ref_std[forI]))
        } else {			    stop(paste("Multiple reference columns identified as:", ref_std[forI])) }
      }
    }
  } else { ref_col <- ref_std }
  
  if(any(is.na(dataset[ , data_col[1]]))) { stop("missing SNP IDs in marker name column") }
  if(any(is.na(ref_set[ ,	ref_col[1]]))) { stop("missing SNP IDs in reference SNP name column") }
  if(any(duplicated(dataset[ , data_col[1]]))) { stop("duplicate SNPs in marker name column") }
  if(any(duplicated(ref_set[ , ref_col[1]]))) { stop("duplicate SNPs in reference SNP name column") }
  
  # The awkward temp_order column is necessary because merge does not maintain the original order of
  # x when there are SNPs that have no match in y.
  order_col <- length(data_std) + 1L
  SNPset <- merge(x = cbind(dataset[ , data_col], 1:nrow(dataset)), y = ref_set[ , ref_col], by.x = 1, by.y = 1, all.x = TRUE, all.y = FALSE, sort = FALSE)
  if(check_header | check_ref) { colnames(SNPset) <- c(data_std, "temp_order", ref_std[-1])
  } else {				 colnames(SNPset)[order_col] <- "temp_order" }
  if(unmatched_data) { SNPset <- SNPset[order(SNPset$temp_order), ] }
  
  missing_list <- is.na(SNPset$EFFECT_ALL) | is.na(SNPset$OTHER_ALL) | is.na(SNPset$MINOR) | is.na(SNPset$MAJOR)
  SNPn_missing <- sum(missing_list)
  if(SNPn_missing > 0L) {
    SNPn_missing_data <- sum(is.na(SNPset$EFFECT_ALL) | is.na(SNPset$OTHER_ALL))
    SNPn_missing_ref	<- sum(is.na(SNPset$MINOR	) | is.na(SNPset$MAJOR	))
    if(SNPn_missing_data > 0L) {
      if(use_log) { save_log(phaseL = 3L, checkL = paste("allele match -", ref_name, sep=" "), typeL = "alleles", SNPL = SNPn_missing_data, allSNPs = log_SNPall, actionL = "excluded from test", noteL = "no allele information: cannot match these SNPs", fileL = paste(save_dir, save_name, sep = "/")) }
      print(paste(" - - missing alleles in", SNPn_missing_data, "SNPs"), quote = FALSE)
    }
    if(SNPn_missing_ref > 0L) {
      if(use_log) { save_log(phaseL = 3L, checkL = paste("allele match -", ref_name, sep=" "), typeL = "incomplete reference", SNPL = SNPn_missing_ref, allSNPs = log_SNPall, actionL = "excluded from test", noteL = "reference does not contain alleles for these SNP", fileL = paste(save_dir, save_name, sep = "/")) }
      print(paste(" - - incomplete reference: no alleles found for", SNPn_missing_ref, "SNPs"), quote = FALSE)
    }
  } else {
    SNPn_missing_data <- 0L
    SNPn_missing_ref	<- 0L
  }
  
  # Part 2: checking strand column for negative-strand SNPs
  if(check_strand) {
    min_list <- SNPset$STRAND == "-" & !is.na(SNPset$STRAND)
    SNPn_min <- sum(min_list)
    SNPn_min_SS <- 0L
    SNPn_min_MM <- 0L
    if(SNPn_min > 0L) { SNPset[min_list, c("EFFECT_ALL", "OTHER_ALL", "STRAND")] <- switch_strand(SNPset[min_list, c("EFFECT_ALL", "OTHER_ALL", "STRAND")], strand_col = TRUE) }
  } else {
    SNPn_min <- NA
    SNPn_min_SS <- NA
    SNPn_min_MM <- NA
  }
  
  # Part 3: checking alleles for strand-mismatches
  switch_list <- !( ( SNPset$EFFECT_ALL == SNPset$MINOR & SNPset$OTHER_ALL == SNPset$MAJOR ) |
                      ( SNPset$OTHER_ALL == SNPset$MINOR & SNPset$EFFECT_ALL == SNPset$MAJOR ) | missing_list )
  SNPn_switch <- sum(switch_list)
  if(SNPn_switch > 0L) {
    if(check_strand & SNPn_min > 0L) { SNPn_min_SS <- sum(min_list[switch_list]) }
    SNPset[switch_list, c("EFFECT_ALL", "OTHER_ALL")] <-
      switch_strand(SNPset[switch_list, c("EFFECT_ALL", "OTHER_ALL")], strand_col = FALSE)
    
    # Part 3b: checking for mismatching SNPs (i.e. SNPs whose alleles do not match even after strand-switch)
    #	Depending on the settings, these are switched back to how they were before entering part 3,
    #	saved to a file, and deleted.
    mismatch_list <- !( ( SNPset$EFFECT_ALL == SNPset$MINOR & SNPset$OTHER_ALL == SNPset$MAJOR ) |
                          ( SNPset$OTHER_ALL == SNPset$MINOR & SNPset$EFFECT_ALL == SNPset$MAJOR ) | missing_list )
    SNPn_mismatch <- sum(mismatch_list)
    
    if(SNPn_mismatch > 0L) {
      missing_list[mismatch_list] <- TRUE
      if(save_mismatches | !delete_mismatches) {
        SNPset[mismatch_list, c("EFFECT_ALL", "OTHER_ALL")] <-
          switch_strand(SNPset[mismatch_list, c("EFFECT_ALL", "OTHER_ALL")], strand_col = FALSE)
      }
      if(save_mismatches) {
        if(check_strand) {
          negative_strand_correction <- min_list[mismatch_list]
          if(SNPn_min > 0L) { SNPn_min_MM <- sum(negative_strand_correction) }
        }
        write.table( SNPset[mismatch_list, -order_col ], paste(save_dir, "/", save_name, "_SNPs_mismatches-", ref_name, ".txt", sep = ""), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
      } else if(check_strand & SNPn_min > 0L) { SNPn_min_MM <- sum(min_list[mismatch_list]) }
      if(delete_mismatches) {				
        SNPset[mismatch_list, "EFFECT_ALL"] <- NA
        if(use_log) { save_log(phaseL = 3L, checkL = paste("allele match -", ref_name, sep=" "), typeL = "allele mismatch", SNPL = SNPn_mismatch, allSNPs = log_SNPall, actionL = "SNPs removed", noteL = "Allele mismatch with reference; strand correction did not correct this", fileL = paste(save_dir, save_name, sep = "/")) }
      } else {
        if(use_log) { save_log(phaseL = 3L, checkL = paste("allele match -", ref_name, sep=" "), typeL = "allele mismatch", SNPL = SNPn_mismatch, allSNPs = log_SNPall, actionL = "-", noteL = "Allele mismatch with reference; strand correction did not correct this", fileL = paste(save_dir, save_name, sep = "/")) }
      }
      if(!use_log) { print(paste(" - - warning:", SNPn_mismatch, "mismatches found - strand correction did not correct this"), quote = FALSE) }
    }
  } else { SNPn_mismatch <- 0L }
  
  # Part 4: matching allele configuration with reference
  flip_list <- SNPset$OTHER_ALL == SNPset$MINOR & SNPset$EFFECT_ALL == SNPset$MAJOR & !missing_list
  SNPn_flip <- sum(flip_list)
  if(SNPn_flip > 0L) {
    temp_al2 <- SNPset$OTHER_ALL
    SNPset$OTHER_ALL <- ifelse(flip_list, SNPset$EFFECT_ALL, SNPset$OTHER_ALL)
    SNPset$EFFECT_ALL <- ifelse(flip_list, temp_al2, SNPset$EFFECT_ALL)
    if(col_effect) { SNPset$EFFECT <- ifelse(flip_list, -SNPset$EFFECT, SNPset$EFFECT) }
    if(check_FRQ) { SNPset$EFF_ALL_FREQ <- ifelse(flip_list, 1 - SNPset$EFF_ALL_FREQ, SNPset$EFF_ALL_FREQ) }
  }
  
  # Part 5: checking for "ambiguous" SNPs, i.e. SNPs with an A/T or C/G allele configuration
  ambiguous_list <- (( SNPset$EFFECT_ALL == "A" & SNPset$OTHER_ALL == "T" ) |
                       ( SNPset$EFFECT_ALL == "T" & SNPset$OTHER_ALL == "A" ) |
                       ( SNPset$EFFECT_ALL == "C" & SNPset$OTHER_ALL == "G" ) |
                       ( SNPset$EFFECT_ALL == "G" & SNPset$OTHER_ALL == "C" ) ) & !missing_list
  SNPn_ambiguous <- sum(ambiguous_list)
  if(check_FRQ & SNPn_ambiguous > 0L) { SNPn_suspect <- sum(ambiguous_list & !is.na(SNPset$EFF_ALL_FREQ) & !is.na(SNPset$MAF) & ( ( SNPset$EFF_ALL_FREQ > 0.65 & SNPset$MAF < 0.35 ) | ( SNPset$EFF_ALL_FREQ < 0.35 & SNPset$MAF > 0.65 ) ))
  } else { SNPn_suspect <- if(check_FRQ) 0L else NA }
  
  # Part 6: correlating allele frequency with reference
  
  if(plot_FRQ) {
    if(missing(HQ_subset)) {
      HQ_subset <- vector(mode = "logical", length = nrow(SNPset))
      plot_HQ <- FALSE
    } else {
      if(is.numeric(HQ_subset)) {
        temp <- vector(mode = "logical", length = nrow(SNPset))
        temp[HQ_subset] <- TRUE
        HQ_subset <- temp
      }
      plot_HQ <- TRUE
    }
  }
  
  FRQcor <- NA
  FRQcor_ambi <- NA
  FRQcor_others <- NA
  SNPn_diffEAF <- NA

  if(check_FRQ) {
    FRQ_all <- SNPset[ , "EFF_ALL_FREQ"]
    if(SNPn_missing + SNPn_mismatch > 0L) { FRQ_all[missing_list] <- NA }
    diffEAF_list <- abs(FRQ_all - SNPset$MAF) > threshold_diffEAF & !is.na(FRQ_all) & !is.na(SNPset$MAF)
    SNPn_diffEAF <- sum(diffEAF_list)
    if(SNPn_diffEAF > 0L) {
      save_log(phaseL = 3L, checkL = paste("allele match -", ref_name, sep=" "), typeL = "allele frequency", SNPL = SNPn_diffEAF, allSNPs = log_SNPall, actionL = if(delete_diffEAF) "Markers removed" else "-", noteL = paste("Allele frequencies deviate from those in", ref_name,"( > ", 100 * threshold_diffEAF, "% )"), fileL = paste(save_dir, save_name, sep = "/"))
      if(!use_log) { print(paste(" - - warning:", SNPn_diffEAF, "SNPs whose allele-frequency differs strongly from reference"), quote = FALSE) }
      if(delete_diffEAF) {
        write.table(data.frame(SNPset[diffEAF_list, -order_col], DIFFERENCE = abs(SNPset$EFF_ALL_FREQ[diffEAF_list] - SNPset$MAF[diffEAF_list]) ), paste(save_dir, "/", save_name, "_SNPs_EAFdifferent-", ref_name, ".txt", sep = ""), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
        SNPset[diffEAF_list, "EFFECT_ALL"] <- NA
      }
    }    
    if(sum( !( is.na(FRQ_all) | is.na(SNPset$MAF) )) > 10L ) {
      FRQcor <- cor(FRQ_all, SNPset$MAF, use = "na.or.complete")
      if(FRQcor < threshold_r) {
        if(use_log) save_log(phaseL = 3L, checkL = paste("allele match -", ref_name, sep=" "), typeL = "allele frequency", SNPL = log_SNPall, allSNPs = log_SNPall, actionL = "-", noteL = paste("Allele frequencies correlate poorly with those in", ref_name,"( r < ", threshold_r, ")"), fileL = paste(save_dir, save_name, sep = "/"))
        print(paste(" - - warning: allele frequency correlates poorly with reference (r = ", round(FRQcor, digits = 3), ")", sep = ""), quote = FALSE)
      }
      if(plot_FRQ & (FRQcor < threshold_r | !plot_if_threshold)) {
        jpeg(paste(save_dir, "/", save_name, "_graph_EAF-", ref_name, ".jpg", sep = ""), width = 720, height = 720, res = 144)
        plot(0.5, 0.5, xlim = c(0,1), ylim = c(0,1), col = "white",
             main=paste(ref_name, "allele-frequency correlation"), sub = dataname, cex.sub = 1.3,
             xlab="Reported allele frequency", ylab="Minor allele frequency")
        points(FRQ_all[!HQ_subset], SNPset$MAF[!HQ_subset], pch = 20, col = ifelse(plot_HQ, "grey", "black"), cex = 0.8)
        points(FRQ_all[HQ_subset], SNPset$MAF[HQ_subset], pch = 20, col = "black", cex = 0.8)
        if(plot_HQ) legend(0, 0.94, c(" Low quality", "High quality"), pch = 20, col = c("grey", "black"))
        text(0.1, 0.98, paste("r =", round(FRQcor, digits = 3)), pos = 4, cex=1.2, col = ifelse(FRQcor < threshold_r, "red", "black") )
        dev.off()
      }
      
      if(check_ambiguous & SNPn_ambiguous > 0L) {
        FRQ_ambi <- ifelse(ambiguous_list, FRQ_all, NA)
        FRQ_rest <- ifelse(ambiguous_list, NA, FRQ_all)
        if(sum( !( is.na(FRQ_ambi) | is.na(SNPset$MAF) )) > 10L ) {
          FRQcor_ambi <- cor(FRQ_ambi, SNPset$MAF, use = "na.or.complete")
          if(FRQcor_ambi < threshold_r) { 
            if(use_log) save_log(phaseL = 3L, checkL = paste("allele match -", ref_name, sep=" "), typeL = "ambiguous SNPs", SNPL = SNPn_ambiguous, allSNPs = SNPn_ambiguous, actionL = "-" , noteL = paste("Ambiguous SNPs' allele frequencies correlate poorly with those in", ref_name,"( r <", threshold_r, ")"), fileL = paste(save_dir, save_name, sep = "/"))
            print(paste(" - - warning: allele frequency of ambiguous SNPs correlates poorly with reference (r = ", round(FRQcor_ambi, digits = 3), ")", sep = ""), quote = FALSE)
          }
          SNPn_diffEAF_ambi <- sum(diffEAF_list[ambiguous_list])
          if(SNPn_diffEAF_ambi > 0L) { save_log(phaseL = 3L, checkL = paste("allele match -", ref_name, sep=" "), typeL = "ambiguous SNPs", SNPL = SNPn_diffEAF_ambi, allSNPs = SNPn_ambiguous, actionL = ifelse( delete_diffEAF, "Markers removed", "-"), noteL = paste("Ambiguous SNPs' allele frequencies deviate from those in", ref_name,"( > ", 100 * threshold_diffEAF, "% )"), fileL = paste(save_dir, save_name, sep = "/")) }
          if(plot_FRQ & (FRQcor_ambi < threshold_r | !plot_if_threshold)) {
            jpeg(paste(save_dir, "/", save_name, "_graph_EAF-", ref_name, "-ambiguous.jpg", sep = ""), width = 720, height = 720, res = 144)
            plot(0.5, 0.5, xlim = c(0,1), ylim = c(0,1), col = "white",
                 main=paste(ref_name, "allele-frequency correlation"), sub = dataname, cex.sub = 1.3,
                 xlab="Reported allele frequency", ylab="Minor allele frequency")
            mtext("ambiguous SNPs", line = 0.5, cex = 1.2)
            points(FRQ_ambi[!HQ_subset], SNPset$MAF[!HQ_subset], pch = 20, col = ifelse(plot_HQ, "grey", "black"), cex = 0.8)
            points(FRQ_ambi[HQ_subset], SNPset$MAF[HQ_subset], pch = 20, col = "black", cex = 0.8)
            if(plot_HQ) legend(0, 0.94, c(" Low quality", "High quality"), pch = 20, col = c("grey", "black"))
            text(0.1, 0.98, paste("r =", round(FRQcor_ambi, digits = 3)), pos = 4, cex=1.2, col = ifelse(FRQcor_ambi < threshold_r, "red", "black"))
            dev.off()
          }
        } else {
          if(use_log) { save_log(phaseL = 3L, checkL = paste("allele match -", ref_name, sep=" "), typeL = "ambiguous SNPs", SNPL = SNPn_ambiguous, allSNPs = log_SNPall, actionL = "-", noteL = "Insufficient frequency data to correlate ambiguous SNPs' allele frequencies with reference", fileL = paste(save_dir, save_name, sep = "/"))
          } else { print(" - - insufficient non-missing allele frequencies to correlate ambiguous SNPs' allele frequencies with reference", quote = FALSE) }
        }
        if(sum( !( is.na(FRQ_rest) | is.na(SNPset$MAF) )) > 10L ) {
          FRQcor_others <- cor(FRQ_rest, SNPset$MAF, use = "na.or.complete")
          if(FRQcor_others < threshold_r) { 
            if(use_log) save_log(phaseL = 3L, checkL = paste("allele match -", ref_name), typeL = "non-ambiguous SNPs", SNPL = nrow(SNPset) - SNPn_missing - SNPn_mismatch - SNPn_ambiguous, allSNPs = nrow(SNPset) - SNPn_missing - SNPn_mismatch - SNPn_ambiguous, actionL = "-", noteL = paste("Non-ambiguous SNPs' allele frequencies correlate poorly with those in", ref_name,"( r <",threshold_r, ")"), fileL = paste(save_dir, save_name, sep = "/") )
            print(paste(" - - warning: allele frequency of non-ambiguous SNPs correlates poorly with reference (r = ", round(FRQcor_others, digits = 3), ")", sep = ""), quote = FALSE)
          }
          SNPn_diffEAF_others <- sum(diffEAF_list[!ambiguous_list])
          if(SNPn_diffEAF_others > 0L) { save_log(phaseL = 3L, checkL = paste("allele match -", ref_name, sep=" "), typeL = "non-ambiguous SNPs", SNPL = SNPn_diffEAF_others, allSNPs = nrow(SNPset) - SNPn_missing - SNPn_mismatch - SNPn_ambiguous, actionL = ifelse( delete_diffEAF, "Markers removed", "-"), noteL = paste("Non-ambiguous SNPs' allele frequencies deviate from those in", ref_name,"( > ", 100 * threshold_diffEAF, "% )"), fileL = paste(save_dir, save_name, sep = "/")) }
          if(plot_FRQ & (FRQcor_others < threshold_r | !plot_if_threshold)) {
            jpeg(paste(save_dir, "/", save_name, "_graph_EAF-", ref_name, "-nonambiguous.jpg", sep = ""),
                 width = 720, height = 720, res = 144)
            plot(0.5, 0.5, xlim = c(0,1), ylim = c(0,1), col = "white",
                 main=paste(ref_name, "allele-frequency correlation"), sub = dataname, cex.sub = 1.3,
                 xlab="Reported allele frequency", ylab="Minor allele frequency")
            mtext("non-ambiguous SNPs", line = 0.5, cex = 1.2)
            points(FRQ_rest[!HQ_subset], SNPset$MAF[!HQ_subset], pch = 20, col = ifelse(plot_HQ, "grey", "black"), cex = 0.8)
            points(FRQ_rest[HQ_subset], SNPset$MAF[HQ_subset], pch = 20, col = "black", cex = 0.8)
            if(plot_HQ) legend(0, 0.94, c(" Low quality", "High quality"), pch = 20, col = c("grey", "black"))
            text(0.1, 0.98, paste("r =", round(FRQcor_others, digits = 3)), pos = 4, cex=1.2, col = ifelse(FRQcor_others < threshold_r, "red", "black"))
            dev.off()
          }
        } else {
          if(use_log) { save_log(phaseL = 3L, checkL = paste("allele match -", ref_name, sep=" "), typeL = "non-ambiguous SNPs", SNPL = nrow(SNPset) - SNPn_missing - SNPn_mismatch - SNPn_ambiguous, allSNPs = log_SNPall, actionL = "-", noteL = "Insufficient frequency data to correlate non-ambiguous SNPs' allele frequencies with reference", fileL = paste(save_dir, save_name, sep = "/"))
          } else { print(" - - insufficient non-missing allele frequencies to correlate non-ambiguous SNPs' allele frequencies with reference", quote = FALSE) }
        }
      }
    } else {
      if(use_log) { save_log(phaseL = 3L, checkL = paste("allele match -", ref_name, sep=" "), typeL = "allele frequency", SNPL = nrow(SNPset), allSNPs = log_SNPall, actionL = "-", noteL = "Insufficient frequency data to correlate with reference", fileL = paste(save_dir, save_name, sep = "/"))
      } else { print(" - - insufficient non-missing allele frequencies to correlate with reference", quote = FALSE) }
    }
  }
  
  match_alleles_list <- list(FRQ_cor = FRQcor, FRQ_cor_ambiguous = FRQcor_ambi, FRQ_cor_nonambi = FRQcor_others,
                             n_SNPs = nrow(SNPset), n_missing = SNPn_missing, n_missing_data = SNPn_missing_data, n_missing_ref = SNPn_missing_ref,
                             n_negative_strand = SNPn_min, n_negative_switch = SNPn_min_SS, n_negative_mismatch = SNPn_min_MM,
                             n_strandswitch = SNPn_switch, n_mismatch = SNPn_mismatch,
                             n_flipped = SNPn_flip, n_ambiguous = SNPn_ambiguous, n_suspect = SNPn_suspect, n_diffEAF = SNPn_diffEAF,
                             
                             MARKER = if(return_SNPs | return_ref_values) { SNPset$MARKER } else { NULL },
                             EFFECT_ALL = if(return_SNPs) { SNPset$EFFECT_ALL } else { NULL },
                             OTHER_ALL = if(return_SNPs) { SNPset$OTHER_ALL } else { NULL },
                             STRAND = if(check_strand & return_SNPs) { SNPset$STRAND } else { NULL },
                             EFFECT = if(col_effect & return_SNPs) { SNPset$EFFECT } else { NULL },
                             EFF_ALL_FREQ = if(check_FRQ & return_SNPs) { SNPset$EFF_ALL_FREQ } else { NULL },
                             
                             ref_MINOR = if(return_ref_values) { SNPset$MINOR } else { NULL },
                             ref_MAJOR = if(return_ref_values) { SNPset$MAJOR } else { NULL },
                             ref_MAF = if(check_FRQ & return_ref_values) { SNPset$MAF } else { NULL }
  )	
  return(match_alleles_list)
}
