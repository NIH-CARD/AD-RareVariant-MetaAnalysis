################################################################################
# STEP 1: TRANS-ANCESTRY META-ANALYSIS
################################################################################
#
# Description:
#   Performs trans-ancestry meta-analysis of gene-based burden test results
#   from multiple cohorts and ancestries using Fisher's and Stouffer's methods.
#
# Methods:
#   - Fisher's combined probability test
#   - Stouffer's Z-score method (weighted by sqrt(N_effective))
#   - Multiple testing correction: Bonferroni and FDR (Benjamini-Hochberg)
#   - Genomic inflation factor (λ and λ₁₀₀₀) calculation
#   - Enhanced robust classification (FDR, Bonferroni, Ultra)
#
# Input:
#   - REGENIE gene burden test results (.regenie.gz files)
#   - Multiple MAF thresholds (0.01, 0.03, 0.05)
#   - Variant types: all variants and coding variants only
#   - Cohorts: UKB (EUR) and ADSP (EUR, AFR, AMR, AAC, AJ, EAS, CAH)
#
# Output:
#   - trans_ancestry_COMPLETE.txt: Complete meta-analysis results
#   - lambda_summary.txt: Genomic inflation factors
#   - robust_methods_comparison.txt: Method comparison summary
#   - FDR/Bonferroni significant genes by method and parameters
#   - ROBUST results: genes passing multiple significance criteria
#
# Requirements:
#   - R packages: data.table, dplyr
#
# Author: Marzieh Khani
# Date: 02/04/2026
################################################################################

library(data.table)
library(dplyr)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

sample_info <- data.frame(
  cohort = c("UKB", "ADSP", "ADSP", "ADSP", "ADSP", "ADSP", "ADSP", "ADSP"),
  ancestry = c("EUR", "EUR", "AFR", "AMR", "AAC", "AJ", "EAS", "CAH"),
  n_cases = c(4394, 6391, 848, 1245, 1158, 575, 923, 1623),
  n_controls = c(54816, 8417, 1993, 2897, 1858, 809, 1222, 2343)
)

sample_info$cohort_ancestry <- paste(sample_info$cohort, sample_info$ancestry, sep="_")
sample_info$n_total <- sample_info$n_cases + sample_info$n_controls

maf_codes <- c("MAF01", "MAF03", "MAF05")
maf_codes_ukb <- c("MAF001", "MAF003", "MAF005")
maf_values <- c(0.01, 0.03, 0.05)
maf_labels <- c("MAF001", "MAF003", "MAF005")
variant_types <- c("all", "coding")

# Set paths to your REGENIE results directories
ADSP_DIR <- "path/to/ADSP/results"  
UKB_DIR <- "path/to/UKB/results"    

cat("Sample sizes:\n")
print(sample_info)
cat(sprintf("\nNote: EUR in 2 independent cohorts (%.1f%% of samples)\n\n",
            100 * sum(sample_info$n_total[sample_info$ancestry == "EUR"]) / sum(sample_info$n_total)))

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

calculate_neff <- function(n_cases, n_controls) {
  (4 * n_cases * n_controls) / (n_cases + n_controls)
}

fishers_method <- function(pvalues) {
  pvalues <- pvalues[!is.na(pvalues)]
  pvalues[pvalues <= 0] <- 1e-300  
  pvalues[pvalues > 1] <- 1
  pvalues <- pvalues[pvalues > 0 & pvalues <= 1]
  if(length(pvalues) < 2) return(list(p = NA, chisq = NA, df = NA))
  chi_sq <- -2 * sum(log(pvalues))
  df <- 2 * length(pvalues)
  p <- pchisq(chi_sq, df, lower.tail = FALSE)
  return(list(p = p, chisq = chi_sq, df = df))
}

stouffers_method <- function(pvalues, weights) {
  valid <- !is.na(pvalues) & !is.na(weights)
  pvalues <- pvalues[valid]
  weights <- weights[valid]
  pvalues[pvalues <= 0] <- 1e-300  
  pvalues[pvalues >= 1] <- 1 - 1e-10
  pvalues <- pvalues[pvalues > 0 & pvalues < 1]
  weights <- weights[1:length(pvalues)]
  if(length(pvalues) < 2) return(list(p = NA, z = NA))
  z_scores <- qnorm(pvalues / 2, lower.tail = FALSE)
  weighted_z <- sum(weights * z_scores) / sqrt(sum(weights^2))
  p <- 2 * pnorm(-abs(weighted_z))
  return(list(p = p, z = weighted_z))
}

calculate_lambda <- function(pvalues) {
  pvalues <- pvalues[!is.na(pvalues) & pvalues > 0 & pvalues < 1]
  if(length(pvalues) < 10) return(NA)
  chisq <- qchisq(1 - pmin(pvalues, 1 - 1e-15), df = 1)
  lambda <- median(chisq) / qchisq(0.5, df = 1)
  return(lambda)
}

calculate_lambda_z <- function(zscores) {
  zscores <- zscores[!is.na(zscores) & is.finite(zscores)]
  if(length(zscores) < 10) return(NA)
  chisq <- zscores^2  
  lambda <- median(chisq) / qchisq(0.5, df = 1)
  return(lambda)
}

calculate_lambda_1000 <- function(pvalues, n_cases, n_controls) {
  lambda <- calculate_lambda(pvalues)
  if(is.na(lambda)) return(NA)
  n_eff_study <- calculate_neff(n_cases, n_controls)
  n_eff_1000 <- 1000
  lambda_1000 <- 1 + (lambda - 1) * (n_eff_1000 / n_eff_study)
  return(lambda_1000)
}

calculate_lambda_1000_z <- function(zscores, n_cases, n_controls) {
  lambda <- calculate_lambda_z(zscores)
  if(is.na(lambda)) return(NA)
  n_eff_study <- calculate_neff(n_cases, n_controls)
  n_eff_1000 <- 1000
  lambda_1000 <- 1 + (lambda - 1) * (n_eff_1000 / n_eff_study)
  return(lambda_1000)
}

cat("Functions defined\n\n")

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("========================================================================\n")
cat("Loading REGENIE Results\n")
cat("========================================================================\n\n")

all_results <- list()
counter <- 1

for(maf_idx in 1:length(maf_codes)) {
  maf_adsp <- maf_codes[maf_idx]
  maf_ukb <- maf_codes_ukb[maf_idx]
  maf_val <- maf_values[maf_idx]
  
  for(vtype in variant_types) {
    cat(sprintf("\n=== MAF %s, %s ===\n", maf_val, vtype))
    
    for(i in seq_len(nrow(sample_info))[-1]) {  
      ancestry <- sample_info$ancestry[i]
      cohort_id <- sample_info$cohort_ancestry[i]
      file_path <- sprintf("%s/new_results_step2_%s_%s_%s_SKATO_meta.regenie.gz",
                          ADSP_DIR, ancestry, maf_adsp, vtype)
      
      if(file.exists(file_path)) {
        cat(sprintf("  [%d] %s\n", counter, basename(file_path)))
        df <- fread(file_path)
        df$GENE <- gsub("\\.M1\\..*", "", df$ID)
        df$PVAL <- 10^(-df$LOG10P)
        df$PVAL[df$PVAL == 0 | is.na(df$PVAL)] <- 1e-300
        df$PVAL[df$PVAL > 1] <- 1
        df$cohort <- "ADSP"
        df$ancestry <- ancestry
        df$cohort_ancestry <- cohort_id
        df$maf_threshold <- maf_val
        df$variant_type <- vtype
        df$n_cases <- sample_info$n_cases[i]
        df$n_controls <- sample_info$n_controls[i]
        df$n_total <- sample_info$n_total[i]
        all_results[[counter]] <- df
        counter <- counter + 1
      }
    }
    
    ukb_file <- sprintf("%s/EUR_%s_%s_SKATO_meta.regenie.gz", UKB_DIR, maf_ukb, vtype)
    if(file.exists(ukb_file)) {
      cat(sprintf("  [%d] %s\n", counter, basename(ukb_file)))
      df <- fread(ukb_file)
      df$GENE <- gsub("\\.M1\\..*", "", df$ID)
      df$PVAL <- 10^(-df$LOG10P)
      df$PVAL[df$PVAL == 0 | is.na(df$PVAL)] <- 1e-300
      df$PVAL[df$PVAL > 1] <- 1
      df$cohort <- "UKB"
      df$ancestry <- "EUR"
      df$cohort_ancestry <- "UKB_EUR"
      df$maf_threshold <- maf_val
      df$variant_type <- vtype
      df$n_cases <- sample_info$n_cases[1]
      df$n_controls <- sample_info$n_controls[1]
      df$n_total <- sample_info$n_total[1]
      all_results[[counter]] <- df
      counter <- counter + 1
    }
  }
}

combined_data <- rbindlist(all_results, fill=TRUE)
cat(sprintf("\nLoaded %d files, %d rows, %d genes\n\n", 
    length(all_results), nrow(combined_data), length(unique(combined_data$GENE))))

combined_data <- combined_data %>%
  filter(!is.na(PVAL), PVAL > 0, PVAL <= 1, !is.na(GENE))

# ==============================================================================
# META-ANALYSIS
# ==============================================================================

cat("Running meta-analysis...\n")
cat("Note: Stouffer's method weighted by sqrt(Neff) per cohort for proper case-control adjustment\n\n")

meta_results <- combined_data %>%
  group_by(GENE, maf_threshold, variant_type) %>%
  summarize(
    n_cohorts = n_distinct(cohort_ancestry),  
    n_ancestries = n_distinct(ancestry),
    total_cases = sum(n_cases),
    total_controls = sum(n_controls),
    total_n = sum(n_total),
    min_p = min(PVAL),
    max_p = max(PVAL),
    median_p = median(PVAL),
    cohorts = paste(cohort_ancestry, collapse=";"),
    fisher_result = list(fishers_method(PVAL)),
    stouffer_result = list(stouffers_method(PVAL, sqrt(calculate_neff(n_cases, n_controls)))),
    .groups = "drop"
  ) %>%
  mutate(
    fisher_p = sapply(fisher_result, function(x) x$p),
    fisher_log10p = -log10(fisher_p),
    stouffer_p = sapply(stouffer_result, function(x) x$p),
    stouffer_log10p = -log10(stouffer_p),
    stouffer_z = sapply(stouffer_result, function(x) x$z)  
  ) %>%
  select(-fisher_result, -stouffer_result)

cat("Complete\n\n")

# ==============================================================================
# MULTIPLE TESTING CORRECTION
# ==============================================================================

cat("Applying corrections...\n")

meta_results <- meta_results %>%
  group_by(maf_threshold, variant_type) %>%
  mutate(
    n_tests = n(),
    fisher_bonf = p.adjust(fisher_p, method = "bonferroni"),
    stouffer_bonf = p.adjust(stouffer_p, method = "bonferroni"),
    fisher_fdr = p.adjust(fisher_p, method = "BH"),
    stouffer_fdr = p.adjust(stouffer_p, method = "BH"),
    fisher_gwas = fisher_p < 5e-8,
    fisher_bonf_sig = fisher_bonf < 0.05,
    fisher_fdr_sig = fisher_fdr < 0.05,
    fisher_nominal = fisher_p < 0.05,
    stouffer_gwas = stouffer_p < 5e-8,
    stouffer_bonf_sig = stouffer_bonf < 0.05,
    stouffer_fdr_sig = stouffer_fdr < 0.05,
    stouffer_nominal = stouffer_p < 0.05,
    fisher_sig_level = case_when(
      fisher_gwas ~ "Genome-wide",
      fisher_bonf_sig ~ "Bonferroni",
      fisher_fdr_sig ~ "FDR",
      fisher_nominal ~ "Nominal",
      TRUE ~ "Not significant"
    ),
    stouffer_sig_level = case_when(
      stouffer_gwas ~ "Genome-wide",
      stouffer_bonf_sig ~ "Bonferroni",
      stouffer_fdr_sig ~ "FDR",
      stouffer_nominal ~ "Nominal",
      TRUE ~ "Not significant"
    )
  ) %>%
  ungroup()

cat("Complete\n\n")

# ==============================================================================
# ENHANCED ROBUST GENE CLASSIFICATION
# ==============================================================================

cat("Defining robust gene categories...\n")

meta_results <- meta_results %>%
  mutate(
    # both methods FDR significant
    robust_fdr = fisher_fdr_sig & stouffer_fdr_sig,
    
    # both methods Bonferroni significant
    robust_bonf = fisher_bonf_sig & stouffer_bonf_sig,
    
    # most stringent - both methods pass BOTH corrections
    robust_ultra = (fisher_fdr_sig & fisher_bonf_sig) & 
                   (stouffer_fdr_sig & stouffer_bonf_sig),
    
    # Classification label for reporting
    robust_category = case_when(
      robust_ultra ~ "Ultra-robust (both corrections)",
      robust_fdr & !robust_bonf ~ "Robust (FDR only)",
      robust_bonf & !robust_fdr ~ "Robust (Bonferroni only)",
      TRUE ~ "Not robust"
    )
  ) %>%
  arrange(maf_threshold, variant_type, fisher_p)

# Summary of robust classifications
robust_summary <- meta_results %>%
  group_by(maf_threshold, variant_type) %>%
  summarize(
    n_genes = n(),
    robust_fdr = sum(robust_fdr, na.rm = TRUE),
    robust_bonf = sum(robust_bonf, na.rm = TRUE),
    robust_ultra = sum(robust_ultra, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nRobust gene summary by correction method:\n")
print(robust_summary)

cat("\nDetailed breakdown by category:\n")
category_breakdown <- meta_results %>%
  count(maf_threshold, variant_type, robust_category) %>%
  tidyr::pivot_wider(names_from = robust_category, values_from = n, values_fill = 0)
print(category_breakdown)
cat("\n")

# ==============================================================================
# LAMBDA CALCULATION
# ==============================================================================

cat("========================================================================\n")
cat("CALCULATING LAMBDA & LAMBDA_1000\n")
cat("========================================================================\n\n")

lambda_summary <- meta_results %>%
  group_by(maf_threshold, variant_type) %>%
  summarize(
    n_genes = n(),
    total_cases = max(total_cases),  
    total_controls = max(total_controls),
    n_eff = calculate_neff(max(total_cases), max(total_controls)),
    lambda_fisher = calculate_lambda(fisher_p),
    lambda_1000_fisher = calculate_lambda_1000(fisher_p, max(total_cases), max(total_controls)),
    lambda_stouffer = calculate_lambda_z(stouffer_z),  
    lambda_1000_stouffer = calculate_lambda_1000_z(stouffer_z, max(total_cases), max(total_controls)),
    .groups = "drop"
  )

cat("Lambda Summary:\n")
print(lambda_summary)

cat("\nInterpretation:\n")
for(i in 1:nrow(lambda_summary)) {
  maf <- lambda_summary$maf_threshold[i]
  vtype <- lambda_summary$variant_type[i]
  lf1000 <- lambda_summary$lambda_1000_fisher[i]
  ls1000 <- lambda_summary$lambda_1000_stouffer[i]
  
  cat(sprintf("MAF < %.2f (%s):\n", maf, vtype))
  cat(sprintf("  Fisher λ₁₀₀₀ = %.3f", lf1000))
  if(lf1000 < 1.05) cat(" [OK]\n") else if(lf1000 <= 1.10) cat(" [WARNING]\n") else cat(" [FAIL]\n")
  cat(sprintf("  Stouffer λ₁₀₀₀ = %.3f", ls1000))
  if(ls1000 < 1.05) cat(" [OK]\n\n") else if(ls1000 <= 1.10) cat(" [WARNING]\n\n") else cat(" [FAIL]\n\n")
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("========================================================================\n")
cat("SUMMARY\n")
cat("========================================================================\n\n")

summary_table <- meta_results %>%
  group_by(maf_threshold, variant_type) %>%
  summarize(
    n_genes = n(),
    fisher_gwas = sum(fisher_gwas),
    fisher_bonf = sum(fisher_bonf_sig),
    fisher_fdr = sum(fisher_fdr_sig),
    stouffer_fdr = sum(stouffer_fdr_sig),
    robust_fdr = sum(robust_fdr),
    robust_bonf = sum(robust_bonf),
    robust_ultra = sum(robust_ultra),
    .groups = "drop"
  )

print(summary_table)

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

cat("\n========================================================================\n")
cat("SAVING RESULTS\n")
cat("========================================================================\n\n")

dir.create("results/step1_trans_ancestry", showWarnings = FALSE, recursive = TRUE)

fwrite(meta_results, "results/step1_trans_ancestry/trans_ancestry_COMPLETE.txt",
       sep = "\t", quote = FALSE)
cat("trans_ancestry_COMPLETE.txt\n")

fwrite(lambda_summary, "results/step1_trans_ancestry/lambda_summary.txt",
       sep = "\t", quote = FALSE)
cat("lambda_summary.txt\n")

fwrite(robust_summary, "results/step1_trans_ancestry/robust_methods_comparison.txt",
       sep = "\t", quote = FALSE)
cat("robust_methods_comparison.txt\n")

# Save robust gene lists
for(maf_val in maf_values) {
  maf_label <- maf_labels[which(maf_values == maf_val)]
  for(vtype in variant_types) {
    
    # FDR Fisher 
    fdr_fisher <- meta_results %>%
      filter(maf_threshold == maf_val, variant_type == vtype, fisher_fdr_sig) %>%
      arrange(fisher_p)
    if(nrow(fdr_fisher) > 0) {
      fwrite(fdr_fisher, sprintf("results/step1_trans_ancestry/FDR_FISHER_%s_%s.txt", maf_label, vtype),
             sep = "\t", quote = FALSE)
      cat(sprintf("FDR_FISHER_%s_%s.txt (%d genes)\n", maf_label, vtype, nrow(fdr_fisher)))
    }
    
    # FDR Stouffer
    fdr_stouffer <- meta_results %>%
      filter(maf_threshold == maf_val, variant_type == vtype, stouffer_fdr_sig) %>%
      arrange(stouffer_p)
    if(nrow(fdr_stouffer) > 0) {
      fwrite(fdr_stouffer, sprintf("results/step1_trans_ancestry/FDR_STOUFFER_%s_%s.txt", maf_label, vtype),
             sep = "\t", quote = FALSE)
      cat(sprintf("FDR_STOUFFER_%s_%s.txt (%d genes)\n", maf_label, vtype, nrow(fdr_stouffer)))
    }
    
    # Bonferroni Fisher
    bonf_fisher <- meta_results %>%
      filter(maf_threshold == maf_val, variant_type == vtype, fisher_bonf_sig) %>%
      arrange(fisher_p)
    if(nrow(bonf_fisher) > 0) {
      fwrite(bonf_fisher, sprintf("results/step1_trans_ancestry/BONF_FISHER_%s_%s.txt", maf_label, vtype),
             sep = "\t", quote = FALSE)
      cat(sprintf("BONF_FISHER_%s_%s.txt (%d genes)\n", maf_label, vtype, nrow(bonf_fisher)))
    }
    
    # Bonferroni Stouffer
    bonf_stouffer <- meta_results %>%
      filter(maf_threshold == maf_val, variant_type == vtype, stouffer_bonf_sig) %>%
      arrange(stouffer_p)
    if(nrow(bonf_stouffer) > 0) {
      fwrite(bonf_stouffer, sprintf("results/step1_trans_ancestry/BONF_STOUFFER_%s_%s.txt", maf_label, vtype),
             sep = "\t", quote = FALSE)
      cat(sprintf("BONF_STOUFFER_%s_%s.txt (%d genes)\n", maf_label, vtype, nrow(bonf_stouffer)))
    }
    
    # FDR robust
    fdr_robust <- meta_results %>%
      filter(maf_threshold == maf_val, variant_type == vtype, robust_fdr) %>%
      arrange(fisher_p)
    if(nrow(fdr_robust) > 0) {
      fwrite(fdr_robust, 
             sprintf("results/step1_trans_ancestry/ROBUST_FDR_%s_%s.txt", maf_label, vtype),
             sep = "\t", quote = FALSE)
      cat(sprintf("ROBUST_FDR_%s_%s.txt (%d genes)\n", maf_label, vtype, nrow(fdr_robust)))
    }
    
    # Bonferroni robust
    bonf_robust <- meta_results %>%
      filter(maf_threshold == maf_val, variant_type == vtype, robust_bonf) %>%
      arrange(fisher_p)
    if(nrow(bonf_robust) > 0) {
      fwrite(bonf_robust, 
             sprintf("results/step1_trans_ancestry/ROBUST_BONF_%s_%s.txt", maf_label, vtype),
             sep = "\t", quote = FALSE)
      cat(sprintf("ROBUST_BONF_%s_%s.txt (%d genes)\n", maf_label, vtype, nrow(bonf_robust)))
    }
    
    # Ultra-robust
    ultra_robust <- meta_results %>%
      filter(maf_threshold == maf_val, variant_type == vtype, robust_ultra) %>%
      arrange(fisher_p)
    if(nrow(ultra_robust) > 0) {
      fwrite(ultra_robust, 
             sprintf("results/step1_trans_ancestry/ROBUST_ULTRA_%s_%s.txt", maf_label, vtype),
             sep = "\t", quote = FALSE)
      cat(sprintf("ROBUST_ULTRA_%s_%s.txt (%d genes)\n", maf_label, vtype, nrow(ultra_robust)))
    }
  }
}

