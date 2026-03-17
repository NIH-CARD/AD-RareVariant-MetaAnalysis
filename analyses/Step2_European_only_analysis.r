################################################################################
# STEP 2: EUR-ONLY META-ANALYSIS
################################################################################
#
# Description:
#   Performs EUR-only meta-analysis combining UKB EUR and ADSP EUR cohorts
#   using Fisher's and Stouffer's methods. Only analyzes genes present in
#   BOTH EUR cohorts for maximum comparability.
#
# Methods:
#   - Fisher's combined probability test
#   - Stouffer's Z-score method (weighted by sqrt(N_effective))
#   - Multiple testing correction: Bonferroni and FDR (Benjamini-Hochberg)
#   - Genomic inflation factor (λ and λ₁₀₀₀) calculation
#   - Enhanced robust classification (FDR, Bonferroni, Ultra)
#
# Input:
#   - REGENIE gene burden test results from UKB EUR and ADSP EUR
#   - Multiple MAF thresholds (0.01, 0.03, 0.05)
#   - Variant types: all variants and coding variants only
#
# Output:
#   - eur_only_COMPLETE.txt: Complete EUR-only meta-analysis results
#   - lambda_summary_EUR.txt: Genomic inflation factors
#   - robust_methods_comparison_EUR.txt: Method comparison summary
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
  cohort = c("UKB", "ADSP"),
  ancestry = c("EUR", "EUR"),
  n_cases = c(4394, 6391),
  n_controls = c(54816, 8417)
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

cat("EUR-only sample sizes:\n")
print(sample_info)
cat(sprintf("\nTotal EUR: %d cases, %d controls\n\n",
            sum(sample_info$n_cases), sum(sample_info$n_controls)))

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
# LOAD EUR DATA
# ==============================================================================

cat("========================================================================\n")
cat("Loading EUR Results\n")
cat("========================================================================\n\n")

all_results <- list()
counter <- 1

for(maf_idx in 1:length(maf_codes)) {
  maf_adsp <- maf_codes[maf_idx]
  maf_ukb <- maf_codes_ukb[maf_idx]
  maf_val <- maf_values[maf_idx]
  
  for(vtype in variant_types) {
    cat(sprintf("\n=== MAF %s, %s ===\n", maf_val, vtype))
    
    # ADSP EUR
    adsp_file <- sprintf("%s/new_results_step2_EUR_%s_%s_SKATO_meta.regenie.gz",
                        ADSP_DIR, maf_adsp, vtype)
    if(file.exists(adsp_file)) {
      cat(sprintf("  [%d] %s\n", counter, basename(adsp_file)))
      df <- fread(adsp_file)
      df$GENE <- gsub("\\.M1\\..*", "", df$ID)
      df$PVAL <- 10^(-df$LOG10P)
      df$PVAL[df$PVAL == 0 | is.na(df$PVAL)] <- 1e-300
      df$PVAL[df$PVAL > 1] <- 1
      df$cohort <- "ADSP"
      df$ancestry <- "EUR"
      df$cohort_ancestry <- "ADSP_EUR"
      df$maf_threshold <- maf_val
      df$variant_type <- vtype
      df$n_cases <- 6391
      df$n_controls <- 8417
      df$n_total <- 14808
      all_results[[counter]] <- df
      counter <- counter + 1
    }
    
    # UKB EUR
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
      df$n_cases <- 4394
      df$n_controls <- 54816
      df$n_total <- 59210
      all_results[[counter]] <- df
      counter <- counter + 1
    }
  }
}

combined_eur <- rbindlist(all_results, fill=TRUE)
cat(sprintf("\nLoaded %d files, %d rows\n\n", length(all_results), nrow(combined_eur)))

combined_eur <- combined_eur %>%
  filter(!is.na(PVAL), PVAL > 0, PVAL <= 1, !is.na(GENE))

# ==============================================================================
# META-ANALYSIS
# ==============================================================================

cat("Running EUR-only meta-analysis...\n")
cat("Note: Stouffer's method weighted by sqrt(Neff) per cohort for proper case-control adjustment\n\n")

eur_meta <- combined_eur %>%
  group_by(GENE, maf_threshold, variant_type) %>%
  summarize(
    n_cohorts = n_distinct(cohort_ancestry),  
    total_cases = sum(n_cases),
    total_controls = sum(n_controls),
    total_n = sum(n_total),
    ukb_p = PVAL[cohort == "UKB"][1],
    adsp_p = PVAL[cohort == "ADSP"][1],
    min_p = min(PVAL),
    max_p = max(PVAL),
    fisher_result = list(fishers_method(PVAL)),
    stouffer_result = list(stouffers_method(PVAL, sqrt(calculate_neff(n_cases, n_controls)))),
    .groups = "drop"
  ) %>%
  filter(n_cohorts == 2) %>%  
  mutate(
    fisher_p = sapply(fisher_result, function(x) x$p),
    fisher_log10p = -log10(fisher_p),
    stouffer_p = sapply(stouffer_result, function(x) x$p),
    stouffer_log10p = -log10(stouffer_p),
    stouffer_z = sapply(stouffer_result, function(x) x$z),  
    ukb_log10p = -log10(ukb_p),
    adsp_log10p = -log10(adsp_p)
  ) %>%
  select(-fisher_result, -stouffer_result)

cat(sprintf("Complete: %d genes (filtered to BOTH EUR cohorts)\n\n", nrow(eur_meta)))

# ==============================================================================
# MULTIPLE TESTING CORRECTION
# ==============================================================================

cat("Applying corrections...\n")

eur_meta <- eur_meta %>%
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
    both_cohorts_sig = ukb_p < 0.05 & adsp_p < 0.05,
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

eur_meta <- eur_meta %>%
  mutate(
    # both methods FDR significant
    robust_fdr = fisher_fdr_sig & stouffer_fdr_sig,
    
    # both methods Bonferroni significant
    robust_bonf = fisher_bonf_sig & stouffer_bonf_sig,
    
    # most stringent - both methods pass BOTH corrections
    robust_ultra = (fisher_fdr_sig & fisher_bonf_sig) & 
                   (stouffer_fdr_sig & stouffer_bonf_sig),
    
    # Classification label
    robust_category = case_when(
      robust_ultra ~ "Ultra-robust (both corrections)",
      robust_fdr & !robust_bonf ~ "Robust (FDR only)",
      robust_bonf & !robust_fdr ~ "Robust (Bonferroni only)",
      TRUE ~ "Not robust"
    )
  ) %>%
  arrange(maf_threshold, variant_type, fisher_p)

# Summary
robust_summary <- eur_meta %>%
  group_by(maf_threshold, variant_type) %>%
  summarize(
    n_genes = n(),
    robust_fdr = sum(robust_fdr, na.rm = TRUE),
    robust_bonf = sum(robust_bonf, na.rm = TRUE),
    robust_ultra = sum(robust_ultra, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nRobust gene summary (EUR-only):\n")
print(robust_summary)

cat("\nDetailed breakdown by category:\n")
category_breakdown <- eur_meta %>%
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

lambda_summary <- eur_meta %>%
  group_by(maf_threshold, variant_type) %>%
  summarize(
    n_genes = n(),
    total_cases = sum(unique(total_cases)),  
    total_controls = sum(unique(total_controls)),
    n_eff = calculate_neff(sum(unique(total_cases)), sum(unique(total_controls))),
    lambda_fisher = calculate_lambda(fisher_p),
    lambda_1000_fisher = calculate_lambda_1000(fisher_p, sum(unique(total_cases)), sum(unique(total_controls))),
    lambda_stouffer = calculate_lambda_z(stouffer_z),  
    lambda_1000_stouffer = calculate_lambda_1000_z(stouffer_z, sum(unique(total_cases)), sum(unique(total_controls))),
    .groups = "drop"
  )

cat("Lambda Summary (EUR):\n")
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
cat("EUR SUMMARY\n")
cat("========================================================================\n\n")

summary_eur <- eur_meta %>%
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
    both_cohorts = sum(both_cohorts_sig),
    .groups = "drop"
  )

print(summary_eur)

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

cat("\n========================================================================\n")
cat("SAVING RESULTS\n")
cat("========================================================================\n\n")

dir.create("results/step2_eur_only", showWarnings = FALSE, recursive = TRUE)

fwrite(eur_meta, "results/step2_eur_only/eur_only_COMPLETE.txt",
       sep = "\t", quote = FALSE)
cat("eur_only_COMPLETE.txt\n")

fwrite(lambda_summary, "results/step2_eur_only/lambda_summary_EUR.txt",
       sep = "\t", quote = FALSE)
cat("lambda_summary_EUR.txt\n")

fwrite(robust_summary, "results/step2_eur_only/robust_methods_comparison_EUR.txt",
       sep = "\t", quote = FALSE)
cat("robust_methods_comparison_EUR.txt\n")

# Save robust lists
for(maf_val in maf_values) {
  maf_label <- maf_labels[which(maf_values == maf_val)]
  for(vtype in variant_types) {
    
    # FDR Fisher
    fdr_fisher <- eur_meta %>%
      filter(maf_threshold == maf_val, variant_type == vtype, fisher_fdr_sig) %>%
      arrange(fisher_p)
    if(nrow(fdr_fisher) > 0) {
      fwrite(fdr_fisher, sprintf("results/step2_eur_only/FDR_FISHER_EUR_%s_%s.txt", maf_label, vtype),
             sep = "\t", quote = FALSE)
      cat(sprintf("FDR_FISHER_EUR_%s_%s.txt (%d genes)\n", maf_label, vtype, nrow(fdr_fisher)))
    }
    
    # FDR Stouffer
    fdr_stouffer <- eur_meta %>%
      filter(maf_threshold == maf_val, variant_type == vtype, stouffer_fdr_sig) %>%
      arrange(stouffer_p)
    if(nrow(fdr_stouffer) > 0) {
      fwrite(fdr_stouffer, sprintf("results/step2_eur_only/FDR_STOUFFER_EUR_%s_%s.txt", maf_label, vtype),
             sep = "\t", quote = FALSE)
      cat(sprintf("FDR_STOUFFER_EUR_%s_%s.txt (%d genes)\n", maf_label, vtype, nrow(fdr_stouffer)))
    }
    
    # Bonferroni Fisher
    bonf_fisher <- eur_meta %>%
      filter(maf_threshold == maf_val, variant_type == vtype, fisher_bonf_sig) %>%
      arrange(fisher_p)
    if(nrow(bonf_fisher) > 0) {
      fwrite(bonf_fisher, sprintf("results/step2_eur_only/BONF_FISHER_EUR_%s_%s.txt", maf_label, vtype),
             sep = "\t", quote = FALSE)
      cat(sprintf("BONF_FISHER_EUR_%s_%s.txt (%d genes)\n", maf_label, vtype, nrow(bonf_fisher)))
    }
    
    # Bonferroni Stouffer
    bonf_stouffer <- eur_meta %>%
      filter(maf_threshold == maf_val, variant_type == vtype, stouffer_bonf_sig) %>%
      arrange(stouffer_p)
    if(nrow(bonf_stouffer) > 0) {
      fwrite(bonf_stouffer, sprintf("results/step2_eur_only/BONF_STOUFFER_EUR_%s_%s.txt", maf_label, vtype),
             sep = "\t", quote = FALSE)
      cat(sprintf("BONF_STOUFFER_EUR_%s_%s.txt (%d genes)\n", maf_label, vtype, nrow(bonf_stouffer)))
    }
    
    # FDR robust
    fdr <- eur_meta %>%
      filter(maf_threshold == maf_val, variant_type == vtype, robust_fdr) %>%
      arrange(fisher_p)
    if(nrow(fdr) > 0) {
      fwrite(fdr, sprintf("results/step2_eur_only/ROBUST_FDR_EUR_%s_%s.txt", maf_label, vtype),
             sep = "\t", quote = FALSE)
      cat(sprintf("ROBUST_FDR_EUR_%s_%s.txt (%d genes)\n", maf_label, vtype, nrow(fdr)))
    }
    
    # Bonferroni robust
    bonf <- eur_meta %>%
      filter(maf_threshold == maf_val, variant_type == vtype, robust_bonf) %>%
      arrange(fisher_p)
    if(nrow(bonf) > 0) {
      fwrite(bonf, sprintf("results/step2_eur_only/ROBUST_BONF_EUR_%s_%s.txt", maf_label, vtype),
             sep = "\t", quote = FALSE)
      cat(sprintf("ROBUST_BONF_EUR_%s_%s.txt (%d genes)\n", maf_label, vtype, nrow(bonf)))
    }
    
    # Ultra-robust
    ultra <- eur_meta %>%
      filter(maf_threshold == maf_val, variant_type == vtype, robust_ultra) %>%
      arrange(fisher_p)
    if(nrow(ultra) > 0) {
      fwrite(ultra, sprintf("results/step2_eur_only/ROBUST_ULTRA_EUR_%s_%s.txt", maf_label, vtype),
             sep = "\t", quote = FALSE)
      cat(sprintf("ROBUST_ULTRA_EUR_%s_%s.txt (%d genes)\n", maf_label, vtype, nrow(ultra)))
    }
  }
}
