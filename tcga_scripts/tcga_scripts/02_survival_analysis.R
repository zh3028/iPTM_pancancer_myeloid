# iPTM Module Survival Analysis Script
library(SummarizedExperiment)
library(survival)
library(survminer)
library(tidyverse)

setwd("/home/zm3028/research_platform/data_5/4_tcga_survival")

iptm_genes <- c("LDHB", "SLC16A1", "HDAC1", "HDAC2", "HDAC3",
                "SLC16A3", "EP300", "SIRT2", "HNMT")

cancer_files <- list.files("data", pattern = "_data.rds$", full.names = TRUE)
cancer_files <- cancer_files[!grepl("all_survival|pan_cancer", cancer_files)]

all_surv <- list()

for (f in cancer_files) {
    project <- gsub("_data.rds", "", basename(f))
    data <- readRDS(f)
    expr <- assay(data, "tpm_unstrand")
    rowData_info <- rowData(data)
    iptm_idx <- which(rowData_info$gene_name %in% iptm_genes)
    if (length(iptm_idx) < 5) next
    iptm_expr <- expr[iptm_idx, , drop = FALSE]
    rownames(iptm_expr) <- rowData_info$gene_name[iptm_idx]
    iptm_expr_log <- log2(iptm_expr + 1)
    iptm_zscore <- t(scale(t(iptm_expr_log)))
    iptm_score <- colMeans(iptm_zscore, na.rm = TRUE)
    clinical <- colData(data)
    surv_df <- data.frame(
        sample = colnames(expr), Cancer = project,
        iPTM_score = iptm_score,
        OS_time = clinical$days_to_death,
        OS_status = ifelse(clinical$vital_status == "Dead", 1, 0),
        stringsAsFactors = FALSE)
    surv_df <- surv_df[!is.na(surv_df$OS_time) & surv_df$OS_time > 0, ]
    if (nrow(surv_df) > 10) all_surv[[project]] <- surv_df
}

combined_surv <- do.call(rbind, all_surv)
combined_surv <- combined_surv %>%
    group_by(Cancer) %>%
    mutate(iPTM_group = ifelse(iPTM_score > median(iPTM_score), "High", "Low")) %>%
    ungroup()
combined_surv$iPTM_group <- factor(combined_surv$iPTM_group, levels = c("Low", "High"))

fit <- survfit(Surv(OS_time, OS_status) ~ iPTM_group, data = combined_surv)
logrank <- survdiff(Surv(OS_time, OS_status) ~ iPTM_group, data = combined_surv)
p_val <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)
cox_fit <- coxph(Surv(OS_time, OS_status) ~ iPTM_score + Cancer, data = combined_surv)

saveRDS(combined_surv, "data/pan_cancer_survival.rds")
cat("Log-rank p =", round(p_val, 4), "\n")
