# Generate Figure 5: iPTM Module and Pan-Cancer Prognosis
library(survival)
library(survminer)
library(tidyverse)

setwd("/home/zm3028/research_platform/data_5/4_tcga_survival")
combined_surv <- readRDS("data/pan_cancer_survival.rds")

# Figure 5A: Pan-cancer KM curve
fit_pan <- survfit(Surv(OS_time, OS_status) ~ iPTM_group, data = combined_surv)
p_pan <- ggsurvplot(
    fit_pan, data = combined_surv,
    pval = TRUE, risk.table = FALSE,
    title = "Pan-Cancer: iPTM Module and Overall Survival (n = 292)",
    xlab = "Time (days)", ylab = "Overall Survival Probability",
    palette = c("#2C3E50", "#E74C3C"),
    legend.title = "iPTM Group", legend.labs = c("Low", "High"))
ggsave("results/figures/Figure5A_pan_cancer_survival.pdf", 
       print(p_pan$plot), width = 5, height = 5)

# Figure 5B: Forest plot
cancer_types <- unique(combined_surv$Cancer)
hr_data <- data.frame()
for (cancer in cancer_types) {
    sub_data <- combined_surv[combined_surv$Cancer == cancer, ]
    if (nrow(sub_data) < 15) next
    cox_fit <- coxph(Surv(OS_time, OS_status) ~ iPTM_score, data = sub_data)
    cox_sum <- summary(cox_fit)
    hr_data <- rbind(hr_data, data.frame(
        Cancer = cancer, HR = cox_sum$coefficients[2],
        lower = cox_sum$conf.int[3], upper = cox_sum$conf.int[4],
        p = cox_sum$coefficients[5]))
}

cox_pan <- coxph(Surv(OS_time, OS_status) ~ iPTM_score + Cancer, 
                 data = combined_surv)
cox_pan_sum <- summary(cox_pan)
hr_data <- rbind(hr_data, data.frame(
    Cancer = "Pan-Cancer (adjusted)",
    HR = cox_pan_sum$coefficients[1, 2],
    lower = cox_pan_sum$conf.int[1, 3],
    upper = cox_pan_sum$conf.int[1, 4],
    p = cox_pan_sum$coefficients[1, 5]))

write.csv(hr_data, "results/tables/Figure5_HR_data.csv", row.names = FALSE)
