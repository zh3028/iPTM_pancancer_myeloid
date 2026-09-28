# ============================================================
# TCGA Data Download Script
# ============================================================

library(TCGAbiolinks)

setwd("/home/zm3028/research_platform/data_5/4_tcga_survival")

target_cancers <- c("TCGA-ESCA", "TCGA-PAAD", "TCGA-THCA", "TCGA-UCEC")

for (project in target_cancers) {
    out_file <- paste0("data/", project, "_data.rds")
    if (file.exists(out_file)) {
        cat("✅ Already exists, skipping:", project, "\n")
        next
    }
    
    tryCatch({
        query <- GDCquery(
            project = project,
            data.category = "Transcriptome Profiling",
            data.type = "Gene Expression Quantification",
            workflow.type = "STAR - Counts"
        )
        GDCdownload(query, directory = "data/GDCdata", method = "api", 
                    files.per.chunk = 50)
        data <- GDCprepare(query, directory = "data/GDCdata")
        saveRDS(data, file = out_file)
        cat("✅ Downloaded:", project, "\n")
    }, error = function(e) {
        cat("❌ Failed:", project, ":", e$message, "\n")
    })
}
